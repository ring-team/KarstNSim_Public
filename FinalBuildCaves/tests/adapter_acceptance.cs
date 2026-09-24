// Acceptance consumer for fbs::caves::to_moocow_csv output.
//
// Compiled ONLY against an unmodified copy of MOOCoW's own Domain and
// Infrastructure projects (git archive of a MOOCoW commit, see
// adapter_acceptance.py). Nothing here parses or validates CSV itself: every
// verdict comes from MOOCoW's CsvImport, CaveDatasetValidator, CsvExport,
// MigrationRunner and CaveDatasetRepository.
//
// usage: adapter_acceptance <fixtures-dir> <scratch-dir>
//   fixtures-dir  *.csv written by caves_adapter_test --write-fixtures
//   scratch-dir   empty temporary directory; a throwaway SQLite file is made here

using System.Security.Cryptography;
using System.Text;
using MOOCoW.Domain;
using MOOCoW.Infrastructure.Persistence;
using MOOCoW.Infrastructure.Persistence.Csv;
using MOOCoW.Infrastructure.Persistence.Migrations;
using MOOCoW.Infrastructure.Persistence.Repositories;
using MOOCoW.Infrastructure.Persistence.SeedData;

if (args.Length != 2)
{
    Console.Error.WriteLine("usage: adapter_acceptance <fixtures-dir> <scratch-dir>");
    return 2;
}

var fixturesDir = args[0];
var scratch = args[1];
var failures = 0;

void Check(bool condition, string message)
{
    Console.WriteLine((condition ? "PASS " : "FAIL ") + message);
    if (!condition) failures++;
}

string Sha(string text) => Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(text))).ToLowerInvariant();

string Query(Database database, string sql)
{
    using var command = database.Connection.CreateCommand();
    command.CommandText = sql;
    using var reader = command.ExecuteReader();
    var sb = new StringBuilder();
    while (reader.Read())
    {
        for (var i = 0; i < reader.FieldCount; i++)
            sb.Append(reader.IsDBNull(i) ? "<null>" : reader.GetValue(i).ToString()).Append('\u001f');
        sb.Append('\n');
    }
    return sb.ToString();
}

var fixtures = Directory.GetFiles(fixturesDir, "*.csv").OrderBy(p => p, StringComparer.Ordinal).ToArray();
Check(fixtures.Length > 0, $"found {fixtures.Length} fixture(s) in {fixturesDir}");

// A throwaway database that already holds MOOCoW's own seed map. It stands in
// for pre-existing records (no user data is read or written anywhere).
var dbPath = Path.Combine(scratch, "acceptance.db");
Check(!File.Exists(dbPath), "scratch database does not exist before the run");
using var database = Database.CreateFile(dbPath);
new MigrationRunner(database).EnsureSchema();
var repository = new CaveDatasetRepository(database);
var existing = SeedDataset.Create();
repository.SaveDataset(existing, "seed", "pre-existing records");
var existingCsv = CsvExport.ExportDataset(repository.LoadDataset(existing.Map.Id)!);
var existingRevision = repository.GetCurrentRevision(existing.Map.Id);
var schemaBefore = Query(database, "SELECT type, name, tbl_name, sql FROM sqlite_master ORDER BY type, name;");
var migrationsBefore = Query(database, "SELECT name, applied_at FROM schema_migrations ORDER BY name;");
Console.WriteLine($"INFO schema sha256 {Sha(schemaBefore)}; migrations: {migrationsBefore.Replace('\u001f', ' ').Trim()}");

var imported = new List<(string File, EntityId MapId, string Csv)>();
foreach (var path in fixtures)
{
    var name = Path.GetFileName(path);
    var text = File.ReadAllText(path, new UTF8Encoding(false, true));
    CaveDataset dataset;
    try
    {
        dataset = CsvImport.ImportDataset(text);
    }
    catch (Exception e)
    {
        Check(false, $"{name}: CsvImport.ImportDataset threw {e.GetType().Name}: {e.Message}");
        continue;
    }
    Check(true, $"{name}: CsvImport.ImportDataset parsed map {dataset.Map.Id} " +
                $"({dataset.Caverns.Count} caverns, {dataset.FloorRegions.Count} floors, {dataset.CliffEdges.Count} cliffs, " +
                $"{dataset.Ports.Count} ports, {dataset.Tunnels.Count} tunnels, {dataset.ElevationBands.Count} bands)");

    var issues = CaveDatasetValidator.Validate(dataset);
    foreach (var issue in issues)
        Console.WriteLine($"     {issue.Severity} {issue.Code}: {issue.Message}");
    Check(issues.Count == 0, $"{name}: CaveDatasetValidator.Validate reports {issues.Count} issue(s)");

    // Byte-identical re-export proves no value was trimmed, merged, normalised,
    // reordered or dropped by MOOCoW's parser and domain constructors.
    var reexport = CsvExport.ExportDataset(dataset);
    Check(reexport == text, $"{name}: CsvExport(CsvImport(csv)) is byte-identical to the adapter output");
    if (reexport != text)
    {
        var a = text.Split('\n'); var b = reexport.Split('\n');
        for (var i = 0; i < Math.Min(a.Length, b.Length); i++)
            if (a[i] != b[i]) { Console.WriteLine($"     first diff line {i + 1}:\n       adapter: {a[i]}\n       moocow:  {b[i]}"); break; }
    }

    // A generated package is a NEW dataset: its map ID and code are unused.
    var fresh = repository.GetCurrentRevision(dataset.Map.Id) == 0 &&
                new MapRepository(database).GetById(dataset.Map.Id) is null &&
                new MapRepository(database).GetByCode(dataset.Map.Code) is null;
    Check(fresh, $"{name}: map ID and code are not present in the target database");
    if (!fresh) continue;

    // Same repository path as `MOOCoW.Cli import` (validate, then transactional save).
    string? backup;
    try
    {
        backup = repository.ImportDatasetWithBackup(dataset, 0, "fbs-caves-acceptance", "generated package");
    }
    catch (Exception e)
    {
        Check(false, $"{name}: CaveDatasetRepository.ImportDatasetWithBackup threw {e.GetType().Name}: {e.Message}");
        continue;
    }
    Check(backup is null, $"{name}: import created revision 1 of a new map (no replacement backup needed)");
    Check(repository.GetCurrentRevision(dataset.Map.Id) == 1, $"{name}: stored at data revision 1");
    var storedHash = Query(database, $"SELECT content_hash FROM revisions WHERE map_id = '{dataset.Map.Id}' AND data_revision = 1;")
        .TrimEnd('\n').TrimEnd('\u001f');
    Check(storedHash == Sha(text + "\nrevision=1"),
          $"{name}: MOOCoW revision content_hash equals sha256(adapter bytes + \"\\nrevision=1\") = {storedHash}");
    var loaded = repository.LoadDataset(dataset.Map.Id)!;
    var stored = CsvExport.ExportDataset(loaded);
    Check(stored == reexport, $"{name}: SQLite round trip exports identical CSV");
    if (stored != reexport)
    {
        var a = reexport.Split('\n'); var b = stored.Split('\n');
        for (var i = 0; i < Math.Min(a.Length, b.Length); i++)
            if (a[i] != b[i]) { Console.WriteLine($"     first diff line {i + 1}:\n       adapter: {a[i]}\n       stored:  {b[i]}"); break; }
    }
    imported.Add((name, dataset.Map.Id, reexport));
}

// Re-importing the same records under another new map ID collides on MOOCoW's
// database-wide primary keys. The save is transactional: nothing changes.
var constructedFixture = Path.Combine(fixturesDir, "constructed.csv");
if (File.Exists(constructedFixture) && imported.Count > 0)
{
    var text = File.ReadAllText(constructedFixture);
    var original = CsvImport.ImportDataset(text);
    var otherMap = "0f0c7d2e-9a41-4b8e-8d3c-2e7f61a0b9ff";
    var clash = CsvImport.ImportDataset(text.Replace(original.Map.Id.ToString(), otherMap));
    var mapsBefore = Query(database, "SELECT id, code, data_revision FROM maps ORDER BY id;");
    var rejected = false;
    try { repository.ImportDatasetWithBackup(clash, 0, "fbs-caves-acceptance", "colliding record IDs"); }
    catch (Exception e) { rejected = true; Console.WriteLine($"     rejected: {e.GetType().Name}: {e.Message}"); }
    Check(rejected, "records whose IDs already exist under another map are rejected (IDs are database-wide keys)");
    Check(Query(database, "SELECT id, code, data_revision FROM maps ORDER BY id;") == mapsBefore &&
          repository.LoadDataset(EntityId.Parse(otherMap)) is null,
          "the rejected import rolled back completely (no map, no records)");
}

// Existing records and the schema are untouched by every import.
Check(CsvExport.ExportDataset(repository.LoadDataset(existing.Map.Id)!) == existingCsv,
      "pre-existing seed map exports byte-identically after all imports");
Check(repository.GetCurrentRevision(existing.Map.Id) == existingRevision, "pre-existing seed map revision unchanged");
var schemaAfter = Query(database, "SELECT type, name, tbl_name, sql FROM sqlite_master ORDER BY type, name;");
Check(schemaAfter == schemaBefore, $"sqlite_master unchanged (sha256 {Sha(schemaAfter)})");
Check(Query(database, "SELECT name, applied_at FROM schema_migrations ORDER BY name;") == migrationsBefore,
      "no migration was added or re-applied");
foreach (var (file, mapId, csv) in imported)
    Check(CsvExport.ExportDataset(repository.LoadDataset(mapId)!) == csv, $"{file}: still intact after later imports");

// Evidence for the adapter's refusals: what the unmodified importer does to
// values the adapter rejects instead of emitting.
var constructedPath = Path.Combine(fixturesDir, "constructed.csv");
if (File.Exists(constructedPath))
{
    var baseline = File.ReadAllText(constructedPath);

    var padded = CsvImport.ImportDataset(baseline.Replace(",Side pocket,", ",Side pocket ,"));
    Check(padded.Caverns.Any(c => c.Name == "Side pocket"),
          "probe: trailing whitespace in a name is silently trimmed by MOOCoW (adapter rejects it)");

    var dupTags = CsvImport.ImportDataset(baseline.Replace("\"[\"\"pocket\"\"]\"", "\"[\"\"pocket\"\",\"\"POCKET\"\"]\""));
    Check(dupTags.Caverns.Any(c => c.Tags.Count == 1 && c.Tags[0].Equals("pocket", StringComparison.OrdinalIgnoreCase)),
          "probe: case-insensitive duplicate tags are silently merged by MOOCoW (adapter rejects them)");

    var facingRow = baseline.Split('\n').First(l => l.Contains(",P-2,West mouth,"));
    var scaled = facingRow.Replace(",-1,0,0,", ",-2,0,0,");
    var rescaled = CsvImport.ImportDataset(baseline.Replace(facingRow, scaled));
    Check(rescaled.Ports.Any(p => p.Code == "P-2" && p.Facing.X == -1),
          "probe: a non-unit facing is rescaled by MOOCoW's UnitVector3 (adapter rejects it)");

    // A generation parameter literally named "parameters": MOOCoW's own
    // domain and exporter accept it, but its import/repository readers unwrap a
    // root "parameters" property. Stored in a separate throwaway database.
    var tunnelLine = baseline.Split('\n').First(l => l.Contains(",T-1,"));
    var routeJson = "\"{\"\"route\"\":\"\"hand, \\u0022made\\u0022\"\"}\"";
    Check(tunnelLine.Contains(routeJson), "probe setup: tunnel row generation_json located");
    var unwrapped = tunnelLine.Replace(routeJson, "\"{\"\"parameters\"\":\"\"literal\"\",\"\"route\"\":\"\"x\"\"}\"");
    var unwrappedImports = true;
    try { CsvImport.ImportDataset(baseline.Replace(tunnelLine, unwrapped)); } catch (InvalidOperationException) { unwrappedImports = false; }
    Check(!unwrappedImports, "probe: CsvImport rejects generation_json {\"parameters\":\"literal\",...}");
    var wrapped = tunnelLine.Replace(routeJson,
        "\"{\"\"parameters\"\":{\"\"parameters\"\":\"\"literal\"\",\"\"route\"\":\"\"x\"\"}}\"");
    var keyed = CsvImport.ImportDataset(baseline.Replace(tunnelLine, wrapped));
    Check(keyed.Tunnels.Any(t => t.Generation.Parameters.TryGetValue("parameters", out var v) && v == "literal"),
          "probe: the legacy wrapped form does import a 'parameters' key");
    var keyedExportImports = true;
    try { CsvImport.ImportDataset(CsvExport.ExportDataset(keyed)); } catch (InvalidOperationException) { keyedExportImports = false; }
    Check(!keyedExportImports, "probe: MOOCoW's own CsvExport output for that record cannot be re-imported");
    using (var probeDb = Database.CreateFile(Path.Combine(scratch, "parameters-probe.db")))
    {
        new MigrationRunner(probeDb).EnsureSchema();
        var probeRepo = new CaveDatasetRepository(probeDb);
        var keyedIssues = CaveDatasetValidator.Validate(keyed).Count;
        probeRepo.SaveDataset(keyed, "probe", "parameters key");
        var reloads = true;
        try { probeRepo.LoadDataset(keyed.Map.Id); } catch (InvalidOperationException) { reloads = false; }
        Check(keyedIssues == 0 && !reloads,
              "probe: a 'parameters' key passes validation and saves, but the map can no longer be loaded (adapter rejects it)");
    }

    var brokenOk = false;
    try
    {
        var broken = CsvImport.ImportDataset(baseline.Replace(",Side pocket,", ",\"Side\npocket\","));
        brokenOk = broken.Caverns.Any(c => c.Name == "Side\npocket");
    }
    catch (Exception) { }
    Check(!brokenOk, "probe: a quoted line break cannot be imported by MOOCoW (adapter rejects it)");
}

Console.WriteLine(failures == 0 ? "ACCEPTANCE PASSED" : $"ACCEPTANCE FAILED ({failures})");
return failures == 0 ? 0 : 1;
