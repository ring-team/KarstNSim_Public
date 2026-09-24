#!/usr/bin/env python3
"""Exercise real generation, independent packages and failure ownership."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary = str(Path(sys.argv[1]).resolve())
map_id = "7e1f4efb-12fb-4c84-bf46-af4c0111e85d"
with tempfile.TemporaryDirectory(prefix="caves cli ") as temporary:
    root = Path(temporary)
    request = root / "request.json"
    request.write_text(json.dumps({"world_seed": 42, "loop_count": 0}))
    output = root / "generated package"

    def run(*extra):
        return subprocess.run([binary, "--request", str(request), "--out", str(output),
                               "--map-id", map_id, *extra], capture_output=True, text=True,
                              timeout=120)

    result = run("--work-limit", "1")
    assert result.returncode != 0 and not output.exists(), result
    request.write_text('{"unknown_field": 1}')
    result = run()
    assert result.returncode != 0 and not output.exists(), result
    request.write_text(json.dumps({"world_seed": 42, "loop_count": 0}))
    result = run()
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["complete"] and not manifest["database_migration_required"]
    assert manifest["map_id"] == map_id and len(manifest["regions"]) == 1
    snapshot = {p.name: p.read_bytes() for p in output.iterdir()}
    assert "dataset.csv" in snapshot and "region-0.json" in snapshot
    result = run()
    assert result.returncode != 0
    assert snapshot == {p.name: p.read_bytes() for p in output.iterdir()}
print("CLI generated a fresh package and preserved existing output on failure")
