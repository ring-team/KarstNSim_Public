#!/usr/bin/env python3
"""Unmodified-MOOCoW acceptance run for fbs::caves::to_moocow_csv.

Everything happens in a fresh temporary directory:
  1. adapter CSV fixtures are produced by the native adapter test
     (--adapter-test BIN, or --standalone which compiles the adapter with the
     test-only region stub, or --fixtures DIR with existing *.csv files);
  2. MOOCoW's Domain + Infrastructure sources are copied with `git archive`
     from a commit of the MOOCoW checkout (read-only; nothing is written to
     the checkout and its global.json is not used);
  3. tests/adapter_acceptance.{cs,csproj} are copied next to it and built with
     the local .NET SDK, restoring only from the local NuGet package cache;
  4. the consumer imports, validates, re-exports and stores every fixture in a
     throwaway SQLite file using MOOCoW's own code.

.NET is optional: without `dotnet` the script prints SKIP and exits 77 (the
CTest/automake "skipped" convention). No global tool or package is installed.
"""

import argparse
import hashlib
import io
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SDK = os.path.dirname(HERE)
DEFAULT_MOOCOW = "/home/max/Documents/Codex/projects/external/most-optimal-omni-crafter-of-worlds"
COPIED = ["Code/src/MOOCoW.Domain", "Code/src/MOOCoW.Infrastructure", "Directory.Build.props"]
KEY_FILES = [
    "Code/src/MOOCoW.Infrastructure/Persistence/Csv/CsvImport.cs",
    "Code/src/MOOCoW.Infrastructure/Persistence/Csv/CsvExport.cs",
    "Code/src/MOOCoW.Infrastructure/Persistence/Csv/CsvText.cs",
    "Code/src/MOOCoW.Domain/DomainValidation.cs",
    "Code/src/MOOCoW.Domain/CaveContracts.cs",
    "Code/src/MOOCoW.Infrastructure/Persistence/Migrations/InitialSchema.cs",
]


def run(cmd, **kw):
    print("+ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--adapter-test", help="built caves_adapter_test executable (real region library)")
    src.add_argument("--fixtures", help="directory of adapter *.csv fixtures")
    src.add_argument("--standalone", action="store_true", help="compile adapter + region stub with $CXX (default)")
    ap.add_argument("--moocow", default=os.environ.get("MOOCOW_CHECKOUT", DEFAULT_MOOCOW), help="MOOCoW git checkout")
    ap.add_argument("--ref", default="HEAD", help="MOOCoW commit to archive (default HEAD)")
    ap.add_argument("--keep", action="store_true", help="keep the temporary directory")
    args = ap.parse_args()

    dotnet = shutil.which("dotnet")
    if not dotnet:
        print("SKIP: dotnet not found; the native SDK does not need .NET")
        return 77
    if not os.path.isdir(os.path.join(args.moocow, ".git")):
        print(f"SKIP: MOOCoW checkout not found at {args.moocow}")
        return 77

    tmp = tempfile.mkdtemp(prefix="fbs-moocow-acceptance-")
    print(f"scratch: {tmp}")
    try:
        fixtures = os.path.join(tmp, "fixtures")
        os.makedirs(fixtures)
        if args.fixtures:
            for name in sorted(os.listdir(args.fixtures)):
                if name.endswith(".csv"):
                    shutil.copy(os.path.join(args.fixtures, name), fixtures)
        elif args.adapter_test:
            run([args.adapter_test, "--write-fixtures", fixtures])
        else:
            exe = os.path.join(tmp, "adapter_stub_test")
            cxx = os.environ.get("CXX", "g++")
            run([cxx, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-I" + os.path.join(SDK, "include"),
                 "-DFBS_ADAPTER_STUB_REGION", os.path.join(SDK, "src", "moocow_csv.cpp"),
                 os.path.join(HERE, "adapter_region_stub.cpp"), os.path.join(HERE, "adapter_tests.cpp"), "-o", exe])
            run([exe, "--write-fixtures", fixtures])
        for name in sorted(os.listdir(fixtures)):
            print(f"fixture {name} sha256 {sha256(os.path.join(fixtures, name))}")

        # Read-only copy of committed MOOCoW sources.
        commit = subprocess.run(["git", "-C", args.moocow, "rev-parse", args.ref], check=True,
                                capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", args.moocow, "status", "--porcelain", "--"] + COPIED, check=True,
                               capture_output=True, text=True).stdout.strip()
        archive = subprocess.run(["git", "-C", args.moocow, "archive", "--format=tar", commit] + COPIED,
                                 check=True, capture_output=True).stdout
        moocow_copy = os.path.join(tmp, "moocow")
        with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
            tar.extractall(moocow_copy)
        print(f"MOOCoW commit {commit} copied via git archive ({len(archive)} bytes)")
        print("MOOCoW working tree changes in copied paths (not used): " + (dirty.replace("\n", "; ") or "none"))
        for rel in KEY_FILES:
            print(f"  {rel} sha256 {sha256(os.path.join(moocow_copy, rel))}")

        consumer = os.path.join(tmp, "consumer")
        os.makedirs(consumer)
        shutil.copy(os.path.join(HERE, "adapter_acceptance.cs"), consumer)
        shutil.copy(os.path.join(HERE, "adapter_acceptance.csproj"), consumer)
        packages = os.environ.get("NUGET_PACKAGES", os.path.expanduser("~/.nuget/packages"))
        with open(os.path.join(tmp, "nuget.config"), "w") as f:
            f.write('<?xml version="1.0" encoding="utf-8"?>\n<configuration>\n  <packageSources>\n    <clear />\n'
                    f'    <add key="local-cache" value="{packages}" />\n  </packageSources>\n</configuration>\n')

        env = dict(os.environ, DOTNET_CLI_TELEMETRY_OPTOUT="1", DOTNET_NOLOGO="1",
                   DOTNET_SKIP_FIRST_TIME_EXPERIENCE="1", MSBUILDDISABLENODEREUSE="1")
        project = os.path.join(consumer, "adapter_acceptance.csproj")
        run([dotnet, "build", project, "-c", "Release", "--nologo", "-v", "q",
             "--configfile", os.path.join(tmp, "nuget.config"),
             "-p:MOOCoWCode=" + os.path.join(moocow_copy, "Code"), "-p:UseSharedCompilation=false"],
            env=env, cwd=tmp)
        work = os.path.join(tmp, "work")
        os.makedirs(work)
        dll = os.path.join(consumer, "bin", "Release", "net10.0", "adapter_acceptance.dll")
        result = subprocess.run([dotnet, dll, fixtures, work], env=env, cwd=tmp)
        return result.returncode
    finally:
        if args.keep:
            print(f"kept {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
