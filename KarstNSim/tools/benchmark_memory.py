#!/usr/bin/env python3
"""Measure the unmodified full-resolution scientific base example on Linux.

The domain, geological inputs, sampling resolution and 100-neighbour graph
remain unchanged. Budgets are whole-process peak RSS in decimal MB, not just
graph-array sizes.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(work, assets):
    work.mkdir(parents=True, exist_ok=False)
    (work / "Input_files").symlink_to(assets, target_is_directory=True)
    base = assets / "1_base"
    text = (base / "instructions.txt").read_text()
    text = re.sub(r"(?m)^main_repository:.*$", "main_repository: .", text)
    (work / "instructions.txt").write_text(text)
    shutil.copy2(base / "connectivity_matrix.txt", work / "connectivity_matrix.txt")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path, help="new output directory; existing directories are refused")
    parser.add_argument("--budget-mb", type=float, default=1000, help="default: 1000; 0 measures a baseline without enforcing a budget")
    parser.add_argument("--compare", type=Path, help="prior metrics.json whose exports must match")
    args = parser.parse_args()
    binary, work = args.binary.resolve(strict=True), args.out.resolve()
    assets = Path(__file__).resolve().parents[2] / "Input_files"
    budget = args.budget_mb
    if budget < 0 or not (budget < float("inf")):
        parser.error("budget must be finite and nonnegative")
    if not Path("/usr/bin/time").is_file():
        parser.error("GNU /usr/bin/time is required for process RSS measurement")
    prepare(work, assets)
    started = time.monotonic()
    with (work / "run.log").open("w") as log:
        result = subprocess.run(["/usr/bin/time", "-v", "-o", str(work / "time.txt"),
            str(binary), "instructions.txt"], cwd=work, stdout=log, stderr=subprocess.STDOUT)
    wall = time.monotonic() - started
    report = (work / "time.txt").read_text()
    match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", report)
    rss = int(match.group(1)) * 1024 if match else None
    exports = {p.name: digest(p) for p in sorted((work / "outputs").glob("*"))
               if p.is_file() and not p.name.startswith("karstnsim_console_")}
    points = work / "outputs/base_0_pts.txt"
    point_count = sum(1 for _ in points.open()) - 1 if points.exists() else None
    same = None
    if args.compare:
        reference = json.loads(args.compare.read_text())
        if reference["case"] != "full":
            parser.error("comparison must use the full scientific base example")
        same = exports == reference["exports"]
    budget_ok = rss is not None and (budget == 0 or rss <= budget * 1_000_000)
    metrics = dict(case="full", binary=str(binary), binary_sha256=digest(binary),
        exit_code=result.returncode, wall_seconds=wall, peak_rss_bytes=rss,
        budget_mb=budget, budget_pass=budget_ok, exports_match=same,
        support_points=point_count, exports=exports)
    (work / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))
    return 0 if result.returncode == 0 and exports and budget_ok and same is not False else 1


if __name__ == "__main__":
    raise SystemExit(main())
