#!/usr/bin/env python3
"""Run the reduced KarstNSim regression corpus against a binary and compare exports.

Typical use:

  # record golden hashes from the upstream baseline binary
  run_regression.py --binary /path/to/baseline/karstnsim --out /tmp/rg-base \
      --write-golden golden/baseline_ff55d39.json

  # check a candidate binary against the recorded golden hashes
  run_regression.py --binary /path/to/optimized/karstnsim --out /tmp/rg-opt \
      --golden golden/baseline_ff55d39.json

Each case in cases/<name>/ holds an instructions.txt (plus connectivity_matrix.txt when
the case uses a user matrix). For every case, a fresh work root <out>/<name>/ is created
containing a copy of those files and three symlinks:

  assets   -> repository Input_files/ (original, unmodified example assets)
  fixtures -> tests/regression/fixtures/ (small committed fixture files)
  gen      -> <out>/_gen/ (derived reduced inputs from make_fixtures.py)

The binary runs with cwd set to the work root and "main_repository: .", so outputs land
in <work>/outputs/. Paths never contain the (space-containing) repository path, because
KarstNSim reads path parameters as single whitespace-delimited tokens.

Comparison: every file in outputs/ is compared byte-for-byte by SHA-256, except the
console log (karstnsim_console_<timestamp>.log), which contains timing and absolute
paths. simulation_times.txt is written to the work root, not outputs/, and is ignored.
None of the compared exports (<name>_<seed>_karst.txt, _pts.txt, _s_pts.txt,
_connectivity_matrix.txt, _nghb_graph.txt) contain timestamps or timing values.
Exit codes and the file set are compared as well. On mismatch a line-level summary and,
for tabular numeric files, the maximum absolute numeric difference are reported.

These are reduced-density fixtures (see README.md); they are NOT the original full-size
examples, and golden hashes only apply to the reduced configurations.
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

CONSOLE_LOG_PREFIX = "karstnsim_console_"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def list_cases(selected):
    root = os.path.join(HERE, "cases")
    names = sorted(d for d in os.listdir(root) if os.path.isfile(os.path.join(root, d, "instructions.txt")))
    if selected:
        unknown = set(selected) - set(names)
        if unknown:
            raise SystemExit("unknown case(s): %s" % ", ".join(sorted(unknown)))
        names = [n for n in names if n in selected]
    return names


def prepare_work(case, out, assets, gen_dir):
    work = os.path.join(out, case)
    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(work)
    src = os.path.join(HERE, "cases", case)
    for f in sorted(os.listdir(src)):
        shutil.copy2(os.path.join(src, f), os.path.join(work, f))
    os.symlink(os.path.abspath(assets), os.path.join(work, "assets"))
    os.symlink(os.path.join(HERE, "fixtures"), os.path.join(work, "fixtures"))
    os.symlink(os.path.abspath(gen_dir), os.path.join(work, "gen"))
    return work


GNU_TIME = "/usr/bin/time"


def run_case(binary, work, timeout):
    stdout_path = os.path.join(work, "stdout.txt")
    time_path = os.path.join(work, "time.txt")
    cmd = [binary, "instructions.txt"]
    if os.access(GNU_TIME, os.X_OK):
        # GNU time execs the binary from a small process, giving a clean peak RSS.
        cmd = [GNU_TIME, "-v", "-o", time_path] + cmd
    t0 = time.monotonic()
    with open(stdout_path, "wb") as log:
        proc = subprocess.Popen(cmd, cwd=work, stdin=subprocess.DEVNULL,
                                stdout=log, stderr=subprocess.STDOUT)
        deadline = t0 + timeout
        status, rusage = 0, None
        while True:
            pid, status, rusage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if time.monotonic() > deadline:
                proc.kill()
                pid, status, rusage = os.wait4(proc.pid, 0)
                status = None
                break
            time.sleep(0.05)
        proc.returncode = 0
    wall = time.monotonic() - t0
    if status is None:
        rc = "timeout"
    elif os.WIFSIGNALED(status):
        rc = -os.WTERMSIG(status)
    else:
        rc = os.WEXITSTATUS(status)
    info = {
        "exit_code": rc,
        "wall_s": round(wall, 3),
        "user_s": round(rusage.ru_utime, 3),
        "sys_s": round(rusage.ru_stime, 3),
        "max_rss_kb": None,
    }
    if os.path.exists(time_path):
        with open(time_path) as f:
            for line in f:
                if "Maximum resident set size" in line:
                    info["max_rss_kb"] = int(line.rsplit(":", 1)[1])
                elif rc != "timeout" and "Command terminated by signal" in line:
                    # GNU time exits 128+N when the binary dies by signal N; report -N.
                    info["exit_code"] = -int(line.rsplit(" ", 1)[1])
    return info


def hash_outputs(work):
    outdir = os.path.join(work, "outputs")
    files = {}
    if os.path.isdir(outdir):
        for f in sorted(os.listdir(outdir)):
            if f.startswith(CONSOLE_LOG_PREFIX):
                continue
            p = os.path.join(outdir, f)
            with open(p, "rb") as fh:
                lines = sum(1 for _ in fh)
            files[f] = {"sha256": sha256(p), "bytes": os.path.getsize(p), "lines": lines}
    return files


def _numeric_rows(path, limit=5_000_000):
    rows = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i > limit:
                return None
            rows.append(line.split())
    return rows


def describe_diff(path_a, path_b):
    """Short human summary of how two text exports differ."""
    try:
        a, b = _numeric_rows(path_a), _numeric_rows(path_b)
    except (OSError, UnicodeDecodeError) as e:
        return "unreadable: %s" % e
    if a is None or b is None:
        return "too large to diff"
    msg = ["lines %d vs %d" % (len(a), len(b))]
    first = next((i for i in range(min(len(a), len(b))) if a[i] != b[i]), None)
    ndiff = sum(1 for i in range(min(len(a), len(b))) if a[i] != b[i]) + abs(len(a) - len(b))
    msg.append("%d differing lines" % ndiff)
    if first is not None:
        msg.append("first at line %d: %r vs %r" % (first + 1, " ".join(a[first])[:120], " ".join(b[first])[:120]))
    if len(a) == len(b):
        maxd, nonnum = 0.0, 0
        for ra, rb in zip(a, b):
            if len(ra) != len(rb):
                nonnum += 1
                continue
            for x, y in zip(ra, rb):
                if x == y:
                    continue
                try:
                    maxd = max(maxd, abs(float(x) - float(y)))
                except ValueError:
                    nonnum += 1
        msg.append("max abs numeric diff %.6g, non-numeric/shape diffs %d" % (maxd, nonnum))
    return "; ".join(msg)


def compare(golden, result, ref_dir=None, out_dir=None):
    problems = []
    if golden.get("generated_inputs") != result.get("generated_inputs"):
        problems.append(("_gen", "derived input hashes differ (fixture generation not reproducible)"))
    for case, g in golden["cases"].items():
        r = result["cases"].get(case)
        if r is None:
            continue
        if g["exit_code"] != r["exit_code"]:
            problems.append((case, "exit code %s (golden %s)" % (r["exit_code"], g["exit_code"])))
        gf, rf = g["outputs"], r["outputs"]
        for f in sorted(set(gf) - set(rf)):
            problems.append((case, "missing output %s" % f))
        for f in sorted(set(rf) - set(gf)):
            problems.append((case, "unexpected output %s" % f))
        for f in sorted(set(gf) & set(rf)):
            if gf[f]["sha256"] != rf[f]["sha256"]:
                detail = "bytes %d vs golden %d, lines %d vs %d" % (rf[f]["bytes"], gf[f]["bytes"], rf[f]["lines"], gf[f]["lines"])
                if ref_dir:
                    pa = os.path.join(ref_dir, case, "outputs", f)
                    pb = os.path.join(out_dir, case, "outputs", f)
                    if os.path.exists(pa) and os.path.exists(pb):
                        detail += "; " + describe_diff(pa, pb)
                problems.append((case, "%s differs: %s" % (f, detail)))
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True, help="karstnsim executable to run")
    ap.add_argument("--out", required=True, help="scratch directory for work roots (recreated per case)")
    ap.add_argument("--assets", default=os.path.normpath(os.path.join(HERE, "..", "..", "..", "Input_files")),
                    help="Input_files directory with the original example assets")
    ap.add_argument("--cases", nargs="*", help="subset of case names")
    ap.add_argument("--timeout", type=float, default=900.0, help="per-case timeout in seconds")
    ap.add_argument("--golden", help="golden JSON to compare against")
    ap.add_argument("--reference-out", help="--out directory of the golden run, for line-level diff summaries")
    ap.add_argument("--write-golden", help="write results as golden JSON to this path")
    a = ap.parse_args()

    binary = os.path.abspath(a.binary)
    if not os.access(binary, os.X_OK):
        raise SystemExit("binary not executable: %s" % binary)
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)

    gen_dir = os.path.join(out, "_gen")
    if os.path.exists(gen_dir):
        shutil.rmtree(gen_dir)
    # Generate in a child process: Linux carries a parent's peak RSS across fork+exec,
    # so a large Python heap here would pollute the binary's ru_maxrss.
    gen = subprocess.run([sys.executable, os.path.join(HERE, "make_fixtures.py"), "--assets", a.assets,
                          "--out", gen_dir], check=True, capture_output=True, text=True)
    gen_hashes = {}
    for line in gen.stdout.splitlines():
        h, name = line.split(None, 1)
        gen_hashes[name] = h

    cases = list_cases(a.cases)
    result = {
        "binary": binary,
        "binary_sha256": sha256(binary),
        "assets": os.path.abspath(a.assets),
        "generated_inputs": gen_hashes,
        "cases": {},
    }
    for case in cases:
        work = prepare_work(case, out, a.assets, gen_dir)
        info = run_case(binary, work, a.timeout)
        info["outputs"] = hash_outputs(work)
        result["cases"][case] = info
        print("%-40s rc=%-7s wall=%7.2fs rss=%8.1fMB files=%d" % (
            case, info["exit_code"], info["wall_s"], (info["max_rss_kb"] or 0) / 1024.0, len(info["outputs"])), flush=True)

    with open(os.path.join(out, "result.json"), "w") as f:
        json.dump(result, f, indent=1, sort_keys=True)

    if a.write_golden:
        golden = {
            "note": "Golden hashes for REDUCED regression fixtures only; not the full-size examples.",
            "binary_sha256": result["binary_sha256"],
            "generated_inputs": gen_hashes,
            "cases": {c: {"exit_code": r["exit_code"], "outputs": r["outputs"]} for c, r in result["cases"].items()},
        }
        with open(a.write_golden, "w") as f:
            json.dump(golden, f, indent=1, sort_keys=True)
            f.write("\n")
        print("golden written: %s" % a.write_golden)

    status = 0
    if a.golden:
        with open(a.golden) as f:
            golden = json.load(f)
        missing_cases = sorted(set(golden["cases"]) - set(result["cases"]))
        if a.cases:
            missing_cases = []
        problems = compare(golden, result, a.reference_out, out)
        for c in missing_cases:
            problems.append((c, "case in golden but not run"))
        if problems:
            status = 1
            print("\nREGRESSION: %d problem(s)" % len(problems))
            for case, p in problems:
                print("  [%s] %s" % (case, p))
        else:
            n = sum(len(r["outputs"]) for r in result["cases"].values())
            print("\nPASS: %d case(s), %d export file(s) byte-identical to golden" % (len(result["cases"]), n))
    return status


if __name__ == "__main__":
    sys.exit(main())
