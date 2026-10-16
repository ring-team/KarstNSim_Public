#!/usr/bin/env python3
"""Semantic check of a corrected *_nghb_graph.txt against the original (buggy) export.

Upstream GraphOperations::save_nghb_graph (ff55d39) appends one segment for EVERY
rectangular adjacency slot, including unused slots with target -1, whose end point is
read from samples[-1] (out of bounds). Cost properties are appended only for valid
targets, so in the original export they are shifted relative to the geometry, and the
tail of the property array is read past its end.

For a corrected export (unused slots skipped) this script verifies that:
  1. every segment end point is a point of the run's sampling cloud (<name>_pts.txt) and
     no segment is degenerate;
  2. the corrected geometry sequence equals the original geometry sequence with the
     segments ending outside the sampling cloud (junk slots) removed;
  3. when cost properties are present, the corrected per-segment cost sequence equals the
     first N per-segment costs of the original export (same values, original order,
     only the misalignment removed);
  4. optionally (--edges), the segment count equals twice the undirected edge count of
     an imported NODES/EDGES graph.

Usage: check_nghb_export.py ORIGINAL CORRECTED PTS [--edges INPUT_GRAPH]
"""

import argparse
import sys


def read_segments(path):
    segs = []
    with open(path) as f:
        header = next(f).split()
        rows = [line.split() for line in f]
    if len(rows) % 2:
        raise SystemExit("%s: odd number of rows" % path)
    for k in range(0, len(rows), 2):
        a, b = rows[k], rows[k + 1]
        if a[0] != b[0]:
            raise SystemExit("%s: rows %d/%d have different segment index" % (path, k + 2, k + 3))
        segs.append((tuple(a[1:4]), tuple(b[1:4]), tuple(a[4:]), tuple(b[4:])))
    return header, segs


def read_pts(path):
    with open(path) as f:
        next(f)
        return {tuple(line.split()[1:4]) for line in f if line.strip()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("original")
    ap.add_argument("corrected")
    ap.add_argument("pts")
    ap.add_argument("--edges", help="imported NODES/EDGES graph used by the run")
    a = ap.parse_args()

    pts = read_pts(a.pts)
    h_orig, orig = read_segments(a.original)
    h_corr, corr = read_segments(a.corrected)
    failures = []
    if h_orig != h_corr:
        failures.append("header differs: %s vs %s" % (h_orig, h_corr))

    bad = [i for i, s in enumerate(corr) if s[0] not in pts or s[1] not in pts or s[0] == s[1]]
    if bad:
        failures.append("corrected: %d segments with an end point outside the sampling cloud or degenerate (first index %d)" % (len(bad), bad[0] + 1))

    junk = [s for s in orig if s[1] not in pts]
    kept = [(s[0], s[1]) for s in orig if s[1] in pts]
    if kept != [(s[0], s[1]) for s in corr]:
        failures.append("corrected geometry != original geometry minus junk segments")

    has_props = len(h_corr) > 4
    if has_props:
        if any(s[2] != s[3] for s in corr):
            failures.append("corrected: start/end rows of a segment carry different properties")
        orig_costs = [s[2] for s in orig]
        corr_costs = [s[2] for s in corr]
        if orig_costs[:len(corr_costs)] != corr_costs:
            failures.append("corrected cost sequence != leading original cost sequence")

    if a.edges:
        with open(a.edges) as f:
            n_edges = next(int(l.split()[1]) for l in f if l.startswith("EDGES"))
        if len(corr) != 2 * n_edges:
            failures.append("segment count %d != 2 x %d imported edges" % (len(corr), n_edges))

    junk_ends = {}
    for s in junk:
        junk_ends[s[1]] = junk_ends.get(s[1], 0) + 1
    print("original segments %d, junk (end outside sampling cloud) %d, distinct junk end points %d%s" % (
        len(orig), len(junk), len(junk_ends),
        (", e.g. " + ", ".join(" ".join(p) for p in list(junk_ends)[:3])) if junk_ends else ""))
    print("corrected segments %d, properties %s" % (len(corr), "yes" if has_props else "no"))
    if has_props:
        shifted = sum(1 for o, c in zip([s for s in orig if s[1] in pts], corr) if o[2] != c[2])
        print("original valid segments carrying a shifted (wrong) cost: %d" % shifted)
    for f in failures:
        print("FAIL:", f)
    print("OK" if not failures else "FAILED")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
