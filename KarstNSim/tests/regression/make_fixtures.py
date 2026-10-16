#!/usr/bin/env python3
"""Deterministically derive reduced-size KarstNSim regression inputs.

All derived files are generated from the repository's Input_files assets; nothing
here is random at run time (fixed-seed Python PRNG only). The derived files are
written to an output directory and are not meant to be committed.

Derived files:
  box_coarse.txt      example_box.txt with every positive density radius multiplied
                      by COARSE_FACTOR (NDV cells and IKP unchanged). The original
                      two-level variable density (0.03 / 0.007) stays variable.
  box_medium.txt      same as box_coarse.txt with MEDIUM_FACTOR (denser, ~3x the points).
  box_gradient.txt    example_box.txt with positive density replaced by a vertical
                      and lateral gradient between GRADIENT_MIN and GRADIENT_MAX,
                      exercising many distinct local Poisson radii.
  graph_pts.txt       jittered lattice of sampling points below the topography and
                      outside the no-karst spheres.
  input_graph.txt     symmetric k-nearest-neighbour graph (NODES/EDGES format read
                      by translate_input_graph) over sinks + springs + waypoints +
                      graph_pts, i.e. exactly the simulation sample set of case
                      c05_import_graph.
"""

import argparse
import hashlib
import math
import os
import random
import sys

COARSE_FACTOR = 3.0
MEDIUM_FACTOR = 2.0
GRADIENT_MIN = 0.05
GRADIENT_MAX = 0.12
LATTICE_STEP = 40.0
LATTICE_JITTER = 12.0
LATTICE_SEED = 20260924
GRAPH_K = 12


def _split_box(path):
    header = []
    with open(path) as f:
        for line in f:
            header.append(line)
            if line.startswith("Index"):
                break
        rows = [line.split() for line in f if line.strip()]
    return header, rows


def _box_dims(header):
    dims = {}
    for line in header:
        p = line.split()
        if p and p[0] in ("nu", "nv", "nw"):
            dims[p[0]] = int(p[1])
    return dims["nu"], dims["nv"], dims["nw"]


def write_box_scaled(src, dst, factor):
    header, rows = _split_box(src)
    with open(dst, "w") as out:
        out.writelines(header)
        for r in rows:
            d = float(r[1])
            if d > 0:
                d = min(d * factor, 0.99)
                r[1] = "%.10f" % d
            out.write("\t".join(r) + "\n")


def write_box_gradient(src, dst):
    header, rows = _split_box(src)
    nu, nv, nw = _box_dims(header)
    with open(dst, "w") as out:
        out.writelines(header)
        for r in rows:
            i = int(r[0])
            d = float(r[1])
            if d > 0:
                u = i % nu
                w = i // (nu * nv)
                t = 0.5 * (w / max(nw - 1, 1)) + 0.5 * (u / max(nu - 1, 1))
                d = GRADIENT_MIN + (GRADIENT_MAX - GRADIENT_MIN) * t
                r[1] = "%.10f" % d
            out.write("\t".join(r) + "\n")


def read_tsurf(path):
    verts, tris = [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == "VRTX":
                verts.append((float(p[2]), float(p[3]), float(p[4])))
            elif p[0] == "TRGL":
                tris.append((int(p[2]), int(p[3]), int(p[4])))
    return verts, tris


class TopoQuery:
    """Point-in-triangle 2D lookup returning interpolated topography z (or None)."""

    def __init__(self, path, cell=50.0):
        self.verts, self.tris = read_tsurf(path)
        self.cell = cell
        self.buckets = {}
        for t in self.tris:
            xs = [self.verts[k][0] for k in t]
            ys = [self.verts[k][1] for k in t]
            for bx in range(int(math.floor(min(xs) / cell)), int(math.floor(max(xs) / cell)) + 1):
                for by in range(int(math.floor(min(ys) / cell)), int(math.floor(max(ys) / cell)) + 1):
                    self.buckets.setdefault((bx, by), []).append(t)

    def z_at(self, x, y):
        key = (int(math.floor(x / self.cell)), int(math.floor(y / self.cell)))
        for t in self.buckets.get(key, ()):
            (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = (self.verts[k] for k in t)
            det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
            if det == 0:
                continue
            a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det
            b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det
            c = 1.0 - a - b
            if a >= -1e-9 and b >= -1e-9 and c >= -1e-9:
                return a * z1 + b * z2 + c * z3
        return None


def read_point_tokens(path):
    """Return raw (x, y, z) string tokens of a KarstNSim point set, in file order."""
    pts = []
    with open(path) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 4:
                pts.append((p[1], p[2], p[3]))
    return pts


def read_spheres(path):
    spheres = []
    with open(path) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 5:
                spheres.append((float(p[1]), float(p[2]), float(p[3]), float(p[4])))
    return spheres


def write_graph_inputs(assets, pts_dst, graph_dst):
    base = os.path.join(assets, "1_base")
    topo = TopoQuery(os.path.join(base, "example_topo_surf.txt"))
    spheres = read_spheres(os.path.join(base, "example_nokarstspheres.txt"))
    rng = random.Random(LATTICE_SEED)

    lattice = []
    nx, ny, nz = int(990 / LATTICE_STEP), int(790 / LATTICE_STEP), int(590 / LATTICE_STEP)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                x = i * LATTICE_STEP + rng.uniform(-LATTICE_JITTER, LATTICE_JITTER)
                y = j * LATTICE_STEP + rng.uniform(-LATTICE_JITTER, LATTICE_JITTER)
                z = k * LATTICE_STEP + rng.uniform(-LATTICE_JITTER, LATTICE_JITTER)
                if not (0 < x < 990 and 0 < y < 790 and 0 < z < 590):
                    continue
                zt = topo.z_at(x, y)
                if zt is None or z >= zt - 2.0:
                    continue
                if any((x - sx) ** 2 + (y - sy) ** 2 + (z - sz) ** 2 <= (sr + 5.0) ** 2
                       for sx, sy, sz, sr in spheres):
                    continue
                lattice.append(("%.3f" % x, "%.3f" % y, "%.3f" % z))

    with open(pts_dst, "w") as f:
        f.write("Index\tX\tY\tZ\n")
        for n, p in enumerate(lattice, 1):
            f.write("%d\t%s\t%s\t%s\n" % (n, p[0], p[1], p[2]))

    keys = []
    for name in ("example_sinks.txt", "example_springs.txt", "example_waypoints.txt"):
        keys.extend(read_point_tokens(os.path.join(base, name)))
    nodes_tok = keys + lattice
    nodes = [tuple(float(c) for c in p) for p in nodes_tok]

    cell = 3.0 * LATTICE_STEP
    buckets = {}
    for idx, (x, y, z) in enumerate(nodes):
        buckets.setdefault((int(x // cell), int(y // cell), int(z // cell)), []).append(idx)
    edges = set()
    for idx, (x, y, z) in enumerate(nodes):
        bx, by, bz = int(x // cell), int(y // cell), int(z // cell)
        cand = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for o in buckets.get((bx + dx, by + dy, bz + dz), ()):
                        if o != idx:
                            ox, oy, oz = nodes[o]
                            cand.append(((ox - x) ** 2 + (oy - y) ** 2 + (oz - z) ** 2, o))
        cand.sort()
        for _, o in cand[:GRAPH_K]:
            edges.add((min(idx, o), max(idx, o)))

    with open(graph_dst, "w") as f:
        f.write("# KarstNSim regression input graph: symmetric %d-NN over keypoints + graph_pts.txt\n" % GRAPH_K)
        f.write("NODES %d\n" % len(nodes_tok))
        for n, p in enumerate(nodes_tok):
            f.write("%d %s %s %s\n" % (n, p[0], p[1], p[2]))
        f.write("EDGES %d\n" % len(edges))
        for a, b in sorted(edges):
            f.write("%d %d\n" % (a, b))
    return len(lattice), len(nodes_tok), len(edges)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def generate(assets, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    box = os.path.join(assets, "1_base", "example_box.txt")
    files = {
        "box_coarse.txt": lambda p: write_box_scaled(box, p, COARSE_FACTOR),
        "box_medium.txt": lambda p: write_box_scaled(box, p, MEDIUM_FACTOR),
        "box_gradient.txt": lambda p: write_box_gradient(box, p),
    }
    for name, fn in files.items():
        fn(os.path.join(out_dir, name))
    write_graph_inputs(assets, os.path.join(out_dir, "graph_pts.txt"), os.path.join(out_dir, "input_graph.txt"))
    names = sorted(os.listdir(out_dir))
    return {n: sha256(os.path.join(out_dir, n)) for n in names}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assets", default=os.path.normpath(os.path.join(here, "..", "..", "..", "Input_files")))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    for n, h in generate(a.assets, a.out).items():
        print(h, n)


if __name__ == "__main__":
    sys.exit(main())
