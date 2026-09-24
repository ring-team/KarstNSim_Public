#!/usr/bin/env python3
"""Plot actual region JSON geometry. Optional matplotlib dependency."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon, Rectangle

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("regions", nargs="+", type=Path)
parser.add_argument("--out", required=True, type=Path)
args = parser.parse_args()
regions = [json.loads(path.read_text()) for path in args.regions]
fig, ax = plt.subplots(figsize=(14, 8), layout="constrained")
fig.patch.set_facecolor("#101922")
ax.set_facecolor("#101922")
ax.tick_params(colors="#acbac7")
for spine in ax.spines.values():
    spine.set_color("#526473")
segments, heights = [], []
for region in regions:
    ox, oy, oz = region["origin"]
    sx, sy, _ = region["size"]
    ax.add_patch(Rectangle((ox, oy), sx, sy, fill=False, edgecolor="#61768a", linewidth=1))
    ax.text(ox + sx/2, oy + sy + 8, f"Region {region['key']}", color="#d0dce8", ha="center")
    for cavern in region["caverns"]:
        if cavern["boundary"]:
            ax.add_patch(Polygon([(x+ox, y+oy) for x, y in cavern["boundary"]],
                                 facecolor="#426675", edgecolor="#b3dde3", linewidth=0.7))
    for tunnel in region["tunnels"]:
        centers = [profile["center"] for profile in tunnel["centerline"]]
        for a, b in zip(centers, centers[1:]):
            segments.append([(a[0]+ox, a[1]+oy), (b[0]+ox, b[1]+oy)])
            heights.append((a[2]+b[2])/2+oz)
    for seam in region["boundaries"]:
        port = next(p for p in region["ports"] if p["id"] == seam["port_id"])
        a, b = port["position"], seam["position"]
        ax.plot([a[0]+ox,b[0]+ox], [a[1]+oy,b[1]+oy], color="#f5db85", linewidth=1.8)
        ax.scatter([b[0]+ox], [b[1]+oy], s=18, color="#f5db85", zorder=5)
if segments:
    lines = LineCollection(segments, cmap="YlOrRd", norm=Normalize(min(heights), max(heights)), linewidths=1.8)
    lines.set_array(heights)
    ax.add_collection(lines)
    bar = fig.colorbar(lines, ax=ax, fraction=.022, pad=.02)
    bar.set_label("Tunnel center elevation (m)", color="#d0dce8")
    bar.ax.tick_params(colors="#acbac7")
ax.autoscale_view()
ax.set_aspect("equal")
ax.set_xlabel("World X (m)", color="#d0dce8")
ax.set_ylabel("World Y (m)", color="#d0dce8")
ax.set_title("FinalBuildCaves: generated logical cave geometry", color="#eef4f8", fontsize=16, pad=28)
fig.text(.5, .012, "Actual generated outlines and tunnel profiles. Gold marks region connections. No 3D wall mesh shown.",
         color="#acbac7", fontsize=9, ha="center")
args.out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.out, dpi=160, facecolor=fig.get_facecolor())
print(args.out.resolve())
