"""
build_crossing_design.py
Hub-and-spoke crossing design figure.

Frederick Henry & Storm Bays sits at the centre; each of the four other
genetic clusters was crossed reciprocally with it (36 crosses in each
direction). No crosses were made between two non-FH clusters, so every
between-cluster cross has FH as exactly one parent. Within-cluster diallels
(FH and East Coast Victoria) sit outside the hub structure and are shown
separately.

Fst values are pairwise against FH, from Scharfenstein et al.
Set SHOW_FST = False to omit them.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

SHOW_FST = True

HUB_FILL, HUB_EDGE = "#F5C4B3", "#993C1D"
SPOKE_FILL, SPOKE_EDGE = "#9FE1CB", "#0F6E56"
NOTE_FILL, NOTE_EDGE = "#F1EFE8", "#5F5E5A"
LINE = "#5F5E5A"

BW, BH = 2.35, 0.95          # spoke box width/height
HW, HH = 2.25, 1.05          # hub box

hub = (5.0, 3.3)
spokes = [
    # (cx, cy, line1, line2, fst)
    (1.70, 4.85, "Actaeon Island", "& Southport", 0.18),
    (8.30, 4.85, "St Helens Point", "& Falmouth", 0.43),
    (1.70, 1.75, "North Coast", "Tasmania", 0.51),
    (8.30, 1.75, "East Coast", "Victoria", 0.41),
]


def box(ax, cx, cy, w, h, fill, edge, lines, sub=None, lw=1.0):
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0,rounding_size=0.10",
        facecolor=fill, edgecolor=edge, linewidth=lw, zorder=3))
    n = len(lines)
    top = cy + (0.155 if n == 2 else 0.075) if sub else cy + 0.115
    for i, t in enumerate(lines):
        ax.text(cx, top - i * 0.235, t, ha="center", va="center",
                fontsize=9.5, color=edge, zorder=4)
    if sub:
        ax.text(cx, cy - h / 2 + 0.20, sub, ha="center", va="center",
                fontsize=7.8, color=edge, alpha=0.85, zorder=4)


fig, ax = plt.subplots(figsize=(8.2, 5.4))

# spokes drawn first so boxes sit on top
for cx, cy, _, _, fst in spokes:
    dx, dy = hub[0] - cx, hub[1] - cy
    d = (dx ** 2 + dy ** 2) ** 0.5
    ux, uy = dx / d, dy / d
    pad_out, pad_hub = 1.28, 1.30
    x1, y1 = cx + ux * pad_out, cy + uy * pad_out
    x2, y2 = hub[0] - ux * pad_hub, hub[1] - uy * pad_hub
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="<|-|>", mutation_scale=11,
        color=LINE, linewidth=1.0, shrinkA=0, shrinkB=0, zorder=2))
    if SHOW_FST:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        ax.text(mx, my, f"$F_{{ST}}$ = {fst:.2f}", ha="center", va="center",
                fontsize=8, color=LINE,
                bbox=dict(boxstyle="round,pad=0.22", facecolor="white",
                          edgecolor="none"), zorder=3)

for cx, cy, l1, l2, _ in spokes:
    box(ax, cx, cy, BW, BH, SPOKE_FILL, SPOKE_EDGE, [l1, l2],
        "36 as mother, 36 as father")

box(ax, hub[0], hub[1], HW, HH, HUB_FILL, HUB_EDGE,
    ["Frederick Henry", "& Storm Bays"], lw=1.4)

box(ax, 5.0, 0.80, 3.0, 0.85, NOTE_FILL, NOTE_EDGE,
    ["Within-cluster crosses"], "36 FH, 36 East Coast Vic")

ax.set_xlim(0, 10)
ax.set_ylim(0.15, 6.0)
ax.axis("off")
fig.tight_layout()
fig.savefig("crossing_design.png", dpi=300, bbox_inches="tight",
            facecolor="white")
print("wrote crossing_design.png")