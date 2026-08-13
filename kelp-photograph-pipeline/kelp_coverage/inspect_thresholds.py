"""
inspect_thresholds.py
---------------------
Randomly sample kelp photos and visually check what the coverage threshold is
actually classifying, ROI by ROI.

Coverage in measure_kelp_coverage.py is:
    valid pixels = green >= DARK_CUTOFF        (excludes dark background)
    kelp pixels  = (red / green) > RG_THRESHOLD
    coverage     = kelp / valid, within each ROI

This tool overlays that classification on the image and gives you live sliders
for both constants, so you can see whether 1.05 / 30 is catching the fronds and
nothing else. It never changes your data: it only writes a log of what you
inspected and any photo you flag.

Keys:  n = next photo   p = previous   f = flag this photo   r = reset sliders
       q = quit
Sliders: R:G threshold, dark cutoff.

Usage:
    python inspect_thresholds.py <root> --n 25
    python inspect_thresholds.py <root> --n 25 --seed 1 --date 20260507
"""
import argparse
import csv
import random
import sys
from pathlib import Path

import numpy as np
import matplotlib
# macosx delivers key events but not reliable widget interaction
for _bk in ("QtAgg", "TkAgg"):
    try:
        matplotlib.use(_bk, force=True)
        break
    except Exception:
        continue
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

import measure_kelp_coverage as mk

LOG_COLS = ["photo", "folder", "roi", "coverage", "valid_px", "kelp_px",
            "rg_threshold", "dark_cutoff", "flagged"]


def find_zip_upwards(folder: Path, stop_at: Path):
    """Nearest ROI zip in `folder` or any parent up to and including stop_at.
    A folder-local zip wins over one higher up, matching the pipeline."""
    f = folder
    while True:
        zips = [z for z in f.glob("*.zip") if not z.name.startswith("._")]
        if len(zips) == 1:
            return zips[0]
        if len(zips) > 1:
            print(f"  several zips in {f}, using {sorted(z.name for z in zips)[0]}")
            return sorted(zips)[0]
        if f == stop_at or f.parent == f:
            return None
        f = f.parent


def gather(root: Path, date_filter=None):
    """Every (tif, roi_zip) pair under root."""
    out = []
    folders = [root] + [p for p in root.rglob("*") if p.is_dir()]
    for f in sorted(set(folders)):
        if date_filter and date_filter not in f.name:
            continue
        try:
            tifs = mk.list_tifs(f)
        except Exception:
            continue
        if not tifs:
            continue
        zp = find_zip_upwards(f, root)
        if zp is None:
            print(f"  no ROI zip for {f} (searched up to {root}); "
                  f"showing photos without ROI outlines")
        out.append((f, zp, tifs))
    flat = []
    for folder, zp, tifs in out:
        print(f"  {folder}: {len(tifs)} tifs, roi zip = "
              f"{zp.name if zp else 'none'}")
        for t in tifs:
            flat.append((t, zp))
    return flat


def classify(rgb, rg_thr, dark):
    """Return per-pixel valid and kelp boolean arrays for the whole frame."""
    r = rgb[:, :, 0].astype(float)
    g = rgb[:, :, 1].astype(float)
    valid = g >= dark
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(g > 0, r / g, 0.0)
    return valid, valid & (ratio > rg_thr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--n", type=int, default=20, help="how many photos")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--date", default=None,
                    help="only sample folders whose name contains this")
    ap.add_argument("--log", default="threshold_inspection.csv")
    a = ap.parse_args()

    root = Path(a.root).resolve()
    if not root.exists():
        sys.exit(f"Path does not exist: {root}")
    print(f"scanning {root}")
    pairs = gather(root, a.date)
    if not pairs:
        subs = [p.name for p in sorted(root.iterdir())][:25] if root.is_dir() \
            else []
        print(f"No TIFFs found under {root}")
        print(f"  contents: {subs}")
        print("  (looked for *.tif/*.tiff in this folder and every subfolder)")
        if a.date:
            print(f"  note: --date {a.date} filters to folders whose NAME "
                  f"contains that string; drop it to scan everything")
        sys.exit(1)
    if a.seed is not None:
        random.seed(a.seed)
    sample = random.sample(pairs, min(a.n, len(pairs)))
    print(f"{len(pairs)} photos found; inspecting {len(sample)} at random.")
    print("opening window (keys: n next, p prev, f flag, r reset, q quit)")

    log = open(a.log, "a", newline="")
    writer = csv.DictWriter(log, fieldnames=LOG_COLS)
    if log.tell() == 0:
        writer.writeheader()

    state = {"i": 0, "rgb": None, "masks": None, "flag": False}

    fig = plt.figure(figsize=(13, 7))
    ax_img = fig.add_axes([0.03, 0.22, 0.45, 0.72])
    ax_ovl = fig.add_axes([0.52, 0.22, 0.45, 0.72])
    ax_rg = fig.add_axes([0.15, 0.10, 0.70, 0.03])
    ax_dk = fig.add_axes([0.15, 0.05, 0.70, 0.03])
    s_rg = Slider(ax_rg, "R:G threshold", 0.80, 1.60,
                  valinit=mk.RG_THRESHOLD, valstep=0.01)
    s_dk = Slider(ax_dk, "dark cutoff", 0, 120, valinit=mk.DARK_CUTOFF,
                  valstep=1)

    def load(i):
        tif, zp = sample[i]
        rgb = mk.load_rgb(tif)
        masks = {}
        if zp is not None:
            try:
                masks = mk.load_roi_masks(zp, rgb.shape[:2])
            except Exception as e:
                print(f"ROI load failed for {tif.name}: {e}")
        state.update(rgb=rgb, masks=masks, flag=False)

    def redraw():
        tif, zp = sample[state["i"]]
        rgb, masks = state["rgb"], state["masks"]
        rg_thr, dark = s_rg.val, s_dk.val
        valid, kelp = classify(rgb, rg_thr, dark)

        ax_img.clear(); ax_ovl.clear()
        disp = rgb.astype(float)
        disp = disp / max(disp.max(), 1.0)
        ax_img.imshow(disp)
        ax_ovl.imshow(disp)

        lines = []
        for name, m in sorted(masks.items()):
            if m.shape != rgb.shape[:2]:
                continue
            v = int((valid & m).sum())
            k = int((kelp & m).sum())
            cov = k / v if v else float("nan")
            lines.append((name, cov, v, k))
            ys, xs = np.where(m)
            if len(xs):
                ax_img.add_patch(plt.Rectangle(
                    (xs.min(), ys.min()), xs.max() - xs.min(),
                    ys.max() - ys.min(), fill=False, edgecolor="yellow",
                    lw=1.2))
                ax_img.text(xs.min(), ys.min() - 6, name, color="yellow",
                            fontsize=8)
        # red tint = classified kelp, inside any ROI
        anyroi = np.zeros(rgb.shape[:2], bool)
        for m in masks.values():
            if m.shape == rgb.shape[:2]:
                anyroi |= m
        show = kelp & anyroi if masks else kelp
        ovl = np.zeros((*rgb.shape[:2], 4))
        ovl[show] = [1, 0, 0, 0.45]
        ax_ovl.imshow(ovl)

        for ax in (ax_img, ax_ovl):
            ax.set_xticks([]); ax.set_yticks([])
        ax_img.set_title(f"[{state['i']+1}/{len(sample)}]  {tif.name}\n"
                         f"{tif.parent.name}", fontsize=9)
        txt = "  |  ".join(f"{n}: {c:.3f}" for n, c, _, _ in lines) or "no ROI"
        ax_ovl.set_title(f"classified kelp (red)   coverage -> {txt}\n"
                         f"R:G>{rg_thr:.2f}, green>={int(dark)}"
                         + ("   [FLAGGED]" if state["flag"] else ""),
                         fontsize=9)
        state["lines"] = lines
        fig.canvas.draw_idle()

    def write_rows():
        tif, _ = sample[state["i"]]
        for name, cov, v, k in state.get("lines", []):
            writer.writerow({"photo": tif.name, "folder": tif.parent.name,
                             "roi": name, "coverage": round(cov, 4),
                             "valid_px": v, "kelp_px": k,
                             "rg_threshold": round(s_rg.val, 3),
                             "dark_cutoff": int(s_dk.val),
                             "flagged": state["flag"]})
        log.flush()

    def on_key(e):
        if e.key == "n":
            write_rows()
            if state["i"] < len(sample) - 1:
                state["i"] += 1
                load(state["i"]); redraw()
            else:
                print("last photo")
        elif e.key == "p" and state["i"] > 0:
            state["i"] -= 1
            load(state["i"]); redraw()
        elif e.key == "f":
            state["flag"] = not state["flag"]
            redraw()
        elif e.key == "r":
            s_rg.reset(); s_dk.reset()
        elif e.key == "q":
            write_rows()
            plt.close(fig)

    s_rg.on_changed(lambda _: redraw())
    s_dk.on_changed(lambda _: redraw())
    fig.canvas.mpl_connect("key_press_event", on_key)

    load(0); redraw()
    plt.show()
    log.close()
    print(f"wrote {a.log}")


if __name__ == "__main__":
    main()