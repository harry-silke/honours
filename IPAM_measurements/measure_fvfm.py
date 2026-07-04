"""
measure_fvfm.py
---------------
Interactive Fv/Fm + necrosis measurement for two-beaker Imaging-PAM TIFFs
(frame 0 = Fo, frame 1 = Fm), produced by pim2tiff.py.

Per image you see the contrast-stretched Fm frame with two draggable rectangles
carried forward from the previous image (left = id1, right = id2). Adjust by
dragging the box handles, then accept. One row per beaker per image is written
and appended immediately, so an interrupted session keeps its progress and can
resume (already-measured files are skipped).

Background vs tissue is decided by the Fo floor, so a rectangle may include
dark corners of water without harm; only pixels with Fo >= FO_FLOOR count as
tissue. Necrosis is tracked as its own column: the share of tissue whose
Fv/Fm <= NECROSIS_CUT. Fv/Fm is the mean over all tissue pixels (necrotic
included), matching the ImageJ macro.

Filename convention: ID1_ID2_DDMMYY.tif   (e.g. 5C_4C_080526.tif)

Keys:  drag box handles to adjust   |   enter/space = accept
       r = reset boxes to carried    |   s = skip image (no row)   |   q = quit

Requirements:  pip install numpy tifffile matplotlib

Usage:
    python measure_fvfm.py <folder-or-root> [--out fvfm_results.csv]
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import tifffile

# ---- tuning knobs (nudge per run if the instrument/exposure changes) ----
FO_FLOOR      = 40      # Fo below this = background water, not tissue
NECROSIS_CUT  = 0.10    # tissue with Fv/Fm <= this counts as necrotic
FVFM_MAX      = 0.85    # clamp ceiling for noise; theoretical max ~0.83
FO_FRAME      = 0       # frame index of Fo
FM_FRAME      = 1       # frame index of Fm

COLS = ["photo_filename", "pair_key", "date", "side", "individual",
        "rect_px", "tissue_px", "Fo", "Fm", "FvFm", "necrosis_pct",
        "px_below_0.10", "px_below_0.15", "px_below_0.20",
        "fo_floor", "necrosis_cut"]


def list_tifs(root: Path):
    tifs = sorted(p for p in root.rglob("*.tif")
                  if not p.name.startswith("._"))
    return [p for p in tifs if len(p.stem.split("_")) >= 3]


def load_fo_fm(path: Path):
    arr = tifffile.imread(str(path))
    if arr.ndim != 3 or arr.shape[0] <= FM_FRAME:
        raise ValueError(f"{path.name}: expected a multi-frame PAM TIFF")
    return arr[FO_FRAME].astype(float), arr[FM_FRAME].astype(float)


def compute(fo, fm, rect):
    x0, y0, x1, y1 = rect
    fo_c = fo[y0:y1, x0:x1]
    fm_c = fm[y0:y1, x0:x1]
    rect_px = fo_c.size
    tissue = fo_c >= FO_FLOOR
    n_t = int(tissue.sum())
    if n_t == 0:
        return dict(rect_px=rect_px, tissue_px=0, Fo=np.nan, Fm=np.nan,
                    FvFm=np.nan, necrosis_pct=np.nan,
                    b10=0, b15=0, b20=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        fvfm = np.where(fm_c > 0, (fm_c - fo_c) / fm_c, 0.0)
    fvfm = np.clip(fvfm, 0.0, FVFM_MAX)
    t = fvfm[tissue]
    return dict(
        rect_px=rect_px, tissue_px=n_t,
        Fo=round(float(fo_c[tissue].mean()), 2),
        Fm=round(float(fm_c[tissue].mean()), 2),
        FvFm=round(float(t.mean()), 4),
        necrosis_pct=round(float((t <= NECROSIS_CUT).mean() * 100), 2),
        b10=int((t <= 0.10).sum()), b15=int((t <= 0.15).sum()),
        b20=int((t <= 0.20).sum()))


def parse_name(stem: str):
    p = stem.split("_")
    return p[0], p[1], p[2]          # id1, id2, date (DDMMYY)


def confirm_image(fo, fm, boxes, title):
    """Show Fm with two draggable boxes; return ('accept', boxes) / ('skip',)
    / ('quit',). boxes is {'left': (x0,y0,x1,y1), 'right': (...)}."""
    import matplotlib.pyplot as plt
    from matplotlib.widgets import RectangleSelector

    H, W = fm.shape
    disp = fm.copy()
    pos = disp[disp > 0]
    lo, hi = (np.percentile(pos, [1, 99]) if pos.size else (0, 1))

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.imshow(disp, cmap="gray", vmin=lo, vmax=hi)
    ax.set_title(title, fontsize=10)
    ax.set_axis_off()

    state = {"action": None}

    def mk(side, colour):
        s = RectangleSelector(ax, lambda *a: None, useblit=True,
                              button=[1], interactive=True,
                              props=dict(edgecolor=colour, facecolor="none",
                                         linewidth=1.6))
        x0, y0, x1, y1 = boxes[side]
        s.extents = (x0, x1, y0, y1)
        return s
    sels = {"left": mk("left", "red"), "right": mk("right", "cyan")}

    def read():
        out = {}
        for side, s in sels.items():
            x0, x1, y0, y1 = s.extents
            x0, x1 = sorted((int(round(x0)), int(round(x1))))
            y0, y1 = sorted((int(round(y0)), int(round(y1))))
            x0, x1 = max(0, x0), min(W, x1)
            y0, y1 = max(0, y0), min(H, y1)
            out[side] = (x0, y0, x1, y1)
        return out

    def on_key(e):
        if e.key in (" ", "enter"):
            state["action"] = ("accept", read())
        elif e.key == "s":
            state["action"] = ("skip",)
        elif e.key == "q":
            state["action"] = ("quit",)
        elif e.key == "r":
            for side, s in sels.items():
                x0, y0, x1, y1 = boxes[side]
                s.extents = (x0, x1, y0, y1)
            fig.canvas.draw_idle()
            return
        else:
            return
        plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()
    return state["action"] or ("skip",)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        sys.exit(f"Not a directory: {root}")
    tifs = list_tifs(root)
    if not tifs:
        sys.exit(f"No ID1_ID2_DATE .tif files under {root}")

    out_path = Path(args.out) if args.out else root / "fvfm_results.csv"
    done = set()
    if out_path.exists():
        with open(out_path) as f:
            done = {r["photo_filename"] for r in csv.DictReader(f)}
    new_file = not out_path.exists()
    fout = open(out_path, "a", newline="")
    writer = csv.DictWriter(fout, fieldnames=COLS)
    if new_file:
        writer.writeheader()

    # carried boxes: start as left/right halves of the first image
    H, W = load_fo_fm(tifs[0])[0].shape
    boxes = {"left": (0, 0, W // 2, H), "right": (W // 2, 0, W, H)}

    todo = [t for t in tifs if t.name not in done]
    print(f"{len(tifs)} images, {len(done)} already done, {len(todo)} to do.")

    for i, tif in enumerate(todo, 1):
        try:
            fo, fm = load_fo_fm(tif)
        except Exception as e:
            print(f"skip {tif.name}: {e}")
            continue
        id1, id2, date = parse_name(tif.stem)
        title = (f"[{i}/{len(todo)}]  {tif.name}\n"
                 f"drag handles | enter=accept  r=reset  s=skip  q=quit")
        res = confirm_image(fo, fm, boxes, title)

        if res[0] == "quit":
            print("Quit. Progress saved.")
            break
        if res[0] == "skip":
            print(f"skipped {tif.name}")
            continue

        boxes = res[1]
        pair_key = f"{id1}_{id2}_{date}"
        for side, individual in (("left", id1), ("right", id2)):
            m = compute(fo, fm, boxes[side])
            writer.writerow({
                "photo_filename": tif.name, "pair_key": pair_key, "date": date,
                "side": side, "individual": individual,
                "rect_px": m["rect_px"], "tissue_px": m["tissue_px"],
                "Fo": m["Fo"], "Fm": m["Fm"], "FvFm": m["FvFm"],
                "necrosis_pct": m["necrosis_pct"],
                "px_below_0.10": m["b10"], "px_below_0.15": m["b15"],
                "px_below_0.20": m["b20"],
                "fo_floor": FO_FLOOR, "necrosis_cut": NECROSIS_CUT})
        fout.flush()
        print(f"{tif.name}: L FvFm={compute(fo,fm,boxes['left'])['FvFm']} "
              f"R FvFm={compute(fo,fm,boxes['right'])['FvFm']}")

    fout.close()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()