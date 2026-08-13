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

Keys:  1 = edit left box, 2 = edit right box (drag to draw the active one)
       enter/space = accept  |  r = reset to carried  |  s = skip  |  q = quit

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

import matplotlib
# The macosx backend delivers key events but not the mouse-drag events this
# tool needs. Prefer Qt or Tk, which do. Falls back to the default otherwise.
for _bk in ("QtAgg", "TkAgg"):
    try:
        matplotlib.use(_bk, force=True)
        break
    except Exception:
        continue

# ---- tuning knobs (nudge per run if the instrument/exposure changes) ----
FO_FLOOR      = 40      # Fo below this = background water, not tissue
NECROSIS_CUT  = 0.10    # tissue with Fv/Fm <= this counts as necrotic
FVFM_MAX      = 0.85    # clamp ceiling for noise; theoretical max ~0.83
FO_FRAME      = 0       # frame index of Fo
FM_FRAME      = 1       # frame index of Fm

COLS = ["photo_filename", "pair_key", "date", "side", "individual",
        "roi_px", "tissue_px", "Fo", "Fm", "FvFm", "necrosis_pct",
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


def compute(fo, fm, box):
    """Metrics over the ellipse inscribed in the bounding box `box`
    (x0,y0,x1,y1), so corner leak outside the beaker is excluded."""
    x0, y0, x1, y1 = box
    fo_c = fo[y0:y1, x0:x1]
    fm_c = fm[y0:y1, x0:x1]
    h, w = fo_c.shape
    if h == 0 or w == 0:
        return dict(roi_px=0, tissue_px=0, Fo=np.nan, Fm=np.nan, FvFm=np.nan,
                    necrosis_pct=np.nan, b10=0, b15=0, b20=0)
    yy, xx = np.ogrid[0:h, 0:w]
    cx, cy = (w - 1) / 2.0, (h - 1) / 2.0
    rx, ry = max(w / 2.0, 0.5), max(h / 2.0, 0.5)
    ellipse = ((xx - cx) / rx) ** 2 + ((yy - cy) / ry) ** 2 <= 1.0
    tissue = ellipse & (fo_c >= FO_FLOOR)
    n_t = int(tissue.sum())
    if n_t == 0:
        return dict(roi_px=int(ellipse.sum()), tissue_px=0, Fo=np.nan,
                    Fm=np.nan, FvFm=np.nan, necrosis_pct=np.nan,
                    b10=0, b15=0, b20=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        fvfm = np.where(fm_c > 0, (fm_c - fo_c) / fm_c, 0.0)
    fvfm = np.clip(fvfm, 0.0, FVFM_MAX)
    t = fvfm[tissue]
    return dict(
        roi_px=int(ellipse.sum()), tissue_px=n_t,
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
    """Show Fm; pick a side with 1 (left) or 2 (right), then click-drag to
    draw that box. Raw mouse events, no RectangleSelector. Returns
    ('accept', boxes) / ('skip',) / ('quit',)."""
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    print(f"matplotlib backend: {matplotlib.get_backend()}")

    H, W = fm.shape
    disp = fm.copy()
    pos = disp[disp > 0]
    lo, hi = (np.percentile(pos, [1, 99]) if pos.size else (0, 1))

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.imshow(disp, cmap="gray", vmin=lo, vmax=hi)
    ax.set_axis_off()

    ext = {s: (boxes[s][0], boxes[s][2], boxes[s][1], boxes[s][3])  # x0,x1,y0,y1
           for s in ("left", "right")}
    colour = {"left": "red", "right": "cyan"}
    patch = {}
    for side in ("left", "right"):
        x0, x1, y0, y1 = ext[side]
        e = Ellipse(((x0 + x1) / 2, (y0 + y1) / 2), abs(x1 - x0), abs(y1 - y0),
                    fill=False, edgecolor=colour[side], linewidth=1.6)
        ax.add_patch(e)
        patch[side] = e

    st = {"action": None, "active": "left", "pressed": False,
          "x0": 0, "y0": 0, "lastx": 0, "lasty": 0}

    def set_title():
        a = st["active"]
        ax.set_title(title + f"\nACTIVE: {a} ({colour[a]})", fontsize=10)
        fig.canvas.draw_idle()
    set_title()

    def draw_box(side, x0, x1, y0, y1):
        e = patch[side]
        e.set_center(((x0 + x1) / 2, (y0 + y1) / 2))
        e.width = abs(x1 - x0)
        e.height = abs(y1 - y0)
        fig.canvas.draw_idle()

    def on_press(e):
        if e.inaxes != ax or e.xdata is None:
            return
        st["pressed"] = True
        st["x0"], st["y0"] = e.xdata, e.ydata
        st["lastx"], st["lasty"] = e.xdata, e.ydata
        print(f"press  {st['active']} at ({e.xdata:.0f},{e.ydata:.0f})")

    def on_motion(e):
        if not st["pressed"] or e.inaxes != ax or e.xdata is None:
            return
        st["lastx"], st["lasty"] = e.xdata, e.ydata
        draw_box(st["active"], st["x0"], e.xdata, st["y0"], e.ydata)

    def on_release(e):
        if not st["pressed"]:
            return
        st["pressed"] = False
        x1 = e.xdata if e.xdata is not None else st["lastx"]
        y1 = e.ydata if e.ydata is not None else st["lasty"]
        ext[st["active"]] = (st["x0"], x1, st["y0"], y1)
        draw_box(st["active"], st["x0"], x1, st["y0"], y1)
        print(f"release {st['active']} box -> "
              f"({st['x0']:.0f},{st['y0']:.0f})-({x1:.0f},{y1:.0f})")

    def read():
        out = {}
        for side in ("left", "right"):
            x0, x1, y0, y1 = ext[side]
            x0, x1 = sorted((int(round(x0)), int(round(x1))))
            y0, y1 = sorted((int(round(y0)), int(round(y1))))
            out[side] = (max(0, x0), max(0, y0), min(W, x1), min(H, y1))
        return out

    def on_key(e):
        if e.key == "1":
            st["active"] = "left"; set_title(); return
        if e.key == "2":
            st["active"] = "right"; set_title(); return
        if e.key == "r":
            for side in ("left", "right"):
                ext[side] = (boxes[side][0], boxes[side][2],
                             boxes[side][1], boxes[side][3])
                draw_box(side, *ext[side])
            return
        if e.key in (" ", "enter"):
            st["action"] = ("accept", read())
        elif e.key == "s":
            st["action"] = ("skip",)
        elif e.key == "q":
            st["action"] = ("quit",)
        else:
            return
        plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", on_press)
    fig.canvas.mpl_connect("motion_notify_event", on_motion)
    fig.canvas.mpl_connect("button_release_event", on_release)
    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()
    return st["action"] or ("skip",)


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
                 f"1=left 2=right, drag to draw | enter=accept r=reset "
                 f"s=skip q=quit")
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
                "roi_px": m["roi_px"], "tissue_px": m["tissue_px"],
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