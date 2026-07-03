"""
measure_kelp_coverage.py
------------------------
Measures kelp coverage, physical area (cm2), and reflectance from a folder
of TIFF images using ROI masks from a single zip file.

Reflectance
-----------
For each sample ROI the script reports mean channel intensity (R, G, B) and
their brightness, measured over the KELP-classified pixels by default
(REFLECTANCE_REGION), so it tracks the colour/pigmentation of the tissue
rather than being diluted by how much kelp is present. It also reports mean
R:G, but note that kelp is *defined* by R:G > threshold, so that ratio is
partly circular; use brightness / per-channel values as the bleaching signal.

Lighting drift between sessions will masquerade as a reflectance change over
time. If the ROI zip contains a white-reference ROI (name in WHITE_ROI_NAMES),
channel values are divided by it to give relative reflectance (0-1) that is
comparable across sessions, and white_ref_used is set True. With no white ROI
the raw 0-255 intensities are emitted and white_ref_used is False.

Grouping, sentinels, TIFF layouts, and output-location behaviour are unchanged
from the coverage-only version.

Usage
-----
    python measure_kelp_coverage.py <folder>

Requirements
------------
    pip install tifffile numpy pandas roifile tqdm Pillow

Output
------
    Writes kelp_coverage_results.csv into the image folder.
"""

import sys
import json
import zipfile
import tempfile
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import roifile
from tqdm import tqdm

warnings.filterwarnings("ignore")

# --- Classifier parameters (keep in sync with main pipeline) ---
RG_THRESHOLD  = 1.05   # R:G ratio threshold for kelp classification
DARK_CUTOFF   = 30     # exclude pixels where green channel < this value

# --- Reflectance parameters ---
# "kelp" -> measure reflectance over kelp-classified pixels (tissue colour).
# "roi"  -> measure over the whole valid ROI (confounds with coverage).
REFLECTANCE_REGION = "kelp"
# ROI names treated as a white reference for normalisation (optional).
WHITE_ROI_NAMES = ("white", "White", "reference", "ref", "whitepatch")

# --- Spatial calibration ---
MM_PER_PIXEL  = 1 / 30
CM2_PER_PIXEL = (MM_PER_PIXEL ** 2) / 100   # mm^2 -> cm^2

# --- Grouping / filename convention ---
EXPECTED_COUNTS   = (2, 3, 4)
SENTINEL_ID       = "1000"
ROI_FOR_LEFT_ID   = "left"
ROI_FOR_RIGHT_ID  = "right"


def list_tifs(folder: Path):
    files = (sorted(folder.glob("*.tif"))  + sorted(folder.glob("*.TIF"))
           + sorted(folder.glob("*.tiff")) + sorted(folder.glob("*.TIFF")))
    return [f for f in files if not f.name.startswith("._")]


def resolve_tiffs_folder(folder: Path):
    """Return the folder that holds the TIFFs (itself or a tiffs/ subfolder),
    or None if neither has any. Non-exiting version of find_tiffs_folder."""
    if list_tifs(folder):
        return folder
    for name in ("tiffs", "TIFFs", "TIFFS", "Tiffs"):
        sub = folder / name
        if sub.is_dir() and list_tifs(sub):
            return sub
    return None


def find_tiffs_folder(folder: Path) -> Path:
    r = resolve_tiffs_folder(folder)
    if r is None:
        sys.exit(f"No TIFFs in {folder} and no 'tiffs/' subfolder with TIFFs.")
    return r


def _zips_in(folder: Path):
    return [f for f in folder.glob("*.zip") if not f.name.startswith("._")]


def single_zip(folder: Path):
    """The one ROI zip directly in `folder`, or None. Raises if several."""
    zips = _zips_in(folder)
    if len(zips) > 1:
        raise ValueError(
            f"Multiple zip files in {folder}: {[z.name for z in zips]}\n"
            "Leave exactly one ROI zip in the folder.")
    return zips[0] if zips else None


def find_roi_zip(tiffs_folder: Path) -> Path:
    searched = []
    for folder in (tiffs_folder, tiffs_folder.parent):
        if folder in searched:
            continue
        searched.append(folder)
        zips = [f for f in folder.glob("*.zip") if not f.name.startswith("._")]
        if len(zips) == 1:
            return zips[0]
        if len(zips) > 1:
            raise ValueError(
                f"Multiple zip files in {folder}: {[z.name for z in zips]}\n"
                "Place exactly one ROI zip in the folder."
            )
    raise FileNotFoundError(
        f"No ROI zip found in {tiffs_folder} or {tiffs_folder.parent}"
    )


def load_roi_masks(roi_zip: Path, image_shape: tuple) -> dict:
    from PIL import Image, ImageDraw

    masks = {}
    with zipfile.ZipFile(roi_zip) as zf:
        roi_names = [n for n in zf.namelist() if n.endswith(".roi")]
        with tempfile.TemporaryDirectory() as tmp:
            zf.extractall(tmp)
            for name in roi_names:
                roi_path = Path(tmp) / name
                roi = roifile.roiread(str(roi_path))
                coords = roi.coordinates()
                if coords is None or len(coords) == 0:
                    continue
                # roifile returns (col, row) for ellipse ROIs -> swap to (row, col)
                rows = coords[:, 1].astype(float)
                cols = coords[:, 0].astype(float)
                row_idx = np.clip(rows.astype(int), 0, image_shape[0] - 1)
                col_idx = np.clip(cols.astype(int), 0, image_shape[1] - 1)
                img = Image.new("L", (image_shape[1], image_shape[0]), 0)
                poly = list(zip(col_idx.tolist(), row_idx.tolist()))
                if len(poly) >= 3:
                    ImageDraw.Draw(img).polygon(poly, fill=1)
                masks[Path(name).stem] = np.array(img, dtype=bool)
    return masks


def white_means(rgb: np.ndarray, mask: np.ndarray):
    """Mean R,G,B over a white-reference ROI. Returns None if empty."""
    if mask.shape != rgb.shape[:2]:
        return None
    sel = rgb[mask].astype(float)
    if len(sel) == 0:
        return None
    m = sel.mean(axis=0)
    # guard against zeros so division is safe
    m = np.where(m <= 0, np.nan, m)
    return m  # array [R, G, B]


def measure_roi(rgb: np.ndarray, mask: np.ndarray, white=None, calib=None) -> dict:
    """Coverage, kelp pixel count, and reflectance for one ROI.

    Reflectance precedence:
      calib  -> per-channel (dn -> reference) interpolation, actual reflectance
      white  -> channel means divided by the white reference (relative 0-1)
      neither-> raw 0-255 channel means
    Coverage always uses raw DN, so its R:G threshold is unaffected.
    """
    nan = float("nan")
    empty = dict(coverage=nan, kelp_px=nan,
                 refl_R=nan, refl_G=nan, refl_B=nan,
                 refl_bright=nan, refl_rg=nan)
    if mask.shape != rgb.shape[:2]:
        return empty

    r = rgb[:, :, 0][mask].astype(float)
    g = rgb[:, :, 1][mask].astype(float)
    b = rgb[:, :, 2][mask].astype(float)

    valid = g >= DARK_CUTOFF
    r, g, b = r[valid], g[valid], b[valid]
    if len(g) == 0:
        return empty

    kelp = (r / g) > RG_THRESHOLD
    coverage = float(kelp.mean())
    kelp_px = float(kelp.sum())

    sel = kelp if REFLECTANCE_REGION == "kelp" else np.ones_like(kelp, bool)
    if sel.sum() == 0:
        return dict(coverage=coverage, kelp_px=kelp_px,
                    refl_R=nan, refl_G=nan, refl_B=nan,
                    refl_bright=nan, refl_rg=nan)

    rs, gs, bs = r[sel], g[sel], b[sel]
    if calib is not None:
        rc = np.interp(rs, calib["R"][0], calib["R"][1])
        gc = np.interp(gs, calib["G"][0], calib["G"][1])
        bc = np.interp(bs, calib["B"][0], calib["B"][1])
        R, G, B = rc.mean(), gc.mean(), bc.mean()
        rg = float((rc / gc).mean())
    elif white is not None:
        R, G, B = rs.mean() / white[0], gs.mean() / white[1], bs.mean() / white[2]
        rg = float((rs / gs).mean())
    else:
        R, G, B = rs.mean(), gs.mean(), bs.mean()
        rg = float((rs / gs).mean())

    bright = float(np.nanmean([R, G, B]))
    return dict(coverage=coverage, kelp_px=kelp_px,
                refl_R=float(R), refl_G=float(G), refl_B=float(B),
                refl_bright=bright, refl_rg=rg)


def parse_stem(stem: str):
    parts = stem.split("_")
    if len(parts) == 3:
        id1, id2, date = parts
        n = 1
    elif len(parts) == 4 and parts[3].isdigit():
        id1, id2, date, n_str = parts
        n = int(n_str)
    else:
        return None
    return {"id1": id1, "id2": id2, "date": date, "n": n,
            "pair_key": f"{id1}_{id2}_{date}"}


def load_rgb(path: Path) -> np.ndarray:
    arr = tifffile.imread(str(path))
    if arr.ndim == 4 and arr.shape[-1] == 3:
        arr = arr[0]
    if arr.ndim == 3:
        if arr.shape[-1] == 3:
            return arr
        if arr.shape[0] == 3:
            return np.moveaxis(arr, 0, -1)
    if arr.ndim == 2:
        raise ValueError(
            f"{path.name} is single-channel (shape {arr.shape}); "
            "an RGB image is required."
        )
    raise ValueError(
        f"{path.name}: unrecognised TIFF shape {arr.shape}. Expected one of "
        "(H,W,3), (3,H,W), (1,H,W,3), (N,H,W,3)."
    )


def process_image(tif_path: Path, sample_masks: dict, white_mask,
                  global_white=None, calib=None) -> tuple:
    """Returns (per-ROI dict of measurements, refl_calibrated bool)."""
    rgb = load_rgb(tif_path)
    if calib is not None:
        out = {name: measure_roi(rgb, mask, calib=calib)
               for name, mask in sample_masks.items()}
        return out, True
    if global_white is not None:
        white = global_white
    elif white_mask is not None:
        white = white_means(rgb, white_mask)
    else:
        white = None
    used = white is not None and not np.any(np.isnan(white))
    out = {name: measure_roi(rgb, mask, white if used else None)
           for name, mask in sample_masks.items()}
    return out, used


def find_white_reference(folder: Path):
    """Search folder and up to two parents for white_reference.json.
    Returns (array[R,G,B], path) or (None, None). One reference is enough
    when every session was shot under the same light."""
    for f in (folder, folder.parent, folder.parent.parent):
        p = f / "white_reference.json"
        if p.is_file():
            d = json.loads(p.read_text())
            return np.array([d["R"], d["G"], d["B"]], float), p
    return None, None


def load_calibration(path: Path):
    """Load calibration.json into {ch: (dn_asc, ref_asc)} for np.interp."""
    d = json.loads(Path(path).read_text())
    return {ch: (np.asarray(d[ch]["dn"], float),
                 np.asarray(d[ch]["ref"], float)) for ch in ("R", "G", "B")}


def find_calibration(folder: Path):
    """Search folder and up to two parents for calibration.json."""
    for f in (folder, folder.parent, folder.parent.parent):
        p = f / "calibration.json"
        if p.is_file():
            return load_calibration(p), p
    return None, None


def emit_white_reference(folder: Path, out_path: Path, roi_name: str = None):
    """Measure the neutral patch once (averaged over this folder's images via
    a ROI in its zip) and write it as the global reference. `roi_name` picks
    the ROI explicitly; otherwise the first name in WHITE_ROI_NAMES is used."""
    tif_files = list_tifs(folder)
    if not tif_files:
        sys.exit(f"No TIFFs in {folder} to measure a white patch from.")
    masks = load_roi_masks(find_roi_zip(folder),
                           load_rgb(tif_files[0]).shape[:2])
    if roi_name:
        if roi_name not in masks:
            sys.exit(f"ROI {roi_name!r} not in zip. Available: "
                     f"{sorted(masks)}")
        chosen = roi_name
    else:
        chosen = next((n for n in masks if n in WHITE_ROI_NAMES), None)
        if not chosen:
            sys.exit(f"No white ROI (names {WHITE_ROI_NAMES}) in the zip. "
                     f"Pass --white-roi with one of: {sorted(masks)}")
    wm = masks[chosen]
    vals = []
    for tif in tif_files:
        m = white_means(load_rgb(tif), wm)
        if m is not None and not np.any(np.isnan(m)):
            vals.append(m)
    if not vals:
        sys.exit("White ROI produced no valid pixels.")
    R, G, B = np.mean(vals, axis=0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(
        {"R": float(R), "G": float(G), "B": float(B),
         "n_images": len(vals), "roi": chosen,
         "source_folder": str(folder)}, indent=2))
    print(f"White reference (R,G,B) = {R:.1f}, {G:.1f}, {B:.1f} "
          f"from ROI {chosen!r} over {len(vals)} image(s)")
    print(f"Wrote {out_path}")


def process_folder(target: Path, args, shared_roi_zip: Path = None) -> None:
    folder = find_tiffs_folder(target)

    tif_files = list_tifs(folder)
    if not tif_files:
        print(f"No TIFF files found in {folder}")
        return

    # Folder-local zip wins; otherwise use the shared zip from the top folder.
    roi_zip = single_zip(folder)
    if roi_zip is None and folder != target:
        roi_zip = single_zip(target)
    if roi_zip is None:
        roi_zip = shared_roi_zip
    if roi_zip is None:
        print(f"No ROI zip for {target.name} (none here or in the top folder)")
        return
    print(f"Target:      {target}")
    if folder != target:
        print(f"TIFFs:       {folder}")
    print(f"ROI zip:     {roi_zip}")
    print(f"Images:      {len(tif_files)}")
    print(f"Reflectance: over {REFLECTANCE_REGION} pixels")

    # Resolve reflectance calibration: explicit calibration.json > white ref.
    global_white = None
    calib = None
    if not args.no_white:
        if args.calibration:
            calib = load_calibration(Path(args.calibration))
            print(f"Calibration: {args.calibration} (explicit)")
        else:
            calib, cal_path = find_calibration(folder)
            if calib is not None:
                print(f"Calibration: {cal_path} (actual reflectance)")
        if calib is None:
            if args.white_ref:
                d = json.loads(Path(args.white_ref).read_text())
                global_white = np.array([d["R"], d["G"], d["B"]], float)
                print(f"White ref:   {args.white_ref} (explicit)")
            else:
                global_white, gw_path = find_white_reference(folder)
                if global_white is not None:
                    print(f"White ref:   {gw_path} (global)")

    print()
    print("Loading ROI masks ...")
    img0 = load_rgb(tif_files[0])
    image_shape = img0.shape[:2]
    del img0

    masks = load_roi_masks(roi_zip, image_shape)

    # A per-folder white ROI is only a fallback when no global reference exists.
    white_name = next((n for n in masks if n in WHITE_ROI_NAMES), None)
    white_mask = masks.pop(white_name) if white_name else None
    sample_masks = masks
    roi_names = list(sample_masks.keys())
    print(f"Sample ROIs: {roi_names}")
    if calib is None and global_white is None and white_mask is None \
            and not args.no_white:
        print("Reflectance: raw intensities (no calibration or white ref)")
    print()

    groups = defaultdict(list)
    unparseable, errors = [], []
    any_white_used = False
    for tif in tqdm(tif_files, desc="Processing"):
        parsed = parse_stem(tif.stem)
        if parsed is None:
            unparseable.append(tif.name)
            continue
        try:
            rois, used = process_image(tif, sample_masks, white_mask,
                                       global_white=global_white, calib=calib)
            any_white_used = any_white_used or used
        except Exception as e:
            errors.append((tif.name, str(e)))
            rois = {name: {"coverage": np.nan, "kelp_px": np.nan,
                           "refl_R": np.nan, "refl_G": np.nan, "refl_B": np.nan,
                           "refl_bright": np.nan, "refl_rg": np.nan}
                    for name in roi_names}
            used = False
        groups[parsed["pair_key"]].append({
            "n": parsed["n"], "filename": tif.name, "rois": rois,
            "white_used": used,
            "id1": parsed["id1"], "id2": parsed["id2"], "date": parsed["date"],
        })

    if unparseable:
        print(f"\nWARNING: {len(unparseable)} file(s) did not match "
              f"ID1_ID2_DATE[_N]:")
        for name in unparseable:
            print(f"  {name}")
    if errors:
        print(f"\nWARNING: {len(errors)} file(s) errored during processing:")
        for name, err in errors:
            print(f"  {name}: {err}")

    refl_keys = ["refl_R", "refl_G", "refl_B", "refl_bright", "refl_rg"]

    rows = []
    for pair_key in sorted(groups.keys()):
        photos = sorted(groups[pair_key], key=lambda d: d["n"])
        n_photos = len(photos)
        id1, id2, date = photos[0]["id1"], photos[0]["id2"], photos[0]["date"]
        white_used = all(p["white_used"] for p in photos)

        drop = set()
        if id1 == SENTINEL_ID:
            drop.add(ROI_FOR_LEFT_ID)
        if id2 == SENTINEL_ID:
            drop.add(ROI_FOR_RIGHT_ID)

        flags = []
        if n_photos not in EXPECTED_COUNTS:
            flags.append(f"unexpected_count:{n_photos}")
        if drop:
            flags.append("sentinel_dropped:" + ",".join(sorted(drop)))
        flag = ";".join(flags) if flags else "ok"

        row = {"pair_key": pair_key, "id1": id1, "id2": id2, "date": date,
               "n_photos": n_photos, "flag": flag,
               "white_ref_used": white_used}

        cov_all, area_all, bright_all, rg_all = [], [], [], []
        for name in roi_names:
            if name in drop:
                row[f"coverage_{name}"] = np.nan
                row[f"area_cm2_{name}"] = np.nan
                for k in refl_keys:
                    row[f"{k}_{name}"] = np.nan
                continue

            def mean_of(key):
                vals = [p["rois"][name][key] for p in photos
                        if not np.isnan(p["rois"][name][key])]
                return float(np.mean(vals)) if vals else np.nan

            cov = mean_of("coverage")
            kelp_px = mean_of("kelp_px")
            area = kelp_px * CM2_PER_PIXEL if not np.isnan(kelp_px) else np.nan
            row[f"coverage_{name}"] = cov
            row[f"area_cm2_{name}"] = area
            for k in refl_keys:
                row[f"{k}_{name}"] = mean_of(k)

            if not np.isnan(cov):
                cov_all.append(cov)
            if not np.isnan(area):
                area_all.append(area)
            if not np.isnan(row[f"refl_bright_{name}"]):
                bright_all.append(row[f"refl_bright_{name}"])
            if not np.isnan(row[f"refl_rg_{name}"]):
                rg_all.append(row[f"refl_rg_{name}"])

        row["mean_coverage"] = float(np.mean(cov_all)) if cov_all else np.nan
        row["mean_area_cm2"] = float(np.mean(area_all)) if area_all else np.nan
        row["mean_refl_bright"] = float(np.mean(bright_all)) if bright_all else np.nan
        row["mean_refl_rg"] = float(np.mean(rg_all)) if rg_all else np.nan
        rows.append(row)

    df = pd.DataFrame(rows)
    front = ["pair_key", "id1", "id2", "date", "n_photos", "flag",
             "white_ref_used", "mean_coverage", "mean_area_cm2",
             "mean_refl_bright", "mean_refl_rg"]
    rest = [c for c in df.columns if c not in front]
    df = df[front + rest]

    print("\nPer-group results (summary columns):")
    print(df[front].to_string(index=False))

    print(f"\nGroups: {len(df)}  (flagged: {(df['flag'] != 'ok').sum()})")
    if calib is not None:
        print("Reflectance: calibrated to reference values")
    elif any_white_used:
        print("Reflectance: white-referenced (relative)")
    else:
        print("Reflectance: raw 0-255, NOT corrected; add calibration.json "
              "or a white reference to compare across sessions.")
    print(f"Folder mean coverage:    {df['mean_coverage'].mean()*100:.2f}%")
    print(f"Folder mean reflectance: {df['mean_refl_bright'].mean():.4f}")

    out_csv = target / "kelp_coverage_results.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nResults saved to: {out_csv}")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", help="a shoot folder, OR a directory of dated "
                                    "subfolders each holding TIFFs")
    ap.add_argument("--emit-white-ref", metavar="OUT.json",
                    help="measure the white patch here and write it as the "
                         "global reference, then exit (one-time)")
    ap.add_argument("--white-roi", metavar="NAME",
                    help="name of the neutral/white ROI in the calibration "
                         "zip (used with --emit-white-ref)")
    ap.add_argument("--calibration", metavar="CAL.json",
                    help="explicit calibration.json (measured->reference) to "
                         "apply to reflectance")
    ap.add_argument("--white-ref", metavar="REF.json",
                    help="explicit white_reference.json to apply")
    ap.add_argument("--no-white", action="store_true",
                    help="ignore any white reference; emit raw intensities")
    args = ap.parse_args()

    target = Path(args.folder).resolve()
    if not target.is_dir():
        print(f"Not a directory: {target}")
        sys.exit(1)

    if args.emit_white_ref:
        emit_white_reference(find_tiffs_folder(target),
                             Path(args.emit_white_ref).resolve(),
                             roi_name=args.white_roi)
        return

    # One shoot folder, or a directory of dated subfolders each holding TIFFs.
    if resolve_tiffs_folder(target) is not None:
        shoots = [target]
        shared_roi_zip = None
    else:
        shoots = [s for s in sorted(p for p in target.iterdir() if p.is_dir())
                  if resolve_tiffs_folder(s) is not None]
        if not shoots:
            sys.exit(f"No TIFFs in {target}, and no subfolder contains TIFFs.")
        shared_roi_zip = single_zip(target)   # the top-level shared ROI zip
        print(f"Found {len(shoots)} date folder(s) with TIFFs.")
        if shared_roi_zip is not None:
            print(f"Shared ROI zip: {shared_roi_zip.name}")
        print()

    failed = []
    for s in shoots:
        if len(shoots) > 1:
            print("=" * 60)
        try:
            process_folder(s, args, shared_roi_zip=shared_roi_zip)
        except Exception as e:
            failed.append((s.name, str(e)))
            print(f"[ERROR] {s.name}: {e}")
        print()

    if failed:
        print("Folders that failed (others still processed):")
        for name, err in failed:
            print(f"  {name}: {err}")


if __name__ == "__main__":
    main()