#!/usr/bin/env python3
"""
Standalone repeatability analysis — multi-shot kelp coverage
=============================================================
Point this at a folder of TIFs with repeated shots and it will:
  1. Save a ROI placement QC image for each timepoint (like test_roi_placement.py)
  2. Analyse how much kelp coverage varies between shots of the same sample
  3. Show how precision improves by averaging more shots together

Expects filenames like: C1_C2_100326_1.tif, C1_C2_100326_2.tif etc.

Usage:
    python test_repeatability.py
"""

import re
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import tifffile
import roifile
from skimage.draw import polygon

# ── Point these at your data ──────────────────────────────────────────────────
IMAGES_DIR = "/Users/harrysilke/Desktop/UNI/Honours/imas_hons/kelp-photograph-pipeline/test_coverage/photos"

# One ROI zip per timepoint
ROI_ZIPS = {
    "230426": "/Users/harrysilke/Desktop/UNI/Honours/imas_hons/kelp-photograph-pipeline/test_coverage/ROI_maps/ROI_Map_230426.zip"
}

# One representative image per timepoint for the ROI placement check
ROI_CHECK_IMAGES = {
    "230426": "/Users/harrysilke/Desktop/UNI/Honours/imas_hons/kelp-photograph-pipeline/test_coverage/photos/C336_C323_230426_1.tif"
}


THRESHOLD       = 1.05
MIN_GREEN_VALUE = 30
MIN_SHOTS       = 2
# ─────────────────────────────────────────────────────────────────────────────

FILENAME_RE = re.compile(
    r"^(?P<id1>[A-Za-z0-9]+)_(?P<id2>[A-Za-z0-9]+)_(?P<date>\d{6})_(?P<shot>\d+)$",
    re.IGNORECASE,
)


# ── Image loading ─────────────────────────────────────────────────────────────

def load_image(path: Path) -> np.ndarray:
    img = tifffile.imread(str(path))
    if img.ndim == 4:
        img = img[0]
    if img.dtype != np.uint8:
        img = (img / 256).astype(np.uint8)
    if img.ndim == 3 and img.shape[0] == 3:
        img = np.moveaxis(img, 0, -1)
    return img


# ── ROI masks ─────────────────────────────────────────────────────────────────

def load_masks(roi_zip: str, shape: tuple) -> dict:
    h, w = shape[:2]
    rois = roifile.roiread(roi_zip)
    if not isinstance(rois, list):
        rois = [rois]
    masks = {}
    for roi in rois:
        if roi.roitype == roifile.ROI_TYPE.RECT:
            mask = np.zeros((h, w), dtype=bool)
            mask[max(0, int(roi.top)):min(h, int(roi.bottom)),
                 max(0, int(roi.left)):min(w, int(roi.right))] = True
        else:
            coords = roi.coordinates()
            if coords is None or len(coords) < 3:
                print(f"  [WARN] ROI '{roi.name}' skipped — no coordinates")
                continue
            # roifile returns (col, row) for ellipse ROIs — swap to (row, col)
            rows = coords[:, 1].astype(float)
            cols = coords[:, 0].astype(float)
            rr, cc = polygon(rows, cols, shape=(h, w))
            mask = np.zeros((h, w), dtype=bool)
            mask[rr, cc] = True
        masks[roi.name] = mask
    return masks


# ── Coverage measurement ──────────────────────────────────────────────────────

def measure_coverage(img: np.ndarray, mask: np.ndarray) -> float:
    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rg = np.where(g > MIN_GREEN_VALUE, r / g, 0.0)
    kelp = mask & (rg > THRESHOLD)
    return float(kelp.sum() / mask.sum()) if mask.sum() > 0 else 0.0


# ── ROI placement QC ──────────────────────────────────────────────────────────

def save_roi_check(date_short: str, roi_zip: str) -> None:
    """
    Save a two-panel QC image for one timepoint showing ROI placement
    and kelp classification, identical to test_roi_placement.py.
    """
    img_path = ROI_CHECK_IMAGES.get(date_short)
    if img_path is None:
        print(f"  [ROI CHECK] No representative image configured for {date_short}")
        return
    if not Path(img_path).exists():
        print(f"  [ROI CHECK] Image not found: {img_path}")
        return

    img   = load_image(Path(img_path))
    h, w  = img.shape[:2]
    masks = load_masks(roi_zip, img.shape)

    if not masks:
        print(f"  [ROI CHECK] No usable ROIs for {date_short}")
        return

    r = img[:, :, 0].astype(float)
    g = img[:, :, 1].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rg = np.where(g > MIN_GREEN_VALUE, r / g, np.nan)

    combined_mask = np.zeros((h, w), dtype=bool)
    for mask in masks.values():
        combined_mask |= mask

    kelp_mask = combined_mask & (rg > THRESHOLD)
    overlay   = img.copy()
    overlay[kelp_mask] = [255, 80, 0]

    for name, mask in masks.items():
        n_roi  = int(mask.sum())
        n_kelp = int((mask & (rg > THRESHOLD)).sum())
        print(f"    ROI '{name}': {n_kelp/n_roi:.1%} kelp ({n_kelp}/{n_roi} px)")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f"ROI placement check — {date_short}  |  {Path(img_path).name}",
                 fontsize=11)

    axes[0].imshow(img)
    axes[0].set_title("Original + ROI outlines")
    for name, mask in masks.items():
        rows = np.where(mask.any(axis=1))[0]
        cols = np.where(mask.any(axis=0))[0]
        if len(rows) and len(cols):
            axes[0].add_patch(plt.Circle(
                ((cols[0] + cols[-1]) / 2, (rows[0] + rows[-1]) / 2),
                (cols[-1] - cols[0]) / 2,
                color="cyan", fill=False, linewidth=2
            ))
            axes[0].text(
                (cols[0] + cols[-1]) / 2,
                rows[0] - 20,
                name, color="cyan", fontsize=10, ha="center",
                bbox=dict(facecolor="black", alpha=0.4, pad=2, edgecolor="none")
            )
    axes[0].axis("off")

    n_kelp_total = int(kelp_mask.sum())
    n_roi_total  = int(combined_mask.sum())
    axes[1].imshow(overlay)
    axes[1].set_title(
        f"Kelp classification (orange)\n"
        f"{n_kelp_total} px, {n_kelp_total/n_roi_total:.1%} of ROI"
    )
    axes[1].axis("off")

    plt.tight_layout()
    out = Path(f"~/Desktop/roi_check_{date_short}.png").expanduser()
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close()
    print(f"    Saved: {out}")


# ── Load all images and group by sample ──────────────────────────────────────

def load_groups() -> dict:
    """
    Scan IMAGES_DIR for TIFs matching the _N suffix pattern.
    Returns {(id1, id2, date): {shot_n: mean_coverage_across_rois}}
    """
    groups    = {}
    roi_cache = {}

    tif_files = sorted(Path(IMAGES_DIR).glob("*.tif"))
    if not tif_files:
        print(f"No TIF files found in {IMAGES_DIR}")
        return groups

    print(f"Found {len(tif_files)} TIF files\n")

    for tif_path in tif_files:
        m = FILENAME_RE.match(tif_path.stem)
        if m is None:
            print(f"  [SKIP] {tif_path.name} — doesn't match expected pattern (need _1, _2 suffix)")
            continue

        date_short = m.group("date")
        roi_zip    = ROI_ZIPS.get(date_short)
        if roi_zip is None:
            print(f"  [SKIP] {tif_path.name} — no ROI zip configured for {date_short}")
            continue
        if not Path(roi_zip).exists():
            print(f"  [SKIP] {tif_path.name} — ROI zip not found: {roi_zip}")
            continue

        try:
            img = load_image(tif_path)
            cache_key = (roi_zip, img.shape[:2])
            if cache_key not in roi_cache:
                roi_cache[cache_key] = load_masks(roi_zip, img.shape)
            masks = roi_cache[cache_key]

            if not masks:
                print(f"  [SKIP] {tif_path.name} — no usable ROIs")
                continue

            coverages = [measure_coverage(img, mask) for mask in masks.values()]
            mean_cov  = float(np.mean(coverages))

            group_key = (
                m.group("id1").upper(),
                m.group("id2").upper(),
                date_short,
            )
            shot = int(m.group("shot"))
            groups.setdefault(group_key, {})[shot] = {
                "mean":    mean_cov,
                "per_roi": {name: cov for name, cov in zip(masks.keys(), coverages)},
            }
            print(f"  [OK]  {tif_path.name}  coverage={mean_cov:.1%}")

        except Exception as e:
            warnings.warn(f"Could not process {tif_path.name}: {e}")

    return groups


# ── Statistics ────────────────────────────────────────────────────────────────

def summarise(shots: dict) -> dict:
    values = np.array([v["mean"] for v in shots.values()])
    if len(values) < 2:
        return {"n": len(values), "mean": float(np.mean(values))}
    return {
        "n":      len(values),
        "mean":   float(np.mean(values)),
        "std":    float(np.std(values, ddof=1)),
        "cv":     float(np.std(values, ddof=1) / np.mean(values)) if np.mean(values) > 0 else np.nan,
        "range":  float(np.max(values) - np.min(values)),
        "values": values,
    }


def sem_curve(values: np.ndarray) -> list[tuple[int, float]]:
    return [
        (n, float(np.std(values[:n], ddof=1) / np.sqrt(n)))
        for n in range(2, len(values) + 1)
    ]


# ── Repeatability plots ───────────────────────────────────────────────────────

def plot_repeatability(groups: dict) -> None:
    stats = {k: summarise(v) for k, v in groups.items() if len(v) >= MIN_SHOTS}

    if not stats:
        print(f"\nNo samples with {MIN_SHOTS}+ repeated shots found.")
        print("Make sure filenames end with _1, _2, _3 etc.")
        return

    cvs    = [s["cv"] * 100 for s in stats.values() if "cv" in s and not np.isnan(s["cv"])]
    ranges = [s["range"] * 100 for s in stats.values() if "range" in s]

    by_n = {}
    for s in stats.values():
        if "values" not in s:
            continue
        for n, sem in sem_curve(s["values"]):
            by_n.setdefault(n, []).append(sem * 100)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Kelp coverage repeatability  |  {len(stats)} samples  |  "
        f"threshold R:G > {THRESHOLD}",
        fontsize=12
    )

    # Panel 1: CV distribution
    if cvs:
        axes[0, 0].hist(cvs, bins=15, color="steelblue", edgecolor="none")
        axes[0, 0].axvline(np.median(cvs), color="red", linewidth=1.5,
                           label=f"median = {np.median(cvs):.1f}%")
        axes[0, 0].set_xlabel("Coefficient of variation (%)")
        axes[0, 0].set_ylabel("Number of samples")
        axes[0, 0].set_title("Shot-to-shot variability (CV)\nlower = more repeatable")
        axes[0, 0].legend()

    # Panel 2: range distribution
    if ranges:
        axes[0, 1].hist(ranges, bins=15, color="coral", edgecolor="none")
        axes[0, 1].axvline(np.median(ranges), color="red", linewidth=1.5,
                           label=f"median = {np.median(ranges):.1f} pp")
        axes[0, 1].set_xlabel("Range between shots (percentage points)")
        axes[0, 1].set_ylabel("Number of samples")
        axes[0, 1].set_title("Max - min coverage across shots\nlower = more repeatable")
        axes[0, 1].legend()

    # Panel 3: SEM vs n shots
    if by_n:
        ns   = sorted(by_n.keys())
        meds = [np.median(by_n[n]) for n in ns]
        p25  = [np.percentile(by_n[n], 25) for n in ns]
        p75  = [np.percentile(by_n[n], 75) for n in ns]
        axes[1, 0].plot(ns, meds, "o-", color="steelblue", linewidth=2,
                        label="median SEM")
        axes[1, 0].fill_between(ns, p25, p75, alpha=0.2, color="steelblue",
                                label="IQR")
        axes[1, 0].set_xlabel("Number of shots averaged")
        axes[1, 0].set_ylabel("Standard error of mean (percentage points)")
        axes[1, 0].set_title("Precision gained by averaging more shots\nlower = more precise")
        axes[1, 0].set_xticks(ns)
        axes[1, 0].legend()

    # Panel 4: per-sample bar chart
    labels = [f"{k[0]}/{k[1]}\n{k[2]}" for k in sorted(stats.keys())]
    means  = [stats[k]["mean"] * 100 for k in sorted(stats.keys())]
    stds   = [stats[k].get("std", 0) * 100 for k in sorted(stats.keys())]
    x      = np.arange(len(labels))
    axes[1, 1].bar(x, means, yerr=stds, color="steelblue", alpha=0.7,
                   capsize=3, error_kw={"linewidth": 1})
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
    axes[1, 1].set_ylabel("Mean kelp coverage (%)")
    axes[1, 1].set_title("Per-sample mean ± SD across shots")

    plt.tight_layout()
    out = Path("~/Desktop/repeatability.png").expanduser()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nSaved: {out}")

    # Console summary
    print(f"\n{'Sample':<20} {'Date':<8} {'N':>3} {'Mean%':>7} {'SD%':>7} "
          f"{'CV%':>7} {'Range pp':>9}")
    print("-" * 60)
    for key in sorted(stats.keys()):
        s = stats[key]
        print(
            f"{key[0]+'/'+key[1]:<20} {key[2]:<8} {s['n']:>3} "
            f"{s['mean']*100:>7.1f} "
            f"{s.get('std', 0)*100:>7.1f} "
            f"{s.get('cv', 0)*100:>7.1f} "
            f"{s.get('range', 0)*100:>9.1f}"
        )
    if cvs:
        print(f"\nMedian CV    : {np.median(cvs):.1f}%")
    if ranges:
        print(f"Median range : {np.median(ranges):.1f} percentage points")
    print(f"Samples      : {len(stats)}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    # Step 1: ROI placement QC — one image per timepoint
    print("=" * 60)
    print("Step 1: ROI placement check")
    print("=" * 60)
    for date_short, roi_zip in ROI_ZIPS.items():
        print(f"\n── {date_short}")
        if not Path(roi_zip).exists():
            print(f"  [SKIP] ROI zip not found: {roi_zip}")
            continue
        save_roi_check(date_short, roi_zip)

    # Step 2: repeatability analysis
    print("\n" + "=" * 60)
    print("Step 2: Repeatability analysis")
    print("=" * 60 + "\n")
    groups = load_groups()
    if groups:
        plot_repeatability(groups)
        plt.show()