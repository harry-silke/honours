#!/usr/bin/env python3
"""
Stage 2: Image Registration
============================
Aligns every sample photo to a fixed reference image so the shared ROI map
lands correctly regardless of minor camera shifts between sessions.

For each time point, sample photos are warped to match the reference frame
and saved to a `registered/` subfolder inside the date folder. Stage 3
(run_imagej.py) reads from those registered subfolders.

Calibration images are NOT registered — they are left untouched in the
original date folders and processed separately.

Low-confidence registrations (inlier ratio below the config threshold) are
written to flagged_for_review.txt in the pipeline directory for manual check.

Usage:
    python register_images.py
    python register_images.py --timepoint 10.02.2026
    python register_images.py --dry-run
    python register_images.py --force   # re-register already-done images
"""

import argparse
import shutil
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Run: pip install pyyaml")
    sys.exit(1)

try:
    import cv2
    import numpy as np
except ImportError:
    print("ERROR: OpenCV not installed. Run: pip install opencv-python")
    sys.exit(1)


# ── Config ─────────────────────────────────────────────────────────────────────

def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# ── Image helpers ───────────────────────────────────────────────────────────────

def load_image(path: Path) -> np.ndarray:
    img = cv2.imread(str(path))
    if img is None:
        raise IOError(f"Could not read image: {path}")
    return img


def save_image(img: np.ndarray, path: Path, quality: int = 95) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ext = path.suffix.lower()
    if ext in (".jpg", ".jpeg"):
        cv2.imwrite(str(path), img, [cv2.IMWRITE_JPEG_QUALITY, quality])
    else:
        cv2.imwrite(str(path), img)


def is_calibration_image(filename: str, prefix: str) -> bool:
    return filename.lower().startswith(prefix.lower())


# ── ORB registration ────────────────────────────────────────────────────────────

def build_orb_reference(ref_img: np.ndarray, n_features: int):
    """Detect and compute ORB keypoints/descriptors on the reference image."""
    orb = cv2.ORB_create(n_features)
    ref_gray = cv2.cvtColor(ref_img, cv2.COLOR_BGR2GRAY)
    kp, des = orb.detectAndCompute(ref_gray, None)
    if des is None or len(kp) < 4:
        raise RuntimeError(
            "Could not detect enough keypoints in the reference image. "
            "Check that the reference photo is clear and well-lit."
        )
    return orb, ref_gray, kp, des


def register_orb(
    ref_img: np.ndarray,
    ref_kp,
    ref_des,
    orb,
    target_img: np.ndarray,
    lowe_ratio: float,
    ransac_threshold: float,
) -> tuple[np.ndarray | None, float]:
    """
    Register target_img to the reference frame using ORB + RANSAC homography.

    Returns (warped_image, confidence).
    confidence = inlier_count / good_match_count  (0.0 – 1.0)
    Returns (None, 0.0) if registration fails (too few matches).
    """
    target_gray = cv2.cvtColor(target_img, cv2.COLOR_BGR2GRAY)
    kp, des = orb.detectAndCompute(target_gray, None)

    if des is None or len(kp) < 4:
        return None, 0.0

    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=False)
    raw_matches = bf.knnMatch(ref_des, des, k=2)

    # Lowe's ratio test — keep only unambiguous matches
    good = [m for m, n in raw_matches if m.distance < lowe_ratio * n.distance]

    if len(good) < 4:
        return None, 0.0

    # Points in reference (source) and target (destination) image
    src_pts = np.float32([ref_kp[m.queryIdx].pt for m in good]).reshape(-1, 1, 2)
    dst_pts = np.float32([kp[m.trainIdx].pt for m in good]).reshape(-1, 1, 2)

    # Homography: maps target points → reference frame
    H, mask = cv2.findHomography(dst_pts, src_pts, cv2.RANSAC, ransac_threshold)

    if H is None:
        return None, 0.0

    inliers    = int(mask.ravel().sum())
    confidence = inliers / len(good)

    h, w   = ref_img.shape[:2]
    warped = cv2.warpPerspective(target_img, H, (w, h))

    return warped, confidence


# ── Phase correlation registration ──────────────────────────────────────────────

def register_phase(
    ref_img: np.ndarray,
    target_img: np.ndarray,
) -> tuple[np.ndarray, float]:
    """
    Register target_img to the reference using FFT phase correlation.
    Corrects pure translation only — no rotation or scale.

    Returns (warped_image, confidence).
    """
    ref_gray    = cv2.cvtColor(ref_img,    cv2.COLOR_BGR2GRAY).astype(np.float32)
    target_gray = cv2.cvtColor(target_img, cv2.COLOR_BGR2GRAY).astype(np.float32)

    (dx, dy), response = cv2.phaseCorrelate(ref_gray, target_gray)

    M      = np.float32([[1, 0, dx], [0, 1, dy]])
    h, w   = ref_img.shape[:2]
    warped = cv2.warpAffine(target_img, M, (w, h))

    confidence = float(np.clip(response, 0.0, 1.0))
    return warped, confidence


# ── File discovery ──────────────────────────────────────────────────────────────

def find_sample_photos(date_folder: Path, cal_prefix: str) -> list[Path]:
    """Return all JPEG/TIFF sample photos, excluding calibration images."""
    photos = []
    for ext in ("*.jpg", "*.jpeg", "*.tif", "*.tiff"):
        for p in sorted(date_folder.glob(ext)):
            if not is_calibration_image(p.name, cal_prefix):
                photos.append(p)
    return photos


# ── Flagged review log ──────────────────────────────────────────────────────────

def append_flagged(flag_path: Path, entries: list[str]) -> None:
    with open(flag_path, "a") as f:
        for line in entries:
            f.write(line + "\n")


# ── Main ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 2: register sample photos to a fixed reference frame"
    )
    parser.add_argument("--config",    default="config.yaml")
    parser.add_argument("--timepoint", default=None,
                        help="Process only this time point (e.g. 10.02.2026)")
    parser.add_argument("--dry-run",   action="store_true",
                        help="Report what would be done without writing any files")
    parser.add_argument("--force",     action="store_true",
                        help="Re-register images that already have a registered/ output")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: Config not found: {config_path}")
        sys.exit(1)

    config      = load_config(config_path)
    images_dir  = Path(config["paths"]["images_dir"])
    reg_cfg     = config.get("registration", {})
    cal_prefix  = config.get("calibration_filename_prefix", "calibration_")

    method           = reg_cfg.get("method", "orb").lower()
    threshold        = float(reg_cfg.get("confidence_threshold", 0.75))
    orb_n_features   = int(reg_cfg.get("orb_n_features", 2000))
    lowe_ratio       = float(reg_cfg.get("lowe_ratio", 0.75))
    ransac_threshold = float(reg_cfg.get("ransac_threshold", 5.0))
    reference_rel    = reg_cfg.get("reference_image")

    if not reference_rel:
        print("ERROR: registration.reference_image is not set in config.yaml.")
        print("       Set it to the path of the photo used to draw ROI_Map.zip,")
        print("       relative to images_dir (e.g. 10.03.2026/C25_C17_100326.jpg)")
        sys.exit(1)

    reference_path = images_dir / reference_rel
    if not reference_path.exists():
        print(f"ERROR: Reference image not found: {reference_path}")
        sys.exit(1)

    # -- Load reference image and prepare feature detector -----------------
    print(f"Loading reference image: {reference_path.name}")
    ref_img = load_image(reference_path)

    if method == "orb":
        orb, ref_gray, ref_kp, ref_des = build_orb_reference(ref_img, orb_n_features)
        print(f"  ORB: {len(ref_kp)} keypoints detected in reference image")
    elif method == "phase":
        orb = ref_kp = ref_des = ref_gray = None
        print("  Method: phase correlation (translation-only)")
    else:
        print(f"ERROR: Unknown registration method '{method}'. Use 'orb' or 'phase'.")
        sys.exit(1)

    print(f"  Confidence threshold: {threshold}")
    print()

    # -- Select time points ------------------------------------------------
    time_points = config.get("time_points", [])
    if args.timepoint:
        time_points = [tp for tp in time_points if tp["label"] == args.timepoint]
        if not time_points:
            print(f"ERROR: Time point '{args.timepoint}' not found in config.")
            sys.exit(1)

    flag_path   = config_path.parent / "flagged_for_review.txt"
    flagged     = []
    total_done  = 0
    total_skip  = 0
    total_flag  = 0
    total_fail  = 0

    for tp in time_points:
        label       = tp["label"]
        date_folder = images_dir / label

        if not date_folder.exists():
            print(f"[SKIP] {label} — folder not found: {date_folder}")
            continue

        photos = find_sample_photos(date_folder, cal_prefix)
        if not photos:
            print(f"[SKIP] {label} — no sample photos found")
            continue

        registered_dir = date_folder / "registered"
        print(f"── {label}  ({len(photos)} photos → {registered_dir})")

        for photo in photos:
            out_path = registered_dir / photo.name

            # Skip already-registered images unless --force
            if out_path.exists() and not args.force:
                print(f"  [SKIP] {photo.name}  (already registered; use --force to redo)")
                total_skip += 1
                continue

            # The reference image is copied as-is (confidence = 1.0)
            if photo.resolve() == reference_path.resolve():
                if args.dry_run:
                    print(f"  [REF ] {photo.name}  → copy (reference image, confidence=1.00)")
                else:
                    out_path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(photo, out_path)
                    print(f"  [REF ] {photo.name}  confidence=1.00")
                total_done += 1
                continue

            # Register
            try:
                target_img = load_image(photo)

                if method == "orb":
                    warped, confidence = register_orb(
                        ref_img, ref_kp, ref_des, orb,
                        target_img, lowe_ratio, ransac_threshold,
                    )
                else:
                    warped, confidence = register_phase(ref_img, target_img)

            except Exception as e:
                print(f"  [FAIL] {photo.name}  — {e}")
                flagged.append(f"FAIL\t{label}/{photo.name}\t{e}")
                total_fail += 1
                continue

            if warped is None:
                print(f"  [FAIL] {photo.name}  — too few matches for homography")
                flagged.append(f"FAIL\t{label}/{photo.name}\ttoo few matches")
                total_fail += 1
                continue

            conf_str = f"confidence={confidence:.2f}"

            if confidence < threshold:
                status = "FLAG"
                flagged.append(f"LOW_CONFIDENCE\t{label}/{photo.name}\t{conf_str}")
                total_flag += 1
            else:
                status = "OK  "

            if args.dry_run:
                print(f"  [{status}] {photo.name}  {conf_str}")
            else:
                save_image(warped, out_path)
                print(f"  [{status}] {photo.name}  {conf_str}")

            total_done += 1

        print()

    # -- Write flagged_for_review.txt --------------------------------------
    if flagged and not args.dry_run:
        append_flagged(flag_path, flagged)
        print(f"Flagged images written to: {flag_path}")
    elif flagged and args.dry_run:
        print(f"(dry-run) Would append {len(flagged)} entries to: {flag_path}")

    # -- Summary -----------------------------------------------------------
    print(f"{'='*50}")
    if args.dry_run:
        print("Dry run complete — no files written.")
    else:
        print("Registration complete.")
    print(f"  Registered : {total_done}")
    print(f"  Skipped    : {total_skip}  (already done)")
    print(f"  Flagged    : {total_flag}  (low confidence — check flagged_for_review.txt)")
    print(f"  Failed     : {total_fail}  (could not register)")

    if not args.dry_run and total_done > 0:
        print("\nRegistered images are in: {date_folder}/registered/")
        print("Next: python run_imagej.py")


if __name__ == "__main__":
    main()
