#!/usr/bin/env python3
"""
Stage 2b: Decompress TIFFs for Fiji compatibility
===================================================
Converts compressed 48-bit TIFFs in registered/ to uncompressed equivalents
in registered_uncompressed/, which vanilla ImageJ can open without Bio-Formats.

Run this between Stage 2 (register_images.py) and Stage 3 (run_imagej.py).

Usage:
    python decompress_tifs.py
    python decompress_tifs.py --timepoint 10.03.2026
    python decompress_tifs.py --force   # re-process even if output exists
"""

import argparse
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Run: pip install pyyaml")
    sys.exit(1)

try:
    import tifffile
except ImportError:
    print("ERROR: tifffile not installed. Run: pip install tifffile")
    sys.exit(1)


def decompress_tif(src: Path, dst: Path) -> None:
    img = tifffile.imread(src)
    tifffile.imwrite(dst, img, compression=None, photometric="rgb")


def process_timepoint(label: str, images_dir: Path, force: bool) -> tuple[int, int, int]:
    registered = images_dir / label / "registered"
    out_dir    = images_dir / label / "registered_uncompressed"

    if not registered.exists():
        print(f"  [SKIP] registered/ not found — run Stage 2 first: {registered}")
        return 0, 1, 0

    tifs = [
        f for f in registered.iterdir()
        if f.suffix.lower() in (".tif", ".tiff")
    ]

    if not tifs:
        print(f"  [SKIP] No TIFFs found in {registered}")
        return 0, 1, 0

    out_dir.mkdir(parents=True, exist_ok=True)

    ok = fail = skipped = 0
    for src in sorted(tifs):
        dst = out_dir / src.name
        if dst.exists() and not force:
            print(f"  [SKIP] Already exists: {dst.name}")
            skipped += 1
            continue
        try:
            decompress_tif(src, dst)
            print(f"  [OK]   {src.name}")
            ok += 1
        except Exception as e:
            print(f"  [FAIL] {src.name}: {e}")
            fail += 1

    return ok, skipped, fail


def main():
    parser = argparse.ArgumentParser(
        description="Stage 2b: decompress TIFFs for Fiji headless compatibility"
    )
    parser.add_argument("--config",     default="config.yaml")
    parser.add_argument("--timepoint",  default=None,
                        help="Process only this time point (e.g. 10.03.2026)")
    parser.add_argument("--force",      action="store_true",
                        help="Re-process files even if output already exists")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: Config not found: {config_path}")
        sys.exit(1)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    images_dir  = Path(config["paths"]["images_dir"])
    time_points = config.get("time_points", [])

    if args.timepoint:
        time_points = [tp for tp in time_points if tp["label"] == args.timepoint]
        if not time_points:
            print(f"ERROR: Time point '{args.timepoint}' not found in config.")
            sys.exit(1)

    total_ok = total_skip = total_fail = 0

    for tp in time_points:
        label = tp["label"]
        print(f"── {label}")
        ok, skip, fail = process_timepoint(label, images_dir, args.force)
        total_ok   += ok
        total_skip += skip
        total_fail += fail
        print()

    print("=" * 50)
    print(f"Stage 2b complete.")
    print(f"  Converted : {total_ok}")
    print(f"  Skipped   : {total_skip}")
    print(f"  Failed    : {total_fail}")

    if total_ok > 0:
        print("\nNext: python run_imagej.py")

    if total_fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
