#!/usr/bin/env python3
"""
Add a date segment to already-renamed TIFFs.

    python add_date.py /Volumes/Expansion/Heat_stress_experiment_photos/29-05-2026/tiffs

Renames  127_128_1.tif  ->  127_128_20260529_1.tif
The date is parsed from the parent folder name (DD-MM-YYYY).

Pass --dry-run to preview without changing anything.

python add_date.py /Volumes/Expansion/Heat_stress_experiment_photos/DD-MM-YYYY/tiffs
"""

import re
import sys
from datetime import datetime
from pathlib import Path

PATTERN = re.compile(r'^(\w+_\w+)__(\d+)(\.tiff?)$', re.IGNORECASE)

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry_run = "--dry-run" in sys.argv

    if not args:
        sys.exit("Usage: python add_date.py <tiffs_folder> [--dry-run]")

    tiff_dir = Path(args[0]).expanduser().resolve()
    if not tiff_dir.is_dir():
        sys.exit(f"Not a directory: {tiff_dir}")

    try:
        date_str = datetime.strptime(tiff_dir.parent.name, "%d-%m-%Y").strftime("%Y%m%d")
    except ValueError:
        sys.exit(f"Could not parse date from folder name: {tiff_dir.parent.name!r} (expected DD-MM-YYYY)")

    tiffs = sorted(tiff_dir.glob("*.tif")) + sorted(tiff_dir.glob("*.tiff"))
    if not tiffs:
        sys.exit(f"No TIF files found in {tiff_dir}")

    renamed = skipped = 0
    for path in tiffs:
        m = PATTERN.match(path.name)
        if not m:
            print(f"skip (unexpected name): {path.name}")
            skipped += 1
            continue

        pair, photo, ext = m.group(1), m.group(2), m.group(3)
        new_name = f"{pair}_{date_str}_{photo}{ext}"
        target = path.with_name(new_name)

        if target.exists():
            print(f"skip (already exists): {new_name}")
            skipped += 1
            continue

        print(f"{path.name}  ->  {new_name}")
        if not dry_run:
            path.rename(target)
        renamed += 1

    print(f"\n{'Would rename' if dry_run else 'Renamed'} {renamed}, skipped {skipped}.")
    if dry_run:
        print("Re-run without --dry-run to apply.")


if __name__ == "__main__":
    main()