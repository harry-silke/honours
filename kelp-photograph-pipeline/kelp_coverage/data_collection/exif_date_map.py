"""exif_date_map.py — build the full date map from EXIF alone, and compare."""

import csv
import json
import re
import subprocess
from collections import Counter
from datetime import datetime
from pathlib import Path

RAW_ROOT = Path("/Volumes/Expansion/Heat_stress_experiment_photos/07-05-2026")
MAP_CSV = Path("/Volumes/Expansion/Heat_stress_experiment_photos/07-05-2026/rename_map.csv")
INFERRED_CSV = Path("/Volumes/Expansion/Heat_stress_experiment_photos/07-05-2026/rename_map_dated.csv")


def exif_map(root=RAW_ROOT):
    out = subprocess.run(
        ["exiftool", "-r", "-ext", "ARW", "-T", "-FileName", "-DateTimeOriginal", str(root)],
        capture_output=True, text=True, check=True,
    )
    dates, conflicts, unreadable = {}, [], []
    for line in out.stdout.splitlines():
        name, _, raw = line.partition("\t")
        name, raw = name.strip(), raw.strip()
        if name.startswith("._"):
            continue
        stem = Path(name).stem
        try:
            dt = datetime.strptime(raw, "%Y:%m:%d %H:%M:%S")
        except ValueError:
            unreadable.append(stem)
            continue
        if stem in dates and dates[stem] != dt:
            conflicts.append(stem)
        dates[stem] = dt
    return dates, conflicts, unreadable


def main():
    dates, conflicts, unreadable = exif_map()
    print(f"{len(dates)} frames read, {len(conflicts)} stem conflicts, "
          f"{len(unreadable)} unreadable")

    for day, n in sorted(Counter(d.strftime("%Y%m%d") for d in dates.values()).items()):
        print(f"  {day}: {n}")

    Path("exif_dates.json").write_text(
        json.dumps({k: v.isoformat() for k, v in dates.items()}, indent=2)
    )

    if not INFERRED_CSV.exists():
        return

    agree = disagree = absent = 0
    for row in csv.DictReader(INFERRED_CSV.open()):
        dt = dates.get(row["old_stem"].strip())
        if dt is None:
            absent += 1
        elif dt.strftime("%Y%m%d") == row["date"]:
            agree += 1
        else:
            disagree += 1
            if disagree <= 10:
                print(f"  disagree: {row['old_stem']} inferred {row['date']} "
                      f"exif {dt:%Y%m%d} ({row['source']})")

    print(f"\nvs inferred: {agree} agree, {disagree} disagree, {absent} no EXIF")


if __name__ == "__main__":
    main()