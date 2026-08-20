"""
Backfill blank date/dt columns in a rename_map.csv from raw file EXIF.

Only rows with a blank date are touched; rows that already have a date are
left exactly as they are. A backup of the original CSV is written alongside
it as rename_map.csv.bak before anything changes.

Usage:
    python backfill_dates.py /path/to/shoot/rename_map.csv /path/to/raws

You can pass more than one raws folder if the ARWs are spread across
several shoots:

    python backfill_dates.py rename_map.csv /shoot1/raw /shoot2/raw

If a raw can't be found or has no EXIF DateTimeOriginal, the row is left
blank and reported at the end. As a last resort, pass --folder-date to use
the shoot folder's DD-MM-YYYY name for any row still blank:

    python backfill_dates.py rename_map.csv /shoot/raw --folder-date
"""

import sys
import csv
import shutil
from pathlib import Path

RAW_EXTS = (".ARW", ".arw")


def exif_datetime(path):
    import exifread
    try:
        with open(path, "rb") as f:
            tags = exifread.process_file(f, stop_tag="DateTimeOriginal")
        return str(tags.get("EXIF DateTimeOriginal", ""))  # "YYYY:MM:DD HH:MM:SS"
    except Exception:
        return ""


def find_raw(stem, raw_dirs):
    for d in raw_dirs:
        for ext in RAW_EXTS:
            candidate = d / f"{stem}{ext}"
            if candidate.exists():
                return candidate
    return None


def folder_date_of(raw_dir):
    """Parse DD-MM-YYYY from the raws folder's parent name, like the pipeline does."""
    from datetime import datetime
    try:
        return datetime.strptime(raw_dir.parent.name, "%d-%m-%Y").strftime("%Y%m%d")
    except ValueError:
        return ""


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    use_folder_date = "--folder-date" in sys.argv
    if len(args) < 2:
        sys.exit(__doc__)

    csv_path = Path(args[0]).expanduser().resolve()
    raw_dirs = [Path(a).expanduser().resolve() for a in args[1:]]
    if not csv_path.exists():
        sys.exit(f"CSV not found: {csv_path}")
    for d in raw_dirs:
        if not d.is_dir():
            sys.exit(f"Raws folder not found: {d}")

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    blank = [r for r in rows if not r.get("date", "").strip()]
    print(f"{len(rows)} rows total, {len(blank)} with a blank date.")
    if not blank:
        print("Nothing to do.")
        return

    backup = csv_path.with_suffix(".csv.bak")
    shutil.copy2(csv_path, backup)
    print(f"Backup written: {backup.name}")

    filled, no_raw, no_exif = 0, [], []
    for r in blank:
        raw = find_raw(r["old_stem"], raw_dirs)
        if raw is None:
            no_raw.append(r["old_stem"])
            continue
        dt = exif_datetime(raw)
        if dt:
            r["dt"] = dt
            r["date"] = dt.split()[0].replace(":", "")
            filled += 1
        else:
            # fall back to the shoot folder's date if allowed
            fd = folder_date_of(raw.parent) if use_folder_date else ""
            if fd:
                r["date"] = fd
                filled += 1
            else:
                no_exif.append(r["old_stem"])

    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"\nFilled {filled} of {len(blank)} blank rows. CSV updated: {csv_path}")
    if no_raw:
        print(f"\n{len(no_raw)} stems had no matching ARW in the given folder(s):")
        for s in no_raw[:20]:
            print(f"  {s}")
        if len(no_raw) > 20:
            print(f"  ... and {len(no_raw) - 20} more")
    if no_exif:
        print(f"\n{len(no_exif)} raws had no EXIF DateTimeOriginal "
              f"(re-run with --folder-date to use the folder name):")
        for s in no_exif[:20]:
            print(f"  {s}")
        if len(no_exif) > 20:
            print(f"  ... and {len(no_exif) - 20} more")


if __name__ == "__main__":
    main()