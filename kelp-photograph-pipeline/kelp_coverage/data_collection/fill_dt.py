"""fill_dt.py — fill the dt column in rename_map.csv from a name/created CSV.

Output has exactly the same columns and row order as rename_map.csv.
"""

import csv
from collections import Counter
from datetime import datetime
from pathlib import Path

BASE = Path("/Volumes/Expansion/Heat_stress_experiment_photos/07-05-2026")
DATE_CSV = BASE / "correct_dt_070526_dates_names.csv"
IN_CSV = BASE / "rename_map.csv"
OUT_CSV = BASE / "rename_map_dated.csv"

# True overwrites the existing date column where a timestamp is found.
# False leaves date untouched and only fills dt.
UPDATE_DATE = True


def build_date_map(csv_path=DATE_CSV):
    dates, conflicts, unparsed = {}, [], 0
    with csv_path.open(newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            name = row["name"].strip()
            if name.startswith("._"):
                continue
            try:
                dt = datetime.fromisoformat(row["created"].strip())
            except ValueError:
                unparsed += 1
                continue
            stem = Path(name).stem
            if stem in dates and dates[stem] != dt:
                conflicts.append(stem)
            dates[stem] = dt
    return dates, conflicts, unparsed


def main():
    dates, conflicts, unparsed = build_date_map()
    print(f"{len(dates)} dates loaded, {unparsed} unparsed, {len(conflicts)} stem conflicts")

    with IN_CSV.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if not rows:
        raise SystemExit("rename_map.csv is empty")

    missing, changed_date = [], Counter()
    for row in rows:
        stem = row["old_stem"].strip()
        dt = dates.get(stem)
        if dt is None:
            missing.append(stem)
            continue
        row["dt"] = dt.strftime("%Y%m%d_%H%M%S")
        new_date = dt.strftime("%Y%m%d")
        if new_date != row.get("date", "").strip():
            changed_date[(row.get("date", "").strip(), new_date)] += 1
        if UPDATE_DATE:
            row["date"] = new_date

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    filled = len(rows) - len(missing)
    print(f"{filled} filled, {len(missing)} missing of {len(rows)} rows")
    for (old, new), n in sorted(changed_date.items()):
        print(f"  date {old or '(blank)'} -> {new} on {n} rows")
    for stem in missing[:10]:
        print("  no date found:", stem)

    dup = [k for k, v in Counter(r["dt"] for r in rows if r["dt"].strip()).items() if v > 1]
    print(f"{len(dup)} timestamps shared by more than one row")

    unused = set(dates) - {r["old_stem"].strip() for r in rows}
    print(f"{len(unused)} names in the date CSV with no matching rename_map row")


if __name__ == "__main__":
    main()