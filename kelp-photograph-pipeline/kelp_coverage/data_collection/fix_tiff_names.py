"""
fix_tiff_names.py
-----------------
Repairs TIFF names after rename_map.csv dates were corrected (EXIF backfill).

The problem: renaming_pipeline.py apply renames old_stem.tif -> new_stem.tif.
Once run, the DSC*.tif files are gone, so a corrected CSV can't be re-applied:
apply just prints "no tiff for DSC...". This script bridges that gap by
reconstructing, for every row, what the TIFF is called NOW and what it SHOULD
be called, then renaming current -> corrected.

"Current name" is derived from the OLD CSV (rename_map.csv.bak - the backup
written by backfill_dates.py before it changed anything, holding the dates the
first apply run used). "Corrected name" comes from the fixed rename_map.csv.
TIFFs never renamed at all (still DSC*.tif, e.g. the manual block) are also
handled: their current name is just old_stem.tif.

Nothing is renamed unless every source file resolves and there are no target
collisions; the plan is printed first and --dry-run shows it without touching
anything.

Usage:
    python fix_tiff_names.py /path/to/shoot/rename_map.csv \
        /path/to/shoot/rename_map.csv.bak /path/to/shoot/tiffs --dry-run
    # review, then rerun without --dry-run

After this, rerun the measurement stage for the timepoint and rebuild
master_photos.csv.
"""

import sys
import csv
from pathlib import Path


def numbered_stems(rows):
    """old_stem -> new_stem using the pipeline's exact grouping/sort rules."""
    rs = [dict(r) for r in rows]
    for r in rs:
        r["pk"] = f"{r['left']}_{r['right']}_{r['date']}"
    rs.sort(key=lambda r: (r["pk"], r["dt"], r["old_stem"]))
    counters, out = {}, {}
    for r in rs:
        n = counters.get(r["pk"], 0) + 1
        counters[r["pk"]] = n
        out[r["old_stem"]] = f"{r['pk']}_{n}"
    return out


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry-run" in sys.argv
    if len(args) != 3:
        sys.exit(__doc__)

    new_csv, old_csv, tiff_dir = (Path(a).expanduser().resolve() for a in args)
    for p in (new_csv, old_csv):
        if not p.exists():
            sys.exit(f"Not found: {p}")
    if not tiff_dir.is_dir():
        sys.exit(f"Not a directory: {tiff_dir}")

    new_map = numbered_stems(read_rows(new_csv))   # old_stem -> corrected stem
    old_map = numbered_stems(read_rows(old_csv))   # old_stem -> current stem
                                                   # (as first apply named them)

    plan, already, missing = [], [], []
    for old_stem, target_stem in new_map.items():
        target = tiff_dir / f"{target_stem}.tif"
        # where is this photo on disk right now?
        candidates = [tiff_dir / f"{old_map.get(old_stem, '')}.tif",
                      tiff_dir / f"{old_stem}.tif"]
        current = next((c for c in candidates if c.name != ".tif" and c.exists()),
                       None)
        if current is None:
            missing.append(old_stem)
        elif current == target:
            already.append(old_stem)
        else:
            plan.append((current, target))

    # collision checks on the plan
    targets = [t for _, t in plan]
    dup = {t for t in targets if targets.count(t) > 1}
    clobber = [t for _, t in plan
               if t.exists() and t not in [c for c, _ in plan]]

    print(f"{len(new_map)} rows in corrected CSV")
    print(f"  already correctly named : {len(already)}")
    print(f"  to rename               : {len(plan)}")
    print(f"  source tiff not found   : {len(missing)}")
    if missing:
        print("  missing stems (first 15):")
        for s in missing[:15]:
            print(f"    {s}  (current guess: {old_map.get(s, s)}.tif)")
    if dup:
        sys.exit(f"STOP: duplicate targets in plan: {sorted(p.name for p in dup)}")
    if clobber:
        sys.exit("STOP: these targets already exist and are not themselves "
                 f"being renamed away: {sorted(p.name for p in clobber)}")

    if not plan:
        print("Nothing to do.")
        return

    print("\nPlan (first 20):")
    for cur, tgt in plan[:20]:
        print(f"  {cur.name}  ->  {tgt.name}")
    if len(plan) > 20:
        print(f"  ... and {len(plan) - 20} more")

    if dry:
        print("\nDRY RUN: nothing renamed. Re-run without --dry-run to apply.")
        return

    # two-phase rename via temp names so old->new chains can't collide midway
    tmp_pairs = []
    for i, (cur, tgt) in enumerate(plan):
        tmp = cur.with_name(f"__fixtmp_{i}.tif")
        cur.rename(tmp)
        tmp_pairs.append((tmp, tgt))
    renamed = 0
    for tmp, tgt in tmp_pairs:
        tmp.rename(tgt)
        renamed += 1
    print(f"\nRenamed {renamed} TIFFs.")
    print("Now rerun the measurement stage for this timepoint and rebuild "
          "master_photos.csv.")


if __name__ == "__main__":
    main()