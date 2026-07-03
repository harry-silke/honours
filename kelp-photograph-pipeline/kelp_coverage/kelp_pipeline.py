"""
kelp_pipeline.py
----------------
One command for the kelp photo pipeline. Runs, in order:

  1. build_metadata_from_master.py  workbook -> <root>/crosses.csv
  2. measure_kelp_coverage.py       every date folder -> per-folder
                                     kelp_coverage_results.csv (averaged) and
                                     kelp_photo_results.csv (per photo per side)
  3. stitch_photos.py               -> <root>/master_photos.csv   (per-photo story)
  4. stitch_results.py              -> <root>/master_results.csv  (averaged, for R)

The per-photo master_photos.csv is the main output: one row per photo per side,
with treatment, female population, and parent cultures joined on.

Graphing (analyse_kelp.R) is a separate step because it needs R, not Python;
run it after this from the folder holding master_results.csv and crosses.csv.

All five scripts must sit in the same folder as this one. calibration.json (and
optionally white_reference.json) should already be in <root> so reflectance is
calibrated; without it reflectance is raw and not comparable across sessions.

Usage
-----
    python kelp_pipeline.py <root> --workbook Heat_stress_experiment-3.xlsx

    <root>            the top photos folder, e.g.
                      /Volumes/Expansion/Heat_stress_experiment_photos
    --workbook PATH   the master .xlsx (needed unless --skip-build)
    --skip-build      reuse an existing crosses.csv
    --skip-measure    reuse existing per-folder result CSVs
    --skip-averaged   do not build master_results.csv (per-photo only)
"""

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def run(step, cmd):
    print(f"\n{'=' * 60}\n{step}\n{'=' * 60}")
    print("$ " + " ".join(str(c) for c in cmd))
    r = subprocess.run([sys.executable, *[str(c) for c in cmd]])
    if r.returncode != 0:
        sys.exit(f"\n[STOP] {step} failed (exit {r.returncode}).")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--workbook")
    ap.add_argument("--skip-build", action="store_true")
    ap.add_argument("--skip-measure", action="store_true")
    ap.add_argument("--skip-averaged", action="store_true")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        sys.exit(f"Not a directory: {root}")
    crosses = root / "crosses.csv"

    if not args.skip_build:
        if not args.workbook:
            sys.exit("--workbook is required unless --skip-build.")
        if not Path(args.workbook).is_file():
            sys.exit(f"Workbook not found: {args.workbook}")
        run("1/4  build population metadata",
            [HERE / "build_metadata_from_master.py", args.workbook,
             "--out", root])
    elif not crosses.is_file():
        sys.exit(f"--skip-build but {crosses} is missing.")

    if not (root / "calibration.json").is_file():
        print(f"\n[warn] no calibration.json in {root}; reflectance will be "
              "raw, not comparable across sessions.")

    if not args.skip_measure:
        run("2/4  measure photos", [HERE / "measure_kelp_coverage.py", root])

    run("3/4  stitch per-photo master (master_photos.csv)",
        [HERE / "stitch_photos.py", root, "--crosses", crosses])

    if not args.skip_averaged:
        run("4/4  stitch averaged master (master_results.csv)",
            [HERE / "stitch_results.py", root])

    print(f"\n{'=' * 60}\nDone.")
    print(f"  per-photo story : {root / 'master_photos.csv'}")
    if not args.skip_averaged:
        print(f"  averaged (for R): {root / 'master_results.csv'}")
    print("\nNext, to draw the graphs, from the folder with master_results.csv "
          "and crosses.csv:\n  Rscript analyse_kelp.R")


if __name__ == "__main__":
    main()