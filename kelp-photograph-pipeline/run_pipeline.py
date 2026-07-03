#!/usr/bin/env python3
"""
Master Pipeline Runner
=======================
Runs the kelp reflectance pipeline:

  Stage 1 — run_imagej.py  Measure reflectance and kelp coverage, average
                           over each pair group, and write per-timepoint
                           CSVs + a master CSV across all timepoints.

(Stage 1 is currently the only stage. Earlier versions had a separate
aggregation stage; that's now folded into run_imagej.py.)

Usage:
    python run_pipeline.py
    python run_pipeline.py --timepoint 10.03.2026
    python run_pipeline.py --dry-run
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path


# ── Stage definitions ────────────────────────────────────────────────────────────

STAGES = [
    {
        "number":  1,
        "name":    "Measure reflectance and aggregate",
        "script":  "run_imagej.py",
        "tp_flag": True,
    },
]


# ── Runner ───────────────────────────────────────────────────────────────────────

def run_stage(stage: dict, extra_args: list[str]) -> bool:
    script = Path(stage["script"])
    if not script.exists():
        print(f"  [FAIL] Script not found: {script}")
        return False

    cmd    = [sys.executable, str(script)] + extra_args
    print(f"  $ {' '.join(cmd)}")

    t0      = time.time()
    result  = subprocess.run(cmd, text=True)
    elapsed = time.time() - t0

    if result.returncode != 0:
        print(f"  [FAIL] Stage {stage['number']} exited with code {result.returncode} "
              f"(after {elapsed:.1f}s)")
        return False

    print(f"  [OK]  Completed in {elapsed:.1f}s")
    return True


# ── Main ─────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run the kelp reflectance pipeline"
    )
    parser.add_argument("--config",     default="config.yaml")
    parser.add_argument("--timepoint",  default=None,
                        help="Process only this time point (e.g. 10.03.2026). "
                             "Updates only that date's rows in master_results.csv.")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Pass --dry-run to each stage (no files written)")
    args = parser.parse_args()

    if not Path(args.config).exists():
        print(f"ERROR: Config not found: {args.config}")
        sys.exit(1)

    print("=" * 60)
    print("Kelp Reflectance Pipeline")
    if args.timepoint:
        print(f"Time point : {args.timepoint}")
    print(f"Stages     : {', '.join(str(s['number']) for s in STAGES)}")
    if args.dry_run:
        print("Mode       : DRY RUN — no files will be written")
    print("=" * 60)
    print()

    results        = {}
    pipeline_start = time.time()

    for stage in STAGES:
        n    = stage["number"]
        name = stage["name"]
        print(f"── Stage {n}: {name}")

        extra = ["--config", args.config]
        if args.timepoint and stage["tp_flag"]:
            extra += ["--timepoint", args.timepoint]
        if args.dry_run:
            extra.append("--dry-run")

        success    = run_stage(stage, extra)
        results[n] = success
        print()

        if not success:
            print(f"Pipeline stopped at Stage {n}: {name}")
            print(f"Fix the error above, then re-run with --only-stage {n}")
            break

    elapsed = time.time() - pipeline_start
    print("=" * 60)
    print(f"Pipeline finished in {elapsed:.1f}s")
    print()
    for stage in STAGES:
        n = stage["number"]
        if n in results:
            status = "OK  " if results[n] else "FAIL"
            print(f"  [{status}] Stage {n}: {stage['name']}")
        else:
            print(f"  [----] Stage {n}: {stage['name']}  (not reached)")
    print()

    all_ok = all(results.values()) if results else False
    if all_ok:
        print("All stages completed successfully.")
        print("Per-timepoint CSVs and master_results.csv are in results/.")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
