#!/usr/bin/env python3
"""
Master Pipeline Runner
=======================
Runs the kelp reflectance pipeline in sequence:

  Stage 1 — run_imagej.py         Measure reflectance (pure Python, no Fiji)
  Stage 2 — aggregate_results.py  Consolidate CSVs into final dataframe

Usage:
    python run_pipeline.py
    python run_pipeline.py --timepoint 10.03.2026
    python run_pipeline.py --only-stage 1
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
        "name":    "Measure reflectance",
        "script":  "run_imagej.py",
        "tp_flag": True,
    },
    {
        "number":  2,
        "name":    "Aggregate results",
        "script":  "aggregate_results.py",
        "tp_flag": False,
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
        description="Run the kelp reflectance pipeline (Stages 1–2)"
    )
    parser.add_argument("--config",     default="config.yaml")
    parser.add_argument("--timepoint",  default=None,
                        help="Process only this time point (e.g. 10.03.2026). "
                             "Stage 4 always aggregates everything.")
    parser.add_argument("--only-stage", type=int, default=None, metavar="N",
                        help="Run only stage N (1 or 2)")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Pass --dry-run to each stage (no files written)")
    args = parser.parse_args()

    if not Path(args.config).exists():
        print(f"ERROR: Config not found: {args.config}")
        sys.exit(1)

    only_stage = args.only_stage
    if only_stage is not None and only_stage not in (1, 2):
        print("ERROR: --only-stage must be 1 or 2")
        sys.exit(1)

    stages_to_run = [
        s for s in STAGES
        if only_stage is None or s["number"] == only_stage
    ]

    print("=" * 60)
    print("Kelp Reflectance Pipeline")
    if args.timepoint:
        print(f"Time point : {args.timepoint}")
    print(f"Stages     : {', '.join(str(s['number']) for s in stages_to_run)}")
    if args.dry_run:
        print("Mode       : DRY RUN — no files will be written")
    print("=" * 60)
    print()

    results        = {}
    pipeline_start = time.time()

    for stage in stages_to_run:
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
    for stage in stages_to_run:
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
        print("Your aggregated CSV is in results/")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()