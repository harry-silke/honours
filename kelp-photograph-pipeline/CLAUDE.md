# Kelp Reflectance Image Analysis Pipeline

## Project Overview
This pipeline automates reflectance image analysis of kelp samples photographed as part of
an Honours experiment at the University of [University]. It wraps ImageJ/Fiji to perform
calibrated reflectance measurement across three colour channels (red, green, blue).

The pipeline processes standard RGB photos (JPEG/TIFF) taken under controlled lighting
with a calibration target. Each photo contains exactly 2 kelp samples arranged side by
side in a fixed, consistent layout.

## Researcher
Harry Silke — Honours student
Photos location: /Users/harrysilke/desktop/OneDrive/Honours/Kelp growth photos/
Pipeline location: /Users/harrysilke/desktop/OneDrive/Honours/kelp-pipeline/

## Experiment Structure
- Each time point = one date folder under "Kelp growth photos/"
- Each date folder contains:
  - ~71 sample photos (pilot) → ~500 photos (full experiment)
  - 1 calibration photo named: calibration_DDMMYY.jpg
- Each sample photo contains exactly 2 kelp samples, side by side, fixed layout
- Filename convention: sampleID1_sampleID2_DDMMYY.jpg
  - Example: C1_C2_100226.jpg → left sample = C1, right sample = C2
  - Treatment is encoded in the sample ID (e.g. C = control, T = treatment)
- Date folder format: DD.MM.YYYY (e.g. 10.02.2026)

## Photo Directory Structure
```
/Users/harrysilke/desktop/OneDrive/Honours/Kelp growth photos/
├── 10.02.2026/
│   ├── C1_C2_100226.jpg
│   ├── C3_C4_100226.jpg
│   ├── ...
│   └── calibration_100226.jpg
├── 18.02.2026/
│   ├── C1_C2_180226.jpg
│   └── calibration_180226.jpg
├── 24.02.2026/
└── 10.03.2026/
```

## Pipeline Directory Structure
```
/Users/harrysilke/desktop/OneDrive/Honours/kelp-pipeline/
├── CLAUDE.md                   # This file
├── config.yaml                 # Master config — only file edited regularly
├── scripts/
│   ├── generate_macros.py      # Stage 1: config → ImageJ macro generator
│   ├── register_images.py      # Stage 2: OpenCV image registration / ROI alignment
│   ├── run_imagej.py           # Stage 3: direct TIF measurement + per-group aggregation
│   └── run_pipeline.py         # Master runner
├── macros/                     # Auto-generated .ijm files — never edit manually
├── calibrations/               # ImageJ calibration .txt files (auto-generated)
├── roi_maps/
│   └── ROI_Map.zip             # ONE map, used for every time point forever
└── results/
    ├── raw/                    # Per-timepoint CSVs: <label>_results.csv
    └── master_results.csv      # Concatenation of all per-timepoint CSVs
```

## ROI Map
- One ROI map used for all time points: ROI_Map.zip
- Contains exactly 2 ROIs:
  - "left"  → corresponds to sampleID1 (first ID in filename)
  - "right" → corresponds to sampleID2 (second ID in filename)
- Side by side, fixed positions, never changes
- Located at: /Users/harrysilke/desktop/OneDrive/Honours/kelp-pipeline/roi_maps/ROI_Map.zip

## Calibration
- 6 known reflectance standards: 99.99, 59.10, 36.20, 22.25, 9.00, 3.10
- Function: 4th Degree Polynomial
- One calibration photo per time point: calibration_DDMMYY.jpg in each date folder
- Pipeline auto-detects calibration image, extracts mean gray values, generates .txt files
- Three calibration files produced per time point (red, green, blue channels)

## Key Dependencies
- Python 3.9+
- PyYAML:      pip install pyyaml
- pandas:      pip install pandas
- OpenCV:      pip install opencv-python
- Fiji (ImageJ) — installed at /Applications/Fiji.app

## Fiji Headless Command (macOS)
```bash
/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --headless --run macro.ijm
```

## Final Weekly Workflow (when pipeline is complete)
1. Add new date folder to "Kelp growth photos/" with sample photos + calibration_DDMMYY.jpg
2. Run: python scripts/run_pipeline.py
3. Check flagged_for_review.txt (images where auto-alignment confidence was low)
4. Clean aggregated CSV is ready for R/Python analysis

## Output Data Structure
One row per pair group (ID1_ID2_DATE), averaged across the 2–4 photos of
that pair. Columns:

| Block | Columns |
|---|---|
| Identifiers | `pair_key`, `id1`, `id2`, `date`, `n_photos`, `flag` |
| Calibrated stats | `{side}_{channel}_{stat}` — side ∈ left/right, channel ∈ red/green/blue, stat ∈ mean/stddev/min/modal |
| Kelp coverage | `left_kelp_coverage`, `right_kelp_coverage` |
| Kelp reflectance | `{side}_kelp_mean_reflectance_{channel}` |
| Kelp area | `left_kelp_area_cm2`, `right_kelp_area_cm2` |

`flag` is `"ok"` unless the group has an unexpected photo count (e.g.
`unexpected_count:1`) or contains a sentinel `1000` ID (e.g.
`sentinel_dropped:right`).

## Reading Output in R
```r
df <- read.csv("results/master_results.csv")
```

## Reading Output in Python
```python
import pandas as pd
df = pd.read_csv("results/master_results.csv")
```

## Current Build Status
- [x] Project structure created
- [x] CLAUDE.md written
- [x] config.yaml written
- [x] Stage 1: generate_macros.py — COMPLETE
- [x] Stage 2: register_images.py — COMPLETE
- [x] Stage 3: run_imagej.py — COMPLETE (also produces master_results.csv)
- [x] Master runner: run_pipeline.py — COMPLETE

## Development Notes
- Always validate calibration files exist before generating macros
- Registration confidence threshold defined in config.yaml — images below it go to flagged_for_review.txt
- Filename parser splits on underscore: parts[0]=sampleID1, parts[1]=sampleID2, parts[2]=date, parts[3]=photo number (optional, 1-4)
- All paths in generated macros must use forward slashes (ImageJ requirement, even on Mac)
- ROI label "left" → sampleID1 and "right" → sampleID2; run_imagej.py applies this when building per-side output columns
- Date folders use DD.MM.YYYY format — pipeline parses this when scanning for new time points
- A full `python run_imagej.py` rebuilds master_results.csv from scratch from that run's per-timepoint CSVs. `--timepoint X` updates only that date's rows in the existing master.
