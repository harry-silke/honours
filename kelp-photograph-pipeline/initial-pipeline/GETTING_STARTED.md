# Getting Started — Kelp Reflectance Pipeline
# Harry Silke, Honours

---

## What you have

```
kelp-pipeline/
├── CLAUDE.md                  ← project context (Claude Code reads this automatically)
├── config.yaml                ← the only file you edit regularly
├── scripts/
│   ├── generate_macros.py     ← Stage 1: generates Fiji macros from config
│   ├── register_images.py     ← Stage 2: image registration
│   ├── run_imagej.py          ← Stage 3: direct TIF measurement + per-group aggregation
│   └── run_pipeline.py        ← Master runner
├── macros/                    ← auto-generated .ijm files land here
├── calibrations/              ← calibration .txt files land here
├── roi_maps/                  ← put your ROI_Map.zip here
└── results/
    ├── raw/                   ← per-timepoint CSVs: <label>_results.csv
    └── master_results.csv     ← concatenated across all timepoints
```

---

## Step 1 — Move the pipeline folder

Move `kelp-pipeline/` to:
```
/Users/harrysilke/desktop/OneDrive/Honours/kelp-pipeline/
```

---

## Step 2 — Install Python dependencies

Open Terminal and run:
```bash
pip install pyyaml pandas opencv-python
```

---

## Step 3 — Create your ROI map in Fiji (one time only)

1. Open a representative sample photo in Fiji
2. Analyze > Tools > ROI Manager
3. ROI Manager > More > Options > tick "Use ROI names as labels"
4. Draw a region around the LEFT sample → Add → Rename to `left`
5. Draw a region around the RIGHT sample → Add → Rename to `right`
6. Select both ROIs (Shift+click) → More > Save → save as `ROI_Map.zip`
7. Put `ROI_Map.zip` in `kelp-pipeline/roi_maps/`

This is the only time you ever do this.

---

## Step 4 — Install Claude Code

```bash
curl -fsSL https://claude.ai/install.sh | sh
```

Verify:
```bash
claude --version
```

You need a Claude Pro, Max, or API account to use it.

---

## Step 5 — Start your first Claude Code session

```bash
cd /Users/harrysilke/desktop/OneDrive/Honours/kelp-pipeline
claude
```

Claude reads `CLAUDE.md` automatically and knows the full project context.

Tell it:
> "Let's build Stage 2 — image registration. Read CLAUDE.md first then plan the approach."

---

## What you can run right now

### Generate macros
After running the calibration step in Fiji and saving the .txt files to `calibrations/`:
```bash
python scripts/generate_macros.py
```

### Measure + aggregate
```bash
python scripts/run_imagej.py                       # all timepoints, rebuilds master CSV
python scripts/run_imagej.py --timepoint 10.03.2026 # one date only, updates that date's rows in master
```

---

## Your naming conventions

| Thing | Format | Example |
|---|---|---|
| Sample photos | sampleID1_sampleID2_DATE_N.tif (N optional, 1-4) | C1_C2_100226_2.tif |
| Calibration photo | calibration_DDMMYY.jpg | calibration_100226.jpg |
| Date folders | DD.MM.YYYY | 10.02.2026 |
| Per-timepoint CSV | <label>_results.csv | 10.03.2026_results.csv |
| Master CSV | master_results.csv | (in results/) |

---

## Adding a new time point

Open `config.yaml` and add two lines at the bottom of the `time_points` list:

```yaml
  - label: 17.03.2026
    date_short: 170326
```

That's it.

---

## Final weekly workflow

```bash
# 1. Add new date folder with photos + calibration_DDMMYY.jpg
# 2. Add time point to config.yaml (2 lines)
# 3. Run:
python scripts/run_pipeline.py

# 4. Check flagged_for_review.txt if it exists
# 5. Your master CSV is at results/master_results.csv
```

---

## Useful Claude Code prompts for building out the pipeline

> "Build Stage 2 — image registration using OpenCV ORB. Read CLAUDE.md first."

> "Build Stage 3 — headless Fiji batch runner. Each image should be processed silently."

> "Wire everything into run_pipeline.py so I can run the whole thing with one command."

> "My aggregated CSV has duplicate rows — here's a sample, can you debug it?"

> "Add a QC script that plots mean reflectance per channel over time as a line graph."
