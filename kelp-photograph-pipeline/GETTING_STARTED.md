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
│   ├── aggregate_results.py   ← Stage 4: collates all CSVs into one clean file
│   ├── register_images.py     ← Stage 2: NOT YET BUILT
│   ├── run_imagej.py          ← Stage 3: NOT YET BUILT
│   └── run_pipeline.py        ← Master runner: NOT YET BUILT
├── macros/                    ← auto-generated .ijm files land here
├── calibrations/              ← calibration .txt files land here
├── roi_maps/                  ← put your ROI_Map.zip here
└── results/
    └── raw/                   ← ImageJ CSVs land here
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

## What you can run right now (before Stages 2 & 3 are built)

### Generate macros
After running the calibration step in Fiji and saving the .txt files to `calibrations/`:
```bash
python scripts/generate_macros.py
```

### Aggregate results
After manually processing images in Fiji and saving CSVs to `results/raw/`:
```bash
python scripts/aggregate_results.py
```

Wide format (one row per sample per timepoint):
```bash
python scripts/aggregate_results.py --wide
```

---

## Your naming conventions

| Thing | Format | Example |
|---|---|---|
| Sample photos | sampleID1_sampleID2_DDMMYY.jpg | C1_C2_100226.jpg |
| Calibration photo | calibration_DDMMYY.jpg | calibration_100226.jpg |
| Date folders | DD.MM.YYYY | 10.02.2026 |
| Result CSVs (auto-named) | sampleID1_sampleID2_DDMMYY_channel.csv | C1_C2_100226_red.csv |

---

## Adding a new time point

Open `config.yaml` and add two lines at the bottom of the `time_points` list:

```yaml
  - label: 17.03.2026
    date_short: 170326
```

That's it.

---

## Final weekly workflow (when all stages are built)

```bash
# 1. Add new date folder with photos + calibration_DDMMYY.jpg
# 2. Add time point to config.yaml (2 lines)
# 3. Run:
python scripts/run_pipeline.py

# 4. Check flagged_for_review.txt if it exists
# 5. Your CSV is ready in results/
```

---

## Useful Claude Code prompts for building out the pipeline

> "Build Stage 2 — image registration using OpenCV ORB. Read CLAUDE.md first."

> "Build Stage 3 — headless Fiji batch runner. Each image should be processed silently."

> "Wire everything into run_pipeline.py so I can run the whole thing with one command."

> "My aggregated CSV has duplicate rows — here's a sample, can you debug it?"

> "Add a QC script that plots mean reflectance per channel over time as a line graph."
