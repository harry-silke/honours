# IMAS Honours — Giant Kelp Heat Tolerance & Outbreeding

Research code, pipelines and figures for an Honours project at the **Institute for
Marine and Antarctic Studies (IMAS), University of Tasmania**, studying heat
tolerance and outbreeding success in **giant kelp (*Macrocystis pyrifera*)** across
south-east Australia.

**Author:** Harry Silke · hisilke@utas.edu.au · IMAS, UTas

---

## About the project

Giant kelp forests along south-east Australia sit at the leading edge of the
warming **East Australian Current (EAC)** and are declining with marine
heatwaves. This project asks whether crossing kelp between and within populations
can produce offspring that are more heat tolerant, and whether that tolerance has
a detectable genetic basis.

The core experiment is a **crossing design** (a hub-and-spoke set of within- and
between-population crosses, including selfed lines such as Malacoota and LSI)
whose gametophyte/sporophyte offspring are grown and then put through a **heat
stress ramp** (control 14 °C vs. elevated 17 °C and 20 °C treatments, Mar–May
2026). Throughout the experiment each individual is tracked by:

- **Growth & health** from calibrated photographs (surface-area coverage and
  colour/reflectance as a proxy for condition), and
- **Photosynthetic performance** from **IPAM chlorophyll-fluorescence** (Fv/Fm).

These phenotypes are then linked back to genetics through a **GWAS** (survival as a
binary phenotype) and a **genotype–environment association (GEA)** analysis against
oceanographic conditions at each source site.

This repo holds the analysis code, helper scripts, and the figures/maps that go
into the thesis. It does **not** contain the raw image/sequencing data or most
generated outputs — those are gitignored (see [Data & outputs](#data--outputs)).

---

## Repository layout

| Path | What it is |
|------|------------|
| `kelp-photograph-pipeline/` | Photo-based growth & health pipeline (the main data-processing engine). Calibrates photos, registers them, measures per-individual kelp coverage and reflectance over time, and stitches everything into master metadata tables. |
| `IPAM_measurements/` | Analysis of IPAM Fv/Fm chlorophyll-fluorescence data — heat-stress time-courses and a resistance ranking per genotype. |
| `GWAS/` | Genome-wide association study pipeline for heat-stress survival (QC → variant calling → filtering → population structure → association → annotation → genomic prediction). Has its own detailed [README](GWAS/README.md). |
| `GEA/` | Genotype–environment association support: downloads oceanographic/environmental data (SST, physics, biogeochemistry, waves) for sequenced sites and builds site climatologies. |
| `ImageJ/` | Standalone ImageJ/Fiji macro (`Kelp_Growth_Health_Tracker.ijm`) plus user guide for HSV-based growth/health scoring of beaker-grid photos. |
| `RandFigures/` | R Markdown + Python for thesis figures: population maps, EAC/SST maps, crossing-design diagram, kelp lifecycle, nutrient profiles, and HOBO temperature-logger reports. |
| `PilotProject/` | Pilot-run master sheet (`pilot_master_sheet.xlsx`). |
| `experiment_timeline.py`, `build_simple_timeline.py` | Render the experiment timeline figure (temperature history per treatment, phases, imaging/IPAM/sampling timepoints). |
| `notes`, `todos.todo`, `HONS_methods_v2.docx` | Lab notes, task list, and the methods write-up. |
| `Experimental design.png`, `experimental_design_complete.png`, `kelp_EAC_SST_map*.png` | Key rendered figures kept at the repo root. |

---

## Components in detail

### `kelp-photograph-pipeline/` — growth & health from photos

The data engine of the project. Photos of kelp (two samples per photo, named
`sampleID1_sampleID2_DDMMYY`) are taken each session alongside a ColorChecker
calibration target.

- **`config.yaml`** — the single file you edit regularly: paths, calibration
  standards, registration settings, measurement options, and one block per
  imaging **time point**.
- **`kelp_coverage/kelp_pipeline.py`** — one command to run the whole chain:
  build crossing metadata from the master workbook → measure coverage &
  reflectance for every date folder → stitch into `master_photos.csv` (one row
  per photo per side) and `master_results.csv` (averaged, for R).
- **`kelp_coverage/measure_kelp_coverage.py`** — per-ROI kelp coverage via a
  red/green-ratio threshold, plus calibrated reflectance.
- **`kelp_coverage/data_collection/`** — collection-time utilities: EXIF/date
  backfilling and repair, TIFF naming fixes, PAM→TIFF conversion (`pim2tiff.py`),
  Fv/Fm stitching, and manual labelling helpers.
- **`initial-pipeline/`** — the earlier reflectance-only pipeline (macro
  generation → image registration → Fiji measurement), with its own
  `GETTING_STARTED.md`.
- **`calibrations/`, `macros/`, `roi_maps/`** — calibration `.txt` files,
  auto-generated Fiji `.ijm` macros, and ROI maps.
- **`test_coverage/`, `inspect_thresholds.py`, `test_*thresholds.py`** — scripts
  for tuning and sanity-checking the coverage thresholds.

Requires **Fiji/ImageJ** (path set in `config.yaml`) plus Python packages
`pyyaml pandas opencv-python` (and friends).

### `IPAM_measurements/` — photosynthetic performance

`analyse_fvfm.py` processes stitched Fv/Fm measurements into heat-stress
time-courses by population and a per-genotype **resistance ranking**. Outputs
live in `figures_fvfm/` (`fvfm_timecourse_by_population.png`,
`fvfm_resistance_ranking.png`, `resistance_ranking.csv`).

### `GWAS/` — genome-wide association

A full, self-contained GWAS pipeline for binary heat-stress survival against the
*Macrocystis pyrifera* Macpyr2 reference. Python + R modules orchestrated by
`run_gwas_pipeline.py` and configured via `GWAS/config.yaml`. See
**[GWAS/README.md](GWAS/README.md)** for installation, per-module usage, inputs,
outputs, and troubleshooting.

### `GEA/` — environmental data for genotype–environment association

`download_env_data_kelp_gea.py` pulls oceanographic variables for the sequenced
sites and produces:
- `env_climatology.csv` — long-term site means (for RDA / GEA), checked into the
  repo under `GEA/kelp_env_data/`.
- `env_at_sampling.csv` — conditions in a ±30-day window around each sampling
  event.

Sources include Copernicus Marine (OSTIA SST, GLORYS12 physics &
biogeochemistry, CMEMS waves) and BRAN2020 via NCI OPeNDAP. Needs a free
Copernicus Marine account (`copernicusmarine login`) and
`copernicusmarine xarray netCDF4 h5py numpy pandas tqdm`.

### `ImageJ/` — standalone macro

`Kelp_Growth_Health_Tracker.ijm` scores growth and health from 3×4 beaker-grid
photos using HSV colour analysis and ColorChecker calibration. See
`Kelp_Tracking_User_Guide.md` and `Quick_Start_Checklist.md`.

### `RandFigures/` — thesis figures & maps

R Markdown documents (rendered to `.html`) and Python for publication figures:
population distribution maps (`kelp_map*`, `kelp_populations_pub`), EAC/SST and
nutrient maps (`kelp_EAC_*`), the crossing/hub-and-spoke design
(`cross_design_map.Rmd`, `hub_and_spoke.py`), kelp lifecycle, and the HOBO
temperature-logger report (`hobo_logger_report.rmd`). Supporting data:
`oisst_annual_mean_seaus.csv`, `HONS_GKHE_hobo_logger_data.csv`.

---

## Getting started

Different components have different toolchains — set up only what you need.

```bash
# Photo pipeline
pip install pyyaml pandas opencv-python
# then edit kelp-photograph-pipeline/config.yaml and run:
python kelp-photograph-pipeline/kelp_coverage/kelp_pipeline.py <photos_root> \
    --workbook Heat_stress_experiment.xlsx

# GWAS (see GWAS/README.md for the full setup)
cd GWAS && bash install_dependencies.sh --method conda

# GEA env data
pip install copernicusmarine xarray netCDF4 h5py numpy pandas tqdm
copernicusmarine login
python GEA/download_env_data_kelp_gea.py

# Figures
# open the .Rmd files in RStudio and knit, or render the timeline:
python build_simple_timeline.py
```

---

## Data & outputs

The `.gitignore` deliberately excludes large/derived files: `results/`, `*.csv`
(except `GEA/kelp_env_data/env_climatology.csv`), `*.png`, `*.nc`, `*.tif`, and
`.DS_Store`. Raw photographs, sequencing data, and most generated tables/figures
therefore live outside the repo (e.g. on OneDrive / external drives referenced in
`config.yaml`), not in version control. A handful of final figures are committed
at the repo root for convenience.

---

## Status & notes

Active Honours work. `todos.todo` tracks outstanding tasks (e.g. running the GEA,
further growth-rate and PAM statistics, the Materials & Methods write-up) and
`notes` holds running lab notes. The GWAS module is the most polished and is
documented in full in its own README.
