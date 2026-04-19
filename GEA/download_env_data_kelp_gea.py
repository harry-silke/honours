"""
=============================================================================
Environmental Data Download Script for Giant Kelp GEA – SE Australia
=============================================================================
Sequenced sites only (24 unique sites, 27 sampling rows)
Outputs TWO CSV files:
  1. env_climatology.csv      – long-term means per unique site (for RDA)
  2. env_at_sampling.csv      – conditions in a ±30-day window around each
                                sampling event (one row per sampling visit)

Data sources:
  [A] OSTIA SST L4            – sea surface temperature (Copernicus)
  [B] GLORYS12 physics        – temp, salinity, MLD, currents, bottom temp (Copernicus)
  [C] GLORYS12 biogeochemistry– nitrate, phosphate, chlorophyll, O2 (Copernicus)
  [D] CMEMS wave reanalysis   – significant wave height, period, direction (Copernicus)
  [E] BRAN2020                – temp, salinity, MLD, currents (NCI OPeNDAP, free)

REQUIRED ACCOUNT
-----------------
  Copernicus Marine (free): https://data.marine.copernicus.eu
  Run once:  copernicusmarine login

REQUIRED PACKAGES
-----------------
  pip install copernicusmarine xarray netCDF4 h5py numpy pandas tqdm

USAGE
-----
  1. copernicusmarine login
  2. python download_env_data_kelp_gea.py
=============================================================================
"""

import warnings
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from datetime import timedelta
from tqdm import tqdm

warnings.filterwarnings("ignore")

# =============================================================================
# SEQUENCED SAMPLING ROWS — 27 rows across 24 unique sites
# =============================================================================

SAMPLING_ROWS = [
    {"site": "SIS", "location": "Sisters Bay",              "lat": -43.4260, "lon": 147.0233, "date": "2024-06-25"},
    {"site": "PHK", "location": "Stable Bay",               "lat": -37.7979, "lon": 149.2816, "date": "2024-08-06"},
    {"site": "DIP", "location": "Petrel Point",             "lat": -37.7809, "lon": 149.4323, "date": "2024-08-06"},
    {"site": "LRA", "location": "Little Rame Head",         "lat": -37.6960, "lon": 149.6374, "date": "2024-08-07"},
    {"site": "SDP", "location": "Sandpatch Point reef",     "lat": -37.7096, "lon": 149.6060, "date": "2024-08-07"},
    {"site": "ULV", "location": "Black Jack Rock",          "lat": -41.1361, "lon": 146.1649, "date": "2024-09-24"},
    {"site": "DEV", "location": "Devonport breakwater",     "lat": -41.1670, "lon": 146.3702, "date": "2024-09-24"},
    {"site": "BBB", "location": "Beerbarrel Beach",         "lat": -41.2835, "lon": 148.3599, "date": "2024-10-01"},
    {"site": "GRG", "location": "George Rocks",             "lat": -40.9246, "lon": 148.3276, "date": "2024-10-02"},
    {"site": "COB", "location": "Cobler Rocks",             "lat": -40.8908, "lon": 148.2631, "date": "2024-10-02"},
    {"site": "BCR", "location": "Ben's Creek",              "lat": -40.9409, "lon": 148.3037, "date": "2024-10-02"},
    {"site": "BBB", "location": "Beerbarrel Beach",         "lat": -41.2835, "lon": 148.3599, "date": "2024-10-03"},
    {"site": "SHL", "location": "Shelly Point reef",        "lat": -41.4353, "lon": 148.2785, "date": "2024-10-03"},
    {"site": "BPT", "location": "Bridport – Bentley St",    "lat": -40.9816, "lon": 147.3872, "date": "2024-10-14"},
    {"site": "LSI", "location": "Little Spectacle Island",  "lat": -42.8674, "lon": 147.5995, "date": "2024-10-17"},
    {"site": "IOC", "location": "Isle of Caves",            "lat": -42.9083, "lon": 147.6505, "date": "2024-10-17"},
    {"site": "9PN", "location": "Ninepin Point",            "lat": -43.2839, "lon": 147.1721, "date": "2024-11-04"},
    {"site": "PTI", "location": "Partridge Island",         "lat": -43.3911, "lon": 147.1026, "date": "2024-11-04"},
    {"site": "SIS", "location": "Sisters Bay",              "lat": -43.4260, "lon": 147.0233, "date": "2024-11-04"},
    {"site": "FSH", "location": "Fishers Point",            "lat": -43.5742, "lon": 146.9266, "date": "2024-11-05"},
    {"site": "ACT", "location": "Actaeon Island",           "lat": -43.5258, "lon": 146.9943, "date": "2024-11-05"},
    {"site": "TIN", "location": "Tinderbox – Pierson's Pt", "lat": -43.0482, "lon": 147.3448, "date": "2024-11-05"},
    {"site": "ROB", "location": "Robe – Caravan Park",      "lat": -37.1762, "lon": 139.7361, "date": "2024-11-24"},
    {"site": "ROB", "location": "Robe – Obelisk",           "lat": -37.1762, "lon": 139.7361, "date": "2024-11-24"},
    {"site": "BCH", "location": "Beachport – Penguin Is",   "lat": -37.4839, "lon": 140.0181, "date": "2024-11-24"},
    {"site": "BLK", "location": "Blackfellows Caves Beach", "lat": -37.9506, "lon": 140.4677, "date": "2024-11-25"},
    {"site": "MAC", "location": "Port MacDonnell reef",     "lat": -38.0564, "lon": 140.7022, "date": "2024-11-25"},
]

# Unique sites for climatology (one row per unique site ID)
UNIQUE_SITES = []
seen_ids = set()
for r in SAMPLING_ROWS:
    if r["site"] not in seen_ids:
        UNIQUE_SITES.append({"site": r["site"], "location": r["location"],
                              "lat": r["lat"], "lon": r["lon"]})
        seen_ids.add(r["site"])

# =============================================================================
# CONFIG
# =============================================================================

# Climatology window — long-term means
CLIM_START = "1993-01-01"
CLIM_END   = "2019-12-31"   # BRAN2020 ends here; keeps all sources in sync

# Window around each sampling date for "at time of sampling" extraction
SAMPLE_WINDOW_DAYS = 30     # ±30 days around sampling date

# Bounding box covering all sequenced sites with 1° buffer
LAT_MIN, LAT_MAX = -44.6, -36.2
LON_MIN, LON_MAX = 138.7, 150.6

DEPTH_MIN = 0.49
DEPTH_MAX = 10.0

OUT_DIR = Path("./kelp_env_data")
OUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# COPERNICUS DATASET IDs
# =============================================================================

# [A] SST — OSTIA daily L4 reanalysis, 0.05°
SST_DATASET     = "METOFFICE-GLO-SST-L4-REP-OBS-SST"

# [B] Physics — GLORYS12 monthly, 1/12°
PHY_DATASET     = "cmems_mod_glo_phy_my_0.083deg_P1M-m"

# [C] Biogeochemistry — PISCES monthly, 1/4°
# Product: GLOBAL_MULTIYEAR_BGC_001_029
BGC_DATASET     = "cmems_mod_glo_bgc_my_0.25deg_P1M-m"

# [D] Waves — ERA5-based global wave reanalysis, monthly, 0.5°
# Product: GLOBAL_MULTIYEAR_WAV_001_032
WAV_DATASET     = "cmems_mod_glo_wav_my_0.2deg_P1M-m"

# =============================================================================
# HELPERS
# =============================================================================

def nearest_point(ds, lat, lon):
    lat_dim = "latitude" if "latitude" in ds.dims else "lat"
    lon_dim = "longitude" if "longitude" in ds.dims else "lon"
    return ds.sel({lat_dim: lat, lon_dim: lon}, method="nearest")


def summary_stats(da, prefix):
    """
    Climatological summary stats from a DataArray time series.
    Southern Hemisphere: summer = DJF, winter = JJA.
    """
    vals = da.values.flatten()
    vals = vals[~np.isnan(vals)]
    keys = ["mean", "std", "min", "max", "summer_mean", "winter_mean"]

    if len(vals) == 0:
        return {f"{prefix}_{k}": np.nan for k in keys}

    stats = {
        f"{prefix}_mean": float(np.nanmean(vals)),
        f"{prefix}_std":  float(np.nanstd(vals)),
        f"{prefix}_min":  float(np.nanmin(vals)),
        f"{prefix}_max":  float(np.nanmax(vals)),
    }
    try:
        stats[f"{prefix}_summer_mean"] = float(
            da.sel(time=da["time.season"] == "DJF").mean())
        stats[f"{prefix}_winter_mean"] = float(
            da.sel(time=da["time.season"] == "JJA").mean())
    except Exception:
        stats[f"{prefix}_summer_mean"] = np.nan
        stats[f"{prefix}_winter_mean"] = np.nan
    return stats


def window_mean(da, date_str, days=30):
    """
    Mean of a DataArray within ±days of a sampling date.
    Returns a single float.
    """
    centre = pd.Timestamp(date_str)
    t_min  = centre - timedelta(days=days)
    t_max  = centre + timedelta(days=days)
    try:
        subset = da.sel(time=slice(str(t_min.date()), str(t_max.date())))
        vals = subset.values.flatten()
        vals = vals[~np.isnan(vals)]
        return float(np.nanmean(vals)) if len(vals) > 0 else np.nan
    except Exception:
        return np.nan


def cmems_subset(dataset_id, variables, nc_file, extra_kwargs=None):
    """Download a Copernicus Marine subset, skipping if already cached."""
    if nc_file.exists():
        print(f"  Cached: {nc_file}")
        return
    import copernicusmarine
    kwargs = dict(
        dataset_id=dataset_id,
        variables=variables,
        minimum_latitude=LAT_MIN, maximum_latitude=LAT_MAX,
        minimum_longitude=LON_MIN, maximum_longitude=LON_MAX,
        start_datetime=CLIM_START, end_datetime=CLIM_END,
        output_filename=str(nc_file),
    )
    if extra_kwargs:
        kwargs.update(extra_kwargs)
    copernicusmarine.subset(**kwargs)
    print(f"  Saved: {nc_file}")


# =============================================================================
# [A] OSTIA SST
# =============================================================================

def get_sst():
    print("\n[A] OSTIA SST L4 Reanalysis ...")
    nc = OUT_DIR / "sst_seaus.nc"
    cmems_subset(SST_DATASET, ["analysed_sst"], nc)
    ds = xr.open_dataset(nc)
    if ds["analysed_sst"].attrs.get("units","").lower() in ["kelvin","k"]:
        ds["analysed_sst"] = ds["analysed_sst"] - 273.15
    return ds


# =============================================================================
# [B] GLORYS12 PHYSICS — surface + bottom temperature
# =============================================================================

def get_physics():
    print("\n[B] GLORYS12 Physics ...")
    # Surface layer
    nc_surf = OUT_DIR / "glorys_surface.nc"
    cmems_subset(PHY_DATASET,
                 ["thetao", "so", "mlotst", "uo", "vo"],
                 nc_surf,
                 {"minimum_depth": DEPTH_MIN, "maximum_depth": DEPTH_MAX})
    ds_surf = xr.open_dataset(nc_surf)
    if "depth" in ds_surf.dims:
        ds_surf = ds_surf.isel(depth=0, drop=True)

    # Bottom temperature (deepest available layer ~200m for shelf sites)
    nc_bot = OUT_DIR / "glorys_bottom.nc"
    cmems_subset(PHY_DATASET,
                 ["thetao"],
                 nc_bot,
                 {"minimum_depth": 100.0, "maximum_depth": 200.0})
    ds_bot = xr.open_dataset(nc_bot)
    if "depth" in ds_bot.dims:
        ds_bot = ds_bot.isel(depth=-1, drop=True)
    ds_bot = ds_bot.rename({"thetao": "bottom_temp"})

    return ds_surf, ds_bot


# =============================================================================
# [C] GLORYS12 BIOGEOCHEMISTRY — nutrients, chlorophyll, O2
# =============================================================================

def get_bgc():
    print("\n[C] GLORYS12 Biogeochemistry (PISCES) ...")
    nc = OUT_DIR / "glorys_bgc.nc"
    cmems_subset(BGC_DATASET,
                 ["no3", "po4", "chl", "o2", "phyc"],
                 nc,
                 {"minimum_depth": DEPTH_MIN, "maximum_depth": DEPTH_MAX})
    ds = xr.open_dataset(nc)
    if "depth" in ds.dims:
        ds = ds.isel(depth=0, drop=True)
    return ds


# =============================================================================
# [D] WAVE REANALYSIS
# =============================================================================

def get_waves():
    print("\n[D] CMEMS Wave Reanalysis ...")
    nc = OUT_DIR / "waves_seaus.nc"
    cmems_subset(WAV_DATASET,
                 ["VHM0", "VTPK", "VMDR"],   # sig wave height, peak period, mean direction
                 nc)
    return xr.open_dataset(nc)


# =============================================================================
# [E] BRAN2020 — via NCI OPeNDAP (no account needed)
# =============================================================================

BRAN_BASE = "https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily"
BRAN_VARS = {
    "ocean_temp": ("temp", "bran_temp"),
    "ocean_salt": ("salt", "bran_salinity"),
    "ocean_mld":  ("mld",  "bran_mld"),
    "ocean_u":    ("u",    "bran_u"),
    "ocean_v":    ("v",    "bran_v"),
}


def get_bran(sites):
    """Stream BRAN2020 monthly and accumulate values per site."""
    print("\n[E] BRAN2020 via NCI OPeNDAP ...")
    months = pd.date_range(start=CLIM_START, end="2019-12-01", freq="MS")
    accum = {s["site"]: {p: [] for _, p in BRAN_VARS.values()} for s in sites}

    for dt in tqdm(months, desc="  BRAN chunks"):
        y, m = dt.year, dt.month
        for nc_var, (field, prefix) in BRAN_VARS.items():
            url = f"{BRAN_BASE}/{nc_var}_{y}_{m:02d}.nc"
            try:
                ds = xr.open_dataset(url, engine="netcdf4")
                rn = {old: new for old, new in
                      [("xt_ocean","longitude"),("yt_ocean","latitude"),
                       ("st_ocean","depth"),("Time","time")] if old in ds.dims}
                if rn:
                    ds = ds.rename(rn)
                ds = ds.sel(latitude=slice(LAT_MIN, LAT_MAX),
                            longitude=slice(LON_MIN, LON_MAX))
                if "depth" in ds.dims:
                    ds = ds.isel(depth=0, drop=True)
                for site in sites:
                    pt = nearest_point(ds, site["lat"], site["lon"])
                    if field in pt:
                        v = pt[field].values.flatten()
                        accum[site["site"]][prefix].extend(v[~np.isnan(v)].tolist())
                ds.close()
            except Exception:
                pass
    return accum


# =============================================================================
# BUILD CLIMATOLOGY CSV — one row per unique site
# =============================================================================

def build_climatology(ds_sst, ds_surf, ds_bot, ds_bgc, ds_wav, bran_accum):
    print("\n Building climatology table ...")
    rows = []
    for site in tqdm(UNIQUE_SITES, desc="  Sites"):
        row = {"site": site["site"], "location": site["location"],
               "lat": site["lat"], "lon": site["lon"]}

        # [A] SST
        pt = nearest_point(ds_sst, site["lat"], site["lon"])
        row.update(summary_stats(pt["analysed_sst"], "sst"))

        # [B] Physics surface
        pt = nearest_point(ds_surf, site["lat"], site["lon"])
        for var, pfx in [("thetao","glorys_temp"),("so","glorys_salinity"),
                          ("mlotst","glorys_mld"),("uo","glorys_u"),("vo","glorys_v")]:
            if var in pt:
                row.update(summary_stats(pt[var], pfx))

        # [B] Bottom temperature
        pt = nearest_point(ds_bot, site["lat"], site["lon"])
        if "bottom_temp" in pt:
            row.update(summary_stats(pt["bottom_temp"], "glorys_bottom_temp"))

        # [C] BGC
        pt = nearest_point(ds_bgc, site["lat"], site["lon"])
        for var, pfx in [("no3","nitrate"),("po4","phosphate"),
                          ("chl","chlorophyll"),("o2","oxygen"),("phyc","phytoplankton")]:
            if var in pt:
                row.update(summary_stats(pt[var], pfx))

        # [D] Waves
        pt = nearest_point(ds_wav, site["lat"], site["lon"])
        for var, pfx in [("VHM0","wave_height"),("VTPK","wave_period"),
                          ("VMDR","wave_direction")]:
            if var in pt:
                row.update(summary_stats(pt[var], pfx))

        # [E] BRAN
        for _, prefix in BRAN_VARS.values():
            arr = np.array(bran_accum[site["site"]][prefix])
            if len(arr) > 0:
                row[f"{prefix}_mean"] = float(np.nanmean(arr))
                row[f"{prefix}_std"]  = float(np.nanstd(arr))
                row[f"{prefix}_min"]  = float(np.nanmin(arr))
                row[f"{prefix}_max"]  = float(np.nanmax(arr))
            else:
                for k in ["mean","std","min","max"]:
                    row[f"{prefix}_{k}"] = np.nan

        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# BUILD AT-SAMPLING CSV — one row per sampling visit
# =============================================================================

def build_at_sampling(ds_sst, ds_surf, ds_bot, ds_bgc, ds_wav):
    print("\n Building at-sampling-date table ...")
    rows = []
    for s in tqdm(SAMPLING_ROWS, desc="  Rows"):
        row = {"site": s["site"], "location": s["location"],
               "sampling_date": s["date"], "lat": s["lat"], "lon": s["lon"]}

        # [A] SST
        pt = nearest_point(ds_sst, s["lat"], s["lon"])
        row["sst_sample_mean"] = window_mean(pt["analysed_sst"], s["date"])

        # [B] Physics surface
        pt = nearest_point(ds_surf, s["lat"], s["lon"])
        for var, pfx in [("thetao","glorys_temp"),("so","glorys_salinity"),
                          ("mlotst","glorys_mld"),("uo","glorys_u"),("vo","glorys_v")]:
            if var in pt:
                row[f"{pfx}_sample_mean"] = window_mean(pt[var], s["date"])

        # [B] Bottom temp
        pt = nearest_point(ds_bot, s["lat"], s["lon"])
        if "bottom_temp" in pt:
            row["glorys_bottom_temp_sample_mean"] = window_mean(
                pt["bottom_temp"], s["date"])

        # [C] BGC
        pt = nearest_point(ds_bgc, s["lat"], s["lon"])
        for var, pfx in [("no3","nitrate"),("po4","phosphate"),
                          ("chl","chlorophyll"),("o2","oxygen"),("phyc","phytoplankton")]:
            if var in pt:
                row[f"{pfx}_sample_mean"] = window_mean(pt[var], s["date"])

        # [D] Waves
        pt = nearest_point(ds_wav, s["lat"], s["lon"])
        for var, pfx in [("VHM0","wave_height"),("VTPK","wave_period"),
                          ("VMDR","wave_direction")]:
            if var in pt:
                row[f"{pfx}_sample_mean"] = window_mean(pt[var], s["date"])

        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 65)
    print(" Giant Kelp GEA – Full Environmental Data Collation")
    print(" SE Australia: SA, VIC, TAS  |  Sequenced sites only")
    print("=" * 65)
    print(f"  Unique sites:    {len(UNIQUE_SITES)}")
    print(f"  Sampling rows:   {len(SAMPLING_ROWS)}")
    print(f"  Climatology:     {CLIM_START} → {CLIM_END}")
    print(f"  Sample window:   ±{SAMPLE_WINDOW_DAYS} days around each visit")
    print("=" * 65)
    print()
    print("  Variables collected:")
    print("    SST (OSTIA)              – sea surface temperature")
    print("    Physics (GLORYS12)       – temp, salinity, MLD, u/v currents,")
    print("                               bottom temperature")
    print("    Biogeochemistry (PISCES) – nitrate, phosphate, chlorophyll,")
    print("                               dissolved oxygen, phytoplankton")
    print("    Waves (CMEMS)            – significant wave height, peak period,")
    print("                               mean direction")
    print("    BRAN2020 (NCI OPeNDAP)  – temp, salinity, MLD, u/v currents")
    print()

    # --- Download all datasets ---
    ds_sst          = get_sst()
    ds_surf, ds_bot = get_physics()
    ds_bgc          = get_bgc()
    ds_wav          = get_waves()
    bran_accum      = get_bran(UNIQUE_SITES)

    # --- Build output tables ---
    df_clim   = build_climatology(ds_sst, ds_surf, ds_bot, ds_bgc, ds_wav, bran_accum)
    df_sample = build_at_sampling(ds_sst, ds_surf, ds_bot, ds_bgc, ds_wav)

    # --- Save ---
    clim_csv   = OUT_DIR / "env_climatology.csv"
    sample_csv = OUT_DIR / "env_at_sampling.csv"
    df_clim.to_csv(clim_csv,   index=False)
    df_sample.to_csv(sample_csv, index=False)

    print(f"\n Done!")
    print(f"  Climatology CSV:  {clim_csv}")
    print(f"    Shape: {df_clim.shape[0]} sites × {df_clim.shape[1]} variables")
    print(f"  At-sampling CSV:  {sample_csv}")
    print(f"    Shape: {df_sample.shape[0]} rows × {df_sample.shape[1]} variables")
    print()
    print("  Next: use env_climatology.csv as your X matrix in RDA/GEA")


if _name_ == "_main_":
    main()