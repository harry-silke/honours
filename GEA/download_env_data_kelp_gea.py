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
  [B] GLORYS12 physics        – temp, salinity, MLD, currents (Copernicus)
  [C] GLORYS12 biogeochemistry– nitrate, phosphate, chlorophyll, O2 (Copernicus)
  [D] CMEMS wave reanalysis   – significant wave height, period, direction (Copernicus)
  [E] BRAN2020                – temp, salinity, MLD (NCI OPeNDAP, free)

CHANGES IN THIS REVISION
------------------------
  * Wave direction (VMDR) is a CIRCULAR variable. It is now summarised with a
    circular mean / circular std (vector averaging), and min/max are dropped
    because they are undefined on a compass bearing. Previously the arithmetic
    mean produced physically wrong bearings (e.g. mean of 350 and 10 -> 180).
  * BRAN now uses nearest_wet_point (land-mask aware) like the CMEMS sources,
    so coastal sites no longer silently return all-NaN when the nearest BRAN
    cell is land.
  * Added glorys_speed = sqrt(uo^2 + vo^2), a scalar current speed. This is a
    cleaner RDA predictor than the signed u/v components (which are retained).

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

# [D] Waves — global wave reanalysis, monthly, 0.2°
# Product: GLOBAL_MULTIYEAR_WAV_001_032
WAV_DATASET     = "cmems_mod_glo_wav_my_0.2deg_P1M-m"

# =============================================================================
# HELPERS
# =============================================================================

def nearest_point(ds, lat, lon):
    lat_dim = "latitude" if "latitude" in ds.dims else "lat"
    lon_dim = "longitude" if "longitude" in ds.dims else "lon"
    return ds.sel({lat_dim: lat, lon_dim: lon}, method="nearest")


def nearest_wet_point(ds, lat, lon, ref_var=None, window=8):
    """
    Select the nearest OCEAN grid cell to (lat, lon).

    Coastal kelp sites frequently fall on (or beside) a masked land cell,
    especially on the coarser BGC (0.25 deg) and wave (0.2 deg) grids — the
    single nearest cell is then NaN and the site drops out. This searches a
    window of cells around the target, keeps those where `ref_var` has any
    valid (non-NaN) data, and returns the NEAREST such wet cell. All variables
    in a dataset share the same land mask, so one reference variable suffices.

    Falls back to the plain nearest cell if no wet cell is found in the window
    (the result will be NaN, flagging a genuinely land-locked lookup).
    """
    lat_dim = "latitude" if "latitude" in ds.dims else "lat"
    lon_dim = "longitude" if "longitude" in ds.dims else "lon"

    if ref_var is None:
        ref_var = list(ds.data_vars)[0]

    lats = ds[lat_dim].values
    lons = ds[lon_dim].values
    ilat = int(np.abs(lats - lat).argmin())
    ilon = int(np.abs(lons - lon).argmin())

    lat_sl = slice(max(0, ilat - window), ilat + window + 1)
    lon_sl = slice(max(0, ilon - window), ilon + window + 1)
    sub = ds.isel({lat_dim: lat_sl, lon_dim: lon_sl})

    da = sub[ref_var]
    reduce_dims = [d for d in da.dims if d not in (lat_dim, lon_dim)]
    valid = da.notnull().any(dim=reduce_dims) if reduce_dims else da.notnull()
    valid = valid.transpose(lat_dim, lon_dim)
    valid_np = valid.values

    if not valid_np.any():
        return ds.sel({lat_dim: lat, lon_dim: lon}, method="nearest")

    sub_lats = sub[lat_dim].values
    sub_lons = sub[lon_dim].values
    LON, LAT = np.meshgrid(sub_lons, sub_lats)   # shape (nlat, nlon)
    # approximate great-circle distance; cos-correct longitude for latitude
    dist = np.sqrt((LAT - lat) ** 2 +
                   ((LON - lon) * np.cos(np.deg2rad(lat))) ** 2)
    dist = np.where(valid_np, dist, np.inf)
    j, k = np.unravel_index(int(np.argmin(dist)), dist.shape)
    return sub.isel({lat_dim: j, lon_dim: k})


# --- circular statistics for directional variables (e.g. wave direction) ----

def circular_mean_deg(deg):
    """Vector (circular) mean of bearings in degrees. Returns 0–360 or NaN."""
    deg = np.asarray(deg, dtype=float)
    deg = deg[~np.isnan(deg)]
    if deg.size == 0:
        return np.nan
    r = np.deg2rad(deg)
    ang = np.rad2deg(np.arctan2(np.mean(np.sin(r)), np.mean(np.cos(r))))
    return float(ang % 360.0)


def circular_std_deg(deg):
    """
    Circular standard deviation of bearings in degrees.
    sd = sqrt(-2 ln R), where R is the mean resultant length. Small sd => tight
    directional spread; large sd => near-uniform / no preferred direction.
    """
    deg = np.asarray(deg, dtype=float)
    deg = deg[~np.isnan(deg)]
    if deg.size == 0:
        return np.nan
    r = np.deg2rad(deg)
    R = np.hypot(np.mean(np.sin(r)), np.mean(np.cos(r)))
    R = min(max(R, 1e-12), 1.0)
    return float(np.rad2deg(np.sqrt(-2.0 * np.log(R))))


def summary_stats(da, prefix):
    """
    Climatological summary stats from a linear DataArray time series.
    Southern Hemisphere: summer = DJF, winter = JJA.
    Do NOT use this for directional variables — see circular_summary.
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


def circular_summary(da, prefix):
    """
    Circular summary for a directional DataArray (degrees).
    Emits mean/std and seasonal circular means. No min/max: they are undefined
    on a compass bearing (0 and 359 are adjacent, not extremes).
    """
    out = {
        f"{prefix}_mean": circular_mean_deg(da.values.flatten()),
        f"{prefix}_std":  circular_std_deg(da.values.flatten()),
    }
    try:
        out[f"{prefix}_summer_mean"] = circular_mean_deg(
            da.sel(time=da["time.season"] == "DJF").values.flatten())
        out[f"{prefix}_winter_mean"] = circular_mean_deg(
            da.sel(time=da["time.season"] == "JJA").values.flatten())
    except Exception:
        out[f"{prefix}_summer_mean"] = np.nan
        out[f"{prefix}_winter_mean"] = np.nan
    return out


def window_mean(da, date_str, days=30):
    """Mean of a linear DataArray within ±days of a sampling date."""
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


def window_circular_mean(da, date_str, days=30):
    """Circular mean of a directional DataArray within ±days of a date."""
    centre = pd.Timestamp(date_str)
    t_min  = centre - timedelta(days=days)
    t_max  = centre + timedelta(days=days)
    try:
        subset = da.sel(time=slice(str(t_min.date()), str(t_max.date())))
        return circular_mean_deg(subset.values.flatten())
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
# [B] GLORYS12 PHYSICS — surface
# =============================================================================
# NOTE: "bottom temperature" was dropped. The only depth layers available in
# the requested 100–200 m range (110/131/156/186 m) sit below the seafloor at
# these shallow coastal kelp sites, so every site returned NaN. The variable is
# ill-defined here; use surface/bulk temperature instead.

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

    return ds_surf


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
# NOTE: BRAN u/v (currents) were dropped. Those files live on the staggered
# velocity grid (xu_ocean/yu_ocean) rather than the tracer grid handled below,
# so every fetch failed and the columns were always empty. Surface currents are
# already covered by GLORYS uo/vo, so BRAN currents are redundant here.
BRAN_VARS = {
    "ocean_temp": ("temp", "bran_temp"),
    "ocean_salt": ("salt", "bran_salinity"),
    "ocean_mld":  ("mld",  "bran_mld"),
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
                    # wet-point aware, matching the CMEMS sources; ref_var is
                    # pinned to the field so we don't pick up a bounds variable.
                    pt = nearest_wet_point(ds, site["lat"], site["lon"],
                                           ref_var=field)
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

def build_climatology(ds_sst, ds_surf, ds_bgc, ds_wav, bran_accum):
    print("\n Building climatology table ...")
    rows = []
    for site in tqdm(UNIQUE_SITES, desc="  Sites"):
        row = {"site": site["site"], "location": site["location"],
               "lat": site["lat"], "lon": site["lon"]}

        # [A] SST
        pt = nearest_wet_point(ds_sst, site["lat"], site["lon"])
        row.update(summary_stats(pt["analysed_sst"], "sst"))

        # [B] Physics surface
        pt = nearest_wet_point(ds_surf, site["lat"], site["lon"])
        for var, pfx in [("thetao","glorys_temp"),("so","glorys_salinity"),
                          ("mlotst","glorys_mld"),("uo","glorys_u"),("vo","glorys_v")]:
            if var in pt:
                row.update(summary_stats(pt[var], pfx))
        # scalar current speed — cleaner RDA predictor than signed u/v
        if "uo" in pt and "vo" in pt:
            speed = np.sqrt(pt["uo"] ** 2 + pt["vo"] ** 2)
            row.update(summary_stats(speed, "glorys_speed"))

        # [C] BGC
        pt = nearest_wet_point(ds_bgc, site["lat"], site["lon"])
        for var, pfx in [("no3","nitrate"),("po4","phosphate"),
                          ("chl","chlorophyll"),("o2","oxygen"),("phyc","phytoplankton")]:
            if var in pt:
                row.update(summary_stats(pt[var], pfx))

        # [D] Waves — height/period are linear, direction is circular
        pt = nearest_wet_point(ds_wav, site["lat"], site["lon"])
        for var, pfx in [("VHM0","wave_height"),("VTPK","wave_period")]:
            if var in pt:
                row.update(summary_stats(pt[var], pfx))
        if "VMDR" in pt:
            row.update(circular_summary(pt["VMDR"], "wave_direction"))

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

def build_at_sampling(ds_sst, ds_surf, ds_bgc, ds_wav):
    print("\n Building at-sampling-date table ...")
    rows = []
    for s in tqdm(SAMPLING_ROWS, desc="  Rows"):
        row = {"site": s["site"], "location": s["location"],
               "sampling_date": s["date"], "lat": s["lat"], "lon": s["lon"]}

        # [A] SST
        pt = nearest_wet_point(ds_sst, s["lat"], s["lon"])
        row["sst_sample_mean"] = window_mean(pt["analysed_sst"], s["date"])

        # [B] Physics surface
        pt = nearest_wet_point(ds_surf, s["lat"], s["lon"])
        for var, pfx in [("thetao","glorys_temp"),("so","glorys_salinity"),
                          ("mlotst","glorys_mld"),("uo","glorys_u"),("vo","glorys_v")]:
            if var in pt:
                row[f"{pfx}_sample_mean"] = window_mean(pt[var], s["date"])
        if "uo" in pt and "vo" in pt:
            speed = np.sqrt(pt["uo"] ** 2 + pt["vo"] ** 2)
            row["glorys_speed_sample_mean"] = window_mean(speed, s["date"])

        # [C] BGC
        pt = nearest_wet_point(ds_bgc, s["lat"], s["lon"])
        for var, pfx in [("no3","nitrate"),("po4","phosphate"),
                          ("chl","chlorophyll"),("o2","oxygen"),("phyc","phytoplankton")]:
            if var in pt:
                row[f"{pfx}_sample_mean"] = window_mean(pt[var], s["date"])

        # [D] Waves
        pt = nearest_wet_point(ds_wav, s["lat"], s["lon"])
        for var, pfx in [("VHM0","wave_height"),("VTPK","wave_period")]:
            if var in pt:
                row[f"{pfx}_sample_mean"] = window_mean(pt[var], s["date"])
        if "VMDR" in pt:
            row["wave_direction_sample_mean"] = window_circular_mean(pt["VMDR"], s["date"])

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
    print("    Physics (GLORYS12)       – temp, salinity, MLD, u/v + speed")
    print("    Biogeochemistry (PISCES) – nitrate, phosphate, chlorophyll,")
    print("                               dissolved oxygen, phytoplankton")
    print("    Waves (CMEMS)            – sig wave height, peak period,")
    print("                               mean direction (circular)")
    print("    BRAN2020 (NCI OPeNDAP)  – temp, salinity, MLD")
    print()

    # --- Download all datasets ---
    ds_sst          = get_sst()
    ds_surf         = get_physics()
    ds_bgc          = get_bgc()
    ds_wav          = get_waves()
    bran_accum      = get_bran(UNIQUE_SITES)

    # --- Build output tables ---
    df_clim   = build_climatology(ds_sst, ds_surf, ds_bgc, ds_wav, bran_accum)
    df_sample = build_at_sampling(ds_sst, ds_surf, ds_bgc, ds_wav)

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
    print("  Next: prune to a low-collinearity subset before RDA (see notes).")


if __name__ == "__main__":
    main()