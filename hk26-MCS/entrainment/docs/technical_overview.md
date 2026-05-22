# Technical Overview

## Data sources

### Model data
UM global simulations accessed via an intake catalog at:
`https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml`

Data is on a HEALPix grid at zoom=9 (~13 km resolution, 3,145,728 global cells). Two temporal resolutions are used:
- **PT3H** — temperature (`ta`), relative humidity (`hur`), vertical velocity (`wa`), zonal wind (`ua`) on pressure levels; used for CAPE/LNB/shear/RH700
- **PT1H** — OLR (`rlut`), precipitation (`pr`), precipitable water (`prw`); sampled at 3-hourly times to match PT3H

### MCS tracking data (PyFLEXTRKR)
Hosted on S3 at `hackathon-o.s3-ext.jc.rl.ac.uk`. Two files per model:
- **Pixel mask** (`mcs_mask_hp9_*.zarr`) — hourly, global HEALPix; integer values where each pixel is labelled with a track number (0 = background)
- **Track statistics** (`mcs_tracks_final_*.nc`) — NetCDF with dims `(tracks=43189, times=650)`; contains centroid position, duration, timing, precipitation feature properties etc.

The mask pixel value N corresponds to track index N−1 in the statistics file.

---

## Pipeline steps

### Step 1: `calc_entrainment.py` — per-cell variables

Computes 11 variables for every WAM cell at every 3-hourly timestep, writing to a zarr store.

**WAM region mask** (`lon > 340° or lon < 20°`, `2° < lat < 15°N`) selects ~39,000 of the ~3.1M global HEALPix cells.

**CAPE/CIN/LNB** are computed per vertical profile using MetPy:
1. Convert relative humidity to dewpoint (`dewpoint_from_relative_humidity`)
2. Compute surface parcel profile (`parcel_profile`)
3. `cape_cin` returns CAPE and CIN in J/kg; CAPE is clipped to ≥ 0
4. `el` (equilibrium level) returns LNB pressure (hPa) and temperature T_LNB (K)
5. Profiles with fewer than 15 valid dewpoint levels return NaN

This is embarrassingly parallel over cells. A `multiprocessing.Pool` of 10 workers (one per SLURM CPU) processes each timestep; 10 timesteps per SLURM array task (CHUNK_SIZE=10).

**Brightness temperature** from OLR via the Yang & Slingo (2001) / Minnis & Harrison (1984) formula (sourced from PyFLEXTRKR `ftfunctions.py`):

```
Tf = (OLR / σ)^0.25
Tb = (−a + √(a² + 4·b·Tf)) / (2·b)
```

where a = 1.228, b = −1.106×10⁻³ K⁻¹, σ = 5.67×10⁻⁸ W m⁻² K⁻⁴.

**Other variables:**
- `w_eff = w(500 hPa) / √CAPE` — normalised updraft proxy (dimensionless); NaN where CAPE = 0
- `tb_diff = Tb − T_LNB` — cloud-top overshooting proxy; negative = cloud top above LNB
- `shear = u(600 hPa) − u(850 hPa)` — low-to-mid level zonal wind shear (650 hPa not available)
- `hur700` — relative humidity at 700 hPa
- `pr`, `prw` — precipitation flux and precipitable water from PT1H, sampled at 3H times

**Zarr store layout:** `data/<model>/entrainment_<region>.zarr`, shape `(3249, ~39130)`, chunked `(10, 39130)` — one chunk per SLURM task.

**SLURM submission** via `submit.py`:
1. Checks `donefiles/<model>/init_<region>.done`; runs `--init` if absent
2. Scans `donefiles/<model>/chunk_NNN.done` to find pending chunks
3. Writes a JSON task list (`slurm/tasks/`) and a SLURM array script (`slurm/scripts/`)
4. `sbatch` submits only the pending chunks — safe to re-run at any time

---

### Step 2: `calc_mcs_entrainment.py` — per-track aggregation

Links the per-cell entrainment zarr with the PyFLEXTRKR track mask to produce per-track, per-lifecycle-step statistics with dims `(tracks, times_3h=217)`.

**Time alignment**: the 3-hourly entrainment times are matched to the hourly mask using a `pd.Timestamp` dict lookup. 3152 timesteps overlap (2020-02-01 to 2021-02-28).

**Spatial mapping**: WAM cell HEALPix indices are located in the global mask array via `np.searchsorted`, giving positional indices for fast slicing.

**Track pre-filtering**: tracks whose centroid never enters the buffered WAM bounding box are dropped before any computation (2871 of 43189 pass). An optional `--surface` filter keeps only land (`pf_landfrac > 0.8`) or ocean (`< 0.2`) MCS.

**Vectorised aggregation** (per timestep, all tracks simultaneously):
1. Load mask at the matching hour → select WAM cells → integer array `mask_wam` (0 = background)
2. For each variable, `np.bincount` over `mask_wam` with value weights gives per-track sums and sum-of-squares in O(n_cells) — no Python loop over tracks
3. Mean and std computed as `mean = sum/count`, `std = √(E[x²] − E[x]²)`
4. Each active track's lifecycle step index `li = (global_time − start_basetime) / 3h` determines the output column

**Output** `mcs_entrainment_wam.nc`:

| Variable | Dims | Description |
|----------|------|-------------|
| `{var}_mean`, `{var}_std` | (tracks, times_3h) | Per-track mean/std for each of the 11 entrainment variables |
| `n_wam_cells` | (tracks, times_3h) | Number of MCS mask cells in WAM at each step |
| `base_time` | (tracks, times_3h) | UTC timestamp of each lifecycle step (NaT if inactive) |
| `track_duration_3h` | (tracks,) | Track duration in 3-hourly steps |
| `start_basetime`, `end_basetime` | (tracks,) | Track start/end times |
| `meanlat`, `meanlon` | (tracks, times_3h) | MCS centroid, downsampled from hourly to 3-hourly |

---

### Step 3: `plot_mcs_entrainment.py` — MCS statistics plots

Three figures, all split by JJA/DJF using track start time. Track counts per season are shown in legends.

**Lifecycle composite** — per-track interpolation approach:
- Each track's valid (lifecycle_fraction, value) pairs are linearly interpolated onto 20 fixed bin midpoints spanning [0, 1]
- `np.nan` outside each track's observed range (short tracks don't reach the extremes)
- Median ± IQR taken across tracks — equal weight per MCS regardless of duration

**MCS diurnal cycle** — groups `base_time` hour across all active (track, step) pairs; mean ± std per UTC hour.

**Distributions** — histograms over all active (track, step) entries, binned to p1–p99 to exclude outliers; 240 K reference line on Tb, zero-line annotations on tb_diff and w_eff.

For `w_eff`, the y-axis is clamped to the min/max of the median lines (not the IQR), since near-zero CAPE can produce extreme outliers in individual cells.

---

### Step 4: `plot_diurnal_cycle.py` — background diurnal cycle

Loads the full entrainment zarr (all WAM cells, not just MCS cells) and selects a sub-box over Burkina Faso (9–15°N, 5°W–5°E). Three figures:
- **Map** — WAM region and Burkina Faso box overlaid on coastlines
- **Diurnal cycle** — CAPE and LNB pressure, JJA/DJF, mean ± 1 std
- **Entrainment proxies** — w_eff and Tb_diff diurnal cycles

---

## Key physical variables

| Variable | Interpretation |
|----------|---------------|
| CAPE | Convective instability; higher = more energetic convection |
| CIN | Convective inhibition; more negative = harder to initiate convection |
| LNB | Level of neutral buoyancy pressure; lower pressure = higher cloud tops |
| T_LNB | Temperature at LNB; proxy for cloud-top temperature |
| Tb | IR brightness temperature from OLR; lower = colder/higher cloud tops |
| tb_diff | Tb − T_LNB; negative = cloud top overshoots LNB (deep convection) |
| w_eff | w(500 hPa)/√CAPE; normalised updraft strength, dimensionless |
| shear | u(600) − u(850 hPa); low-level zonal wind shear |
| hur700 | Mid-level relative humidity; modulates entrainment drying |
| pr | Surface precipitation rate |
| prw | Total column water vapour (precipitable water) |

---

## Configuration (`models.py`)

All model keys, analysis region bounds, S3 URLs, and filesystem paths are defined in `models.py`. To add a new model, add an entry to `MODELS`. To add a new region, add an entry to `REGIONS` with lon/lat bounds and buffered pre-filter bounds.

Data is written to scratch (`/work/scratch-nopw2/mmuetz/hk26/hk26-MCS/entrainment/data/`) via `models.data_dir()`. Done files and SLURM artefacts stay under the working directory.
