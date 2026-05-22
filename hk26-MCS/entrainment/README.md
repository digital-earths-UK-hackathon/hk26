# MCS Entrainment Analysis — WAM Region

Analysis of convective entrainment proxies for Mesoscale Convective Systems (MCS) over the West African Monsoon (WAM) region, using UM global model output from the 2026 Digital Earths UK Hackathon.

Four global UM simulations are supported: N2560 RAL3.3, N2560 CoMA9, N1280 GAL9, N1280 CoMA9.

For technical details see [docs/technical_overview.md](docs/technical_overview.md).

---

## Quick start

All commands are run from `hk26-MCS/entrainment/`. Activate the hackathon environment first:

```bash
conda activate hk26_env
```

### 1. Compute entrainment variables (SLURM array)

```bash
python submit.py --model um_glm_n2560_RAL3p3_tuned_hk26
```

This initialises the zarr store if needed, then submits only the chunks whose done files are absent. Re-run at any time to pick up failed or missing chunks. Use `--dry-run` to preview without submitting.

### 2. Link entrainment with MCS tracks

```bash
python calc_mcs_entrainment.py --model um_glm_n2560_RAL3p3_tuned_hk26
```

Optionally filter to land or ocean MCS:

```bash
python calc_mcs_entrainment.py --model um_glm_n2560_RAL3p3_tuned_hk26 --surface land
python calc_mcs_entrainment.py --model um_glm_n2560_RAL3p3_tuned_hk26 --surface ocean
```

Output: `data/<model>/mcs_entrainment_wam[_land|_ocean].nc`

### 3. Plot MCS entrainment statistics

```bash
python plot_mcs_entrainment.py --model um_glm_n2560_RAL3p3_tuned_hk26
```

Produces three figures in `figs/<model>/`:
- `*.mcs_lifecycle_entrainment.png` — per-variable lifecycle composite (JJA/DJF)
- `*.mcs_dc_entrainment.png` — MCS entrainment diurnal cycle
- `*.mcs_distributions.png` — variable distributions (JJA/DJF)

### 4. Plot background diurnal cycle

```bash
python plot_diurnal_cycle.py --model um_glm_n2560_RAL3p3_tuned_hk26
```

Produces map, diurnal cycle, and entrainment proxy figures over a Burkina Faso sub-box.

---

## Directory layout

```
entrainment/
├── models.py                 # Central config: models, regions, URLs, paths
├── submit.py                 # Smart SLURM submitter
├── calc_entrainment.py       # Per-cell CAPE/LNB/Tb/shear/precip computation
├── calc_mcs_entrainment.py   # Link entrainment with PyFLEXTRKR MCS tracks
├── plot_mcs_entrainment.py   # MCS lifecycle / diurnal cycle / distribution plots
├── plot_diurnal_cycle.py     # Background diurnal cycle plots
├── donefiles/<model>/        # Sentinel files tracking completed SLURM chunks
├── slurm/
│   ├── scripts/              # Generated SLURM batch scripts
│   ├── tasks/                # JSON task lists (one per submission)
│   └── output/               # SLURM stdout/stderr
├── data/<model>/             # Zarr stores and output NetCDF files (on scratch)
└── figs/<model>/             # Output figures
```

---

## Models

| Key | Resolution | Config |
|-----|-----------|--------|
| `um_glm_n2560_RAL3p3_tuned_hk26` | N2560 (~13 km) | RAL3.3 |
| `um_glm_n2560_CoMA9_hk26` | N2560 (~13 km) | CoMA9 |
| `um_glm_n1280_GAL9_v2_hk26` | N1280 (~25 km) | GAL9 |
| `um_glm_n1280_CoMA9_hk26` | N1280 (~25 km) | CoMA9 |

All use HEALPix zoom=9, covering 2020-02-01 to 2021-03-01.
