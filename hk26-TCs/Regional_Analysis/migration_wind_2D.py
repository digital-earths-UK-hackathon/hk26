"""
Program: Analysis of North Atlantic Cyclone Activity and
         Sensible Heat Transport

Zexi Sun May 2026
         
Description:
    This script connects to the Digital Earths Cloud Catalog and
    analyzes a single climate model to investigate seasonal
    relationships between cyclone activity and atmospheric
    sensible heat transport over the North Atlantic region.

Methodology:
    1. Connect to the Digital Earths Cloud Catalog and load the
       specified ensemble members.

    2. Reconstruct the global HEALPix geometry using the first
       ensemble member and extract the North Atlantic region:
           Longitude: -60° to 20°
           Latitude: 20°N to 80°N

    3. Divide the selected region into 1° latitude bins for
       zonal averaging.
    
    4. Data Analysis

        A. Extract seasonal datasets:
            - Winter: January 2020

        B. Compute mean-state diagnostics:
            - Sea level pressure variance (PSL variance) as a
              proxy for storm-track intensity
            - Mean sensible heat transport proxy:
                    vT = vas × tas

        C. Compute transient eddy diagnostics:
            - Resample fields to daily means
            - Remove seasonal means to obtain anomalies
            - Calculate:
                    psl'²
                    v'T'
            - Compute Pearson correlation between zonally averaged
              transient cyclone activity and transient heat transport
              within each latitude band

"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import intake
import healpy as hp
import warnings

warnings.filterwarnings('ignore')

print("Accessing to Digital Earths Cloud Catalog...")
url = "https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml"
cat = intake.open_catalog(url)['online']
dataset_key = 'um_glm_n1280_CoMA9_hk26'
ds = cat[dataset_key].to_dask()

# =====================================================================
# 1. Geometry & Regional Masking
# =====================================================================
print("Reconstructing geometry & regional mask...")
n_pixels = ds.sizes['healpix_index']
hp_indices = np.arange(n_pixels)
nside = hp.npix2nside(n_pixels)
lon_vals, lat_vals = hp.pix2ang(nside, hp_indices, nest=True, lonlat=True)
lon_vals = np.where(lon_vals > 180, lon_vals - 360, lon_vals)

# Large North Atlantic & European Sector
lon_min, lon_max = -80, 40
lat_min, lat_max = 30, 80 
region_mask = (lat_vals >= lat_min) & (lat_vals <= lat_max) & (lon_vals >= lon_min) & (lon_vals <= lon_max)
valid_indices = np.where(region_mask)[0]

lat_clean = lat_vals[valid_indices]
lon_clean = lon_vals[valid_indices]

# =====================================================================
# 2. Daily Physics & Vectorized Correlation (Xarray/Dask)
# =====================================================================
print("Slicing Winter data and resampling to daily means...")
# Grab PSL (Sea Level Pressure), UAS (Zonal wind), VAS (Meridional wind), TAS (Temperature)
winter_ds = ds[['psl', 'uas', 'vas', 'tas']].sel(time=slice('2020-01-01', '2020-01-31')).isel(healpix_index=valid_indices)

# Resample to daily means to capture individual storm passages
daily = winter_ds.resample(time='1D').mean()

# Calculate anomalies (deviation from the 30-day mean state)
mean_state = daily.mean(dim='time')
anomalies = daily - mean_state

print("Calculating 2D transient eddy components...")
# 1. Storm Track Intensity (Daily PSL variance)
eddy_intensity = anomalies['psl']**2

# 2. Transient Heat Flux Vectors (u'T' and v'T')
u_flux = anomalies['uas'] * anomalies['tas']
v_flux = anomalies['vas'] * anomalies['tas']

# Total Heat Flux Magnitude
total_flux_mag = np.sqrt(u_flux**2 + v_flux**2)

print("Computing pixel-wise temporal correlation...")
# Vectorized Pearson correlation across the time dimension using Xarray
# This correlates local storm intensity with local heat transport across 30 days
correlation_2d = xr.corr(eddy_intensity, total_flux_mag, dim='time')

# Bring the final 2D arrays into memory for plotting
psl_var_map = eddy_intensity.mean(dim='time').compute().values
flux_mag_map = total_flux_mag.mean(dim='time').compute().values
u_flux_mean = u_flux.mean(dim='time').compute().values
v_flux_mean = v_flux.mean(dim='time').compute().values
corr_map = correlation_2d.compute().values

# =====================================================================
# 3. 2D Contour Plotting
# =====================================================================
fig = plt.figure(figsize=(20, 6))
fig.suptitle('2D Spatial Relationship: Winter Cyclones & Sensible Heat Transport', fontsize=18, y=1.05)
proj = ccrs.PlateCarree()

# Helper function to format map panels
def setup_map(ax, title):
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=proj)
    ax.coastlines(color='black', linewidth=1)
    ax.gridlines(draw_labels=True, linestyle='--', alpha=0.4)
    ax.set_title(title, fontsize=14, pad=10)

# Panel 1: Storm Track Intensity Map
ax1 = plt.subplot(1, 3, 1, projection=proj)
setup_map(ax1, 'Storm Track Intensity\n(Mean $psl\'^2$)')
c1 = ax1.tricontourf(lon_clean, lat_clean, psl_var_map, levels=15, cmap='Purples', transform=proj)
plt.colorbar(c1, ax=ax1, orientation='horizontal', pad=0.08, label='Variance (Pa²)')

# Panel 2: Total Eddy Heat Flux with Directional Vectors
ax2 = plt.subplot(1, 3, 2, projection=proj)
setup_map(ax2, 'Eddy Heat Flux Magnitude & Direction\n($\sqrt{(u\'T\')^2 + (v\'T\')^2}$)')
c2 = ax2.tricontourf(lon_clean, lat_clean, flux_mag_map, levels=15, cmap='YlOrRd', transform=proj)
plt.colorbar(c2, ax=ax2, orientation='horizontal', pad=0.08, label='Flux (m/s * K)')

# Overlay Quiver Vectors (Stride every 80 points to prevent clutter)
stride = 80
ax2.quiver(lon_clean[::stride], lat_clean[::stride], u_flux_mean[::stride], v_flux_mean[::stride], 
           transform=proj, color='black', alpha=0.6, scale=400, headwidth=4)

# Panel 3: Temporal Correlation Contour Map
ax3 = plt.subplot(1, 3, 3, projection=proj)
setup_map(ax3, 'Local Temporal Correlation\n(Intensity vs. Transport)')
# Using a diverging colormap (RdBu) centered on zero for correlation
c3 = ax3.tricontourf(lon_clean, lat_clean, corr_map, levels=np.linspace(-1, 1, 21), cmap='RdBu_r', transform=proj)
plt.colorbar(c3, ax=ax3, orientation='horizontal', pad=0.08, label='Pearson Correlation (r)')

plt.tight_layout()

print("Saving figure to '2d_heat_transport_correlation.png'...")
plt.savefig('2d_heat_transport_correlation.png', dpi=300, bbox_inches='tight')
plt.show()