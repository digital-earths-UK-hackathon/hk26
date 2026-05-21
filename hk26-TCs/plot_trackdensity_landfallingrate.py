import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Point

datasets = {
    "n1280CoMA9v2": "./syclops/n1280CoMA9v2/n1280CoMA9v2_track_20200201-20210301_6h_tc_psl.csv",
    "n1280GAL9v2": "./syclops/n1280GAL9v2/n1280GAL9v2_track_20200201-20210301_6h_tc_psl.csv",
    "n2560CoMA9": "./syclops/n2560CoMA9/n2560CoMA9_track_20200201-20210301_6h_tc_psl.csv",
    "n2560RAL3p3regridv2": "./syclops/n2560RAL3p3regridv2/n2560RAL3p3regridv2_track_20200201-20210301_6h_tc_psl.csv",
}

ibtracs_file = "IBTrACS.ALL.v04r01.nc"
land_file = "ne_110m_land.zip"

start_date = np.datetime64("2020-02-01")
end_date = np.datetime64("2021-03-01")

track_col = "track_id"
lon_col = "lon"
lat_col = "lat"

lon_edges = np.arange(-180, 185, 5)
lat_edges = np.arange(-60, 65, 5)

land = gpd.read_file(land_file)
land_geom = land.union_all()


def wrap_lon_to_180(lon):
    return ((np.asarray(lon) + 180) % 360) - 180


def has_landfall(lon, lat):
    lon = wrap_lon_to_180(lon)
    lat = np.asarray(lat)

    valid = np.isfinite(lon) & np.isfinite(lat)
    lon = lon[valid]
    lat = lat[valid]

    if len(lon) < 2:
        return False

    points = gpd.GeoSeries(
        [Point(x, y) for x, y in zip(lon, lat)],
        crs="EPSG:4326"
    )

    is_land = points.within(land_geom).values

    # Landfall = ocean-to-land transition
    return np.any((~is_land[:-1]) & (is_land[1:]))


def get_csv_density_and_landfall(csv_file):
    df = pd.read_csv(csv_file)
    df.columns = df.columns.str.strip()

    df["time"] = pd.to_datetime(
        dict(
            year=df["year"],
            month=df["month"],
            day=df["day"],
            hour=df["hour"]
        )
    )

    duration = df.groupby(track_col)["time"].agg(
        lambda x: (x.max() - x.min()).total_seconds() / 3600
    )

    valid_tracks = duration[duration > 24].index
    df = df[df[track_col].isin(valid_tracks)].copy()

    total_storms = 0
    landfall_storms = 0

    all_lon = []
    all_lat = []

    for tid, storm in df.groupby(track_col):
        storm = storm.sort_values("time")

        lon = wrap_lon_to_180(storm[lon_col].values)
        lat = storm[lat_col].values

        valid_full = np.isfinite(lon) & np.isfinite(lat)

        lon_full = lon[valid_full]
        lat_full = lat[valid_full]

        if len(lon_full) == 0:
            continue

        total_storms += 1

        if has_landfall(lon_full, lat_full):
            landfall_storms += 1

        valid_plot = valid_full & (lat >= -60) & (lat <= 60)

        all_lon.append(lon[valid_plot])
        all_lat.append(lat[valid_plot])

    if len(all_lon) > 0:
        all_lon = np.concatenate(all_lon)
        all_lat = np.concatenate(all_lat)
    else:
        all_lon = np.array([])
        all_lat = np.array([])

    density, _, _ = np.histogram2d(
        all_lat,
        all_lon,
        bins=[lat_edges, lon_edges]
    )

    landfall_rate = (
        landfall_storms / total_storms * 100
        if total_storms > 0 else np.nan
    )

    return density, total_storms, landfall_storms, landfall_rate


def get_ibtracs_density_and_landfall(ibtracs_file, start_date, end_date):
    ds = xr.open_dataset(ibtracs_file)

    lon = ds["lon"]
    lat = ds["lat"]
    time = ds["time"]

    mask_time = (time >= start_date) & (time < end_date)

    # Data1-4 are 6-hourly, so use only synoptic 6-hourly IBTrACS records
    mask_6h = time.dt.hour.isin([0, 6, 12, 18])
    mask = mask_time & mask_6h

    lon_period = lon.where(mask)
    lat_period = lat.where(mask)

    storm_dim = lon.dims[0]

    total_storms = 0
    landfall_storms = 0

    all_lon = []
    all_lat = []

    for i in range(ds.sizes[storm_dim]):

        lon_i = wrap_lon_to_180(lon_period.isel({storm_dim: i}).values)
        lat_i = lat_period.isel({storm_dim: i}).values

        valid = np.isfinite(lon_i) & np.isfinite(lat_i)

        lon_i = lon_i[valid]
        lat_i = lat_i[valid]

        # >24 h using 6-hourly data means at least 6 points:
        # 0, 6, 12, 18, 24, 30 h
        if len(lon_i) < 6:
            continue

        total_storms += 1

        if has_landfall(lon_i, lat_i):
            landfall_storms += 1

        valid_plot = (lat_i >= -60) & (lat_i <= 60)

        all_lon.append(lon_i[valid_plot])
        all_lat.append(lat_i[valid_plot])

    if len(all_lon) > 0:
        all_lon = np.concatenate(all_lon)
        all_lat = np.concatenate(all_lat)
    else:
        all_lon = np.array([])
        all_lat = np.array([])

    density, _, _ = np.histogram2d(
        all_lat,
        all_lon,
        bins=[lat_edges, lon_edges]
    )

    landfall_rate = (
        landfall_storms / total_storms * 100
        if total_storms > 0 else np.nan
    )

    return density, total_storms, landfall_storms, landfall_rate


ib_density, ib_total, ib_landfall, ib_rate = get_ibtracs_density_and_landfall(
    ibtracs_file,
    start_date,
    end_date
)

print(f"IBTrACS: total={ib_total}, landfall={ib_landfall}, rate={ib_rate:.2f}%")

density_orig = {}
density_diff = {}
landfall_stats = {}

for name, csv_file in datasets.items():

    model_density, total, landfall, rate = get_csv_density_and_landfall(csv_file)

    density_orig[name] = model_density
    density_diff[name] = model_density - ib_density

    landfall_stats[name] = {
        "total": total,
        "landfall": landfall,
        "rate": rate
    }

    print(f"{name}: total={total}, landfall={landfall}, rate={rate:.2f}%")


diff_absmax = max(np.nanmax(np.abs(v)) for v in density_diff.values())
diff_vmin = -diff_absmax
diff_vmax = diff_absmax

density_vmax = max(
    max(np.nanmax(v) for v in density_orig.values()),
    np.nanmax(ib_density)
)


# ======================================================
# Plot IBTrACS original track density
# ======================================================

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())

ax.coastlines(linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linewidth=0.4)
ax.add_feature(cfeature.LAND, alpha=0.35)

pcm = ax.pcolor(
    lon_edges,
    lat_edges,
    ib_density,
    cmap="viridis",
    vmin=0,
    vmax=density_vmax,
    transform=ccrs.PlateCarree()
)

cbar = plt.colorbar(
    pcm,
    ax=ax,
    orientation="horizontal",
    pad=0.07,
    aspect=40
)

cbar.set_label("Track density per 5° × 5° grid box")

ax.text(
    0.02,
    0.03,
    f"IBTrACS: {ib_rate:.1f}% ({ib_landfall}/{ib_total})",
    transform=ax.transAxes,
    ha="left",
    va="bottom",
    fontsize=10,
    bbox=dict(facecolor="white", edgecolor="black", alpha=0.8)
)

gl = ax.gridlines(draw_labels=True, linewidth=0.3, linestyle="--")
gl.top_labels = False
gl.right_labels = False

ax.set_title("IBTrACS Track Density")

plt.tight_layout()
plt.savefig(
    "IBTrACS_track_density_5deg.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()


# ======================================================
# Plot model original density and model minus IBTrACS
# ======================================================

for name in datasets.keys():

    out_name = name.replace(" ", "_")

    rate = landfall_stats[name]["rate"]
    total = landfall_stats[name]["total"]
    landfall = landfall_stats[name]["landfall"]

    text = (
        f"{name}: {rate:.1f}% ({landfall}/{total})\n"
        f"IBTrACS: {ib_rate:.1f}% ({ib_landfall}/{ib_total})"
    )

    # --------------------------------------------------
    # Original model track density
    # --------------------------------------------------

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.4)
    ax.add_feature(cfeature.LAND, alpha=0.35)

    pcm = ax.pcolor(
        lon_edges,
        lat_edges,
        density_orig[name],
        cmap="viridis",
        vmin=0,
        vmax=density_vmax,
        transform=ccrs.PlateCarree()
    )

    cbar = plt.colorbar(
        pcm,
        ax=ax,
        orientation="horizontal",
        pad=0.07,
        aspect=40
    )

    cbar.set_label("Track density per 5° × 5° grid box")

    ax.text(
        0.02,
        0.03,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox=dict(facecolor="white", edgecolor="black", alpha=0.8)
    )

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False

    ax.set_title(f"{name} Track Density")

    plt.tight_layout()
    plt.savefig(
        f"{out_name}_track_density_5deg.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # --------------------------------------------------
    # Model minus IBTrACS track density
    # --------------------------------------------------

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.4)
    ax.add_feature(cfeature.LAND, alpha=0.35)

    pcm = ax.pcolor(
        lon_edges,
        lat_edges,
        density_diff[name],
        cmap="RdBu_r",
        vmin=diff_vmin,
        vmax=diff_vmax,
        transform=ccrs.PlateCarree()
    )

    cbar = plt.colorbar(
        pcm,
        ax=ax,
        orientation="horizontal",
        pad=0.07,
        aspect=40
    )

    cbar.set_label("Track density difference per 5° × 5° grid box")

    ax.text(
        0.02,
        0.03,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox=dict(facecolor="white", edgecolor="black", alpha=0.8)
    )

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False

    ax.set_title(f"{name} minus IBTrACS Track Density")

    plt.tight_layout()
    plt.savefig(
        f"{out_name}_track_density_minus_IBTrACS_5deg.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()


rows = []

for name, stats in landfall_stats.items():
    rows.append({
        "Dataset": name,
        "Total Storms": stats["total"],
        "Landfalling Storms": stats["landfall"],
        "Landfall Rate (%)": stats["rate"]
    })

rows.append({
    "Dataset": "IBTrACS",
    "Total Storms": ib_total,
    "Landfalling Storms": ib_landfall,
    "Landfall Rate (%)": ib_rate
})

results_df = pd.DataFrame(rows)
results_df.to_csv("landfall_statistics_with_IBTrACS.csv", index=False)

print("Saved:")
print("IBTrACS_track_density_5deg.png")
print("Data_*_track_density_5deg.png")
print("Data_*_track_density_minus_IBTrACS_5deg.png")
print("landfall_statistics_with_IBTrACS.csv")
