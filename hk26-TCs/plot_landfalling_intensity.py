import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point

# ======================================================
# USER SETTINGS
# ======================================================
datasets = {
    "10km CoMA9": "./syclops/n1280CoMA9v2/n1280CoMA9v2_track_20200201-20210301_6h_tc_psl.csv",
    "10km GAL": "./syclops/n1280GAL9v2/n1280GAL9v2_track_20200201-20210301_6h_tc_psl.csv",
    "5km CoMA9": "./syclops/n2560CoMA9/n2560CoMA9_track_20200201-20210301_6h_tc_psl.csv",
    "5km RAL": "./syclops/n2560RAL3p3regridv2/n2560RAL3p3regridv2_track_20200201-20210301_6h_tc_psl.csv",
}

ibtracs_file = "IBTrACS.ALL.v04r01.nc"
land_file = "ne_110m_land.zip"

start_date = np.datetime64("2020-02-01")
end_date = np.datetime64("2021-03-01")

track_col = "track_id"
lon_col = "lon"
lat_col = "lat"
wind_col = "sfcWind_max"   # m/s for Data1-4

# IBTrACS is 3-hourly, but Data1-4 are 6-hourly.
# For landfall intensity comparison, use only 6-hourly IBTrACS points.
USE_IBTRACS_6H_ONLY = True

colors = {
    "10km CoMA9": "violet",
    "10km GAL": "b",
    "5km CoMA9": "darkviolet",
    "5km RAL": "darkorange",
    "IBTrACS": "black"
}

# ======================================================
# READ LAND
# ======================================================

land = gpd.read_file(land_file)
land_geom = land.union_all()


def wrap_lon_to_180(lon):
    return ((np.asarray(lon) + 180) % 360) - 180


def landfall_intensity_from_track(lon, lat, intensity):
    """
    Return intensity at the point immediately before first landfall.
    Landfall is defined as first ocean-to-land transition.
    """

    lon = wrap_lon_to_180(lon)
    lat = np.asarray(lat)
    intensity = np.asarray(intensity)

    valid = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(intensity)

    lon = lon[valid]
    lat = lat[valid]
    intensity = intensity[valid]

    if len(lon) < 2:
        return np.nan

    points = gpd.GeoSeries(
        [Point(x, y) for x, y in zip(lon, lat)],
        crs="EPSG:4326"
    )

    is_land = points.within(land_geom).values

    transition = (~is_land[:-1]) & (is_land[1:])

    if np.any(transition):
        # transition from ocean point i to land point i+1
        # use intensity at point i, before landfall
        idx = np.where(transition)[0][0]
        return intensity[idx]

    return np.nan


def read_csv_landfall_intensity_and_lmi(csv_file):
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

    landfall_intensity = []
    landfall_lmi = []

    for tid, storm in df.groupby(track_col):
        storm = storm.sort_values("time")

        lf_int = landfall_intensity_from_track(
            storm[lon_col].values,
            storm[lat_col].values,
            storm[wind_col].values
        )

        if np.isfinite(lf_int):
            landfall_intensity.append(lf_int)
            landfall_lmi.append(storm[wind_col].max())

    return np.asarray(landfall_intensity), np.asarray(landfall_lmi)


def read_ibtracs_landfall_intensity_and_lmi(ibtracs_file, start_date, end_date):
    ds = xr.open_dataset(ibtracs_file)

    lon = ds["lon"]
    lat = ds["lat"]
    time = ds["time"]

    if "usa_wind" in ds:
        wind = ds["usa_wind"]   # kt
    elif "wmo_wind" in ds:
        wind = ds["wmo_wind"]   # kt
    else:
        raise ValueError("Cannot find usa_wind or wmo_wind in IBTrACS file.")

    mask_time = (time >= start_date) & (time < end_date)

    if USE_IBTRACS_6H_ONLY:
        mask_6h = time.dt.hour.isin([0, 6, 12, 18])
        mask = mask_time & mask_6h
    else:
        mask = mask_time

    lon_period = lon.where(mask)
    lat_period = lat.where(mask)
    wind_period = wind.where(mask)

    storm_dim = lon.dims[0]

    landfall_intensity = []
    landfall_lmi = []

    for i in range(ds.sizes[storm_dim]):

        lon_i = lon_period.isel({storm_dim: i}).values
        lat_i = lat_period.isel({storm_dim: i}).values
        wind_i_kt = wind_period.isel({storm_dim: i}).values

        valid = (
            np.isfinite(lon_i)
            & np.isfinite(lat_i)
            & np.isfinite(wind_i_kt)
        )

        lon_i = lon_i[valid]
        lat_i = lat_i[valid]
        wind_i_kt = wind_i_kt[valid]

        if USE_IBTRACS_6H_ONLY:
            if len(lon_i) < 6:
                continue
        else:
            if len(lon_i) < 10:
                continue

        wind_i_ms = wind_i_kt * 0.514444

        lf_int = landfall_intensity_from_track(
            lon_i,
            lat_i,
            wind_i_ms
        )

        if np.isfinite(lf_int):
            landfall_intensity.append(lf_int)
            landfall_lmi.append(np.nanmax(wind_i_ms))

    return np.asarray(landfall_intensity), np.asarray(landfall_lmi)


# ======================================================
# CALCULATE
# ======================================================

all_landfall_intensity = {}
all_landfall_lmi = {}

for name, csv_file in datasets.items():

    lf_int, lf_lmi = read_csv_landfall_intensity_and_lmi(csv_file)

    all_landfall_intensity[name] = lf_int
    all_landfall_lmi[name] = lf_lmi

    print(
        f"{name}: landfalling storms = {len(lf_int)}, "
        f"mean landfall intensity = {np.nanmean(lf_int):.2f} m/s, "
        f"mean landfalling-TC LMI = {np.nanmean(lf_lmi):.2f} m/s"
    )

ib_lf_int, ib_lf_lmi = read_ibtracs_landfall_intensity_and_lmi(
    ibtracs_file,
    start_date,
    end_date
)

all_landfall_intensity["IBTrACS"] = ib_lf_int
all_landfall_lmi["IBTrACS"] = ib_lf_lmi

print(
    f"IBTrACS: landfalling storms = {len(ib_lf_int)}, "
    f"mean landfall intensity = {np.nanmean(ib_lf_int):.2f} m/s, "
    f"mean landfalling-TC LMI = {np.nanmean(ib_lf_lmi):.2f} m/s"
)

# ======================================================
# SAVE TABLE
# ======================================================

rows = []

for name in all_landfall_intensity.keys():
    for lf_int, lf_lmi in zip(
        all_landfall_intensity[name],
        all_landfall_lmi[name]
    ):
        rows.append({
            "Dataset": name,
            "Landfall Intensity Before Landfall (m/s)": lf_int,
            "LMI of Landfalling TC (m/s)": lf_lmi
        })

out_df = pd.DataFrame(rows)
out_df.to_csv(
    "landfalling_intensity_and_LMI_4datasets_IBTrACS.csv",
    index=False
)

# ======================================================
# CATEGORY SETTINGS
# ======================================================

cat_edges = [0, 18, 33, 43, 50, 58, 70, np.inf]
cat_labels = ["TD", "TS", "C1", "C2", "C3", "C4", "C5"]

x = np.arange(len(cat_labels))
n_data = len(all_landfall_intensity)
bar_width = 0.8 / n_data

# ======================================================
# 1. LANDFALLING INTENSITY CATEGORY: PERCENTAGE
# ======================================================

fig, ax = plt.subplots(figsize=(10, 6))

for i, (name, values) in enumerate(all_landfall_intensity.items()):

    if len(values) == 0:
        continue

    counts, _ = np.histogram(values, bins=cat_edges)
    percent = counts / counts.sum() * 100

    ax.bar(
        x + i * bar_width,
        percent,
        width=bar_width,
        color=colors[name],
        edgecolor="black",
        linewidth=0.4,
        alpha=0.8,
        label=f"{name}, N={len(values)}"
    )

ax.set_xticks(x + bar_width * (n_data - 1) / 2)
ax.set_xticklabels(cat_labels)

ax.set_xlabel("Landfalling Intensity Category")
ax.set_ylabel("Percentage (%)")
ax.set_title("Landfalling Intensity Distribution")
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig(
    "landfalling_intensity_category_percentage_4datasets_IBTrACS.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ======================================================
# 2. LANDFALLING INTENSITY CATEGORY: NUMBER
# ======================================================

fig, ax = plt.subplots(figsize=(10, 6))

for i, (name, values) in enumerate(all_landfall_intensity.items()):

    if len(values) == 0:
        continue

    counts, _ = np.histogram(values, bins=cat_edges)

    ax.bar(
        x + i * bar_width,
        counts,
        width=bar_width,
        color=colors[name],
        edgecolor="black",
        linewidth=0.4,
        alpha=0.8,
        label=f"{name}, N={len(values)}"
    )

ax.set_xticks(x + bar_width * (n_data - 1) / 2)
ax.set_xticklabels(cat_labels)

ax.set_xlabel("Landfalling Intensity Category")
ax.set_ylabel("Number of Storms")
ax.set_title("Number of Landfalling Storms by Intensity Category")
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig(
    "landfalling_intensity_category_counts_4datasets_IBTrACS.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ======================================================
# 3. LMI OF LANDFALLING TCs: PERCENTAGE
# ======================================================

fig, ax = plt.subplots(figsize=(10, 6))

for i, (name, values) in enumerate(all_landfall_lmi.items()):

    if len(values) == 0:
        continue

    counts, _ = np.histogram(values, bins=cat_edges)
    percent = counts / counts.sum() * 100

    ax.bar(
        x + i * bar_width,
        percent,
        width=bar_width,
        color=colors[name],
        edgecolor="black",
        linewidth=0.4,
        alpha=0.8,
        label=f"{name}, N={len(values)}"
    )

ax.set_xticks(x + bar_width * (n_data - 1) / 2)
ax.set_xticklabels(cat_labels)

ax.set_xlabel("LMI Category")
ax.set_ylabel("Percentage (%)")
ax.set_title("LMI Distribution of Landfalling TCs")
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig(
    "landfalling_TC_LMI_category_percentage_4datasets_IBTrACS.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ======================================================
# 4. LMI OF LANDFALLING TCs: NUMBER
# ======================================================

fig, ax = plt.subplots(figsize=(10, 6))

for i, (name, values) in enumerate(all_landfall_lmi.items()):

    if len(values) == 0:
        continue

    counts, _ = np.histogram(values, bins=cat_edges)

    ax.bar(
        x + i * bar_width,
        counts,
        width=bar_width,
        color=colors[name],
        edgecolor="black",
        linewidth=0.4,
        alpha=0.8,
        label=f"{name}, N={len(values)}"
    )

ax.set_xticks(x + bar_width * (n_data - 1) / 2)
ax.set_xticklabels(cat_labels)

ax.set_xlabel("LMI Category")
ax.set_ylabel("Number of Storms")
ax.set_title("Number of Landfalling TCs by LMI Category")
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig(
    "landfalling_TC_LMI_category_counts_4datasets_IBTrACS.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Saved:")
print("landfalling_intensity_and_LMI_4datasets_IBTrACS.csv")
print("landfalling_intensity_category_percentage_4datasets_IBTrACS.png")
print("landfalling_intensity_category_counts_4datasets_IBTrACS.png")
print("landfalling_TC_LMI_category_percentage_4datasets_IBTrACS.png")
print("landfalling_TC_LMI_category_counts_4datasets_IBTrACS.png")
