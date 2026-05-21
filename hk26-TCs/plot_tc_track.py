import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

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

start_date = np.datetime64("2020-02-01")
end_date = np.datetime64("2021-03-01")

track_col = "track_id"
wind_col = "sfcWind_max"   # m/s for Data1-4
lon_col = "lon"
lat_col = "lat"

colors = {
    "10km CoMA9": "violet",
    "10km GAL": "b",
    "5km CoMA9": "darkviolet",
    "5km RAL": "darkorange",
    "IBTrACS": "black"
}
# ======================================================
# FUNCTIONS
# ======================================================

def wrap_lon_to_180(lon):
    return ((np.asarray(lon) + 180) % 360) - 180


def split_track_for_plot(lon, lat):
    lon_plot = wrap_lon_to_180(lon)
    lat = np.asarray(lat)

    segments = []
    start = 0

    for i in range(1, len(lon_plot)):
        if abs(lon_plot[i] - lon_plot[i - 1]) > 180:
            segments.append((lon_plot[start:i], lat[start:i]))
            start = i

    segments.append((lon_plot[start:], lat[start:]))
    return segments


def read_ibtracs_lmi_ace(ibtracs_file, start_date, end_date):
    ds = xr.open_dataset(ibtracs_file)

    if "usa_wind" in ds:
        wind = ds["usa_wind"]   # kt
    elif "wmo_wind" in ds:
        wind = ds["wmo_wind"]   # kt
    else:
        raise ValueError("Cannot find usa_wind or wmo_wind in IBTrACS file.")

    time = ds["time"]

    mask_time = (time >= start_date) & (time < end_date)
    wind_period = wind.where(mask_time)

    # LMI, kt -> m/s
    lmi_kt = wind_period.max(dim="date_time", skipna=True)
    lmi_kt = lmi_kt.where(np.isfinite(lmi_kt), drop=True)
    lmi_ms = lmi_kt.values * 0.514444

    # IBTrACS is 3-hourly: ACE weighted by dt / 6 h
    dt_hours = time.diff(dim="date_time") / np.timedelta64(1, "h")
    dt_hours = xr.concat(
        [dt_hours.isel(date_time=0), dt_hours],
        dim="date_time"
    )
    dt_hours = dt_hours.assign_coords(date_time=time["date_time"])

    ace_each = (
        (wind_period.where(wind_period >= 34) ** 2)
        * (dt_hours / 6.0)
    ).sum(dim="date_time", skipna=True) * 1e-4

    ace_each = ace_each.where(ace_each > 0, drop=True)
    total_ace = float(ace_each.sum(skipna=True).values)

    return pd.Series(lmi_ms), total_ace


# ======================================================
# READ DATA1-4, CALCULATE LMI/ACE, PLOT TRACKS
# ======================================================

all_lmi = {}
all_ace = {}
all_genesis_lat = {}

for name, csv_file in datasets.items():

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

    df["wind_kt"] = df[wind_col] * 1.94384

    lmi = df.groupby(track_col)[wind_col].max()

    ace_each = df.groupby(track_col)["wind_kt"].apply(
        lambda x: np.sum((x[x >= 34]) ** 2) * 1e-4
    )

    all_lmi[name] = lmi
    all_ace[name] = ace_each.sum()

    genesis = df.sort_values("time").groupby(track_col).first()
    all_genesis_lat[name] = genesis[lat_col].values

    # ==================================================
    # Track figure
    # ==================================================

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.4)
    ax.add_feature(cfeature.LAND, alpha=0.3)
    ax.add_feature(cfeature.OCEAN, alpha=0.2)

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False

    for j, (tid, storm) in enumerate(df.groupby(track_col)):

        storm = storm.sort_values("time")

        lon = storm[lon_col].values
        lat = storm[lat_col].values

        color = f"C{j % 10}"

        for lon_seg, lat_seg in split_track_for_plot(lon, lat):
            ax.plot(
                lon_seg,
                lat_seg,
                linewidth=1.0,
                alpha=0.7,
                color=color,
                transform=ccrs.PlateCarree()
            )

        ax.scatter(
            wrap_lon_to_180(lon[0]),
            lat[0],
            s=12,
            color=color,
            marker="o",
            transform=ccrs.PlateCarree()
        )

    ax.set_title(f"{name}: Storm Tracks Longer Than 24 h")

    out_name = name.replace(" ", "_")
    plt.tight_layout()
    plt.savefig(
        f"{out_name}_tracks_gt24h.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()


# ======================================================
# READ IBTrACS LMI/ACE
# ======================================================

ibtracs_lmi, ibtracs_ace = read_ibtracs_lmi_ace(
    ibtracs_file,
    start_date,
    end_date
)

all_lmi["IBTrACS"] = ibtracs_lmi
all_ace["IBTrACS"] = ibtracs_ace

# ======================================================
# LMI CATEGORY: PERCENTAGE
# ======================================================

fig, ax = plt.subplots(figsize=(10, 6))

cat_edges = [0, 18, 33, 43, 50, 58, 70, np.inf]
cat_labels = ["TD", "TS", "C1", "C2", "C3", "C4", "C5"]

x = np.arange(len(cat_labels))
n_data = len(all_lmi)
bar_width = 0.8 / n_data

for i, (name, lmi) in enumerate(all_lmi.items()):

    counts, _ = np.histogram(lmi, bins=cat_edges)
    percent = counts / counts.sum() * 100

    ax.bar(
        x + i * bar_width,
        percent,
        width=bar_width,
        color=colors[name],
        edgecolor="black",
        linewidth=0.4,
        alpha=0.8,
        label=f"{name}, ACE={all_ace[name]:.1f}"
    )

ax.set_xticks(x + bar_width * (n_data - 1) / 2)
ax.set_xticklabels(cat_labels)

ax.set_xlabel("LMI Category")
ax.set_ylabel("Percentage (%)")
ax.set_title("LMI Category Distribution")
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig(
    "LMI_category_percentage_4datasets_IBTrACS.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ======================================================
# LMI CATEGORY: STORM NUMBER
# ======================================================

fig, ax = plt.subplots(figsize=(10, 6))

for i, (name, lmi) in enumerate(all_lmi.items()):

    counts, _ = np.histogram(lmi, bins=cat_edges)

    ax.bar(
        x + i * bar_width,
        counts,
        width=bar_width,
        color=colors[name],
        edgecolor="black",
        linewidth=0.4,
        alpha=0.8,
        label=f"{name}, N={len(lmi)}"
    )

ax.set_xticks(x + bar_width * (n_data - 1) / 2)
ax.set_xticklabels(cat_labels)

ax.set_xlabel("LMI Category")
ax.set_ylabel("Number of Storms")
ax.set_title("Number of Storms by LMI Category")
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig(
    "LMI_category_counts_4datasets_IBTrACS.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ======================================================
# GENESIS LATITUDE PDF
# ======================================================

fig, ax = plt.subplots(figsize=(9, 6))

bins = np.linspace(-60, 60, 40)

for name, genesis_lat in all_genesis_lat.items():

    hist, edges = np.histogram(
        genesis_lat,
        bins=bins,
        density=True
    )

    centers = 0.5 * (edges[:-1] + edges[1:])
    pdf_percent = hist * 100

    ax.plot(
        centers,
        pdf_percent,
        linewidth=2,
        color=colors[name],
        label=name
    )

ax.set_xlabel("Genesis Latitude")
ax.set_ylabel("PDF (%)")
ax.set_title("PDF of TC Genesis Latitude")
ax.legend(frameon=False)

plt.tight_layout()
plt.savefig(
    "Genesis_Latitude_PDF_4datasets.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ======================================================
# SAVE SUMMARY
# ======================================================

summary_rows = []

for name in all_lmi.keys():
    summary_rows.append({
        "Dataset": name,
        "Storm Number": len(all_lmi[name]),
        "ACE": all_ace[name],
        "Mean LMI (m/s)": np.nanmean(all_lmi[name]),
        "Max LMI (m/s)": np.nanmax(all_lmi[name])
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("LMI_ACE_summary_4datasets_IBTrACS.csv", index=False)

print("Saved:")
print("Data_*_tracks_gt24h.png")
print("LMI_category_percentage_4datasets_IBTrACS.png")
print("LMI_category_counts_4datasets_IBTrACS.png")
print("Genesis_Latitude_PDF_4datasets.png")
print("LMI_ACE_summary_4datasets_IBTrACS.csv")

print("\nACE values:")
for name, ace in all_ace.items():
    print(f"{name}: {ace:.2f}")
