"""
Paleoclimate Global Temperature Map Series (0 to 50,000 BP)
With Overlaid Radiocarbon Date Occurrences by Material Class:
- Bone Collagen (Purified Collagen Dates)
- Non-Collagen Controls & Inorganics (Charcoal, Wood, Plant Remains, Shell, Bone Apatite)

Capturing the Geographic Disappearance of Bone Collagen across Climate Transitions:
- Panel A (0 - 9 ka BP): Late Holocene Baseline (0 BP)
- Panel B (9 - 12 ka BP): Early Holocene Optimum (9 ka BP)
- Panel C (12 - 14 ka BP): Younger Dryas Cold Reversal (12 ka BP proxy for 12.3 ka)
- Panel D (14 - 20 ka BP): Bølling–Allerød Deglacial Warming (14 ka BP)
- Panel E (20 - 30 ka BP): Last Glacial Maximum (20 ka BP)
- Panel F (30 - 40 ka BP): Late MIS 3 Interstadial (30 ka BP)
- Panel G (40 - 50 ka BP): Mid MIS 3 / AMH Dispersal (40 ka BP)
- Panel H (>= 50 ka BP / Instrumental Boundary): Early MIS 3 / Radiocarbon Wall (50 ka BP)
"""

import os
import shutil
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

# File Paths
output_dir = r"D:\26 Modelling Collagen Hydrolysis\outputs"
artifact_dir = r"C:\Users\matth\.gemini\antigravity-ide\brain\cd905880-511d-460a-944c-c857596f9afc"
nc_path = r"D:\26 Modelling Collagen Hydrolysis\pastclim\Beyer2020_annual_vars_v1.2.2.nc"
coastline_path = r"D:\26 Modelling Collagen Hydrolysis\pastclim\coastline\ne_110m_coastline.shp"
c14_path = r"D:\26 Modelling Collagen Hydrolysis\data\harmonized_c14_master.parquet"

print("Loading datasets...")
ds = xr.open_dataset(nc_path, decode_times=False)
coastline = gpd.read_file(coastline_path)
df_c14 = pd.read_parquet(c14_path)

# Filter valid lat/lon and age
df_c14 = df_c14.dropna(subset=["latitude", "longitude", "c14_age"]).copy()
# Longitude normalization to [-180, 180]
df_c14["longitude"] = np.where(df_c14["longitude"] > 180.0, df_c14["longitude"] - 360.0, df_c14["longitude"])

# Grid coordinates
lats = ds["latitude"].values
lons = ds["longitude"].values

# Definition of the 8 time intervals & corresponding climate model snapshots
slices_info = [
    {
        "panel": "A",
        "time_bp_model": 0,
        "age_min": 0,
        "age_max": 9000,
        "title": "A. 0 – 9 ka BP (Climate: 0 ka Modern Baseline)",
        "badge": "0 – 9 ka BP",
        "context": "Global bone collagen widespread across temperate & southern zones",
        "badge_color": "#1B7837"
    },
    {
        "panel": "B",
        "time_bp_model": -9000,
        "age_min": 9000,
        "age_max": 12000,
        "title": "B. 9 – 12 ka BP (Climate: 9 ka Early Holocene Optimum)",
        "badge": "9 – 12 ka BP",
        "context": "Post-glacial thermal optimum; collagen extends into southern latitudes",
        "badge_color": "#313695"
    },
    {
        "panel": "C",
        "time_bp_model": -12000,
        "age_min": 12000,
        "age_max": 14000,
        "title": "C. 12 – 14 ka BP (Climate: 12 ka Younger Dryas)",
        "badge": "12 – 14 ka BP",
        "context": "Abrupt Younger Dryas cooling; Southern Hemisphere collagen contracts",
        "badge_color": "#053061"
    },
    {
        "panel": "D",
        "time_bp_model": -14000,
        "age_min": 14000,
        "age_max": 20000,
        "title": "D. 14 – 20 ka BP (Climate: 14 ka Bølling–Allerød Warm Pulse)",
        "badge": "14 – 20 ka BP",
        "context": "Deglacial transition; collagen restricted to latitudes > 30°N!",
        "badge_color": "#B2182B"
    },
    {
        "panel": "E",
        "time_bp_model": -20000,
        "age_min": 20000,
        "age_max": 30000,
        "title": "E. 20 – 30 ka BP (Climate: 20 ka Last Glacial Maximum)",
        "badge": "20 – 30 ka BP",
        "context": "LGM peak cold; zero collagen in tropics; permafrost preservation in north",
        "badge_color": "#2166AC"
    },
    {
        "panel": "F",
        "time_bp_model": -30000,
        "age_min": 30000,
        "age_max": 40000,
        "title": "F. 30 – 40 ka BP (Climate: 30 ka Late MIS 3)",
        "badge": "30 – 40 ka BP",
        "context": "Transitional stadial; dense European cave collagen; African/tropical bone absent",
        "badge_color": "#4393C3"
    },
    {
        "panel": "G",
        "time_bp_model": -40000,
        "age_min": 40000,
        "age_max": 50000,
        "title": "G. 40 – 50 ka BP (Climate: 40 ka Mid MIS 3)",
        "badge": "40 – 50 ka BP",
        "context": "Neanderthal-AMH horizon; collagen strictly confined to high-latitude Eurasia",
        "badge_color": "#762A83"
    },
    {
        "panel": "H",
        "time_bp_model": -50000,
        "age_min": 50000,
        "age_max": 65000,
        "title": "H. >= 50 ka BP (Climate: 50 ka Early MIS 3 / 14C Blank Ceiling)",
        "badge": "≥ 50 ka BP",
        "context": "Instrumental 14C background limit; rare infinite/near-blank dates only",
        "badge_color": "#4D4D4D"
    }
]

# Configure styling
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.9

fig, axes = plt.subplots(4, 2, figsize=(20, 21), dpi=300)
plt.subplots_adjust(left=0.045, right=0.89, top=0.94, bottom=0.06, hspace=0.30, wspace=0.10)

cmap = plt.colormaps["RdYlBu_r"].copy()
cmap.set_bad(color="#EBF2F7") # Oceans
norm = mcolors.Normalize(vmin=-30.0, vmax=30.0)

for idx, (ax, info) in enumerate(zip(axes.flatten(), slices_info)):
    t_model = info["time_bp_model"]
    b1 = ds["bio01"].sel(time=t_model).values
    
    # Background temperature field
    ax.set_facecolor("#EBF2F7")
    mesh = ax.pcolormesh(lons, lats, b1, cmap=cmap, norm=norm, shading="auto", zorder=1)
    
    # Coastlines
    coastline.plot(ax=ax, color="#444444", linewidth=0.55, zorder=2)
    
    # Extract 14C dates in this age window
    sub_c14 = df_c14[(df_c14["c14_age"] >= info["age_min"]) & (df_c14["c14_age"] < info["age_max"])]
    
    # Separate into Bone Collagen vs Inorganic / Organic Controls
    is_col = sub_c14["material_category"] == "COLLAGEN"
    is_ctrl = sub_c14["material_category"].isin(["NON_COLLAGEN_ORGANIC_CONTROL", "BONE_APATITE", "SHELL_CONTROL"])
    
    col_dates = sub_c14[is_col]
    ctrl_dates = sub_c14[is_ctrl]
    
    # Plot Controls first (black/dark grey circles with thin border)
    if len(ctrl_dates) > 0:
        ax.scatter(ctrl_dates["longitude"], ctrl_dates["latitude"],
                   c="#1A1A1A", s=8.0, alpha=0.55, marker="o", edgecolors="none", zorder=3,
                   label=f"Inorganic & Controls (Charcoal/Wood/Apatite, N={len(ctrl_dates):,})")
                   
    # Plot Bone Collagen on top (vibrant deep red/magenta diamonds with dark edge)
    if len(col_dates) > 0:
        ax.scatter(col_dates["longitude"], col_dates["latitude"],
                   c="#D9002C", s=14.0, alpha=0.80, marker="D", edgecolors="#40000B", linewidths=0.4, zorder=4,
                   label=f"Bone Collagen (Purified, N={len(col_dates):,})")
    
    # Gridlines and aesthetics
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-60, 86, 30))
    ax.set_xticklabels([r"$180^\circ$", r"$120^\circ\mathrm{W}$", r"$60^\circ\mathrm{W}$", r"$0^\circ$",
                        r"$60^\circ\mathrm{E}$", r"$120^\circ\mathrm{E}$", r"$180^\circ$"], fontsize=8.0, color="#555555")
    ax.set_yticklabels([r"$60^\circ\mathrm{S}$", r"$30^\circ\mathrm{S}$", r"$0^\circ$",
                        r"$30^\circ\mathrm{N}$", r"$60^\circ\mathrm{N}$"], fontsize=8.0, color="#555555")
    ax.grid(True, linestyle=":", color="#777777", alpha=0.40, zorder=2)
    
    # Subplot Header
    ax.set_title(info["title"], fontsize=11.5, fontweight="bold", pad=8, loc="left", color="#111111")
    
    # Badge Box (Top Right)
    ax.text(0.985, 0.94, info["badge"], transform=ax.transAxes,
            ha="right", va="top", fontsize=8.8, fontweight="bold", color="#FFFFFF",
            bbox=dict(boxstyle="round,pad=0.32", facecolor=info["badge_color"], edgecolor="none", alpha=0.92),
            zorder=6)
            
    # Sample Count & Disappearance Metric Box (Bottom Left)
    col_min_lat = col_dates["latitude"].min() if len(col_dates) > 0 else np.nan
    stat_str = f"Bone Col: N={len(col_dates):,} (Min Lat: {col_min_lat:+.1f}°) | Controls: N={len(ctrl_dates):,}" if len(col_dates) > 0 else f"Bone Col: N=0 (Extinct in window) | Controls: N={len(ctrl_dates):,}"
    ax.text(0.015, 0.045, stat_str, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=8.0, fontweight="bold", color="#111111",
            bbox=dict(boxstyle="square,pad=0.25", facecolor="#FFFFFF", edgecolor="#888888", alpha=0.92),
            zorder=6)
            
    # Commentary (Bottom Right)
    ax.text(0.985, 0.045, info["context"], transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.6, style="italic", color="#222222",
            bbox=dict(boxstyle="square,pad=0.25", facecolor="#F8F8F8", edgecolor="#CCCCCC", alpha=0.88),
            zorder=6)
            
    # Mini Legend inside panel (Upper Left)
    ax.legend(loc="upper left", frameon=True, framealpha=0.90, edgecolor="#BBBBBB", fontsize=7.2, markerscale=1.4)

# Shared Colorbar on the right
cbar_ax = fig.add_axes([0.915, 0.12, 0.016, 0.74])
cbar = fig.colorbar(mesh, cax=cbar_ax, orientation="vertical", extend="both")
cbar.set_label("Mean Annual Surface Air Temperature (°C)", fontsize=12.0, fontweight="bold", labelpad=10)
cbar.set_ticks(np.arange(-30, 31, 5))
cbar.ax.tick_params(labelsize=9.0)

# Freezing line marker on colorbar
cbar.ax.axhline(0.0, color="#000000", linewidth=2.0, linestyle="--")
cbar.ax.text(3.5, 0.0, r"$\mathbf{0\,^\circ C}$ (Permafrost / Ice Limit)", transform=cbar.ax.get_yaxis_transform(),
             va="center", ha="left", fontsize=8.5, fontweight="bold", color="#000080")

# Super-Title
fig.suptitle("Global Paleotemperature Evolution & Spatiotemporal Disappearance of Bone Collagen (0 to 50,000 BP)\n"
             "Overlaid Radiocarbon Occurrences: Purified Bone Collagen (Red Diamonds) vs Charcoal/Inorganics (Black Circles)",
             fontsize=15.0, fontweight="bold", y=0.985)

# Save
out_file = os.path.join(output_dir, "Figure_Earth_Paleotemperature_with_C14_Dates_Panel.png")
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.close()
print(f"Successfully generated 14C paleotemperature panel map: {out_file}")

# Mirror to Brain Artifact directory
artifact_file = os.path.join(artifact_dir, "Figure_Earth_Paleotemperature_with_C14_Dates_Panel.png")
shutil.copy2(out_file, artifact_file)
print(f"Mirrored map to artifact directory: {artifact_file}")
