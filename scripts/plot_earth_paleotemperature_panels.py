"""
Paleoclimate Global Temperature Map Series (0 to 50,000 BP)
Capturing Key Climate Transitions & Shifted Epoch Markers:
- 0 BP: Modern / Pre-Industrial Baseline
- 9,000 BP: Early Holocene Optimum (9–8 ka BP)
- 12,000 BP (proxy for 12,300 BP): Younger Dryas Abrupt Cold Reversal (12.9–11.7 ka BP)
- 14,000 BP: Bølling–Allerød Interstadial Rapid Warming (14.7–12.9 ka BP)
- 20,000 BP: Last Glacial Maximum (Peak Global Ice Volume, 26.5–19 ka BP)
- 30,000 BP: Late Marine Isotope Stage 3 (MIS 3) Interstadial Cooling
- 40,000 BP: Middle MIS 3 Interstadial / Anatomically Modern Human Dispersal
- 50,000 BP: Early MIS 3 / Radiocarbon Methodological Limit (~50 ka BP)
"""

import os
import shutil
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Directories
output_dir = r"D:\26 Modelling Collagen Hydrolysis\outputs"
artifact_dir = r"C:\Users\matth\.gemini\antigravity-ide\brain\cd905880-511d-460a-944c-c857596f9afc"
nc_path = r"D:\26 Modelling Collagen Hydrolysis\pastclim\Beyer2020_annual_vars_v1.2.2.nc"
coastline_path = r"D:\26 Modelling Collagen Hydrolysis\pastclim\coastline\ne_110m_coastline.shp"

# Load NetCDF and Coastline
ds = xr.open_dataset(nc_path, decode_times=False)
coastline = gpd.read_file(coastline_path)

# Map grid definition
lats = ds["latitude"].values
lons = ds["longitude"].values

# Shifted panels definition:
# Title, Time BP in NetCDF, Subtitle/Context
slices_info = [
    {
        "panel": "A",
        "time_bp": 0,
        "title": "0 ka BP: Late Holocene / Pre-Industrial Baseline",
        "badge": "0 BP",
        "context": "Modern baseline climate; global mean terrestrial T ≈ 9.3 °C",
        "color": "#1B7837"
    },
    {
        "panel": "B",
        "time_bp": -9000,
        "title": "9 ka BP: Early Holocene Climate Optimum",
        "badge": "9,000 BP (7,050 BCE)",
        "context": "Post-glacial thermal maximum; enhanced African/Asian monsoons",
        "color": "#313695"
    },
    {
        "panel": "C",
        "time_bp": -12000,
        "title": "12 ka BP: Younger Dryas (YD) Cold Reversal",
        "badge": "12,300 BP (10,350 BCE)",
        "context": "Abrupt 1,200-yr return to near-glacial cooling across N. Hemisphere",
        "color": "#053061"
    },
    {
        "panel": "D",
        "time_bp": -14000,
        "title": "14 ka BP: Bølling–Allerød Warm Interstadial",
        "badge": "14,000 BP (12,050 BCE)",
        "context": "Rapid deglacial warming pulse following Heinrich Event 1",
        "color": "#B2182B"
    },
    {
        "panel": "E",
        "time_bp": -20000,
        "title": "20 ka BP: Last Glacial Maximum (LGM)",
        "badge": "20,000 BP (18,050 BCE)",
        "context": "Peak global ice volume, sea level -125 m; global cooling ~4–6 °C",
        "color": "#2166AC"
    },
    {
        "panel": "F",
        "time_bp": -30000,
        "title": "30 ka BP: Late Marine Isotope Stage 3 (MIS 3)",
        "badge": "30,000 BP",
        "context": "Cool transitional stadial prior to LGM ice-sheet expansion",
        "color": "#4393C3"
    },
    {
        "panel": "G",
        "time_bp": -40000,
        "title": "40 ka BP: Mid MIS 3 Interstadial / AMH Dispersal",
        "badge": "40,000 BP",
        "context": "Variable climate, Neanderthal-AMH transition in Eurasia",
        "color": "#762A83"
    },
    {
        "panel": "H",
        "time_bp": -50000,
        "title": "50 ka BP: Early MIS 3 / 14C Methodological Horizon",
        "badge": "50,000 BP",
        "context": "Pre-LGM interstadial warmth; absolute instrumental radiocarbon limit",
        "color": "#4D4D4D"
    }
]

# Style settings
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.9

# Create Figure: 4 rows x 2 cols
fig, axes = plt.subplots(4, 2, figsize=(20, 21), dpi=300)
plt.subplots_adjust(left=0.045, right=0.89, top=0.94, bottom=0.06, hspace=0.30, wspace=0.10)

# Colormap for Temperature: Smooth diverging / scientific thermal palette
cmap = plt.colormaps["RdYlBu_r"].copy()
cmap.set_bad(color="#EBF2F7") # Light pale blue/grey for oceans
norm = mcolors.Normalize(vmin=-30.0, vmax=30.0)

for idx, (ax, info) in enumerate(zip(axes.flatten(), slices_info)):
    t_val = info["time_bp"]
    b1 = ds["bio01"].sel(time=t_val).values
    
    # Plot temperature field
    ax.set_facecolor("#EBF2F7") # Oceans
    mesh = ax.pcolormesh(lons, lats, b1, cmap=cmap, norm=norm, shading="auto", zorder=1)
    
    # Overlay global coastlines
    coastline.plot(ax=ax, color="#222222", linewidth=0.55, zorder=3)
    
    # Gridlines and aesthetics
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-60, 86, 30))
    ax.set_xticklabels([r"$180^\circ$", r"$120^\circ\mathrm{W}$", r"$60^\circ\mathrm{W}$", r"$0^\circ$",
                        r"$60^\circ\mathrm{E}$", r"$120^\circ\mathrm{E}$", r"$180^\circ$"], fontsize=8.2, color="#444444")
    ax.set_yticklabels([r"$60^\circ\mathrm{S}$", r"$30^\circ\mathrm{S}$", r"$0^\circ$",
                        r"$30^\circ\mathrm{N}$", r"$60^\circ\mathrm{N}$"], fontsize=8.2, color="#444444")
    ax.grid(True, linestyle=":", color="#777777", alpha=0.45, zorder=2)
    
    # Subplot Header Banner
    ax.set_title(f"{info['panel']}. {info['title']}", fontsize=11.5, fontweight="bold", pad=8, loc="left", color="#111111")
    
    # Event Badge Box (Top Right inside axis)
    ax.text(0.985, 0.94, info["badge"], transform=ax.transAxes,
            ha="right", va="top", fontsize=9.0, fontweight="bold", color="#FFFFFF",
            bbox=dict(boxstyle="round,pad=0.35", facecolor=info["color"], edgecolor="none", alpha=0.92),
            zorder=5)
    
    # Terrestrial Summary Stats (Bottom Left inside axis)
    valid_t = b1[~np.isnan(b1)]
    mean_t = np.mean(valid_t)
    min_t = np.min(valid_t)
    max_t = np.max(valid_t)
    stat_str = f"Terrestrial Mean: {mean_t:+.1f} °C  (Min: {min_t:+.1f} °C | Max: {max_t:+.1f} °C)"
    ax.text(0.015, 0.045, stat_str, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=8.0, fontweight="bold", color="#111111",
            bbox=dict(boxstyle="square,pad=0.25", facecolor="#FFFFFF", edgecolor="#888888", alpha=0.88),
            zorder=5)
            
    # Context annotation (Bottom Right inside axis)
    ax.text(0.985, 0.045, info["context"], transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.8, style="italic", color="#222222",
            bbox=dict(boxstyle="square,pad=0.25", facecolor="#F5F5F5", edgecolor="#CCCCCC", alpha=0.85),
            zorder=5)

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

# Figure Super-Title
fig.suptitle("Global Paleotemperature Evolution: 0 to 50,000 BP (Beyer et al. 2020 HadCM3 Snapshot Series)\n"
             "Tracking Key Climatic Horizons: Holocene Optimum, Younger Dryas, Bølling–Allerød, LGM & MIS 3",
             fontsize=15.5, fontweight="bold", y=0.985)

# Save figure
out_file = os.path.join(output_dir, "Figure_Earth_Paleotemperature_0_to_50ka_Panel.png")
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.close()
print(f"Successfully generated paleotemperature panel map: {out_file}")

# Mirror to Brain Artifact directory
artifact_file = os.path.join(artifact_dir, "Figure_Earth_Paleotemperature_0_to_50ka_Panel.png")
shutil.copy2(out_file, artifact_file)
print(f"Mirrored map to artifact directory: {artifact_file}")
