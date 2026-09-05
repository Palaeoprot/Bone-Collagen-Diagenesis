# Bone Collagen Diagenesis

**Date & Time:** 2026-09-05 10:32:00 (+02:00)

Repository companion for:
**"An Estimate of Bone Collagen Hydrolysis Rates over Deep Time using Experimental, Palaeoproteomic and Radiocarbon data"**

Remote: [https://github.com/Palaeoprot/Bone-Collagen-Diagenesis](https://github.com/Palaeoprot/Bone-Collagen-Diagenesis)

---

## Overview

For two decades, bone collagen survival in archaeological and paleontological contexts has been modelled using an Arrhenius activation energy of $E_a = 173\text{ kJ/mol}$ derived from high-temperature laboratory heating experiments ($55\text{--}95^\circ\text{C}$, Smith 2002; Collins et al. 2002). However, this model produces an "Arctic Paradox," predicting that collagen should routinely survive for $>50\text{--}70\text{ Ma}$ in permafrost regimes (such as Ellesmere Island, $-10.5^\circ\text{C}$), when empirical survival on Earth never exceeds $\sim 4\text{--}7\text{ Ma}$.

This repository contains all code, harmonised radiocarbon datasets, paleoclimate integrations, and publication figures demonstrating that:

1. **The Radiocarbon Wall Artefact:** Instrumental radiocarbon blank censoring ($\sim 42\text{--}45\text{ ka BP}$) imposes a flat horizontal cutoff across all paleoclimates with mean integrated temperature $\bar{T} < 13.8^\circ\text{C}$, mathematically suppressing apparent thermal ages toward zero in cold environments.
2. **True Empirical Kinetics:** Inverting the uncensored kinetic regime ($\bar{T} \ge 13.8^\circ\text{C}$) yields a true empirical activation energy of **$E_a = 130.8\text{--}133.4\text{ kJ/mol}$** ($R^2 = 0.88$).
3. **Deep-Time Validation:** Extrapolating this empirical envelope simultaneously unifies modern radiocarbon boundaries, Dmanisi rhino dentin collagen ($1.77\text{ Ma}$ at $+11^\circ\text{C}$), Yukon horse ($735\text{ ka}$ at $-9^\circ\text{C}$), and the Ellesmere Bear / High Arctic Camel ceiling ($3.4\text{--}3.9\text{ Ma}$ at $-10.5^\circ\text{C}$).
4. **Spatiotemporal Disappearance:** Beyond $14\text{ ka BP}$ (post-Bølling–Allerød to LGM and MIS 3), bone collagen completely vanishes from all tropical and subtropical regions south of $+30^\circ\text{N}$, while charcoal and inorganic dates persist globally.

---

## Repository Structure

```text
Bone-Collagen-Diagenesis/
│
├── README.md                                          <- Project overview and reproduction guide
├── docs/
│   └── PROJECT_HANDOVER.md                            <- Comprehensive technical handover document
│
├── data/
│   ├── harmonized_c14_master.parquet                  <- 185,189 harmonized global 14C determinations (4.68 MB)
│   ├── collagen_vs_control_thermal_cohort.parquet     <- 36,202 paired collagen vs organic control records (2.30 MB)
│   └── site_paleotemperature_series_0_to_50ka.parquet <- HadCM3 paleotemperature time series per site (0.81 MB)
│
├── scripts/
│   ├── paleoclimate_engine.py                         <- HadCM3 paleoclimate emulator and integration class
│   ├── harmonize_radiocarbon.py                       <- Raw 14C database ingestion & taxonomy classification
│   ├── run_thermal_integration.py                     <- Thermal age & effective temperature numerical integrator
│   ├── fit_quantile_envelope.py                       <- Quantile regression on 14C upper boundary
│   ├── scan_quantiles.py                              <- Quantile sweep demonstrating horizontal clipping
│   ├── zero_intercept_fit.py                          <- Arrhenius regression on uncensored regime (T >= 13.8 °C)
│   ├── plot_comprehensive_deconvolution_4panel.py      <- Generates Figure 1 (4-Panel Publication Figure)
│   ├── plot_earth_paleotemperature_with_c14_panels.py <- Generates Figure 2 (8-Panel Paleomap Atlas)
│   ├── plot_arrhenius_smith_experimental.py           <- Lab vs field Arrhenius rate comparison
│   ├── plot_c14_wall_deconvolution.py                 <- Radiocarbon wall crossover regression plot
│   └── plot_panel_d_standalone.py                     <- Standalone deep-time extrapolation figure
│
└── figures/
    ├── Figure_Comprehensive_4Panel_Radiocarbon_Deconvolution.png
    ├── Figure_Earth_Paleotemperature_with_C14_Dates_Panel.png
    ├── Figure_Arrhenius_Smith_High_Temp_vs_Geological_Field.png
    └── Figure_Panel_D_True_Kinetic_Extrapolation_Deep_Time.png

```

---

## Data Summary

All datasets are compressed in columnar Apache Parquet format for fast loading and full reproducibility:

1. **`harmonized_c14_master.parquet`** ($N = 185,189$):
   - Standardized columns: `sample_id`, `lab_id`, `site_name`, `country`, `latitude`, `longitude`, `c14_age`, `c14_error`, `material_raw`, `material_category`, `cn_ratio`, `collagen_yield_pct`, `taxa`, `source_db`.
   - Material categories: `COLLAGEN`, `NON_COLLAGEN_ORGANIC_CONTROL` (charcoal, wood, seeds), `BONE_APATITE`, `SHELL_CONTROL`, `BONE_UNDIFFERENTIATED`, `SEDIMENT`, `OTHER`.
2. **`collagen_vs_control_thermal_cohort.parquet`** ($N = 36,202$):
   - Exactly paired cohort ($18,101$ purified bone collagen vs $18,101$ non-collagen organic controls) with paleoclimate-integrated mean effective temperatures ($\bar{T}$) and equivalent thermal ages at $10^\circ\mathrm{C}$ reference under $E_a = 173\ \mathrm{kJ/mol}$ and $E_a = 100\ \mathrm{kJ/mol}$.
3. **`site_paleotemperature_series_0_to_50ka.parquet`** ($N = 406,038$):
   - Continuous site-specific paleotemperature time series from the HadCM3 climate emulator (Beyer et al. 2020) across all 37 time steps from $0$ to $50,000\ \mathrm{BP}$ for all 10,974 unique sample coordinates.

---

## How to Reproduce

### Prerequisites
Install dependencies:
```bash
pip install numpy pandas xarray geopandas matplotlib pyarrow fastparquet
```

### Reproducing Figures
To reproduce the primary publication figures:
```bash
# Generate Figure 1: Comprehensive 4-Panel Deconvolution & Arrhenius Rates
python scripts/plot_comprehensive_deconvolution_4panel.py

# Generate Figure 2: 8-Panel Earth Paleotemperature & C14 Disappearance Atlas
python scripts/plot_earth_paleotemperature_with_c14_panels.py
```

---

## Citation & Authorship
- **Matthew Collins**: [https://github.com/Palaeoprot](https://github.com/Palaeoprot)
- **Repository**: [https://github.com/Palaeoprot/Bone-Collagen-Diagenesis](https://github.com/Palaeoprot/Bone-Collagen-Diagenesis)
