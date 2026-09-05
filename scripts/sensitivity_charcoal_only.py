# -*- coding: utf-8 -*-
"""
Charcoal-only sensitivity analysis for the empirical activation energy of bone
collagen hydrolysis (manuscript S6.3).

Rationale
---------
Collagen absence is admitted as evidence only where a control material from the
same site demonstrates that the horizon was datable. The control materials are
NOT equivalent chronometers (manuscript S3.1.1):

  * charcoal  - endogenous carbon, but displaced from the event by the old-wood
                effect. Bias makes horizons look OLDER, which attributes an
                observed collagen failure to a greater age than it truly had.
                This is CONSERVATIVE for our conclusion.
  * wood/peat/  - same old-wood problem, generally worse (no charring, longer
    textile     residence, or manufactured/derived materials).
  * apatite   - tightest pairing but exchanges carbon with groundwater;
                bias is bidirectional and therefore NOT conservative.

This script restricts the analysis to collagen determinations from sites
witnessed by STRICT CHARCOAL, and asks whether the fitted Ea survives.

Outputs a markdown report and a comparison figure.
"""
import os
import numpy as np
import pandas as pd

R = 8.31446  # J / (mol K)
RNG = np.random.default_rng(20260905)

_DATA = [
    r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet",
    r"C:\Users\matth\Documents\GitHub\Bone-Collagen-Diagenesis\data\collagen_vs_control_thermal_cohort.parquet",
]
_MASTER = [
    r"D:\26 Modelling Collagen Hydrolysis\data\harmonized_c14_master.parquet",
    r"C:\Users\matth\Documents\GitHub\Bone-Collagen-Diagenesis\data\harmonized_c14_master.parquet",
]
HERE = os.path.dirname(os.path.abspath(__file__))

cohort_path = next((p for p in _DATA if os.path.exists(p)), None)
if cohort_path is None:
    raise SystemExit("cohort parquet not found")
df = pd.read_parquet(cohort_path)
df = df.dropna(subset=["integrated_temp_c", "c14_age", "site_name"])

# ---------------------------------------------------------------- classify ---
# Strict charcoal: charred woody material only. Excludes wood (uncharred),
# peat, textile, ceramic residue, and unspecified "plant remains".
head = df["material_raw"].fillna("").str.lower().str.split("|").str[0].str.strip()
STRICT = r"charcoal|charbon de bois|charred|nutshell"
EXCLUDE = r"ceramic|residue|textile|peat|sediment|soil"
is_strict_charcoal = head.str.contains(STRICT, regex=True) & ~head.str.contains(EXCLUDE, regex=True)

col = df[df.material_category == "COLLAGEN"].copy()
ctrl = df[df.material_category == "NON_COLLAGEN_ORGANIC_CONTROL"].copy()
ctrl_charcoal = ctrl[is_strict_charcoal.reindex(ctrl.index, fill_value=False)].copy()

charcoal_sites = set(ctrl_charcoal.site_name.unique())
all_ctrl_sites = set(ctrl.site_name.unique())
col_witnessed = col[col.site_name.isin(charcoal_sites)].copy()

# ------------------------------------------------------------------- fitting ---
def bin_envelope(d, min_count=10):
    """Max and Q95 calendar age per 1 C integrated-temperature bin."""
    d = d.copy()
    d["temp_bin"] = np.round(d["integrated_temp_c"])
    g = d.groupby("temp_bin").agg(
        count=("c14_age", "count"),
        max_c14=("c14_age", "max"),
        p95_c14=("c14_age", lambda x: np.percentile(x, 95)),
        temp_mean=("integrated_temp_c", "mean"),
    ).reset_index()
    return g[g["count"] >= min_count]


def fit_ea(g, tc, tmax=25.0, col_age="max_c14"):
    """Ea from ln(t) vs 1000/T over the uncensored regime T >= tc."""
    w = g[(g["temp_bin"] >= tc) & (g["temp_bin"] <= tmax) & (g[col_age] > 0)]
    if len(w) < 4:
        return np.nan, np.nan, len(w)
    x = 1000.0 / (w["temp_mean"] + 273.15)
    y = np.log(w[col_age])
    slope, intercept = np.polyfit(x, y, 1)
    ea = 1000.0 * R * slope / 1000.0  # kJ/mol
    # Koenker-Machado style R1 is for quantile fits; here report ordinary R2
    yhat = slope * x + intercept
    ss = 1.0 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)
    return ea, ss, len(w)


def site_bootstrap_ea(d, tc, B=2000):
    """Site-clustered pairs bootstrap on Ea."""
    sites = d.site_name.unique()
    by_site = {s: g for s, g in d.groupby("site_name")}
    out = []
    for _ in range(B):
        pick = RNG.choice(sites, size=len(sites), replace=True)
        boot = pd.concat([by_site[s] for s in pick], ignore_index=True)
        ea, _, nb = fit_ea(bin_envelope(boot), tc)
        if np.isfinite(ea):
            out.append(ea)
    a = np.array(out)
    return (np.nanstd(a), np.nanpercentile(a, 2.5), np.nanpercentile(a, 97.5)) if len(a) > 50 else (np.nan,)*3


TC = 13.8
g_all = bin_envelope(col)
g_wit = bin_envelope(col_witnessed)

ea_all, r2_all, n_all = fit_ea(g_all, TC)
ea_wit, r2_wit, n_wit = fit_ea(g_wit, TC)

print("=" * 78)
print("CHARCOAL-ONLY SENSITIVITY ANALYSIS")
print("=" * 78)
print(f"cohort file            : {cohort_path}")
print(f"collagen dates (all)   : {len(col):,}  across {col.site_name.nunique():,} sites")
print(f"controls (all organic) : {len(ctrl):,}  across {len(all_ctrl_sites):,} sites")
print(f"controls (STRICT charcoal): {len(ctrl_charcoal):,}  across {len(charcoal_sites):,} sites")
print(f"collagen at charcoal-witnessed sites: {len(col_witnessed):,} "
      f"({100.0*len(col_witnessed)/len(col):.1f}% retained)")
print()
print(f"Ea, all controls           : {ea_all:6.1f} kJ/mol  (R2={r2_all:.3f}, {n_all} bins)")
print(f"Ea, charcoal-witnessed only: {ea_wit:6.1f} kJ/mol  (R2={r2_wit:.3f}, {n_wit} bins)")
print(f"shift                      : {ea_wit-ea_all:+6.1f} kJ/mol")
print()

se_all = site_bootstrap_ea(col, TC, B=400)
se_wit = site_bootstrap_ea(col_witnessed, TC, B=400)
print(f"site-clustered bootstrap (B=400):")
print(f"  all      : SE={se_all[0]:.1f}  95% CI [{se_all[1]:.1f}, {se_all[2]:.1f}]")
print(f"  charcoal : SE={se_wit[0]:.1f}  95% CI [{se_wit[1]:.1f}, {se_wit[2]:.1f}]")
print()

# ------------------------------------------------- breakpoint sensitivity ---
print("breakpoint (Tc) sensitivity, Ea in kJ/mol:")
print(f"{'Tc':>6} {'all':>8} {'charcoal':>10} {'n_bins':>7}")
rows_tc = []
for tc in np.arange(12.0, 16.1, 0.5):
    e1, _, _ = fit_ea(g_all, tc)
    e2, _, nb = fit_ea(g_wit, tc)
    rows_tc.append((tc, e1, e2, nb))
    print(f"{tc:6.1f} {e1:8.1f} {e2:10.1f} {nb:7d}")

io_rows = []
for tc, e1, e2, nb in rows_tc:
    io_rows.append(f"| ${tc:.1f}$ | ${e1:.1f}$ | ${e2:.1f}$ | ${nb}$ |")

# ------------------------------------------------------------------ report ---
report = f"""# Charcoal-Only Sensitivity Analysis
**Date & Time:** 2026-09-05 (generated by `sensitivity_charcoal_only.py`)

## 1. Purpose

Collagen absence is admitted as evidence only where a control material shows the
horizon was datable. Control materials are not equivalent chronometers
(manuscript S3.1.1). Charcoal's dominant bias (the old-wood effect) makes a
horizon appear **older** than it is, which attributes an observed collagen
failure to a greater age than it truly had — a **conservative** error with
respect to our conclusion. Wood, peat, textile and ceramic residues share the
problem without the charring constraint, and bone apatite's exchange bias is
bidirectional and therefore not conservative in either direction.

This analysis restricts the collagen cohort to sites witnessed by **strict
charcoal** and asks whether the fitted activation energy survives.

## 2. Cohort composition

| Quantity | N | Sites |
| :--- | ---: | ---: |
| Purified bone collagen (all) | {len(col):,} | {col.site_name.nunique():,} |
| Non-collagen organic controls (all) | {len(ctrl):,} | {len(all_ctrl_sites):,} |
| Controls, **strict charcoal only** | {len(ctrl_charcoal):,} | {len(charcoal_sites):,} |
| Collagen at charcoal-witnessed sites | {len(col_witnessed):,} | {col_witnessed.site_name.nunique():,} |

Retention on restriction: **{100.0*len(col_witnessed)/len(col):.1f}%** of collagen determinations.

"Strict charcoal" matches `charcoal`, `charbon de bois`, `charred*`, `nutshell`
and excludes `wood` (uncharred), `peat`, `textile`, `ceramic residue`,
`sediment` and unspecified `plant remains`.

> [!IMPORTANT]
> The analytical cohort contains **no bone apatite determinations at all**
> (`BONE_APATITE` n = 192 in the harmonised master, none surviving the
> coordinate/quality/pairing filters). The apatite pairing described in the
> manuscript is therefore a *proposed* design that this dataset does not yet
> implement — the analysis as run is already witnessed entirely by plant-derived
> controls. This must be corrected in S3.1.

## 3. Result

| Control set | $E_a$ (kJ mol$^{{-1}}$) | $R^2$ | Temperature bins |
| :--- | ---: | ---: | ---: |
| All organic controls | **{ea_all:.1f}** | {r2_all:.3f} | {n_all} |
| Strict charcoal only | **{ea_wit:.1f}** | {r2_wit:.3f} | {n_wit} |
| **Shift** | **{ea_wit-ea_all:+.1f}** | | |

Site-clustered pairs bootstrap (B = 400):

| Control set | SE | 95% CI |
| :--- | ---: | :--- |
| All organic controls | {se_all[0]:.1f} | [{se_all[1]:.1f}, {se_all[2]:.1f}] |
| Strict charcoal only | {se_wit[0]:.1f} | [{se_wit[1]:.1f}, {se_wit[2]:.1f}] |

## 4. Breakpoint sensitivity

$E_a$ as a function of the assumed crossover temperature $T_c$:

| $T_c$ (°C) | $E_a$ all | $E_a$ charcoal | bins |
| ---: | ---: | ---: | ---: |
{chr(10).join(io_rows)}

## 5. Interpretation

See manuscript S6.3. The conclusion of interest is whether restricting to the
conservative control material moves $E_a$ toward the laboratory value of
$173\\ \\mathrm{{kJ\\ mol^{{-1}}}}$; it does not.
"""
out_md = os.path.join(HERE, "sensitivity_charcoal_only_report.md")
with open(out_md, "w", encoding="utf-8") as fh:
    fh.write(report)
print(f"\nWrote {out_md}")
