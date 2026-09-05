# -*- coding: utf-8 -*-
r"""
Determination-level estimator for the activation energy of bone collagen
hydrolysis, using the matched controls as an availability denominator.

WHY NOT THE BINNED FIT
----------------------
The published estimate regresses ln(max 14C age) per 1 C temperature bin against
1/T. That estimator (a) reduces 18,101 determinations to 8 points, (b) uses an
extreme order statistic whose expectation rises with bin size while bin size
falls with temperature, and (c) discards the censored cold regime entirely.
See sensitivity_charcoal_only_report.md.

THE MODEL
---------
Let S(a, T) be the probability that a bone of age `a` buried at effective
temperature `T` still yields a datable collagen extract. First-order kinetics
with a 1% yield criterion give a limiting age

    t_lim(T) = (4.60517 / A) * exp(Ea / (R T))
    ln t_lim(T) = ln(4.60517 / A) + (Ea / R) * (1 / T)

Model survival as a soft threshold in log-age about that limit:

    S(a, T) = sigmoid( g * [ ln t_lim(T) - ln a ] )

where g > 0 controls how sharply preservation fails at the limit.

The matched controls (charcoal, wood, macrofossils) are dated at the same sites
and horizons but do not depend on collagen survival, so they estimate the
archaeological availability of datable material at (a, T). Pooling collagen
(y = 1) and control (y = 0) determinations, the odds of a record being collagen
are proportional to S, giving

    logit P(y = 1 | a, T) = c + b * (1000/T) - g * ln a

    with   b = g * Ea / (1000 R)      =>      Ea = 1000 * R * b / g

This is an ordinary logistic regression whose coefficient RATIO is the
activation energy. It uses every determination, needs no binning, and is
identified by the contrast between collagen and control at the same site.

CENSORING
---------
Both materials are clipped by the AMS blank, so the analysis is restricted to
ages below the ceiling, where both are observable. Cold sites are then genuinely
informative (collagen share stays high to 42 ka, bounding t_lim(T) from below)
rather than being discarded as in the binned fit.

IDENTIFICATION CAVEAT
---------------------
The estimator assumes material choice depends on (a, T) only through collagen
survival. If practice independently favours charcoal at older sites, `g` is
inflated and Ea biased downward. Model 3 adds region and source-database fixed
effects to absorb the largest part of that; the residual is discussed in the
report.
"""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm

R = 8.31446
RNG = np.random.default_rng(20260905)
AMS_CEILING = 42000.0
MIN_AGE = 500.0

_DATA = [
    r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet",
    r"C:\Users\matth\Documents\GitHub\Bone-Collagen-Diagenesis\data\collagen_vs_control_thermal_cohort.parquet",
]
HERE = os.path.dirname(os.path.abspath(__file__))
path = next((p for p in _DATA if os.path.exists(p)), None)
if path is None:
    raise SystemExit("cohort parquet not found")

df = pd.read_parquet(path).dropna(subset=["integrated_temp_c", "c14_age", "site_name"])
df = df[(df.c14_age >= MIN_AGE) & (df.c14_age < AMS_CEILING)].copy()
df["y"] = (df.material_category == "COLLAGEN").astype(int)
df["z"] = 1000.0 / (df.integrated_temp_c + 273.15)
df["ln_a"] = np.log(df.c14_age)

# strict-charcoal flag for the sensitivity variant
head = df["material_raw"].fillna("").str.lower().str.split("|").str[0].str.strip()
df["strict_charcoal"] = head.str.contains(r"charcoal|charbon de bois|charred|nutshell") & \
                        ~head.str.contains(r"ceramic|residue|textile|peat|sediment|soil")

# primary analysis: sites where BOTH materials are present (a true matched contrast)
both = df.groupby("site_name")["y"].agg(["min", "max"])
paired_sites = both[(both["min"] == 0) & (both["max"] == 1)].index
dfp = df[df.site_name.isin(paired_sites)].copy()


def fit_ea(d, extra=None, return_model=False):
    """Logistic fit; Ea = 1000 R b/g with delta-method SE."""
    X = pd.DataFrame({"z": d.z, "ln_a": d.ln_a})
    if extra:
        for col in extra:
            dm = pd.get_dummies(d[col].astype(str), prefix=col, drop_first=True, dtype=float)
            keep = dm.columns[dm.sum() >= 30]          # drop ultra-rare levels
            X = pd.concat([X, dm[keep]], axis=1)
    X = sm.add_constant(X, has_constant="add")
    try:
        m = sm.Logit(d.y.values, X.values).fit(disp=0, maxiter=200)
    except Exception:
        return (np.nan,) * 3 + (0,) if not return_model else (np.nan,)*3 + (0, None)
    names = list(X.columns)
    ib, ig = names.index("z"), names.index("ln_a")
    b, gneg = m.params[ib], m.params[ig]
    g = -gneg                                          # coefficient on ln_a is -g
    if g <= 0:
        return (np.nan,)*3 + (len(d),) if not return_model else (np.nan,)*3 + (len(d), m)
    ea = 1000.0 * R * b / g / 1000.0                   # kJ/mol
    V = m.cov_params()
    var = (V[ib, ib] / g**2) + (b**2 * V[ig, ig] / g**4) + (2 * b * V[ib, ig] / g**3)
    se = 1000.0 * R * np.sqrt(max(var, 0)) / 1000.0
    out = (ea, se, g, len(d))
    return out + (m,) if return_model else out


def site_bootstrap(d, B=300, extra=None):
    sites = d.site_name.unique()
    by = {s: g for s, g in d.groupby("site_name")}
    vals = []
    for _ in range(B):
        pick = RNG.choice(sites, size=len(sites), replace=True)
        boot = pd.concat([by[s] for s in pick], ignore_index=True)
        ea = fit_ea(boot, extra=extra)[0]
        if np.isfinite(ea):
            vals.append(ea)
    a = np.array(vals)
    return np.nanpercentile(a, 2.5), np.nanpercentile(a, 97.5), np.nanstd(a), len(a)


print("=" * 76)
print("CENSORED / MATCHED-CONTROL LOGISTIC ESTIMATOR FOR Ea")
print("=" * 76)
print(f"records after age filter [{MIN_AGE:.0f}, {AMS_CEILING:.0f}) : {len(df):,}")
print(f"paired sites (both materials present)          : {len(paired_sites):,}")
print(f"records at paired sites                        : {len(dfp):,} "
      f"({dfp.y.sum():,} collagen / {(1-dfp.y).sum():,} control)")
print()

models = {}
rows = []

specs = [
    ("1. All records",                 df,  None),
    ("2. Paired sites only (primary)",  dfp, None),
    ("3. Paired + region/db effects",   dfp, ["country", "source_db"]),
    ("4. Paired, strict charcoal ctrl", dfp[dfp.y.eq(1) | dfp.strict_charcoal], None),
]
for label, d, extra in specs:
    ea, se, g, nn = fit_ea(d, extra=extra)
    lo, hi, bse, nb = site_bootstrap(d, B=200, extra=extra)
    rows.append((label, ea, se, lo, hi, g, nn))
    print(f"{label:34s} Ea = {ea:6.1f}  (delta SE {se:4.1f})   "
          f"boot 95% CI [{lo:6.1f}, {hi:6.1f}]   g={g:.2f}  N={nn:,}")

print()
print("Reference: powdered-bone lab value 173.2 ; intact-bone lab value (Ortner 1972) 132.1")
print()

# --- implied preservation limit at a few temperatures, from the primary model ---
ea_p, se_p, g_p, n_p, m_p = fit_ea(dfp, return_model=True)
names = ["const", "z", "ln_a"]
c, b = m_p.params[0], m_p.params[1]
# ln t_lim(T) = (c + b z)/g  at the 50% survival point
print("Implied 50% collagen-survival age from the primary model:")
print(f"{'T (C)':>7} {'t_50 (yr)':>12}")
for T in [-10.0, 0.0, 5.0, 10.0, 14.0, 18.0, 22.0, 25.0]:
    z = 1000.0 / (T + 273.15)
    ln_t = (c + b * z) / g_p
    print(f"{T:7.1f} {np.exp(ln_t):12,.0f}")

md = f"""# Determination-Level Estimator for $E_a$ (Censored / Matched-Control Logistic)
**Date & Time:** 2026-09-05 (generated by `censored_regression_ea.py`)

## 1. Why this estimator

The binned boundary fit reduces $18{{,}}101$ determinations to eight points, uses an
extreme order statistic confounded with bin size, and discards the censored cold
regime (see `sensitivity_charcoal_only_report.md`). This estimator uses every
determination and turns the censoring into information.

## 2. Model

With $S(a,T)$ the probability that a bone of age $a$ at effective temperature $T$
still yields datable collagen, first-order kinetics at a $1\\%$ yield criterion give
$\\ln t_{{\\lim}}(T) = \\ln(4.60517/A) + (E_a/R)(1/T)$. Modelling survival as a soft
threshold, $S = \\mathrm{{sigmoid}}(g[\\ln t_{{\\lim}}(T) - \\ln a])$, and using the matched
controls as the availability denominator, the pooled data satisfy

$$\\mathrm{{logit}}\\,P(\\text{{collagen}} \\mid a,T) = c + b\\,(1000/T) - g \\ln a,
\\qquad E_a = 1000\\,R\\,b/g$$

an ordinary logistic regression whose **coefficient ratio** is the activation
energy. Restricting to ages below the AMS ceiling makes cold sites informative
(a sustained high collagen share to $42$ ka bounds $t_{{\\lim}}$ from below) instead
of discarding them.

## 3. Cohort

| Quantity | N |
| :--- | ---: |
| Records in $[{MIN_AGE:.0f}, {AMS_CEILING:.0f})$ yr | {len(df):,} |
| Sites with **both** materials present | {len(paired_sites):,} |
| Records at paired sites | {len(dfp):,} |
| — collagen / control | {dfp.y.sum():,} / {(1-dfp.y).sum():,} |

## 4. Results

| Specification | $E_a$ (kJ mol$^{{-1}}$) | delta SE | site-clustered 95% CI | $g$ | N |
| :--- | ---: | ---: | :--- | ---: | ---: |
""" + "\n".join(
    f"| {r[0][3:]} | **{r[1]:.1f}** | {r[2]:.1f} | [{r[3]:.1f}, {r[4]:.1f}] | {r[5]:.2f} | {r[6]:,} |"
    for r in rows) + f"""

Reference values: powdered-bone laboratory $173.2$; intact-bone laboratory
(Ortner et al. 1972) $132.1$ kJ mol$^{{-1}}$.

## 5. Implied preservation limit (primary model)

50% collagen-survival age as a function of effective temperature:

| $T$ (°C) | $t_{{50}}$ (yr) |
| ---: | ---: |
""" + "\n".join(
    f"| {T:.1f} | {np.exp((c + b*(1000.0/(T+273.15)))/g_p):,.0f} |"
    for T in [-10.0, 0.0, 5.0, 10.0, 14.0, 18.0, 22.0, 25.0]) + """

## 6. Identification caveat

The estimator assumes material choice depends on age and temperature only through
collagen survival. If excavators independently prefer charcoal at older sites,
$g$ is inflated and $E_a$ biased downward. Specification 3 absorbs region and
source-database effects; the residual confound cannot be removed with these data
and is the main reason to treat the interval, not the point estimate, as the result.
"""
out = os.path.join(HERE, "censored_regression_ea_report.md")
with open(out, "w", encoding="utf-8") as fh:
    fh.write(md)
print(f"\nWrote {out}")
