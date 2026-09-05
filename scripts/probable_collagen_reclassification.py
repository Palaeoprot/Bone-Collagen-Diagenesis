# -*- coding: utf-8 -*-
r"""
Reclassification test: treat undifferentiated bone as PROBABLE COLLAGEN.

RATIONALE (M. Collins, 2026-09-05)
----------------------------------
Two observations motivate this test.

1. Unless a record explicitly states cremated / calcined / burnt, a bone
   radiocarbon date is *probably* a collagen date. Routine practice is to
   extract collagen; the `BONE_UNDIFFERENTIATED` class in our harmonisation is
   largely a metadata-reporting gap, not a chemistry gap. It holds 22,198
   records -- larger than the explicit `COLLAGEN` class itself.

2. Bone apatite dating is rare because it is disliked, technically difficult
   and controversial -- not because our harmonisation failed to find it. Its
   near-absence (192 records, ~70 genuinely bone) is therefore a fact about
   disciplinary practice. Apatite cannot supply the specimen-level control the
   design would ideally use, and no amount of re-parsing will change that.

Cremated/calcined bone IS excluded: those dates are made on the carbonate of
the calcined mineral, not on collagen, and behave as a different material.

The test asks how much the activation energy estimate moves when the probable
collagen dates are admitted -- and, critically, whether they populate the
21-24 C temperature gap that currently leaves the binned fit resting on a
single 5-site bin at 25 C.
"""
import os, time
import numpy as np
import pandas as pd
import statsmodels.api as sm
from paleoclimate_engine import PaleoclimateEngine

R = 8.31446
RNG = np.random.default_rng(20260905)
AMS_CEILING, MIN_AGE = 42000.0, 500.0
EA_INT = 173.0          # Ea used for the temperature integration itself (as published)
HERE = os.path.dirname(os.path.abspath(__file__))

master = pd.read_parquet(r"D:\26 Modelling Collagen Hydrolysis\data\harmonized_c14_master.parquet")
cohort = pd.read_parquet(r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet")

# ------------------------------------------------ select probable collagen ---
bu = master[master.material_category == "BONE_UNDIFFERENTIATED"].copy()
h = bu.material_raw.fillna("").str.lower()
CREMATED = r"cremat|calcin|burn|charred|ivory|antler"
bu["is_cremated"] = h.str.contains(CREMATED, regex=True)
prob = bu[~bu.is_cremated].dropna(subset=["latitude", "longitude", "c14_age"]).copy()
prob = prob[(prob.c14_age >= MIN_AGE) & (prob.c14_age < AMS_CEILING)]
print(f"BONE_UNDIFFERENTIATED      : {len(bu):,}")
print(f"  cremated/calcined/burnt  : {bu.is_cremated.sum():,}  (excluded - apatite chemistry)")
print(f"  probable collagen usable : {len(prob):,}")

# --------------------------------------- paleoclimate integration for those ---
engine = PaleoclimateEngine()
coords = prob[["latitude", "longitude"]].drop_duplicates().values
print(f"caching {len(coords):,} unique coordinates...")
cache, t0 = {}, time.time()
for i, (lat, lon) in enumerate(coords):
    if i % 2000 == 0:
        print(f"  {i:,}/{len(coords):,}  ({time.time()-t0:.0f}s)")
    tbp, tc, tk = engine.get_temperature_series(lat, lon)
    cache[(lat, lon)] = (tbp, tk)
print(f"cached in {time.time()-t0:.0f}s")

r_gas = 8.314462618e-3
t_ref_k = 283.15
temps = []
for _, r in prob.iterrows():
    tbp, tk = cache[(r.latitude, r.longitude)]
    m = tbp <= r.c14_age
    st, sk = tbp[m].tolist(), tk[m].tolist()
    if not st or st[-1] < r.c14_age:
        st.append(r.c14_age); sk.append(float(np.interp(r.c14_age, tbp, tk)))
    st = np.array(st); sk = np.nan_to_num(np.array(sk), nan=t_ref_k)
    arrh = np.exp(-(EA_INT / r_gas) * (1.0 / sk - 1.0 / t_ref_k))
    ta = np.trapezoid(arrh, st) if len(st) >= 2 else arrh[0] * r.c14_age
    ratio = max(ta / max(r.c14_age, 1.0), 1e-12)
    inv_t = 1.0 / t_ref_k - (r_gas / EA_INT) * np.log(ratio)
    temps.append(1.0 / inv_t - 273.15)
prob["integrated_temp_c"] = temps
prob = prob.replace([np.inf, -np.inf], np.nan).dropna(subset=["integrated_temp_c"])
prob = prob[prob.integrated_temp_c.between(-30, 32)]
print(f"probable collagen with integrated temperature: {len(prob):,}")

# ------------------------------------------------------ assemble the cohorts ---
base = cohort.dropna(subset=["integrated_temp_c", "c14_age", "site_name"]).copy()
base = base[(base.c14_age >= MIN_AGE) & (base.c14_age < AMS_CEILING)]
base["y"] = (base.material_category == "COLLAGEN").astype(int)

prob["y"] = 1
prob["material_category"] = "PROBABLE_COLLAGEN"
ext = pd.concat([base, prob[base.columns.intersection(prob.columns)]], ignore_index=True)
for d in (base, ext):
    d["z"] = 1000.0 / (d.integrated_temp_c + 273.15)
    d["ln_a"] = np.log(d.c14_age)

# ------------------------------------------------------------ the estimators ---
def binned_ea(d, tc=13.8, tmax=25.0, stat="max"):
    c = d[d.y == 1].copy(); c["b"] = np.round(c.integrated_temp_c)
    g = c.groupby("b").agg(n=("c14_age", "count"),
                           v=("c14_age", "max" if stat == "max" else (lambda x: np.percentile(x, 95))),
                           tm=("integrated_temp_c", "mean")).reset_index()
    g = g[(g.n >= 10) & (g.b >= tc) & (g.b <= tmax) & (g.v > 0)]
    if len(g) < 4: return np.nan, len(g), g
    s, _ = np.polyfit(1000.0 / (g.tm + 273.15), np.log(g.v), 1)
    return 1000.0 * R * s / 1000.0, len(g), g

def logistic_ea(d):
    X = sm.add_constant(pd.DataFrame({"z": d.z, "ln_a": d.ln_a}), has_constant="add")
    try: m = sm.Logit(d.y.values, X.values).fit(disp=0, maxiter=200)
    except Exception: return np.nan, np.nan, np.nan
    b, g = m.params[1], -m.params[2]
    if g <= 0: return np.nan, np.nan, g
    ea = 1000.0 * R * b / g / 1000.0
    V = m.cov_params()
    var = V[1,1]/g**2 + b**2*V[2,2]/g**4 + 2*b*V[1,2]/g**3
    return ea, 1000.0*R*np.sqrt(max(var,0))/1000.0, g

print("\n" + "="*74)
print("RECLASSIFICATION TEST: undifferentiated bone admitted as probable collagen")
print("="*74)
lines = []
for label, d in [("baseline (explicit COLLAGEN only)", base), ("extended (+ probable collagen)", ext)]:
    eb, nb, g = binned_ea(d)
    ebp, nbp, _ = binned_ea(d, stat="p95")
    el, se, gg = logistic_ea(d)
    ncol = int(d.y.sum())
    print(f"\n{label}")
    print(f"  collagen-class records : {ncol:,}")
    print(f"  binned max   Ea = {eb:6.1f}  ({nb} bins)")
    print(f"  binned Q95   Ea = {ebp:6.1f}  ({nbp} bins)")
    print(f"  logistic     Ea = {el:6.1f}  (SE {se:.1f}, g={gg:.2f})")
    lines.append((label, ncol, eb, nb, ebp, el, se))

_, _, gtab = binned_ea(ext)
_, _, gbase = binned_ea(base)
print("\n--- warm-bin coverage (n>=10), 14-25 C ---")
print("baseline bins:", sorted(gbase.b.astype(int).tolist()))
print("extended bins:", sorted(gtab.b.astype(int).tolist()))
print("\nextended warm bins:")
print(gtab.to_string(index=False))

md = f"""# Reclassification Test: Undifferentiated Bone as Probable Collagen
**Date & Time:** 2026-09-05 (generated by `probable_collagen_reclassification.py`)

## 1. Rationale

Unless a record states *cremated*, *calcined* or *burnt*, a bone radiocarbon
date is probably a collagen date: routine practice is to extract collagen, and
the `BONE_UNDIFFERENTIATED` class is largely a metadata-reporting gap rather
than a chemistry gap. At {len(bu):,} records it is larger than the explicit
`COLLAGEN` class. Cremated and calcined bone are excluded, since those dates are
made on the carbonate of the calcined mineral and are a different material.

Note also that the rarity of bone apatite dates ({(master.material_category=='BONE_APATITE').sum()} records,
~70 genuinely bone) reflects the fact that apatite dating is disliked,
technically difficult and controversial — it is a fact about practice, not a
harmonisation failure, and it cannot be repaired by re-parsing.

| Class | N |
| :--- | ---: |
| `BONE_UNDIFFERENTIATED` total | {len(bu):,} |
| cremated / calcined / burnt (excluded) | {bu.is_cremated.sum():,} |
| probable collagen, usable | {len(prob):,} |

## 2. Effect on the estimate

| Cohort | collagen records | binned max $E_a$ | bins | binned $Q_{{95}}$ $E_a$ | logistic $E_a$ | SE |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
""" + "\n".join(f"| {l} | {c:,} | {a:.1f} | {n} | {p:.1f} | {e:.1f} | {s:.1f} |"
                for l, c, a, n, p, e, s in lines) + f"""

## 3. Warm-bin coverage

The binned estimator previously rested on a single 5-site bin at 25 °C with no
admissible bins at 21–24 °C.

* baseline admissible warm bins: {sorted(gbase.b.astype(int).tolist())}
* extended admissible warm bins: {sorted(gtab.b.astype(int).tolist())}

```
{gtab.to_string(index=False)}
```

## 4. Interpretation

See manuscript §6.3. The relevant question is not whether the point estimate
moves but whether admitting ~{len(prob)//1000}k additional probable-collagen
determinations stabilises the warm end of the fit.
"""
out = os.path.join(HERE, "probable_collagen_reclassification_report.md")
open(out, "w", encoding="utf-8").write(md)
print(f"\nWrote {out}")
