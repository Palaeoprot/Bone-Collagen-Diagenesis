import os, sys, time
import pandas as pd
import numpy as np
from paleoclimate_engine import PaleoclimateEngine

data_dir = r"D:\26 Modelling Collagen Hydrolysis\data"
master_path = os.path.join(data_dir, "harmonized_c14_master.parquet")
print(f"Loading master dataset from {master_path}...")
df = pd.read_parquet(master_path)

# Filter for analysis:
# 1. Purified Bone Collagen vs Non-Collagen Organic Controls
# 2. Deconvolved from instrumental limit: calendar ages < 42,000 BP (to eliminate instrumental blank ceiling bias)
# 3. Ages > 500 BP (historical samples excluded)
df_analysis = df[
    (df["material_category"].isin(["COLLAGEN", "NON_COLLAGEN_ORGANIC_CONTROL"])) &
    (df["c14_age"] >= 500) &
    (df["c14_age"] <= 42000)
].copy()

print(f"Total analysis cohort: {len(df_analysis)} records.")
print(df_analysis["material_category"].value_counts())

# Downsample negative controls to balance sample size for compute efficiency while retaining robust stats
collagen_df = df_analysis[df_analysis["material_category"] == "COLLAGEN"]
control_df = df_analysis[df_analysis["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"].sample(n=len(collagen_df), random_state=42)

cohort = pd.concat([collagen_df, control_df]).reset_index(drop=True)
print(f"Balanced study cohort: {len(cohort)} records ({len(collagen_df)} collagen, {len(control_df)} controls).")

# Initialize paleoclimate engine
engine = PaleoclimateEngine()

# Pre-cache unique grid coordinates to drastically speed up processing
coords = cohort[["latitude", "longitude"]].drop_duplicates().values
print(f"Extracting paleotemperature histories across {len(coords)} unique coordinate locations...")

coord_cache = {}
t0 = time.time()
for i, (lat, lon) in enumerate(coords):
    if i % 2000 == 0:
        print(f"  Processed {i} / {len(coords)} coordinates ({(time.time()-t0):.1f}s)...")
    t_bp, tc, tk = engine.get_temperature_series(lat, lon)
    coord_cache[(lat, lon)] = (t_bp, tk)

print(f"Coordinate caching completed in {(time.time()-t0):.1f}s.")

# Function to compute thermal ages for a specific Ea
def compute_cohort_thermal_ages(cohort_df, ea_kj: float, t_ref_c: float = 10.0, delta_t_micro: float = 0.0):
    t_ref_k = t_ref_c + 273.15
    r_gas = 8.314462618e-3 # kJ / (mol * K)
    
    thermal_ages = []
    for idx, r in cohort_df.iterrows():
        lat = r["latitude"]
        lon = r["longitude"]
        cal_age = r["c14_age"]
        
        times_bp, tk = coord_cache[(lat, lon)]
        sub_temps = tk + delta_t_micro
        
        mask = times_bp <= cal_age
        sub_t = times_bp[mask].tolist()
        sub_k = sub_temps[mask].tolist()
        
        if len(sub_t) == 0 or sub_t[-1] < cal_age:
            t_interp = np.interp(cal_age, times_bp, sub_temps)
            sub_t.append(cal_age)
            sub_k.append(t_interp)
            
        sub_t = np.array(sub_t)
        sub_k = np.nan_to_num(np.array(sub_k), nan=t_ref_k)
        
        arrh = np.exp(-(ea_kj / r_gas) * (1.0 / sub_k - 1.0 / t_ref_k))
        
        if len(sub_t) >= 2:
            th_age = np.trapezoid(arrh, sub_t) if hasattr(np, "trapezoid") else np.trapz(arrh, sub_t)
        else:
            th_age = arrh[0] * cal_age
            
        thermal_ages.append(float(th_age))
        
    return np.array(thermal_ages)

# Compute thermal ages at Ea = 173 kJ/mol (Lab benchmark) and Ea = 100 kJ/mol (Uncatalyzed solution benchmark)
print("\nComputing thermal ages at Ea = 173 kJ/mol (Collins lab estimate)...")
t0 = time.time()
cohort["thermal_age_173"] = compute_cohort_thermal_ages(cohort, ea_kj=173.0)
print(f"Done in {(time.time()-t0):.1f}s.")

print("Computing thermal ages at Ea = 100 kJ/mol (Uncatalyzed solution dipeptide)...")
cohort["thermal_age_100"] = compute_cohort_thermal_ages(cohort, ea_kj=100.0)

out_cohort = os.path.join(data_dir, "collagen_vs_control_thermal_cohort.parquet")
cohort.to_parquet(out_cohort, index=False)
print(f"Saved thermal age cohort to: {out_cohort}")

print("\n=== Thermal Age Distribution Summary (Years @ 10°C) ===")
for cat in ["COLLAGEN", "NON_COLLAGEN_ORGANIC_CONTROL"]:
    sub = cohort[cohort["material_category"] == cat]
    print(f"\nCategory: {cat} (N = {len(sub):,})")
    print(f"  Cal Age Mean: {sub['c14_age'].mean():,.1f} | 95th %: {sub['c14_age'].quantile(0.95):,.1f} | Max: {sub['c14_age'].max():,.1f}")
    print(f"  Therm Age (173 kJ) Mean: {sub['thermal_age_173'].mean():,.1f} | 95th %: {sub['thermal_age_173'].quantile(0.95):,.1f} | Max: {sub['thermal_age_173'].max():,.1f}")
    print(f"  Therm Age (100 kJ) Mean: {sub['thermal_age_100'].mean():,.1f} | 95th %: {sub['thermal_age_100'].quantile(0.95):,.1f} | Max: {sub['thermal_age_100'].max():,.1f}")
