import pandas as pd
import pyreadr
import os, sys, re
from material_classifier import classify_material

data_dir = r"D:\26 Modelling Collagen Hydrolysis\data"

# 1. Ingest MEGA14C (Herrando-Pérez et al. 2026)
print("Ingesting MEGA14C...")
mega_path = os.path.join(data_dir, "MEGA14C_Dataset_figshare.xlsx")
df_mega = pd.read_excel(mega_path, sheet_name="MEGA14C_formatted")

mega_records = []
for _, r in df_mega.iterrows():
    lat = r.get("Latitude")
    lon = r.get("Longitude")
    age = r.get("Date")
    err = r.get("Error")
    mat_raw = str(r.get("Material", ""))
    pretreat = str(r.get("Pretreatment", ""))
    combined_mat = f"{mat_raw} | {pretreat}"
    
    # In MEGA14C, all records are purified collagen by definition of the dataset inclusion criteria
    category = "COLLAGEN"
    
    cn = r.get("CNratio")
    try:
        cn_val = float(cn) if pd.notna(cn) and float(cn) > 0 else None
    except:
        cn_val = None
        
    # Robust error parsing (handle asymmetric +X/-Y or text)
    try:
        err_val = float(err) if pd.notna(err) else None
    except:
        m = re.findall(r'\d+', str(err))
        err_val = float(m[0]) if m else None
        
    mega_records.append({
        "sample_id": str(r.get("Sample_ID") or r.get("Code")),
        "lab_id": str(r.get("Code")),
        "site_name": str(r.get("Site")),
        "country": str(r.get("Country")),
        "latitude": float(lat) if pd.notna(lat) else None,
        "longitude": float(lon) if pd.notna(lon) else None,
        "c14_age": float(age) if pd.notna(age) else None,
        "c14_error": err_val,
        "material_raw": combined_mat,
        "material_category": category,
        "cn_ratio": cn_val,
        "collagen_yield_pct": None,
        "taxa": str(r.get("Genus_Species")),
        "source_db": "MEGA14C"
    })
df_mega_clean = pd.DataFrame(mega_records)
print(f"Loaded {len(df_mega_clean)} MEGA14C records.")

# 2. Ingest p3k14c (Bird et al. 2022 / 2024)
print("Ingesting p3k14c...")
p3k_path = os.path.join(data_dir, "p3k14c_data.rda")
res_p3k = pyreadr.read_r(p3k_path)
df_p3k = res_p3k["p3k14c_data"]

p3k_records = []
for _, r in df_p3k.iterrows():
    lat = r.get("Lat")
    lon = r.get("Long")
    age = r.get("Age")
    err = r.get("Error")
    mat = str(r.get("Material", ""))
    category = classify_material(mat)
    
    p3k_records.append({
        "sample_id": str(r.get("LabID")),
        "lab_id": str(r.get("LabID")),
        "site_name": str(r.get("SiteName")),
        "country": str(r.get("Country")),
        "latitude": float(lat) if pd.notna(lat) else None,
        "longitude": float(lon) if pd.notna(lon) else None,
        "c14_age": float(age) if pd.notna(age) else None,
        "c14_error": float(err) if pd.notna(err) else None,
        "material_raw": mat,
        "material_category": category,
        "cn_ratio": None,
        "collagen_yield_pct": None,
        "taxa": str(r.get("Taxa")),
        "source_db": f"p3k14c_{r.get('Source')}"
    })
df_p3k_clean = pd.DataFrame(p3k_records)
print(f"Loaded {len(df_p3k_clean)} p3k14c records.")

# 3. Combine and Deduplicate
df_all = pd.concat([df_mega_clean, df_p3k_clean], ignore_index=True)
print(f"Total raw records before dedup: {len(df_all)}")

# Filter to records with valid coordinates and valid positive ages
df_valid = df_all.dropna(subset=["latitude", "longitude", "c14_age"]).copy()
df_valid = df_valid[(df_valid["c14_age"] >= 0) & (df_valid["c14_age"] <= 65000)]

# Prefer MEGA14C for duplicate lab_ids as it has curated chemistry
df_valid["priority"] = df_valid["source_db"].apply(lambda s: 1 if "MEGA14C" in s else 2)
df_valid = df_valid.sort_values(by=["priority"]).drop_duplicates(subset=["lab_id", "c14_age"], keep="first")
df_valid = df_valid.drop(columns=["priority"])

out_csv = os.path.join(data_dir, "harmonized_c14_master.csv")
out_parquet = os.path.join(data_dir, "harmonized_c14_master.parquet")
df_valid.to_csv(out_csv, index=False)
df_valid.to_parquet(out_parquet, index=False)

print("\n=== Master Harmonized Dataset Summary ===")
print(f"Total deduplicated valid records: {len(df_valid)}")
print("Records by material category:")
print(df_valid["material_category"].value_counts())
print(f"\nSaved to:\n  {out_csv}\n  {out_parquet}")
