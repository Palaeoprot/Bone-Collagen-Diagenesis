import os
import sys
import re
import pandas as pd
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from material_classifier import classify_material

DATA_DIR = r"D:\26 Modelling Collagen Hydrolysis\data"
RAW_DIR = os.path.join(DATA_DIR, "raw_sources")

def clean_lab_id(val):
    if not val or pd.isna(val):
        return None
    s = str(val).strip()
    s = re.sub(r'\s+', '-', s)
    return s if len(s) > 1 else None

def parse_float(val):
    if pd.isna(val):
        return None
    try:
        f = float(val)
        return f if np.isfinite(f) else None
    except:
        m = re.findall(r'[-+]?\d*\.?\d+', str(val))
        return float(m[0]) if m else None

def parse_mega14c(filepath):
    print("Parsing MEGA14C...")
    df = pd.read_excel(filepath, sheet_name="MEGA14C_formatted")
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("Latitude"))
        lon = parse_float(r.get("Longitude"))
        age = parse_float(r.get("Date"))
        err = parse_float(r.get("Error"))
        mat_raw = str(r.get("Material", ""))
        pretreat = str(r.get("Pretreatment", ""))
        combined_mat = f"{mat_raw} | {pretreat}".strip(" |")
        cn = parse_float(r.get("CNratio"))
        lab = clean_lab_id(r.get("Code"))
        
        records.append({
            "sample_id": str(r.get("Sample_ID") or lab),
            "lab_id": lab,
            "site_name": str(r.get("Site") or ""),
            "country": str(r.get("Country") or ""),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": combined_mat,
            "material_category": "COLLAGEN", # Curated ultrafiltered collagen by definition
            "cn_ratio": cn if (cn and cn > 0) else None,
            "collagen_yield_pct": None,
            "taxa": str(r.get("Genus_Species") or ""),
            "source_db": "MEGA14C",
            "priority": 1
        })
    return pd.DataFrame(records)

def parse_xronos(filepath):
    print("Parsing XRONOS...")
    df = pd.read_csv(filepath, encoding="utf-8", low_memory=False)
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("lat"))
        lon = parse_float(r.get("lng"))
        age = parse_float(r.get("bp"))
        err = parse_float(r.get("std"))
        mat_raw = str(r.get("material") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("labnr"))
        src_db = str(r.get("source_database") or "xronos")
        
        records.append({
            "sample_id": str(r.get("id") or lab),
            "lab_id": lab,
            "site_name": str(r.get("site") or ""),
            "country": str(r.get("country") or ""),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": str(r.get("species") or ""),
            "source_db": f"XRONOS_{src_db}",
            "priority": 3
        })
    return pd.DataFrame(records)

def parse_p3k14c(filepath):
    print("Parsing p3k14c...")
    df = pd.read_csv(filepath, encoding="utf-8", low_memory=False)
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("Lat"))
        lon = parse_float(r.get("Long"))
        age = parse_float(r.get("Age"))
        err = parse_float(r.get("Error"))
        mat_raw = str(r.get("Material") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("LabID"))
        src = str(r.get("Source") or "p3k14c")
        
        records.append({
            "sample_id": lab,
            "lab_id": lab,
            "site_name": str(r.get("SiteName") or ""),
            "country": str(r.get("Country") or ""),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": str(r.get("Taxa") or ""),
            "source_db": f"p3k14c_{src}",
            "priority": 3
        })
    return pd.DataFrame(records)

def parse_14cpalaeolithic(filepath):
    print("Parsing 14C Palaeolithic Europe...")
    df = pd.read_excel(filepath, header=1)
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("Lat"))
        lon = parse_float(r.get("Long"))
        age = parse_float(r.get("Age"))
        err = parse_float(r.get("pm"))
        mat_raw = str(r.get("sample") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("labref"))
        
        records.append({
            "sample_id": str(r.get("id") or lab),
            "lab_id": lab,
            "site_name": str(r.get("sitename") or ""),
            "country": str(r.get("country") or ""),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": None,
            "source_db": "14cpalaeolithic_v33",
            "priority": 2
        })
    return pd.DataFrame(records)

def parse_calpal(filepath):
    print("Parsing CalPal...")
    df = pd.read_csv(filepath, sep="\t", encoding="utf-8", on_bad_lines="skip", low_memory=False)
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("LATITUDE"))
        lon = parse_float(r.get("LONGITUDE"))
        age = parse_float(r.get("C14AGE"))
        err = parse_float(r.get("C14STD"))
        mat_raw = str(r.get("MATERIAL") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("LABNR"))
        
        records.append({
            "sample_id": str(r.get("ID") or lab),
            "lab_id": lab,
            "site_name": str(r.get("SITE") or ""),
            "country": str(r.get("COUNTRY") or ""),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": str(r.get("SPECIES") or ""),
            "source_db": "CalPal_2020",
            "priority": 2
        })
    return pd.DataFrame(records)

def parse_radonb(filepath):
    print("Parsing Radon-B...")
    df = pd.read_csv(filepath, sep="\t", encoding="utf-8", on_bad_lines="skip", low_memory=False)
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("LATITUDE"))
        lon = parse_float(r.get("LONGITUDE"))
        age = parse_float(r.get("C14AGE"))
        err = parse_float(r.get("C14STD"))
        mat_raw = str(r.get("MATERIAL") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("LABNR"))
        
        records.append({
            "sample_id": str(r.get("ID") or lab),
            "lab_id": lab,
            "site_name": str(r.get("SITE") or ""),
            "country": str(r.get("COUNTRY") or ""),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": str(r.get("SPECIES") or ""),
            "source_db": "RadonB",
            "priority": 2
        })
    return pd.DataFrame(records)

def parse_rxpand(filepath):
    print("Parsing RXPAND (South America)...")
    df = pd.read_csv(filepath, encoding="utf-8", low_memory=False)
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("Latitude"))
        lon = parse_float(r.get("Longitude"))
        age = parse_float(r.get("C14Age"))
        err = parse_float(r.get("C14SD"))
        mat_raw = str(r.get("Material") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("LabCode"))
        
        records.append({
            "sample_id": lab,
            "lab_id": lab,
            "site_name": str(r.get("Site") or ""),
            "country": "South America",
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": None,
            "source_db": "RXPAND_SouthAmerica",
            "priority": 2
        })
    return pd.DataFrame(records)

def parse_austarch(filepath):
    print("Parsing AustArch (Australia)...")
    df = pd.read_csv(filepath, encoding="latin1", low_memory=False)
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("LATITUDE"))
        lon = parse_float(r.get("LONGITUDE"))
        age = parse_float(r.get("AGE"))
        err = parse_float(r.get("ERROR"))
        mat_raw = str(r.get("MATERIAL") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("LAB_CODE"))
        
        records.append({
            "sample_id": str(r.get("ADSID") or lab),
            "lab_id": lab,
            "site_name": str(r.get("SITE") or ""),
            "country": "Australia",
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": None,
            "source_db": "AustArch",
            "priority": 2
        })
    return pd.DataFrame(records)

def parse_14sea(filepath):
    print("Parsing 14SEA (Southeast Asia)...")
    df = pd.read_excel(filepath)
    records = []
    for _, r in df.iterrows():
        age = parse_float(r.get("Date BP"))
        err = parse_float(r.get("±"))
        mat_raw = str(r.get("Material") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("Lab. no."))
        
        records.append({
            "sample_id": lab,
            "lab_id": lab,
            "site_name": str(r.get("Site") or ""),
            "country": str(r.get("Country") or "Southeast Asia"),
            "latitude": None, # Coordinates resolved in XRONOS cross-match if available
            "longitude": None,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": None,
            "source_db": "14SEA",
            "priority": 2
        })
    return pd.DataFrame(records)

def parse_sard(filepath):
    print("Parsing SARD (Southern Africa)...")
    df = pd.read_csv(filepath, encoding="utf-8")
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("DecdegS"))
        lon = parse_float(r.get("DecdegE"))
        age = parse_float(r.get("Date"))
        err = parse_float(r.get("Uncertainty"))
        mat_raw = str(r.get("Material dated") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("Lab ID"))
        
        records.append({
            "sample_id": lab,
            "lab_id": lab,
            "site_name": str(r.get(" Site") or r.get("Site") or ""),
            "country": str(r.get("Country") or "South Africa"),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": str(r.get(" Species") or ""),
            "source_db": "SARD",
            "priority": 2
        })
    return pd.DataFrame(records)

def parse_caribbean(filepath):
    print("Parsing Caribbean 14C...")
    df = pd.read_csv(filepath, encoding="utf-8")
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("Lat"))
        lon = parse_float(r.get("Lon"))
        age = parse_float(r.get("Age"))
        err = parse_float(r.get("Error"))
        mat_raw = str(r.get("Material") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("LabNo"))
        
        records.append({
            "sample_id": str(r.get("UniqID") or lab),
            "lab_id": lab,
            "site_name": str(r.get("SiteName") or ""),
            "country": str(r.get("Country.Territory") or "Caribbean"),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": None,
            "source_db": "Caribbean14C",
            "priority": 2
        })
    return pd.DataFrame(records)

def parse_adrac(filepath):
    print("Parsing aDRAC (Central Africa)...")
    df = pd.read_csv(filepath, encoding="utf-8", low_memory=False)
    records = []
    for _, r in df.iterrows():
        lat = parse_float(r.get("LAT"))
        lon = parse_float(r.get("LONG"))
        age = parse_float(r.get("C14AGE"))
        err = parse_float(r.get("C14STD"))
        mat_raw = str(r.get("MATERIAL") or "")
        cat = classify_material(mat_raw)
        lab = clean_lab_id(r.get("LABNR"))
        
        records.append({
            "sample_id": lab,
            "lab_id": lab,
            "site_name": str(r.get("SITE") or ""),
            "country": str(r.get("COUNTRY") or "Central Africa"),
            "latitude": lat,
            "longitude": lon,
            "c14_age": age,
            "c14_error": err,
            "material_raw": mat_raw,
            "material_category": cat,
            "cn_ratio": None,
            "collagen_yield_pct": None,
            "taxa": None,
            "source_db": "aDRAC",
            "priority": 2
        })
    return pd.DataFrame(records)

def run_harmonization():
    print("=== Commencing Global Radiocarbon Database Harmonization ===")
    
    dfs = []
    
    # 1. MEGA14C
    mega_fp = os.path.join(DATA_DIR, "MEGA14C_Dataset_figshare.xlsx")
    if os.path.exists(mega_fp):
        dfs.append(parse_mega14c(mega_fp))
        
    # 2. Primary Regional Repositories
    parsers = [
        ("14cpalaeolithic_v33.xlsx", parse_14cpalaeolithic),
        ("CalPal_2020_08_20.tsv", parse_calpal),
        ("radonb_daily.txt", parse_radonb),
        ("rxpand_south_america.csv", parse_rxpand),
        ("austarch_australia.csv", parse_austarch),
        ("14sea_southeast_asia.xlsx", parse_14sea),
        ("sard_south_africa.csv", parse_sard),
        ("caribbean_14c.csv", parse_caribbean),
        ("adrac_central_africa.csv", parse_adrac),
        ("p3k14c_scrubbed_fuzzed.csv", parse_p3k14c),
        ("xronos_data.csv", parse_xronos)
    ]
    
    for fname, parser_fn in parsers:
        fp = os.path.join(RAW_DIR, fname)
        if os.path.exists(fp):
            try:
                df_parsed = parser_fn(fp)
                print(f"  -> Extracted {len(df_parsed):,} records from {fname}")
                dfs.append(df_parsed)
            except Exception as e:
                print(f"  [ERROR] Failed to parse {fname}: {e}")
                
    print("\nConcatenating all extracted datasets...")
    df_all = pd.concat(dfs, ignore_index=True)
    print(f"Total raw compiled records: {len(df_all):,}")
    
    # Filter to valid determinations
    print("\nFiltering to valid entries with valid coordinates and ages [0, 65,000 BP]...")
    df_valid = df_all.dropna(subset=["latitude", "longitude", "c14_age"]).copy()
    df_valid = df_valid[
        (df_valid["c14_age"] >= 0) & 
        (df_valid["c14_age"] <= 65000) &
        (df_valid["latitude"] >= -90) & (df_valid["latitude"] <= 90) &
        (df_valid["longitude"] >= -180) & (df_valid["longitude"] <= 180)
    ]
    print(f"Valid geolocated records: {len(df_valid):,}")
    
    # Hierarchical Deduplication
    print("\nPerforming Hierarchical Deduplication...")
    # Sort by priority (1 = MEGA14C, 2 = Regional, 3 = Meta-aggregators)
    # Also prioritize non-null lab_id
    df_valid["has_lab"] = df_valid["lab_id"].notna().astype(int)
    df_valid = df_valid.sort_values(by=["priority", "has_lab"], ascending=[True, False])
    
    # Deduplicate on (lab_id, c14_age) where lab_id is valid, and on (latitude, longitude, c14_age, material_category) where lab_id is missing
    df_with_lab = df_valid[df_valid["lab_id"].notna()].drop_duplicates(subset=["lab_id", "c14_age"], keep="first")
    df_no_lab = df_valid[df_valid["lab_id"].isna()].drop_duplicates(subset=["latitude", "longitude", "c14_age", "material_category"], keep="first")
    
    df_master = pd.concat([df_with_lab, df_no_lab], ignore_index=True)
    df_master = df_master.drop(columns=["has_lab", "priority"])
    
    # Output file paths
    out_csv = os.path.join(DATA_DIR, "harmonized_c14_master_expanded.csv")
    out_parquet = os.path.join(DATA_DIR, "harmonized_c14_master_expanded.parquet")
    
    print(f"\nWriting master expanded datasets...")
    df_master.to_csv(out_csv, index=False)
    df_master.to_parquet(out_parquet, index=False)
    
    print("\n================ Master Harmonized Dataset Summary ================")
    print(f"Total Unique Deduplicated Records: {len(df_master):,}")
    print("\nBreakdown by Material Category:")
    cat_counts = df_master["material_category"].value_counts()
    for cat, cnt in cat_counts.items():
        pct = cnt / len(df_master) * 100
        print(f"  {cat:30s}: {cnt:7,} ({pct:5.2f}%)")
        
    print(f"\nSaved successfully to:\n  {out_csv}\n  {out_parquet}")
    return df_master

if __name__ == "__main__":
    run_harmonization()
