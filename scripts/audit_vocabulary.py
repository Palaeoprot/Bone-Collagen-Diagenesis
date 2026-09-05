import pandas as pd
import numpy as np

p = r"D:\26 Modelling Collagen Hydrolysis\data\harmonized_c14_master.parquet"
df = pd.read_parquet(p)

audit_samples = []
for cat in ["COLLAGEN", "NON_COLLAGEN_ORGANIC_CONTROL"]:
    sub = df[df["material_category"] == cat]
    samp = sub.sample(n=50, random_state=42)
    audit_samples.append(samp)

audit_df = pd.concat(audit_samples)
print("=== 100-Record Vocabulary Audit ===")
mismatches = 0
for idx, r in audit_df.iterrows():
    cat = r["material_category"]
    raw = str(r["material_raw"]).lower()
    
    # Validation logic
    is_valid = False
    if cat == "COLLAGEN":
        if any(k in raw for k in ["collagen", "gelatin", "ultrafiltration", "xad", "amino acid"]):
            is_valid = True
    elif cat == "NON_COLLAGEN_ORGANIC_CONTROL":
        if any(k in raw for k in ["charcoal", "wood", "seed", "grain", "plant", "peat"]):
            is_valid = True
            
    if not is_valid:
        mismatches += 1
        print(f"FLAGGED: {cat} <-- '{r['material_raw']}' ({r['source_db']})")

print(f"Total inspected: 100")
print(f"Confirmed valid: {100 - mismatches} / 100 ({100 - mismatches}%)")
print(f"Mismatches: {mismatches}")
