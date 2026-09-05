import os, re

dir_path = r"c:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis\05_Palaeoproteomics_and_Collagen_Screening\cleaned_markdown"
for fn in os.listdir(dir_path):
    if "2026" in fn and "Holarctic" in fn or "Hola" in fn:
        full_p = os.path.join(dir_path, fn)
        with open(full_p, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
        print(f"File: {fn} ({len(content)} chars)")
        links = re.findall(r'https?://[^\s\)\"\'>]+', content)
        for l in links:
            if any(k in l.lower() for k in ["zenodo", "dryad", "figshare", "github", "10.5061", "10.5281"]):
                print("Repo link:", l)
