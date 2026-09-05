import os
import sys
import re
import json
import shutil
import datetime
import email
from email import policy
from bs4 import BeautifulSoup
from markdownify import markdownify as md

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

PROJ_DIR = r"C:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis"
RAW_DIR = os.path.join(PROJ_DIR, "raw_papers")
LIB_BASE = r"C:\Users\matth\Documents\GitHub\Palaeoprot-Publications\library\thermal_age\collagen_hydrolysis_kinetics"
META_PATH = os.path.join(PROJ_DIR, "papers_metadata.json")

SUBFOLDERS = {
    '01_Peptide_Bond_and_Protein_Hydrolysis_Biochemistry': [
        'qian', 'hill', 'rawitscher', 'synge', 'kahne', 'bryant', 'radzicka', 'robert m. smith',
        'sun, yi', 'david, rolf', 'white, robert', 'schreiner', 'akbarian', 'instability challenges',
        'heats of hydrolysis', 'hydrolysis of a peptide bond', 'uncatalyzed rate', 'rates of uncatalyzed',
        'ph dependent mechanisms', 'structural and chemical changes induced', 'heat induced hydrolytic',
        'pan, bin', 'molecular mechanism of hydrolysis', 'ph-rate profile', 'hydrolysis of proteins',
        'hydrolytic stability of biomolecules', '1-s2.0-001670379390540d', '1-s2.0-s0065323308603885',
        '310430a0', 'biochemj00961', 'ja00230a041', 'ja01476a003', 'kinetics of peptide hydrolysis'
    ],
    '02_Collagen_Fibril_Architecture_and_Mechanochemistry': [
        'flynn', 'debnath', 'buhr', 'okada', 'degradation of collagen suture', 'single-fibril erosion',
        'modeling collagen fibril degradation', 'triple helix reduces the susceptibility',
        '1-s2.0-014296129290165k', 's10237-012-0399-2'
    ],
    '03_Bone_Diagenesis_and_Degradation_Modelling': [
        'polymer model', '1-s2.0-0141391094901139',
        'basic mathematical simulation', '1-s2.0-s0305440385700192',
        'survival of organic matter in bone', 'archaeometry - 2002 - collins',
        'predicting protein decomposition', 'rstb.1999.0359',
        'protein diagenesis: a reassessment', 'collins riley 2000',
        'fossil proteins in vertebrate calcified', 'armstrong',
        'glutamine deamidation in collagen', 'bone degradation using glutamine',
        'assessing the distribution of asian palaeolithic',
        'ftir-based model for the diagenetic alteration',
        'diagnosis of archaeological bones', 'archaeological collagen: why worry',
        'survival of immunological epitopes', 'a practical approach to the identification of low temperature',
        'age determination based on amino acid racemization', 'mitterer', 'bf00807707'
    ],
    '04_Radiocarbon_Databases_and_Chronologies': [
        'xronos', 'neotoma', 'intchron', 'canadian archaeological radiocarbon', 'card',
        'development of the intcal', 'the nerd dataset', 'aegean history and archaeology',
        'pacific archaeology radiocarbon', 'new radiocarbon database for the lower 48',
        'archaeological radiocarbon database for southern africa', 'making vertebrate fossil radiocarbon',
        'paleobiology database application'
    ],
    '05_Palaeoproteomics_and_Collagen_Screening': [
        'cappellini', 'wadsworth', 'harvey', 'procopio', 'peters, carli', 'malegori', 'boudin',
        'talamo', 'here we go again', 'buckley', 'grotta guattari',
        'bone need not remain an elephant', 'subtropical australia is preserved',
        'holarctic mammal collagen', 'near-infrared hyperspectral',
        'not just old but old and cold', 'weighing the mass spectrometric', 'mammoth femur',
        'collagen fingerprinting', 'characterization of proteomes extracted'
    ],
    '06_Theory_Quantum_and_Ensemble_Methods': [
        'reiher', 'tesei', 'bulow', 'blow', 'liu, hongbin', 'elucidating reaction mechanisms',
        'conformational ensembles of the human intrinsically', 'phase-separation propensities', 'quantum computing'
    ]
}

def sanitize_filename(name):
    # Remove problematic filesystem characters
    s = re.sub(r'[\\/*?:"<>|]', '', name)
    s = re.sub(r'\s+', ' ', s).strip()
    return s[:150]

def clean_mhtml_to_md(mhtml_path, out_md_path, doc_title, author=None, year=None, doi=None):
    with open(mhtml_path, "rb") as fp:
        msg = email.message_from_binary_file(fp, policy=policy.default)

    html_content = ""
    for part in msg.walk():
        if "html" in part.get_content_type():
            payload = part.get_payload(decode=True)
            charset = part.get_content_charset() or "utf-8"
            html_content = payload.decode(charset, errors="replace")
            break

    if not html_content:
        return False

    soup = BeautifulSoup(html_content, "html.parser")

    nature_main = soup.find("main", class_=re.compile(r"c-article-main-column"))
    pmc_main = soup.find("main", id="main-content")
    plos_main = soup.find("main", id="main-content")
    wiley_main = soup.find("section", class_=re.compile(r"article-section__content"))
    cambridge_main = soup.find("div", class_=re.compile(r"article-content|core-content"))

    if nature_main:
        main_elem = nature_main
        for bad in main_elem.find_all(class_=re.compile(r"(c-article-recommendations|c-context-bar|c-article-share-box|c-site-messages)", re.I)):
            bad.decompose()
    elif pmc_main:
        main_elem = pmc_main
        for bad in main_elem.find_all(class_=re.compile(r"(usa-banner|pmc-footer|share-tool|social-share)", re.I)):
            bad.decompose()
    elif plos_main:
        main_elem = plos_main
        for tab in main_elem.find_all(class_=re.compile(r"(tab-list|article-tabs|nav-tabs)", re.I)):
            tab.decompose()
    elif wiley_main:
        main_elem = wiley_main
    elif cambridge_main:
        main_elem = cambridge_main
    else:
        main_elem = (
            soup.find("article", class_=re.compile(r"article|content|paper", re.I))
            or soup.find("main")
            or soup.find("article")
            or soup.body
        )

    if not main_elem:
        return False

    for tag in main_elem(["script", "style", "nav", "footer", "noscript", "svg", "form", "iframe", "header"]):
        tag.decompose()

    for tag in main_elem.find_all(class_=re.compile(r"(sidebar|navigation|menu|advert|banner|cookie|share|social|metric)", re.I)):
        tag.decompose()

    md_content = md(str(main_elem), heading_style="ATX", strip=["img"])
    md_content = re.sub(r"\n{3,}", "\n\n", md_content).strip()

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S (+02:00)")
    meta_lines = []
    if author: meta_lines.append(f"**Authors:** {author}")
    if year: meta_lines.append(f"**Year:** {year}")
    if doi: meta_lines.append(f"**DOI:** [{doi}](https://doi.org/{doi})")
    
    meta_str = "\n".join(meta_lines) + "\n\n" if meta_lines else ""
    header = f"# {doc_title}\n\n**Date & Time:** {now}\n\n{meta_str}"
    final_md = header + md_content

    os.makedirs(os.path.dirname(os.path.abspath(out_md_path)), exist_ok=True)
    with open(out_md_path, "w", encoding="utf-8") as out_fp:
        out_fp.write(final_md)
    return True

def run_migration():
    print("=== Migrating Literature to Library ===")
    
    with open(META_PATH, "r", encoding="utf-8") as f:
        metadata_list = json.load(f)
        
    meta_by_fn = {p["filename"]: p for p in metadata_list}
    
    # Process all files in raw_papers
    if os.path.exists(RAW_DIR):
        raw_files = os.listdir(RAW_DIR)
        print(f"Found {len(raw_files)} files in raw_papers/.")
        
        for fn in raw_files:
            src_fp = os.path.join(RAW_DIR, fn)
            meta = meta_by_fn.get(fn, {})
            title = meta.get("title") or os.path.splitext(fn)[0]
            author = meta.get("author") or ""
            year = meta.get("year") or ""
            doi = meta.get("doi") or ""
            
            # Determine subfolder
            target_sub = meta.get("folder")
            if not target_sub:
                text_to_match = f"{fn} {title} {author}".lower()
                for folder, keywords in SUBFOLDERS.items():
                    if any(kw.lower() in text_to_match for kw in keywords):
                        target_sub = folder
                        break
            if not target_sub:
                target_sub = "01_Peptide_Bond_and_Protein_Hydrolysis_Biochemistry"
                
            dest_dir = os.path.join(LIB_BASE, target_sub)
            os.makedirs(dest_dir, exist_ok=True)
            
            # Construct standard academic filename
            first_author = author.split(",")[0].strip() if author else "Unknown"
            if first_author and year:
                clean_base = f"{first_author} {year} - {sanitize_filename(title)}"
            else:
                clean_base = sanitize_filename(title)
                
            if fn.lower().endswith(".mhtml"):
                out_md = os.path.join(dest_dir, f"{clean_base}.md")
                out_mhtml = os.path.join(dest_dir, f"{clean_base}.mhtml")
                success = clean_mhtml_to_md(src_fp, out_md, title, author, year, doi)
                if success:
                    shutil.copy2(src_fp, out_mhtml)
                    print(f"[CONVERTED & FILED] {clean_base}.md -> {target_sub}")
                else:
                    shutil.copy2(src_fp, out_mhtml)
                    print(f"[FILED MHTML] {clean_base}.mhtml -> {target_sub}")
            elif fn.lower().endswith(".pdf"):
                out_pdf = os.path.join(dest_dir, f"{clean_base}.pdf")
                shutil.copy2(src_fp, out_pdf)
                print(f"[FILED PDF] {clean_base}.pdf -> {target_sub}")

    # Process Peters et al. 2023 and Colin Smith (2002) thesis
    peters_dir = os.path.join(LIB_BASE, "05_Palaeoproteomics_and_Collagen_Screening", "Peters_2023_Subtropical_Australia_Collagen")
    os.makedirs(peters_dir, exist_ok=True)
    
    peters_patterns = ["43247_2023", "peters_supp", "MOESM"]
    for f in os.listdir(PROJ_DIR):
        if any(pat in f for pat in peters_patterns):
            src_fp = os.path.join(PROJ_DIR, f)
            if os.path.isfile(src_fp):
                dst_fp = os.path.join(peters_dir, f)
                shutil.copy2(src_fp, dst_fp)
                print(f"[FILED PETERS DATA] {f} -> {peters_dir}")
                
    smith_src = os.path.join(PROJ_DIR, "Smith02.pdf")
    if os.path.exists(smith_src):
        smith_dst = os.path.join(r"C:\Users\matth\Documents\GitHub\Palaeoprot-Publications\library\thermal_age", "Smith 2002 - The Rate of Collagen Degradation in Bone (PhD Thesis).pdf")
        shutil.copy2(smith_src, smith_dst)
        print(f"[FILED SMITH THESIS] -> {smith_dst}")
        
    print("\n=== Cleaning Up Loose Files in Project Folder ===")
    # Remove raw_papers folder
    if os.path.exists(RAW_DIR):
        shutil.rmtree(RAW_DIR)
        print(f"Removed redundant raw_papers directory: {RAW_DIR}")
        
    # Remove cleaned_papers if empty
    cleaned_dir = os.path.join(PROJ_DIR, "cleaned_papers")
    if os.path.exists(cleaned_dir):
        shutil.rmtree(cleaned_dir)
        print(f"Removed cleaned_papers directory.")

    # Remove loose ESM/MOESM files from project root
    for f in os.listdir(PROJ_DIR):
        if any(pat in f for pat in peters_patterns) or f == "Smith02.pdf":
            fp = os.path.join(PROJ_DIR, f)
            if os.path.isfile(fp):
                os.remove(fp)
                print(f"Removed loose project file: {f}")

    # Remove temporary subfolders 01_ to 06_ if they were created in project root
    for sub in SUBFOLDERS:
        sub_p = os.path.join(PROJ_DIR, sub)
        if os.path.exists(sub_p):
            shutil.rmtree(sub_p)
            print(f"Removed redundant project subfolder: {sub}")

    print("\n=== Literature Migration to Library Successfully Completed ===")

if __name__ == "__main__":
    run_migration()
