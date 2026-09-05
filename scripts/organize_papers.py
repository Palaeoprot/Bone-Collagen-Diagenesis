import os, json, shutil, re

proj_dir = r'c:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis'
raw_dir = os.path.join(proj_dir, 'raw_papers')

subfolders = {
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

with open(os.path.join(proj_dir, 'papers_metadata.json'), 'r', encoding='utf-8') as f:
    papers = json.load(f)

for sub in subfolders:
    os.makedirs(os.path.join(proj_dir, sub), exist_ok=True)

unassigned = []
assigned_count = 0

for p in papers:
    fn = p['filename']
    title = p.get('title') or ''
    author = p.get('author') or ''
    text_to_match = f"{fn} {title} {author}".lower()
    
    assigned = None
    for folder, keywords in subfolders.items():
        if any(kw.lower() in text_to_match for kw in keywords):
            assigned = folder
            break
            
    if assigned:
        assigned_count += 1
        p['folder'] = assigned
        src = os.path.join(raw_dir, fn)
        dst = os.path.join(proj_dir, assigned, fn)
        shutil.copy2(src, dst)
    else:
        unassigned.append(p)

print(f"Assigned: {assigned_count} / {len(papers)}")
if unassigned:
    print("Unassigned:")
    for u in unassigned:
        print(" ", u['filename'], "-->", u.get('title'))

with open(os.path.join(proj_dir, 'papers_metadata.json'), 'w', encoding='utf-8') as f:
    json.dump(papers, f, indent=2, ensure_ascii=False)
