import re

def classify_material(mat_str):
    if not mat_str or not isinstance(mat_str, str):
        return "UNKNOWN"
    s = mat_str.strip().lower()
    
    # 1. Purified Bone Collagen (Target Signal)
    # Includes English, German (kollagen/gelatine), French (collagène/gélatine), Spanish (colágeno/gelatina), Italian (collagene)
    if re.search(r'(collagen|gelatin|kollagen|gelatine|collag[eè]ne|g[eé]latine|col[aá]geno|hydroxyproline|amino\s*acids?|ultrafilter|xad|amidase|v-filter|longin)', s):
        if not re.search(r'(apatite|bioapatite|carbonate|carbonat|carbonato|enamel|schmelz|émail|esmalte)', s):
            return "COLLAGEN"

    # 2. Bone Apatite / Carbonate (Explicitly Excluded)
    if re.search(r'(apatite|bioapatite|carbonate|carbonat|carbonato|enamel|schmelz|émail|esmalte|fraction\s+inorganique)', s):
        return "BONE_APATITE"
        
    # 3. Whole / Undifferentiated Bone & Teeth (Quarantine)
    if re.search(r'(bone|skeleton|mandible|tibia|femur|humerus|skull|rib|vertebra|antler|horn|tooth|teeth|phalanx|metatarsal|metacarpal|ivory|dentine|knochen|geweih|zahn|zähne|skelett|elfenbein|\bos\b|squelette|bois\s+de\s+cerf|\bdent\b|ivoire|hueso|esqueleto|asta|diente|marfil|molar|premolar|crane|cranium)', s):
        return "BONE_UNDIFFERENTIATED"
        
    # 4. Organic Non-Collagen Controls (Negative Controls)
    if re.search(r'(charcoal|charred|wood|timber|seed|grain|cereal|plant|macrofossil|flora|bark|leaf|leaves|peat|textile|linen|fibre|fiber|pitch|resin|coprolite|dung|straw|mastic|papyrus|holzkohle|holz|samen|korn|getreide|pflanze|makrofossil|torf|rinde|blatt|gewebe|leinen|pech|harz|charbon|bois|graine|c[eé]r[eé]ale|plante|macrofossile|tourbe|[eé]corce|feuille|textile|lin|r[eé]sine|carb[oó]n|madera|semilla|grano|planta|macrof[oó]sil|turba|corteza|hoja|lino|carbonella|legno|seme|tessuto)', s):
        return "NON_COLLAGEN_ORGANIC_CONTROL"
        
    # 5. Shell / Marine Carbonate Controls
    if re.search(r'(shell|mollusc|mollusk|conch|marine|oyster|gastropod|bivalve|coral|foram|muschel|schnecke|schale|koralle|coquille|mollusque|hu[iî]tre|corail|concha|molusco|ostra|valva|conchiglia|ostrica|corallo)', s):
        return "SHELL_CONTROL"
        
    # 6. Sediment / Soil / Humic / Speleothem
    if re.search(r'(sediment|soil|humus|humic|silt|clay|gyttja|stalagmite|speleothem|calcite|travertine|tufa|boden|humins[aä]ure|sédiment|sol|terre|sedi)', s):
        return "SEDIMENT"

    return "OTHER"

if __name__ == "__main__":
    test_cases = [
        ("collagen bone", "COLLAGEN"),
        ("Ultrafiltration(>30kDa)-Gelatin(NaOH/Collagen)", "COLLAGEN"),
        ("bone collagen", "COLLAGEN"),
        ("Kollagen", "COLLAGEN"),
        ("collagène extrait d'os", "COLLAGEN"),
        ("colágeno óseo", "COLLAGEN"),
        ("amino acids from bone", "COLLAGEN"),
        ("hydroxyproline", "COLLAGEN"),
        ("bone", "BONE_UNDIFFERENTIATED"),
        ("Knochen", "BONE_UNDIFFERENTIATED"),
        ("os animal", "BONE_UNDIFFERENTIATED"),
        ("hueso de mamífero", "BONE_UNDIFFERENTIATED"),
        ("mandible", "BONE_UNDIFFERENTIATED"),
        ("antler", "BONE_UNDIFFERENTIATED"),
        ("geweih", "BONE_UNDIFFERENTIATED"),
        ("tooth (molar)", "BONE_UNDIFFERENTIATED"),
        ("bone apatite", "BONE_APATITE"),
        ("tooth enamel carbonate", "BONE_APATITE"),
        ("charcoal", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("CHARCOAL", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("Holzkohle", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("charbon de bois", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("carbón vegetal", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("charred seeds", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("wood", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("Holz", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("grain (charred)", "NON_COLLAGEN_ORGANIC_CONTROL"),
        ("SHELL", "SHELL_CONTROL"),
        ("marine shell", "SHELL_CONTROL"),
        ("Muschelschale", "SHELL_CONTROL"),
        ("coquille marine", "SHELL_CONTROL"),
        ("sediment", "SEDIMENT"),
        ("unknown", "OTHER")
    ]
    for raw, expected in test_cases:
        res = classify_material(raw)
        assert res == expected, f"Failed on '{raw}': got '{res}', expected '{expected}'"
    print(f"All {len(test_cases)} material classification unit test cases passed successfully!")

