import sys
import re
import pandas as pd

# ============================================================
# Nettoyage de la colonne `melting_point` → `melting_point_k`
# Comportement fail-fast : s'arrête au PREMIER format non reconnu.
# ============================================================

# ------------------------------------------------------------
# Formats pris en charge (enrichi au fur et à mesure) :
#   '485 K' / '558.0(1) K' → val K   (incertitude cristallographique (...) optionnelle)
#   '345deg. K' / '548 deg. K (dec.)' → val K  (notation deg. K, annotation optionnelle)
#   '190-192 deg.C'  → milieu + 273.15 (plage °C → K)
#   '327-328 K'      → milieu          (plage K, tiret optionnellement entouré d'espaces)
#   'decomposes …' / 'thermal decomposition' → None  (_QUALITATIVE : préfixe 'decompos', décomp. ≠ fusion)
#   'above …' / 'below …' / 'over …'  → None  (_QUALITATIVE : bornes imprécises)
#   '72deg.C(dec.)' / 'from 175deg.C(dec.)' → val+273.15  (°C ±, préfixe 'from' optionnel, annot. décomp.)
#   '178-181deg.C(dec.)' / '199-200deg.C(lit.)' → milieu + 273.15  (plage °C + annotation quelconque)
#   '519 K(dec.)' / '510K.(dec.)' / '690 K (decomp.)' → val  (K + annotation décomp., (dec.) ou (decomp.))
#   '491.0-493.0 K(dec.)'      → milieu              (plage K avec décomposition)
#   '540 deg.C decomp.'     → val+273.15  (°C + 'decomp.' sans parenthèses)
#   '234deg.C' / '-51deg.C' / 'at 313deg.C' → val+273.15  (°C ±, préfixe 'at' optionnel → K)
#   '618.9-619.4 K.(dec.)' → milieu   (plage K avec point après K et décomposition)
#   '220-222 deg. C'       → milieu + 273.15  (plage °C avec espace entre deg. et C)
#   '313-343 K decomposition' → milieu        (plage K avec décomposition en toutes lettres)
#   'greater than …' / '> …'  → None          (_QUALITATIVE : borne inférieure non précise)
#   'sublimes' / 'sublimation' / '(subl.)' → None  (_QUALITATIVE : préfixe 'subl', toutes variantes)
#   '… (sub.) …'                      → None  (_QUALITATIVE : sublimation, notation (sub.))
#   'approximately …' / 'about …' / 'around …' / 'ca. …' → None  (_QUALITATIVE : valeur approx.)
# ============================================================

SOURCE_FILE = '../csd_all.csv'
COLUMN = 'melting_point'
CLEAN_COLUMN = 'melting_point_k'

# Valeurs qualitatives / non convertibles, attendues à None (à enrichir si besoin) :
_QUALITATIVE = (
    'decompos',
    'above',
    'below',
    'over',
    'greater than',
    '>',
    'subl',
    '(sub.',
    'approximately',
    'about',
    'around',
    'ca.',
    'less than',
    '',
)


def parse_value(raw):
    """Convertit une valeur brute en valeur propre, ou None si le format n'est pas (encore) pris en charge."""
    if pd.isna(raw):
        return None
    s = str(raw).strip()
    if s.lower() in ('nan', 'none', '', '?'):
        return None
    if any(q in s.lower() for q in _QUALITATIVE):
        return None

    # --- AJOUTER ICI une règle par format reconnu ---
    # '528' → val K  (valeur numérique seule sans unité, traitée comme K)
    m = re.match(r'^([\d.]+)$', s)
    if m:
        return float(m.group(1))

    # '485 K' / '558.0(1) K' → val K  (incertitude cristallographique optionnelle)
    m = re.match(r'^([\d.]+)(?:\(\d+\))?\s*K$', s)
    if m:
        return float(m.group(1))

    # '345deg. K' / '548 deg. K (dec.)' → val K
    m = re.match(r'^([\d.]+)\s*deg\.\s*K(\s*\([^)]+\))?$', s, re.IGNORECASE)
    if m:
        return float(m.group(1))

    # '325-327 °deg.C' → midpoint + 273.15  (variante °deg.C)
    m = re.match(r'^([\d.]+)\s*-\s*([\d.]+)\s*°\s*deg\.\s*C\.?$', s, re.IGNORECASE)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2 + 273.15

    # 'at 165-175deg.C' / '190-192 deg.C' / '220-222 deg. C' → midpoint + 273.15
    m = re.match(r'^(?:at\s+)?([\d.]+)\s*-\s*([\d.]+)\s*deg\.\s*C\.?$', s, re.IGNORECASE)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2 + 273.15

    # '491.0-493.0 K(dec.)' → midpoint
    m = re.match(r'^([\d.]+)\s*-\s*([\d.]+)\s*K\s*\(dec\.?\)$', s, re.IGNORECASE)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2

    # '478(dec.) K' → val K  (annotation décomp. avant l'unité)
    m = re.match(r'^([\d.]+)\s*\(dec(?:omp)?\.?\)\s*K$', s, re.IGNORECASE)
    if m:
        return float(m.group(1))

    # '519 K(dec.)' / '510K.(dec.)' / '690 K (decomp.)' → val
    m = re.match(r'^([\d.]+)\s*K\.?\s*\(dec(?:omp)?\.?\)$', s, re.IGNORECASE)
    if m:
        return float(m.group(1))

    # '178-181deg.C(dec.)' / '199-200deg.C(lit.)' → midpoint + 273.15  (plage °C + annotation)
    m = re.match(r'^([\d.]+)\s*-\s*([\d.]+)\s*deg\.\s*C\.?\s*\([^)]+\)$', s, re.IGNORECASE)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2 + 273.15

    # '72deg.C(dec.)' / 'from 175deg.C(dec.)' / 'at 180deg.C(dec.)' → val+273.15
    m = re.match(r'^(?:(?:from|at)\s+)?(-?[\d.]+)\s*deg\.\s*(?:C\.?\s*)?\(dec(?:omp)?\.?\)$', s, re.IGNORECASE)
    if m:
        return float(m.group(1)) + 273.15

    # 'decomp. 330 deg.C' → val + 273.15  (préfixe decomp. + °C)
    m = re.match(r'^decomp\.\s*(-?[\d.]+)\s*deg\.\s*C\.?$', s, re.IGNORECASE)
    if m:
        return float(m.group(1)) + 273.15

    # '540 deg.C decomp.' → val + 273.15  (°C + décomp. sans parenthèses)
    m = re.match(r'^([\d.]+)\s*deg\.\s*C\.?\s*decomp\.?$', s, re.IGNORECASE)
    if m:
        return float(m.group(1)) + 273.15

    # '234deg.C' / '97.5deg. C' / '-51deg.C' / 'at 313deg.C' → val + 273.15
    m = re.match(r'^(?:at\s+)?(-?[\d.]+)\s*deg\.\s*C\.?$', s, re.IGNORECASE)
    if m:
        return float(m.group(1)) + 273.15

    # '618.9-619.4 K.(dec.)' → midpoint  (plage K avec point après K et décomposition)
    m = re.match(r'^([\d.]+)\s*-\s*([\d.]+)\s*K\.\s*\(dec\.?\)$', s, re.IGNORECASE)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2

    # '683 K decomp.' → val  (K + decomp. sans parenthèses)
    m = re.match(r'^([\d.]+)\s*K\s+decomp\.?$', s, re.IGNORECASE)
    if m:
        return float(m.group(1))

    # '313-343 K decomposition' → midpoint  (plage K avec décomposition en toutes lettres)
    m = re.match(r'^([\d.]+)\s*-\s*([\d.]+)\s*K\s+decomposition$', s, re.IGNORECASE)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2

    # '433-435' → midpoint K  (plage sans unité, traitée comme K)
    m = re.match(r'^([\d.]+)\s*-\s*([\d.]+)$', s)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2

    # '419-420 deg. K' → midpoint  (plage K avec notation deg. K)
    m = re.match(r'^([\d.]+)\s*-\s*([\d.]+)\s*deg\.\s*K$', s, re.IGNORECASE)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2

    # '327-328 K' / '441- 443 K' → midpoint
    m = re.match(r'^([\d.]+)\s*-\s*([\d.]+)\s*K$', s)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2

    return None


def _is_expected_none(raw):
    """True si la valeur est légitimement non convertible (vide, manquante, qualitative)."""
    if pd.isna(raw):
        return True
    s = str(raw).strip()
    return s.lower() in ('nan', 'none', '', '?') or any(q in s.lower() for q in _QUALITATIVE)


def main():
    print(f"Chargement de {SOURCE_FILE}...")
    df = pd.read_csv(SOURCE_FILE, low_memory=False)
    col = df[COLUMN]

    print(f"Analyse ligne par ligne de '{COLUMN}'...")
    # Parcours ligne par ligne : on s'ARRÊTE à la PREMIÈRE valeur non reconnue.
    for i, raw in enumerate(col):
        if _is_expected_none(raw):
            continue
        if parse_value(raw) is None:
            print(f"\n⛔ Format non reconnu (ligne {i}) :")
            print(f"   {raw!r}")
            print("\n→ Ajoute UNE règle pour ce format dans parse_value(), documente-le en tête, puis relance.")
            print("   (Si la valeur n'est pas convertible : ajoute-la à _QUALITATIVE.)")
            sys.exit(1)

    # La boucle est allée au bout sans rien rencontrer d'inconnu : tout est reconnu.
    print("✅ Tous les formats sont reconnus. Sauvegarde...")
    results = col.map(parse_value)
    if CLEAN_COLUMN in df.columns:
        df = df.drop(columns=[CLEAN_COLUMN])
    idx = df.columns.get_loc(COLUMN) + 1
    df.insert(idx, CLEAN_COLUMN, results)
    df.to_csv(SOURCE_FILE, index=False)
    print(f"✅ Terminé ! Colonne '{CLEAN_COLUMN}' ajoutée dans {SOURCE_FILE}.")


if __name__ == '__main__':
    main()
