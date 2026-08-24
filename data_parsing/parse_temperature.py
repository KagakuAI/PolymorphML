import os
import sys
import re
import pandas as pd

# ============================================================
# Formats de température pris en charge (enrichi au fur et à mesure) :
#   '293 K'              → 293.0    (ex: '293 K')
#   'at 293 K'           → 293.0    (ex: 'at 398 K')
#   'at -135 deg.C'      → 138.15   (ex: 'at -135 deg.C')
#   'at -60deg.C.'       → 213.15   (ex: 'at -60deg.C.')  — pas d'espace, point final
#   '240-260 K'          → 250.0    (ex: '240-260 K')     — milieu de plage K
#   'at -50 to -70deg.C' → 213.15   (ex: 'at -50 to -70deg.C') — milieu de plage °C
#   'at about 120 K'     → 120.0    (ex: 'at about 120 K')
#   'at below 28 deg.C'  → 301.15   (ex: 'at below 28 deg.C')  — borne supérieure
#   'at above 310 K'     → 310.0    (ex: 'at above 310 K')
#   'at around 130 K'    → 130.0    (ex: 'at around 130 K')
#   'at 115.7K K'        → 115.7    (typo CSD, K collé au nombre)
#   'at 133 deg K'       → 133.0    (variante deg K)
#   'at low temperature' → None     (qualitatif, non convertible)
#   '?'                  → None     (donnée manquante)
# ============================================================

SOURCE_FILE = '../csd_all.csv'
COLUMN = 'temperature'
CLEAN_COLUMN = 'temperature_K'

# Valeurs qualitatives / non convertibles, attendues à None :
_QUALITATIVE = (
    'low temperature', 'high temperature', 'room temperature',
    'ambient temperature', 'liquid nitrogen', 'low temp',
    'not reported', 'no temp reported',
)
# liquid nitrogen ≈ 80 K
# ambient/room temperature ≈ 18/25 °C


def parse_value(raw):
    """Convertit une valeur brute en Kelvin, ou None si le format n'est pas (encore) pris en charge."""
    if pd.isna(raw):
        return None
    s = str(raw).strip()
    if s.lower() in ('nan', 'none', '', '?'):
        return None
    if any(q in s.lower() for q in _QUALITATIVE):
        return None

    # "293 K" / "293.5 K"
    m = re.match(r'^([\d.]+)\s*K$', s)
    if m:
        return float(m.group(1))

    # "at 293 K" / "at about 120 K" / "at above/below/around 310 K" / "at 115.7K K" / "at 133 deg K"
    m = re.match(r'^at\s+(?:about|above|below|around)?\s*([\d.]+)\s*K?\s*(?:deg\s+)?K$', s, re.IGNORECASE)
    if m:
        return float(m.group(1))

    # "at -135 deg.C" / "at -60deg.C." / "at below 28 deg.C"
    m = re.match(r'^(?:at\s+)?(?:about|below|above)?\s*([+-]?[\d.]+)\s*deg\.C\.?$', s, re.IGNORECASE)
    if m:
        return float(m.group(1)) + 273.15

    # "240-260 K" — plage K, milieu
    m = re.match(r'^([\d.]+)-([\d.]+)\s*K$', s)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2

    # "at -50 to -70deg.C" — plage °C, milieu converti
    m = re.match(r'^(?:at\s+)?([+-]?[\d.]+)\s+to\s+([+-]?[\d.]+)\s*deg\.C\.?$', s, re.IGNORECASE)
    if m:
        return ((float(m.group(1)) + float(m.group(2))) / 2) + 273.15

    return None


def _is_expected_none(raw):
    """True si la valeur est légitimement non convertible (vide, manquante, qualitative)."""
    if pd.isna(raw):
        return True
    s = str(raw).strip()
    return s.lower() in ('nan', 'none', '', '?') or any(q in s.lower() for q in _QUALITATIVE)


def _progress_file():
    """Petit fichier d'état qui mémorise la ligne où l'analyse s'est arrêtée."""
    return f'.progress_{COLUMN}.txt'


def main():
    # Passer 'restart' en argument pour ignorer l'état et tout réanalyser depuis le début.
    restart = 'restart' in sys.argv[1:]
    progress = _progress_file()

    start = 0
    if not restart and os.path.exists(progress):
        try:
            start = int(open(progress).read().strip())
        except ValueError:
            start = 0

    # Phase de mise au point : on ne lit QUE la colonne concernée (bien plus rapide).
    print(f"Lecture de la colonne '{COLUMN}'...")
    col = pd.read_csv(SOURCE_FILE, usecols=[COLUMN], low_memory=False)[COLUMN]
    n = len(col)
    if start >= n:
        start = 0

    print(f"Analyse ligne par ligne à partir de la ligne {start:,} / {n:,}...")
    # On s'ARRÊTE à la PREMIÈRE valeur non reconnue, et on mémorise sa position.
    for i in range(start, n):
        raw = col.iat[i]
        if _is_expected_none(raw):
            continue
        if parse_value(raw) is None:
            with open(progress, 'w') as f:
                f.write(str(i))
            print(f"\n⛔ Format non reconnu (ligne {i}) :")
            print(f"   {raw!r}")
            print("\n→ Ajoute UNE règle dans parse_value(), documente-la en tête, puis relance : il reprendra à cette ligne.")
            print("   (Valeur non convertible ? ajoute-la à _QUALITATIVE.)")
            sys.exit(1)

    # Tout est reconnu jusqu'au bout : on efface l'état, on relit le CSV complet et on sauvegarde.
    if os.path.exists(progress):
        os.remove(progress)
    print("✅ Tous les formats sont reconnus. Relecture complète et sauvegarde...")
    df = pd.read_csv(SOURCE_FILE, low_memory=False)
    results = df[COLUMN].map(parse_value)
    if CLEAN_COLUMN in df.columns:
        df = df.drop(columns=[CLEAN_COLUMN])
    idx = df.columns.get_loc(COLUMN) + 1
    df.insert(idx, CLEAN_COLUMN, results)
    df.to_csv(SOURCE_FILE, index=False)
    print(f"✅ Terminé ! Colonne '{CLEAN_COLUMN}' ajoutée dans {SOURCE_FILE}.")


if __name__ == '__main__':
    main()
