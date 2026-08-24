import os
import sys
import re
import pandas as pd

# ============================================================
# GABARIT DE RÉFÉRENCE — comportement fail-fast.
# Copier ce fichier en parse_<colonne>.py, régler COLUMN ci-dessous,
# puis enrichir UNIQUEMENT parse_value() (et _QUALITATIVE) au fil des retours.
# Ne JAMAIS modifier la logique de main() : elle doit s'arrêter au PREMIER
# format non reconnu, sans rien écrire tant que tout n'est pas pris en charge.
# ============================================================

# ------------------------------------------------------------
# Formats pris en charge (enrichi au fur et à mesure) :
#   (aucun pour l'instant — ajouter une ligne par format dès qu'une règle est créée)
# ============================================================

SOURCE_FILE = '../csd_all.csv'
COLUMN = 'À_REMPLIR'              # nom de la colonne source à nettoyer
CLEAN_COLUMN = f'{COLUMN}_clean'  # nom de la nouvelle colonne propre

# Valeurs qualitatives / non convertibles, attendues à None (à enrichir si besoin) :
_QUALITATIVE = (
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
    # Exemple de structure (à adapter / supprimer) :
    # m = re.match(r'^([\d.]+)\s*UNITE$', s)
    # if m:
    #     return float(m.group(1))

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
