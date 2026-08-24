import os
import sys
import re
import pandas as pd

# ============================================================
# Formats pris en charge :
#   - couleurs de base directes          : 'red', 'blue', 'colorless', ...
#   - synonymes                          : 'violet'→purple, 'gold'→yellow, 'grey'→gray, ...
#   - qualificatifs strippés             : 'pale yellow'→yellow, 'dark red'→red, ...
#   - suffixe -ish                       : 'yellowish'→yellow, 'reddish'→red, ...
#   - couleurs composées (tiret/espace)  : 'red-brown'→red/brown, 'blue green'→blue/green
#   - séparateurs 'to'/'and'             : 'colorless to yellow'→colorless/yellow
#   - qualificatifs seuls sans couleur   : 'pale', 'slightly colored'→colorless
# ============================================================

SOURCE_FILE = '../csd_all.csv'
COLUMN      = 'color'
CLEAN_COLUMN = 'color_clean'

_WHITELIST = {
    'colorless', 'white', 'black', 'gray',
    'red', 'orange', 'yellow', 'green', 'blue', 'purple', 'pink', 'brown',
}

_SYNONYMS = {
    # colorless
    'colourless':  'colorless',
    'transparent': 'colorless',
    'clear':       'colorless',
    # gray
    'grey':        'gray',
    'silver':      'gray',
    'silvery':     'gray',
    # purple
    'violet':      'purple',
    'mauve':       'purple',
    'lilac':       'purple',
    'lavender':    'purple',
    'indigo':      'purple',
    # yellow
    'gold':        'yellow',
    'golden':      'yellow',
    'lemon':       'yellow',
    'straw':       'yellow',
    # orange
    'amber':       'orange',
    # brown
    'beige':       'brown',
    'tan':         'brown',
    'khaki':       'brown',
    # white
    'cream':       'white',
    'ivory':       'white',
    # red
    'crimson':     'red',
    'scarlet':     'red',
    'vermilion':   'red',
    'burgundy':    'red',
    'maroon':      'red',
    'cherry':      'red',
    'brick':       'red',
    'copper':      'red',
    'bronze':      'brown',
    'rust':        'red',
    # formes -ish dont le strip naïf ([:-3]) ne donne pas la bonne base
    'reddish':     'red',
    'orangish':    'orange',
    'bluish':      'blue',
    'purplish':    'purple',
    'whitish':     'white',
    # couleurs inhabituelles
    'rufous':      'red',
    'carmine':     'red',
    'wine':        'red',
    'aqua':        'blue',
    'lightyellow': 'yellow',
    'lightblue':   'blue',
    'lightgreen':  'green',
    'lightbrown':  'brown',
    'lightorange': 'orange',
    'lightpink':   'pink',
    'lightpurple': 'purple',
    'darkyellow':  'yellow',
    'darkorange':  'orange',
    'darkpurple':  'purple',
    'darkred':     'red',
    'darkblue':    'blue',
    'darkgreen':   'green',
    'darkbrown':   'brown',
    'darkblack':   'black',
    'darkwhite':   'white',
    'darkpink':    'pink',
    'darkgray':    'gray',
    'darkgrey':    'gray',
    'fawn':        'brown',
    'tangerine':   'orange',
    'flaxen':      'yellow',
    'apricot':     'orange',
    'ochre':       'yellow',
    'ocher':       'yellow',
    'olivine':     'green',
    'amaranthine': 'purple',
    'deepred':     'red',
    'deepblue':    'blue',
    'deepgreen':   'green',
    'deepbrown':   'brown',
    'deeppurple':  'purple',
    'taupe':       'brown',
    'plum':        'purple',
    'caramel':     'brown',
    'glaucous':    'blue',
    'garnet':      'red',
    'mazarine':    'blue',
    'chartreuse':  'yellow/green',
    'jasmine':     'white',
    'claret':      'red',
    'bordeaux':    'purple/red',
    'atrovirens':  'green',
    'mustard':     'yellow',
    'honey':       'yellow',
    'chocolate':   'brown',
    'saffron':     'yellow',
    'flesh':       'pink',
    'colourles':   'colorless',
    'grape':       'purple',
    'tawny':       'brown/orange',
    'wheat':       'yellow',
    # pink
    'magenta':     'pink',
    'rose':        'pink',
    'salmon':      'pink',
    'fuchsia':     'pink',
    # blue
    'navy':        'blue',
    'cobalt':      'blue',
    'cyan':        'blue',
    'turquoise':   'blue',
    'azure':       'blue',
    # green
    'olive':       'green',
    'lime':        'green',
    'teal':        'green',
    'emerald':     'green',
    'jade':        'green',
    # verdure / nature
    'grass':       'green',
    'forest':      'green',
    # canari / jaune nature
    'canary':      'yellow',
    'primrose':    'yellow',
    # ex-qualificatifs qui sont en réalité des couleurs
    'cognac':      'brown',
    'vine':        'green',
    'rusty':       'red',
    'fire':        'red/orange',
    'metallic':    'gray',
    # batch 13
    'ruby':        'red',
    'vermillion':  'red',
    'coloress':    'colorless',
    'apple':       'green',
    'colouless':   'colorless',
    'maple':       'brown/orange',
    'celadon':     'green',
    'jacinth':     'red/orange',
    'sorrel':      'brown/red',
    'kermesinus':  'red',
    'lightred':    'red',
    'peach':       'pink/orange',
    'farblos':     'colorless',
    'palegreen':   'green',
    'nacarat':     'red/orange',
    'claybank':    'brown/yellow',
    'topaz':       'yellow/brown',
    'palebrown':   'brown',
    'grassy':      'green',
    'coffee':      'brown',
    'raspberry':   'red/pink',
    # couleurs composées (mot + couleur de base) — noms de couleurs spécifiques
    'bottle green':  'green',
    'steel blue':    'blue',
    'steel gray':    'gray',
    'steel grey':    'gray',
    'ink blue':      'blue',
    'ink black':     'black',
    'moss green':    'green',
    'laurel green':  'green',
    'slate blue':    'blue',
    'slate gray':    'gray',
    'slate grey':    'gray',
    'royal blue':    'blue',
    'royal purple':  'purple',
    'sea green':     'green',
    'sea blue':      'blue',
    'chrome yellow': 'yellow',
    'chrome green':  'green',
    'sapphire blue': 'blue',
    # mots seuls dont la couleur est non ambiguë sans contexte
    'bottle':    'green',
    'steel':     'gray',
    'ink':       'black',
    'moss':      'green',
    'laurel':    'green',
    'slate':     'gray',
    'royal':     'blue',
    'sapphire':  'blue',
    # 'sea' et 'chrome' seuls restent inconnus (trop ambigus sans la couleur associée)
    # couleurs composées ex-qualificatifs : water, palm, electric, turkish, ice
    'water blue':    'blue',
    'water green':   'green',
    'electric blue': 'blue',
    'electric green':'green',
    'turkish blue':  'blue',
    'turkish green': 'green',
    'ice blue':      'blue',
    'ice white':     'white',
    'palm green':    'green',
    # couleurs composées supplémentaires
    'off-white':   'white',
    'off white':   'white',
    'blood red':   'red',
    'sky blue':    'blue',
    'blood':       'red',
    'sky':         'blue',
    'buff':        'yellow/brown',
    'aquamarine':  'blue/green',
    'paleyellow':  'yellow',
    'sea-green':   'green',
    'sea green':   'green',
    'chrome-yellow': 'yellow',
    'blackberry':  'purple/black',
    'jujube':      'red/brown',
    'jujube red':  'red/brown',
    'sea-blue':    'blue',
    'sea blue':    'blue',
    'celery':      'green',
    'browm':       'brown',
    'brawn':       'brown',
    'reseda':      'green',
    'bule':        'blue',
    'rod':         'red',
    'lutescent':   'yellow',
    'limpid':      'colorless',
    'snow':        'white',
    'snow-white':  'white',
    'snow white':  'white',
    'courless':    'colorless',
    'carnation':   'pink',
    'anthracite':  'gray',
    'flame':       'red/orange',
    'flame-red':   'red',
    'flame red':   'red',
    'colorles':    'colorless',
    'colourlos':   'colorless',
    'coral':       'red/orange',
    'cinnabar':    'red',
    'colorle':     'colorless',
    'multicoloured': None,
    'multicolored':  None,
    'mint':        'green',
    'mint green':  'green',
    'peony':       'pink/red',
    'brass':       'yellow/brown',
    'roux':        'brown/red',
    'ultramarine': 'blue',
    'henna':       'brown/red',
    'redbrown':    'red/brown',
    'carrot':      'orange',
    'malachite':   'green',
    'flavescens':  'yellow',
    'colurless':   'colorless',
    'cuticolor':   None,
    'coulorless':  'colorless',
    'intenseorange': 'orange',
    'lyons':       'blue',
    'yellowgreen': 'yellow/green',
    'nut':         'brown',
    'lemonade':    'yellow',
    'achromatic':  'colorless',
    'haki':        'brown',
    'rosa':        'pink',
    'colerless':   'colorless',
    'colour less': 'colorless',
    'color less':  'colorless',
    'v':           None,
    'redblue':     'red/blue',
    'sand':        'yellow/brown',
    'mulberry':    'purple/red',
    'chrome':      None,
    'pumpkin':     'orange',
    'ginger':      'brown/red',
    'yelllow':     'yellow',
    'harki':       'brown',
    'mossgreen':   'green',
    'fulvous':     'yellow/brown',
    'coloreless':  'colorless',
    'turqouise':   'blue',
    'jaune':       'yellow',
    'sea wave':    'blue/green',
    'organ':       'orange',
    'çolorless':   'colorless',
    'bluegreen':   'blue/green',
    'clourless':   'colorless',
    'yelloq':      'yellow',
    'darkcyan':    'blue',
    'british racing green': 'green',
    'british':     'green',
    'organge':     'orange',
    'redness':     'red',
    'creamy':      'white',
    'brwon':       'brown',
    'sundown':     'orange/pink',
    'sepia':       'brown',
    'jonquil':     'yellow',
    'bordo':       'purple/red',
    'amythyst':    'purple',
    'amethyst':    'purple',
    'lettuce':     'green',
    'bown':        'brown',
    'blck':        'black',
    'less':        'colorless',
    'turkish':     'blue',
    'sangria':     'red/purple',
    'courlorless': 'colorless',
    'greem':       'green',
    'tea':         'yellow/brown',
    'tea yellow':  'yellow/brown',
    'midnight':    'blue',
    'midnight blue': 'blue',
    'midnight-blue': 'blue',
    'organe':      'orange',
    'cerise':      'red',
    'aubergine':   'purple',
    'loden':       'green',
    'loden green': 'green',
    'citron':      'yellow',
    'citron yellow': 'yellow',
    'mahogany':    'brown/red',
    'greeny':      'green',
    'orangey':     'orange',
    'madder':      'red',
    'selsel':      None,
    'metal grey':  'gray',
    'metal gray':  'gray',
    'metal':       'gray',
    'almond':      'brown',
    'sap':         'green',
    'sap green':   'green',
    'pearl blue':  'blue',
    'pearl':       'white',
    'kaki':        'brown',
    'coluorless':  'colorless',
    'orang':       'orange',
    'turquiose':   'blue',
    'jet':         'black',
    'jet-black':   'black',
    'jet black':   'black',
    'puple':       'purple',
    'marsh':       'green/brown',
    'orchid':      'pink/purple',
    'graphite':    'gray',
    'permanganate': 'purple',
    'poppy':       'red',
    'res':         'red',
    'colourness':  'colorless',
    'geeen':       'green',
    'greenblack':  'green/black',
    'vinous':      'red/purple',
    'ebony':       'black',
    'yellowless':  'yellow',
    'melon':       'orange/yellow',
    'sage':        'green',
    'sage green':  'green',
    'cololess':    'colorless',
    'bourdeuxe':   'purple/red',
    'yelowish':    'yellow',
    'corless':     'colorless',
    'rot':         'red',
    'oreange':     'orange',
    'turquise':    'blue',
    'velloy':      'yellow',
    'colouress':   'colorless',
    'balck':       'black',
    'pine':        'green',
    'pine green':  'green',
    'myrtle':      'green',
    'gree':        'green',
    'rubylith':    'red',
    'limegreen':   'green',
    'port':        'red',
    'grenadine':   'red/pink',
    'yello':       'yellow',
    'plue':        'blue',
    'maize':       'yellow',
    'granite':     'gray',
    'rosy':        'pink',
    'colourl':     'colorless',
    'pinky':       'pink',
    'borwn':       'brown',
    'coloueless':  'colorless',
    'bluesh':      'blue',
    'bistre':      'brown',
    'peacock':     'green',
    'olourless':   'colorless',
    'yelow':       'yellow',
    'modena':      'purple/red',
    'puce':        'brown/pink',
    'pista':       'green',
    'metalliclightoran': 'orange',
    'unk':         None,
    'kelly':       'green',
    'kelly green': 'green',
    'liver':       'brown/red',
    'petrol':      'blue',
    'petrol blue': 'blue',
    'cornflower':  'blue',
    'cornflower blue': 'blue',
    'yeloow':      'yellow',
    'willow':      'green',
    'willow green': 'green',
    'sandy':       'yellow/brown',
    'marine':      'blue',
    'marine blue': 'blue',
    'chestnut':    'brown',
    'voilet':      'purple',
    'sulfur':      'yellow',
    'sulphur':     'yellow',
    'cinnamon':    'brown/orange',
}

# Mots à ignorer : uniquement les modificateurs d'intensité/de brillance/de forme
# qui ne donnent AUCUNE information sur la couleur, quelle que soit la couleur à côté.
_QUALIFIERS_RE = re.compile(
    r'\b(pale|light|dark|deep|bright|vivid|faint|intense|dull|faded|slightly|'
    r'very|quite|strong|strongly|weak|brilliant|medium|moderate|nearly|almost|rich|intensely|extremely|baby|hot|electric|vibrant|acid|'
    r'dirty|dusty|muted|powdery|shiny|lustrous|neon|fluorescent|iridescent|sparkling|shining|shifted|streaky|'
    r'colored|coloured|colour|color|slight|faintly|dichroic|dichroism|pleochroism|pleochromism|pleochromic|tribochromic|thermochromic|piezochromic|trichroic|trichoric|dichromic|pleochrome|like|fragment|prism|off|milk|milky|plae|shine|luster|lustre|lusterous|drab|pleochroic|luminescent|translucent|transluscent|tranclucent|transparent|transparant|opaque|sauce|ht|lt|from|smokey|smoky|sligthly|earthy|glow|pales|pallid|palid|disrupted|indian|gloss|tinged|half|inclusion|irregular|right|ligh|'
    r'with|hue|tinge|tint|tinted|shade|'
    r'crystalline|microcrystalline|crystal|cristal|crystall|crystals|powder|solid|amorphous|'
    r'waxy|oily|glassy|needle|needles|plate|plates|plated|block|blocks|prismatic|hexagonal|rectangle|showing|wedge|shard|sheet|ligth|'
    r'under|polarized|polarised|unpolarized)\b'
)

_UNKNOWN = object()


def _resolve_token(token):
    """Résout un token individuel → couleur canonique, None (ignorable) ou _UNKNOWN."""
    t = token.strip()
    if not t:
        return None
    # Ignorer les tokens purement non-alphabétiques (ponctuation résiduelle, chiffres)
    if not any(c.isalpha() for c in t):
        return None
    if t in _WHITELIST:
        return t
    if t in _SYNONYMS:
        return _SYNONYMS[t]
    # suffixe -ish : reddish → red
    if t.endswith('ish'):
        base = t[:-3]
        if base in _WHITELIST:
            return base
        if base in _SYNONYMS:
            return _SYNONYMS[base]
    return _UNKNOWN


def parse_value(raw):
    if pd.isna(raw):
        return None
    s = str(raw).strip().lower()
    if s in ('nan', 'none', '', '?'):
        return None

    # Descriptions textuelles longues → pas une couleur
    if len(s) > 60:
        return None

    # Vérification directe avant tout traitement
    if s in _WHITELIST:
        return s
    if s in _SYNONYMS:
        v = _SYNONYMS[s]
        return '/'.join(sorted(v.split('/'))) if v else v

    # Supprimer les annotations entre parenthèses : 'green(dark)' → 'green'
    s = re.sub(r'\([^)]*\)', '', s).strip()

    # Supprimer les tokens numériques parasites : '0.11 light colorless' → 'light colorless'
    s = re.sub(r'\b\d+[\d.,]*\b', '', s).strip()

    # Supprimer ponctuation parasite : 'green?' → 'green', guillemets typo → ''
    # Ne pas supprimer les virgules ici — elles sont normalisées en '/' plus bas
    s = re.sub(r'[?!.;:`\'\"]+', '', s).strip()

    # Normaliser les séparateurs de couleurs composées : 'to', 'and', 'or', ',' → '/'
    s = re.sub(r'\s*,\s*', '/', s)
    s = re.sub(r'\s+to\s+|\bto\b', '/', s)
    s = re.sub(r'\s+and\s+|\band\b', '/', s)
    s = re.sub(r'\s+or\s+|\bor\b', '/', s)

    # Supprimer les qualificatifs (pale, dark, light, ...)
    s = _QUALIFIERS_RE.sub('', s)
    s = re.sub(r'\s+', ' ', s).strip()

    # Si plus rien après le stripping → était uniquement des qualificatifs → colorless
    if not s:
        return 'colorless'

    # Vérification après stripping qualificatifs
    if s in _WHITELIST:
        return s
    if s in _SYNONYMS:
        v = _SYNONYMS[s]
        return '/'.join(sorted(v.split('/'))) if v else v

    # Séparer les couleurs composées : tirets, espaces, slashes, virgules
    parts = re.split(r'[\-/,\s]+', s)
    colors = []
    seen = set()
    for part in parts:
        part = part.strip()
        if not part:
            continue
        resolved = _resolve_token(part)
        if resolved is _UNKNOWN:
            return _UNKNOWN
        if resolved is not None:
            for c in resolved.split('/'):
                if c not in seen:
                    seen.add(c)
                    colors.append(c)

    if not colors:
        return None

    return '/'.join(sorted(colors))


def _is_expected_none(raw):
    if pd.isna(raw):
        return True
    s = str(raw).strip().lower()
    return s in ('nan', 'none', '', '?')


def _find_unknown_token(raw):
    """Retourne le token non reconnu pour l'afficher dans le message d'erreur."""
    s = str(raw).strip().lower()
    s = re.sub(r'\s+to\s+', '/', s)
    s = re.sub(r'\s+and\s+', '/', s)
    s = _QUALIFIERS_RE.sub('', s).strip()
    s = re.sub(r'\s+', ' ', s)
    for part in re.split(r'[\-/\s]+', s):
        part = part.strip()
        if not part:
            continue
        if _resolve_token(part) is _UNKNOWN:
            return part
    return raw


def _progress_file():
    return f'.progress_{COLUMN}.txt'


def _origin_file():
    return f'.progress_{COLUMN}_origin.txt'


def main():
    restart = 'restart' in sys.argv[1:]
    progress = _progress_file()

    if restart:
        for f in (progress, _origin_file()):
            if os.path.exists(f):
                os.remove(f)

    start = 0
    if not restart and os.path.exists(progress):
        try:
            start = int(open(progress).read().strip())
        except ValueError:
            start = 0

    print(f"Lecture de la colonne '{COLUMN}'...")
    col = pd.read_csv(SOURCE_FILE, usecols=[COLUMN], low_memory=False)[COLUMN]

    # Travailler sur valeurs uniques triées par fréquence : les inconnues
    # les plus impactantes remontent en premier.
    vc = col.dropna().value_counts()
    unique_vals = list(vc.index)
    n = len(unique_vals)
    n_total_rows = col.notna().sum()

    if start >= n:
        start = 0

    print(f"Analyse de {n:,} valeurs uniques à partir de l'index {start}...\n")

    BATCH_SIZE = 5
    errors = []

    for i in range(start, n):
        raw = unique_vals[i]
        if _is_expected_none(raw):
            continue
        result = parse_value(raw)
        if result is _UNKNOWN:
            unknown = _find_unknown_token(raw)
            if unknown and any(e[2] == unknown for e in errors):
                continue
            errors.append((i, raw, unknown))
            if len(errors) >= BATCH_SIZE:
                break

    if errors:
        with open(progress, 'w') as f:
            f.write(str(errors[0][0]))

        # Mémoriser l'index de départ de la zone inconnue (écrit une seule fois)
        origin_file = _origin_file()
        if not os.path.exists(origin_file):
            with open(origin_file, 'w') as f:
                f.write(str(errors[0][0]))
        origin = int(open(origin_file).read().strip())

        done = errors[0][0] - origin
        total_unknown = n - origin
        pct = done / total_unknown * 100 if total_unknown else 0

        print(f"Progression zone inconnue : {done}/{total_unknown} valeurs uniques traitées ({pct:.1f}%)\n")

        for idx, (i, raw, unknown) in enumerate(errors, 1):
            print(f"  ⛔ [{idx}/{len(errors)}]  {raw!r}  ({vc[raw]:,} occurrences)")
            print(f"         token inconnu : {unknown!r}\n")

        print("→ Ajoute dans _SYNONYMS (variante/synonyme) ou _QUALIFIERS_RE (mot à ignorer)")
        print("→ Relance ensuite (reprend ici). Passer 'restart' pour tout réanalyser.")
        sys.exit(1)

    # Tout reconnu → sauvegarder
    if os.path.exists(progress):
        os.remove(progress)
    if os.path.exists(_origin_file()):
        os.remove(_origin_file())
    print("✅ Toutes les valeurs reconnues. Lecture complète et sauvegarde...")
    df = pd.read_csv(SOURCE_FILE, low_memory=False)

    mapping = {v: (None if (r := parse_value(v)) is _UNKNOWN else r) for v in vc.index}
    results = df[COLUMN].map(lambda x: mapping.get(x) if pd.notna(x) else None)

    if CLEAN_COLUMN in df.columns:
        df = df.drop(columns=[CLEAN_COLUMN])
    idx = df.columns.get_loc(COLUMN) + 1
    df.insert(idx, CLEAN_COLUMN, results)
    df.to_csv(SOURCE_FILE, index=False)
    print(f"✅ Terminé ! Colonne '{CLEAN_COLUMN}' ajoutée dans {SOURCE_FILE}.")


if __name__ == '__main__':
    main()
