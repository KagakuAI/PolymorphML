#!/usr/bin/env python3
"""
analyse_whitelist.py — Détection de synonymes dans WHITELIST via PubChem CID

Usage :
    python analyse_whitelist.py solvent_data.py

Sorties :
    whitelist_cache.json   — cache des résultats PubChem (permet de reprendre si coupé)
    whitelist_report.txt   — rapport lisible
    whitelist_data.json    — données complètes (synonymes, renommages, non résolus)
"""

import ast
import json
import time
import sys
import os
import urllib.request
import urllib.parse
import urllib.error
from collections import defaultdict

# ── Paramètres ────────────────────────────────────────────────────────────────
DELAY        = 0.22          # délai entre requêtes (~4.5 req/s, limite PubChem = 5/s)
MAX_RETRIES  = 3             # tentatives par composé avant abandon
CACHE_FILE   = "whitelist_cache.json"
REPORT_FILE  = "whitelist_report.txt"
DATA_FILE    = "whitelist_data.json"


# ── Extraction ─────────────────────────────────────────────────────────────────
def extract_from_script(path):
    """Extrait WHITELIST et SYNONYMS depuis solvent_data.py via ast."""
    with open(path, encoding="utf-8") as f:
        source = f.read()
    tree = ast.parse(source)
    whitelist, synonyms = None, None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    if target.id == "WHITELIST":
                        whitelist = ast.literal_eval(node.value)
                    elif target.id == "SYNONYMS":
                        synonyms = ast.literal_eval(node.value)
    if whitelist is None:
        raise ValueError("WHITELIST introuvable dans le fichier.")
    return sorted(whitelist), synonyms or {}


# ── PubChem ────────────────────────────────────────────────────────────────────
def fetch_pubchem(name):
    """
    Interroge PubChem PUG REST pour un nom de composé.
    Retourne dict {"cid": int|None, "iupac": str|None, "smiles": str|None}.
    """
    encoded = urllib.parse.quote(name, safe="")
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/property/IUPACName,CanonicalSMILES/JSON"
    )
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "whitelist-checker/1.0 (research)"}
            )
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read())
            props = data.get("PropertyTable", {}).get("Properties", [])
            if props:
                p = props[0]
                return {
                    "cid":   p.get("CID"),
                    "iupac": p.get("IUPACName"),
                    "smiles":p.get("CanonicalSMILES"),
                }
            return {"cid": None, "iupac": None, "smiles": None}

        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"cid": None, "iupac": None, "smiles": None}
            if e.code in (429, 503) and attempt < MAX_RETRIES - 1:
                wait = 2 ** (attempt + 1)
                print(f"\n  ⚠ Rate-limit ({e.code}) — attente {wait}s...", end="")
                time.sleep(wait)
                continue
            return {"cid": None, "iupac": None, "smiles": None}

        except Exception as exc:
            if attempt < MAX_RETRIES - 1:
                time.sleep(1)
                continue
            print(f"\n  ✗ Erreur pour '{name}': {exc}")
            return {"cid": None, "iupac": None, "smiles": None}

    return {"cid": None, "iupac": None, "smiles": None}


# ── Barre de progression ───────────────────────────────────────────────────────
def progress(i, total, name, found):
    pct  = i / total
    fill = int(pct * 40)
    bar  = "█" * fill + "░" * (40 - fill)
    tag  = "✓" if found else "✗"
    label = (name[:38] + "..") if len(name) > 40 else name
    print(f"\r[{bar}] {i:>3}/{total} {tag} {label:<42}", end="", flush=True)


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) < 2:
        print("Usage: python analyse_whitelist.py <parse_solvent.py>")
        sys.exit(1)

    source_path = sys.argv[1]
    print(f"\n{'='*65}")
    print(f"  Analyse whitelist — {source_path}")
    print(f"{'='*65}\n")

    # 1. Extraction
    print("① Extraction de _WHITELIST et _SYNONYMS...")
    entries, synonyms = extract_from_script(source_path)
    print(f"   Whitelist : {len(entries)} entrées")
    print(f"   Synonyms  : {len(synonyms)} clés\n")

    # 2. Cross-check statique
    conflicts = [(e, synonyms[e]) for e in entries if e in synonyms]

    # 3. Chargement du cache (reprise si coupé)
    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE) as f:
            cache = json.load(f)
        already = sum(1 for e in entries if e in cache)
        print(f"② Cache trouvé ({already}/{len(entries)} déjà résolus) — reprise.")
    else:
        print("② Aucun cache — démarrage complet.")

    # 4. Requêtes PubChem
    to_fetch = [e for e in entries if e not in cache]
    if to_fetch:
        eta = len(to_fetch) * DELAY
        print(f"\n   {len(to_fetch)} composés à interroger (~{eta:.0f}s)\n")
        for i, name in enumerate(to_fetch, 1):
            result = fetch_pubchem(name)
            cache[name] = result
            progress(i, len(to_fetch), name, result["cid"] is not None)
            # Sauvegarde du cache toutes les 50 entrées
            if i % 50 == 0:
                with open(CACHE_FILE, "w") as f:
                    json.dump(cache, f, ensure_ascii=False)
            time.sleep(DELAY)
        # Sauvegarde finale du cache
        with open(CACHE_FILE, "w") as f:
            json.dump(cache, f, ensure_ascii=False)
        print(f"\n\n   Cache sauvegardé → {CACHE_FILE}\n")
    else:
        print("   Tous les composés sont en cache, pas de requêtes.\n")

    # 5. Analyse des résultats
    by_cid      = defaultdict(list)
    not_found   = []

    for name in entries:
        info = cache.get(name, {"cid": None, "iupac": None, "smiles": None})
        if info["cid"] is not None:
            by_cid[info["cid"]].append(name)
        else:
            not_found.append(name)

    # Groupes synonymes (même CID, ≥2 entrées whitelist)
    synonym_groups = sorted(
        [(cid, names) for cid, names in by_cid.items() if len(names) > 1],
        key=lambda x: -len(x[1]),
    )

    # Suggestions de renommage (nom whitelist ≠ nom IUPAC PubChem)
    rename_suggestions = sorted(
        [
            (name, cache[name]["iupac"], cache[name]["cid"])
            for name in entries
            if name in cache
            and cache[name]["cid"] is not None
            and cache[name]["iupac"]
            and cache[name]["iupac"].lower() != name.lower()
        ],
        key=lambda x: x[0],
    )

    # 6. Rapport
    lines = []
    sep   = "─" * 65

    def h(title):
        lines.append("")
        lines.append(sep)
        lines.append(f"  {title}")
        lines.append(sep)

    lines.append("=" * 65)
    lines.append("  RAPPORT — ANALYSE WHITELIST parse_solvent.py")
    lines.append("=" * 65)
    lines.append(f"  Entrées analysées  : {len(entries)}")
    lines.append(f"  Résolus PubChem    : {len(entries) - len(not_found)}")
    lines.append(f"  Non résolus        : {len(not_found)}")
    lines.append(f"  Synonymes (col. A) : {len(synonym_groups)}")
    lines.append(f"  Renommages (col. B): {len(rename_suggestions)}")
    lines.append(f"  Cross-check        : {len(conflicts)} ({sum(1 for _,t in conflicts if _ != t)} réel(s))")

    # A — Synonymes
    h("COL. A — SYNONYMES  (même CID PubChem → même molécule)")
    if synonym_groups:
        for cid, names in synonym_groups:
            iupac_ref = cache[names[0]].get("iupac", "?")
            lines.append(f"\n  CID {cid}  |  IUPAC : {iupac_ref}")
            for n in names:
                lines.append(f"    {'→' if n != names[0] else '✓'} {n}")
    else:
        lines.append("  Aucun synonyme détecté via CID.")

    # B — Renommages
    h(f"COL. B — RENOMMAGES SUGGÉRÉS  ({len(rename_suggestions)} entrées)")
    lines.append(f"  {'Nom actuel (whitelist)':<42}  →  Nom IUPAC PubChem")
    lines.append(f"  {'─'*42}  {'─'*30}")
    for name, iupac, cid in rename_suggestions:
        lines.append(f"  {name:<42}  →  {iupac}")

    # C — Non résolus
    h(f"NON RÉSOLUS  ({len(not_found)} entrées — vérifier via Wikipedia)")
    for n in not_found:
        lines.append(f"    - {n}")

    # D — Cross-check statique
    h("CROSS-CHECK STATIQUE — entrées whitelist présentes dans _SYNONYMS")
    for entry, target in conflicts:
        tag = "[CONFLIT RÉEL]" if entry != target else "[auto-référence]"
        arrow = f"  →  {target}" if entry != target else ""
        lines.append(f"  {tag:<18}  {entry}{arrow}")

    lines.append("")
    lines.append("=" * 65)

    report = "\n".join(lines)
    print(report)

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"\n  Rapport sauvegardé  → {REPORT_FILE}")

    # 7. JSON complet
    data = {
        "summary": {
            "total": len(entries),
            "resolved": len(entries) - len(not_found),
            "not_found": len(not_found),
            "synonym_groups": len(synonym_groups),
            "rename_suggestions": len(rename_suggestions),
        },
        "synonym_groups": [
            {
                "cid": cid,
                "iupac_ref": cache[names[0]].get("iupac"),
                "smiles_ref": cache[names[0]].get("smiles"),
                "names": names,
            }
            for cid, names in synonym_groups
        ],
        "rename_suggestions": [
            {"current": n, "iupac": i, "cid": c}
            for n, i, c in rename_suggestions
        ],
        "not_found": not_found,
        "cross_check_conflicts": [
            {"entry": e, "target": t, "type": "real" if e != t else "self"}
            for e, t in conflicts
        ],
        "all_results": {
            name: cache.get(name, {"cid": None, "iupac": None, "smiles": None})
            for name in entries
        },
    }
    with open(DATA_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(f"  Données JSON sauvegardées → {DATA_FILE}\n")


if __name__ == "__main__":
    main()
