#!/usr/bin/env python3
"""
analyse_whitelist.py v2 — Détection de synonymes dans WHITELIST via PubChem CID
                           + nom préféré PubChem (Title) pour col. B

Usage :
    python analyse_whitelist.py solvent_data.py

Sorties :
    whitelist_cache.json   — cache (permet de reprendre si coupé)
    whitelist_report.txt   — rapport lisible
    whitelist_data.json    — données complètes
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
DELAY        = 0.22          # ~4.5 req/s (limite PubChem = 5/s)
MAX_RETRIES  = 3
BATCH_SIZE   = 50            # pour la récupération des titres par lots
CACHE_FILE   = "whitelist_cache.json"
REPORT_FILE  = "whitelist_report.txt"
DATA_FILE    = "whitelist_data.json"


# ── Extraction ─────────────────────────────────────────────────────────────────
def extract_from_script(path):
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
        raise ValueError("WHITELIST introuvable.")
    return sorted(whitelist), synonyms or {}


# ── PubChem — CID + IUPACName + SMILES ────────────────────────────────────────
def fetch_pubchem(name):
    """Retourne {"cid", "iupac", "smiles"} pour un nom donné."""
    encoded = urllib.parse.quote(name, safe="")
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/property/IUPACName,CanonicalSMILES/JSON"
    )
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "whitelist-checker/2.0"}
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
                    "preferred": None,  # rempli plus tard par lot
                }
            return {"cid": None, "iupac": None, "smiles": None, "preferred": None}

        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"cid": None, "iupac": None, "smiles": None, "preferred": None}
            if e.code in (429, 503) and attempt < MAX_RETRIES - 1:
                wait = 2 ** (attempt + 1)
                print(f"\n  ⚠ Rate-limit ({e.code}) — attente {wait}s...", end="")
                time.sleep(wait)
                continue
            return {"cid": None, "iupac": None, "smiles": None, "preferred": None}
        except Exception as exc:
            if attempt < MAX_RETRIES - 1:
                time.sleep(1)
                continue
            return {"cid": None, "iupac": None, "smiles": None, "preferred": None}

    return {"cid": None, "iupac": None, "smiles": None, "preferred": None}


# ── PubChem — Title (nom préféré) par lots de CIDs ────────────────────────────
def fetch_titles_batch(cids):
    """
    Retourne un dict {cid: title} pour une liste de CIDs.
    Utilise l'endpoint /description qui donne le 'Title' officiel PubChem
    (ex: 'Water', 'Ethanol', 'Cumene' — pas 'oxidane', 'ethanol', etc.)
    """
    cid_str = ",".join(str(c) for c in cids)
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid_str}/description/JSON"
    )
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "whitelist-checker/2.0"}
            )
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.loads(resp.read())
            result = {}
            for info in data.get("InformationList", {}).get("Information", []):
                cid   = info.get("CID")
                title = info.get("Title")
                if cid and title:
                    result[cid] = title
            return result
        except urllib.error.HTTPError as e:
            if e.code in (429, 503) and attempt < MAX_RETRIES - 1:
                time.sleep(2 ** (attempt + 1))
                continue
            return {}
        except Exception:
            if attempt < MAX_RETRIES - 1:
                time.sleep(1)
                continue
            return {}
    return {}


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

    # 3. Cache
    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE) as f:
            cache = json.load(f)
        already = sum(1 for e in entries if e in cache)
        print(f"② Cache trouvé ({already}/{len(entries)} déjà résolus) — reprise.")
    else:
        print("② Aucun cache — démarrage complet.")

    # 4. Requêtes PubChem — CID + IUPAC + SMILES
    to_fetch = [e for e in entries if e not in cache]
    if to_fetch:
        eta = len(to_fetch) * DELAY
        print(f"\n   {len(to_fetch)} composés à interroger (~{eta:.0f}s)\n")
        for i, name in enumerate(to_fetch, 1):
            result = fetch_pubchem(name)
            cache[name] = result
            progress(i, len(to_fetch), name, result["cid"] is not None)
            if i % 50 == 0:
                with open(CACHE_FILE, "w") as f:
                    json.dump(cache, f, ensure_ascii=False)
            time.sleep(DELAY)
        with open(CACHE_FILE, "w") as f:
            json.dump(cache, f, ensure_ascii=False)
        print(f"\n\n   Cache sauvegardé → {CACHE_FILE}\n")
    else:
        print("   Tous les composés sont en cache.\n")

    # 5. Récupération des noms préférés (Title) par lots
    print("③ Récupération des noms préférés PubChem (Title) par lots...")
    cids_needed = [
        cache[name]["cid"]
        for name in entries
        if name in cache and cache[name]["cid"] is not None
        and cache[name].get("preferred") is None
    ]
    cids_needed = list(set(cids_needed))  # dédoublonne

    if cids_needed:
        title_map = {}
        batches = [cids_needed[i:i+BATCH_SIZE] for i in range(0, len(cids_needed), BATCH_SIZE)]
        for idx, batch in enumerate(batches, 1):
            print(f"\r   Lot {idx}/{len(batches)} ({len(batch)} CIDs)...", end="", flush=True)
            titles = fetch_titles_batch(batch)
            title_map.update(titles)
            time.sleep(DELAY * 2)   # un peu plus de marge pour les batch

        # Injecte dans le cache
        for name in entries:
            info = cache.get(name, {})
            if info.get("cid") and info["cid"] in title_map:
                info["preferred"] = title_map[info["cid"]]
                cache[name] = info

        with open(CACHE_FILE, "w") as f:
            json.dump(cache, f, ensure_ascii=False)
        print(f"\n   Noms préférés récupérés ({len(title_map)}) — cache mis à jour.\n")
    else:
        print("   Noms préférés déjà en cache.\n")

    # 6. Analyse
    by_cid    = defaultdict(list)
    not_found = []
    for name in entries:
        info = cache.get(name, {"cid": None})
        if info["cid"] is not None:
            by_cid[info["cid"]].append(name)
        else:
            not_found.append(name)

    synonym_groups = sorted(
        [(cid, names) for cid, names in by_cid.items() if len(names) > 1],
        key=lambda x: -len(x[1]),
    )

    # Col B : utilise le nom PRÉFÉRÉ (Title) — ignore si identique (insensible à la casse)
    rename_suggestions = []
    for name in entries:
        info = cache.get(name, {})
        if not info.get("cid"):
            continue
        preferred = info.get("preferred")
        if preferred and preferred.lower() != name.lower():
            rename_suggestions.append((name, preferred, info["cid"]))
    rename_suggestions.sort(key=lambda x: x[0])

    # 7. Rapport
    lines = []
    sep   = "─" * 65

    def h(title):
        lines.append("")
        lines.append(sep)
        lines.append(f"  {title}")
        lines.append(sep)

    lines.append("=" * 65)
    lines.append("  RAPPORT — ANALYSE WHITELIST parse_solvent.py  (v2)")
    lines.append("=" * 65)
    lines.append(f"  Entrées analysées  : {len(entries)}")
    lines.append(f"  Résolus PubChem    : {len(entries) - len(not_found)}")
    lines.append(f"  Non résolus        : {len(not_found)}")
    lines.append(f"  Synonymes (col. A) : {len(synonym_groups)}")
    lines.append(f"  Renommages (col. B): {len(rename_suggestions)}")
    lines.append(f"  Cross-check        : {len(conflicts)} ({sum(1 for e,t in conflicts if e!=t)} réel(s))")
    lines.append("")
    lines.append("  Note col. B : utilise le nom PRÉFÉRÉ PubChem (Title),")
    lines.append("  pas le nom IUPAC strict (pas d'oxidane pour water, etc.)")

    # A — Synonymes
    h("COL. A — SYNONYMES  (même CID PubChem → même molécule)")
    if synonym_groups:
        for cid, names in synonym_groups:
            preferred_ref = cache[names[0]].get("preferred") or cache[names[0]].get("iupac", "?")
            lines.append(f"\n  CID {cid}  |  Nom préféré : {preferred_ref}")
            for n in names:
                lines.append(f"    {'✓' if n == names[0] else '→'} {n}")
        lines.append("")
        lines.append("  ⚠  Note : ✓/→ = ordre alphabétique seulement, pas une décision chimique.")
        lines.append("     Choisir le canonique en fonction du nom préféré PubChem.")
    else:
        lines.append("  Aucun synonyme détecté via CID.")

    # B — Renommages
    h(f"COL. B — RENOMMAGES SUGGÉRÉS  (nom whitelist ≠ nom préféré PubChem)")
    lines.append(f"  {'Nom actuel (whitelist)':<42}  →  Nom préféré PubChem")
    lines.append(f"  {'─'*42}  {'─'*30}")
    for name, preferred, cid in rename_suggestions:
        lines.append(f"  {name:<42}  →  {preferred}")

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
    print(f"\n  Rapport → {REPORT_FILE}")

    # JSON complet
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
                "preferred_ref": cache[names[0]].get("preferred"),
                "iupac_ref":     cache[names[0]].get("iupac"),
                "smiles_ref":    cache[names[0]].get("smiles"),
                "names": names,
            }
            for cid, names in synonym_groups
        ],
        "rename_suggestions": [
            {"current": n, "preferred": p, "cid": c}
            for n, p, c in rename_suggestions
        ],
        "not_found": not_found,
        "cross_check_conflicts": [
            {"entry": e, "target": t, "type": "real" if e != t else "self"}
            for e, t in conflicts
        ],
        "all_results": {
            name: cache.get(name, {"cid": None, "iupac": None, "smiles": None, "preferred": None})
            for name in entries
        },
    }
    with open(DATA_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(f"  Données JSON → {DATA_FILE}\n")


if __name__ == "__main__":
    main()
