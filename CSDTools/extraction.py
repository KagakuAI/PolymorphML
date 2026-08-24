import re
from concurrent.futures import ProcessPoolExecutor

from ccdc import io
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import time

RACEMIC_RE = re.compile(
    r'\brac-|\(rac\)|\bDL-|\(DL\)|\bD,L-|\(RS\)|\(R,S\)|\(R/S\)|\(\+/-\)|\bracemic\b',
    re.IGNORECASE,
)

# ─── Configuration ────────────────────────────────────────────────────────────
# To skip a property, remove it from the relevant list/set.
#
# ENTRY_PROPS / MOLECULE_PROPS / CRYSTAL_PROPS : properties that return native
# Python types (bool, int, float, str, None) — extracted via generic loop.
#
# MANUAL_PROPS : properties that return non-native types and require explicit
# handling below. Remove from the set to skip extraction entirely.

ENTRY_PROPS = [
    "identifier",
    "chemical_name",
    "is_organic",
    "has_disorder",
    "disorder_details",
    "color",
    "polymorph",
    "phase_transition",
    "radiation_source",
    "temperature",
    "r_factor",
    "calculated_density",
    "melting_point",
    "melting_point_default_units",
    # "heat_capacity",        # Solubility Platform not available in license
    # "heat_capacity_notes",
    # "heat_of_fusion",
    # "heat_of_fusion_notes",
    # "solubility_data",
    "bioactivity",
    "ccdc_number",
    "doi",
    "solvent",
]

MOLECULE_PROPS = [
    "smiles",
    "molecular_volume",
    "molecular_weight",
]

CRYSTAL_PROPS = [
    "spacegroup_symbol",
    "z_prime",
    "z_value",
    "cell_volume",
]

MANUAL_PROPS = {
    "publications",               # Citation tuple → "pub|pub|..."  +  publication_years
    "deposition_date",            # datetime.date  → "YYYY-MM-DD"
    "spacegroup_number_and_setting",  # (int, int) tuple → spacegroup_number + spacegroup_setting
    "symmetry_operators",         # str tuple → "op|op|..."
    "inchi",                      # generated via mol.generate_inchi() — includes Z'>1 count prefix
    "inchi_single",               # generated via mol.components[0].generate_inchi() — always single-molecule
    "cell_lengths",               # non-native object → cell_length_a/b/c
    "cell_angles",                # non-native object → cell_angle_alpha/beta/gamma
    "reduced_cell",               # Crystal.ReducedCell → reduced_cell_a/b/c/alpha/beta/gamma
    "has_metal",                  # computed from mol.atoms
    "num_components",             # computed from mol.components
    "is_single_species",          # True if all components share the same formula (pure crystal, Z'>1 ok)
    "is_racemic_name",            # derived from chemical_name via RACEMIC_RE
    "smiles_canonical",           # RDKit canonical SMILES of largest fragment
}

COLUMN_ORDER = [
    # ENTRY
    "identifier", "chemical_name", "is_racemic_name",
    "is_organic", "has_metal", "num_components", "is_single_species",
    "has_disorder", "disorder_details", "color",
    "polymorph", "phase_transition", "radiation_source", "temperature", "r_factor",
    "calculated_density", "melting_point", "melting_point_default_units",
    "bioactivity", "ccdc_number", "doi", "publications", "publication_years", "deposition_date", "solvent",
    # MOLECULE
    "smiles", "smiles_canonical", "molecular_volume", "molecular_weight", "inchi", "inchi_single",
    # CRYSTAL
    "spacegroup_number", "spacegroup_setting", "spacegroup_symbol",
    "symmetry_operators", "z_prime", "z_value",
    "cell_length_a", "cell_length_b", "cell_length_c",
    "cell_angle_alpha", "cell_angle_beta", "cell_angle_gamma",
    "cell_volume",
    "reduced_cell_a", "reduced_cell_b", "reduced_cell_c",
    "reduced_cell_alpha", "reduced_cell_beta", "reduced_cell_gamma",
]

# ─── Extraction ────────────────────────────────────────────────────────────────

NATIVE_TYPES = (bool, int, float, str, type(None))


def chunk_extraction(chunk_range):
    start, stop = chunk_range
    results = []
    warnings = set()
    with io.EntryReader('CSD') as reader:
        for i in range(start, stop):
            entry = reader[i]
            mol = entry.molecule
            crys = entry.crystal

            row = {}
            for prop in ENTRY_PROPS:
                try:
                    val = getattr(entry, prop, None)
                    if not isinstance(val, NATIVE_TYPES):
                        warnings.add(f"[ENTRY] {prop} -> {type(val).__name__} : {repr(val)[:60]}")
                        val = None
                    row[prop] = val
                except Exception:
                    row[prop] = None

            for prop in MOLECULE_PROPS:
                try:
                    val = getattr(mol, prop, None)
                    if not isinstance(val, NATIVE_TYPES):
                        warnings.add(f"[MOLECULE] {prop} -> {type(val).__name__} : {repr(val)[:60]}")
                        val = None
                    row[prop] = val
                except Exception:
                    row[prop] = None

            for prop in CRYSTAL_PROPS:
                try:
                    val = getattr(crys, prop, None)
                    if not isinstance(val, NATIVE_TYPES):
                        warnings.add(f"[CRYSTAL] {prop} -> {type(val).__name__} : {repr(val)[:60]}")
                        val = None
                    row[prop] = val
                except Exception:
                    row[prop] = None

            if "spacegroup_number_and_setting" in MANUAL_PROPS:
                try:
                    sg = crys.spacegroup_number_and_setting
                    row["spacegroup_number"] = sg[0] if sg else None
                    row["spacegroup_setting"] = sg[1] if sg else None
                except Exception:
                    row["spacegroup_number"] = None
                    row["spacegroup_setting"] = None

            if "publications" in MANUAL_PROPS:
                try:
                    pubs = entry.publications
                    row["publications"] = "|".join(str(p) for p in pubs)
                    years = [str(p.year) for p in pubs if p.year is not None and p.year >= 1800]
                    row["publication_years"] = "|".join(years) if years else None
                except Exception:
                    row["publications"] = None
                    row["publication_years"] = None

            if "deposition_date" in MANUAL_PROPS:
                try:
                    d = entry.deposition_date
                    row["deposition_date"] = d.isoformat() if d is not None else None
                except Exception:
                    row["deposition_date"] = None

            if "symmetry_operators" in MANUAL_PROPS:
                try:
                    row["symmetry_operators"] = "|".join(crys.symmetry_operators)
                except Exception:
                    row["symmetry_operators"] = None

            if "inchi" in MANUAL_PROPS:
                try:
                    row["inchi"] = mol.generate_inchi().inchi
                except Exception:
                    row["inchi"] = None

            if "inchi_single" in MANUAL_PROPS:
                try:
                    row["inchi_single"] = mol.components[0].generate_inchi().inchi
                except Exception:
                    row["inchi_single"] = None

            if "cell_lengths" in MANUAL_PROPS:
                try:
                    cl = crys.cell_lengths
                    row["cell_length_a"] = cl.a
                    row["cell_length_b"] = cl.b
                    row["cell_length_c"] = cl.c
                except Exception:
                    row["cell_length_a"] = None
                    row["cell_length_b"] = None
                    row["cell_length_c"] = None

            if "cell_angles" in MANUAL_PROPS:
                try:
                    ca = crys.cell_angles
                    row["cell_angle_alpha"] = ca.alpha
                    row["cell_angle_beta"] = ca.beta
                    row["cell_angle_gamma"] = ca.gamma
                except Exception:
                    row["cell_angle_alpha"] = None
                    row["cell_angle_beta"] = None
                    row["cell_angle_gamma"] = None

            if "reduced_cell" in MANUAL_PROPS:
                try:
                    rc = crys.reduced_cell
                    row["reduced_cell_a"]     = rc.cell_lengths.a
                    row["reduced_cell_b"]     = rc.cell_lengths.b
                    row["reduced_cell_c"]     = rc.cell_lengths.c
                    row["reduced_cell_alpha"] = rc.cell_angles.alpha
                    row["reduced_cell_beta"]  = rc.cell_angles.beta
                    row["reduced_cell_gamma"] = rc.cell_angles.gamma
                except Exception:
                    for col in ("reduced_cell_a", "reduced_cell_b", "reduced_cell_c",
                                "reduced_cell_alpha", "reduced_cell_beta", "reduced_cell_gamma"):
                        row[col] = None

            if "has_metal" in MANUAL_PROPS:
                row["has_metal"] = any(atom.is_metal for atom in mol.atoms)

            if "num_components" in MANUAL_PROPS:
                row["num_components"] = len(mol.components)

            if "is_single_species" in MANUAL_PROPS:
                comps = mol.components
                if len(comps) <= 1:
                    row["is_single_species"] = True
                else:
                    formulas = [c.formula for c in comps]
                    row["is_single_species"] = len(set(formulas)) == 1

            if "is_racemic_name" in MANUAL_PROPS:
                name = row.get("chemical_name")
                row["is_racemic_name"] = bool(RACEMIC_RE.search(name)) if name else False

            if "smiles_canonical" in MANUAL_PROPS:
                smiles = row.get("smiles")
                if isinstance(smiles, str):
                    fragment = max(smiles.split("."), key=len)
                    mol = Chem.MolFromSmiles(fragment)
                    row["smiles_canonical"] = Chem.MolToSmiles(mol) if mol else None
                else:
                    row["smiles_canonical"] = None

            results.append(row)
    return results, warnings


def run_extraction(output_path="csd_all.csv", num_workers=24, chunk_size=500, max_entries=None):
    with io.EntryReader('CSD') as reader:
        n = len(reader)
    if max_entries:
        n = min(n, max_entries)

    chunk_ranges = [(i, min(i + chunk_size, n)) for i in range(0, n, chunk_size)]
    pbar = tqdm(total=n, unit="entry", desc=f"Extraction ({num_workers} workers)", smoothing=0.1)

    start = time.time()
    all_results = []
    all_warnings = set()
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for chunk_results, chunk_warnings in executor.map(chunk_extraction, chunk_ranges):
            all_results.extend(chunk_results)
            all_warnings.update(chunk_warnings)
            pbar.update(len(chunk_results))
    pbar.close()

    delta_time = round(time.time() - start, 1)
    df = pd.DataFrame(all_results)
    df = df[[col for col in COLUMN_ORDER if col in df.columns]]
    df.to_csv(output_path, index=False)

    print(f"Processed {len(df)//1000}k entries in {delta_time}s ({num_workers} workers)")
    if all_warnings:
        print("\n⚠️  Non-native types detected (set to None):")
        for w in sorted(all_warnings):
            print(" ", w)

    return df


if __name__ == "__main__":
    run_extraction()
