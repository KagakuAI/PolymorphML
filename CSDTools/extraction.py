from concurrent.futures import ProcessPoolExecutor

from ccdc import io
from tqdm import tqdm
import pandas as pd
import time

from CSDTools.extraction_config import (
    ENTRY_PROPS, CRYSTAL_PROPS, MOLECULE_PROPS, COLUMN_ORDER
)

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

            try:
                sg = crys.spacegroup_number_and_setting
                row["spacegroup_number"] = sg[0] if sg else None
                row["spacegroup_setting"] = sg[1] if sg else None
            except Exception:
                row["spacegroup_number"] = None
                row["spacegroup_setting"] = None

            try:
                pubs = entry.publications
                row["publications"] = "|".join(str(p) for p in pubs)
                years = [str(p.year) for p in pubs if p.year is not None and p.year >= 1800]
                row["publication_years"] = "|".join(years) if years else None
            except Exception:
                row["publications"] = None
                row["publication_years"] = None

            try:
                row["symmetry_operators"] = "|".join(crys.symmetry_operators)
            except Exception:
                row["symmetry_operators"] = None

            try:
                row["inchi"] = mol.generate_inchi().inchi
            except Exception:
                row["inchi"] = None

            try:
                cl = crys.cell_lengths
                row["cell_length_a"] = cl.a
                row["cell_length_b"] = cl.b
                row["cell_length_c"] = cl.c
            except Exception:
                row["cell_length_a"] = None
                row["cell_length_b"] = None
                row["cell_length_c"] = None

            try:
                ca = crys.cell_angles
                row["cell_angle_alpha"] = ca.alpha
                row["cell_angle_beta"] = ca.beta
                row["cell_angle_gamma"] = ca.gamma
            except Exception:
                row["cell_angle_alpha"] = None
                row["cell_angle_beta"] = None
                row["cell_angle_gamma"] = None

            row["has_metal"] = any(atom.is_metal for atom in mol.atoms)
            row["num_components"] = len(mol.components)

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
    df = pd.DataFrame(all_results)[COLUMN_ORDER]
    df.to_csv(output_path, index=False)

    print(f"Processed {len(df)//1000}k entries in {delta_time}s ({num_workers} workers)")
    if all_warnings:
        print("\n⚠️  Non-native types detected (set to None):")
        for w in sorted(all_warnings):
            print(" ", w)

    return df


if __name__ == "__main__":
    run_extraction()
