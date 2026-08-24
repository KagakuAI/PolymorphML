"""patch_smiles_canonical.py — adds 'smiles_canonical' column to csd_all.csv.

RDKit canonical SMILES of the largest fragment (handles Z'>1 multi-component SMILES).
Used by forms_for_smiles() for reliable SMILES-based lookup.
"""
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

TARGET_CSV   = "../csd_cleaned.csv"
INSERT_AFTER = "smiles"


def canonical_smiles(smiles):
    if not isinstance(smiles, str):
        return None
    fragment = max(smiles.split("."), key=len)
    mol = Chem.MolFromSmiles(fragment)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def patch_csv(csv_path=TARGET_CSV):
    df = pd.read_csv(csv_path, low_memory=False)
    if "smiles_canonical" in df.columns:
        df = df.drop(columns=["smiles_canonical"])

    tqdm.pandas(desc="smiles_canonical")
    df["smiles_canonical"] = df["smiles"].progress_apply(canonical_smiles)

    filled = df["smiles_canonical"].notna().sum()
    print(f"  smiles_canonical: {filled:,} / {len(df):,} filled ({100 * filled / len(df):.1f}%)")

    cols = list(df.columns)
    cols.remove("smiles_canonical")
    idx = cols.index(INSERT_AFTER) + 1 if INSERT_AFTER in cols else len(cols)
    cols.insert(idx, "smiles_canonical")
    df = df[cols]

    df.to_csv(csv_path, index=False)
    print(f"\nPatched {csv_path}: 'smiles_canonical' inserted after '{INSERT_AFTER}'.")


if __name__ == "__main__":
    patch_csv()
