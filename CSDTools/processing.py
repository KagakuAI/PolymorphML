import pandas as pd
from rdkit import Chem

LENGTH_PARAMS   = ["reduced_cell_a", "reduced_cell_b", "reduced_cell_c"]
ANGLE_PARAMS    = ["reduced_cell_alpha", "reduced_cell_beta", "reduced_cell_gamma"]
DISCRETE_PARAMS = ["spacegroup_number", "z_prime", "z_value"]
CELL_PARAMS     = LENGTH_PARAMS + ANGLE_PARAMS  # for display
# 3% covers thermal expansion across ~200 K for typical organic crystals
# 1° absolute is sufficient to distinguish genuinely different polymorphs
LENGTH_RTOL = 0.03
ANGLE_ATOL  = 2.0   # degrees — covers thermal expansion (~150K to RT)


def all_forms_identical(group):
    """True if all crystal forms in the group are indistinguishable."""
    for param in DISCRETE_PARAMS:
        vals = group[param].dropna()
        if vals.nunique() > 1:
            return False
    for param in LENGTH_PARAMS:
        vals = group[param].dropna()
        if len(vals) < 2:
            continue
        ref = vals.mean()
        if ref > 0 and (vals.max() - vals.min()) / ref > LENGTH_RTOL:
            return False
    for param in ANGLE_PARAMS:
        vals = group[param].dropna()
        if len(vals) < 2:
            continue
        if vals.max() - vals.min() > ANGLE_ATOL:
            return False
    return True


def find_duplicate_clusters(group, length_rtol=LENGTH_RTOL, angle_atol=ANGLE_ATOL):
    """Groups a molecule's crystal forms into clusters of probably-identical structures.

    Pairwise union-find: two forms land in the same cluster if their discrete params
    match and their cell params agree within tolerance, so a single outlier no longer
    prevents detecting duplicates among the rest.
    Lengths use a relative tolerance (LENGTH_RTOL), angles an absolute one (ANGLE_ATOL).
    """
    rows = group.reset_index(drop=True)
    n = len(rows)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    def pair_identical(a, b):
        for param in DISCRETE_PARAMS:
            if pd.notna(a[param]) and pd.notna(b[param]) and a[param] != b[param]:
                return False
        for param in LENGTH_PARAMS:
            va, vb = a[param], b[param]
            if pd.isna(va) or pd.isna(vb):
                continue
            ref = (va + vb) / 2
            if ref > 0 and abs(va - vb) / ref > length_rtol:
                return False
        for param in ANGLE_PARAMS:
            va, vb = a[param], b[param]
            if pd.isna(va) or pd.isna(vb):
                continue
            if abs(va - vb) > angle_atol:
                return False
        return True

    for i in range(n):
        for j in range(i + 1, n):
            if pair_identical(rows.iloc[i], rows.iloc[j]):
                union(i, j)

    clusters = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(rows.iloc[i]["identifier"])
    return list(clusters.values())


def forms_for_smiles(smiles, df):
    """Return all crystal forms in df for the molecule identified by SMILES."""
    fragments = smiles.split(".")
    smiles_single = max(fragments, key=len) if len(fragments) > 1 else smiles
    mol = Chem.MolFromSmiles(smiles_single)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    canonical = Chem.MolToSmiles(mol)
    return df[df["smiles_canonical"] == canonical].copy()


def filter_entries_basic(df):
    """Minimal quality filter: organic, no metal, single species, valid SMILES."""
    mask = (
        (df["is_organic"] == True) &
        (df["has_metal"] == False) &
        (df["is_single_species"] == True) &
        (df["smiles_canonical"].notna())
    )
    return df[mask].copy()


def filter_entries(df):
    mask = (
        (df["is_organic"] == True) &
        (df["has_metal"] == False) &
        (df["is_single_species"] == True) &
        (df["is_racemic_name"] == False) &
        (df["r_factor"].notna()) &
        (df["r_factor"] <= 10.0) &
        (df["has_disorder"] == False) &
        (df["calculated_density"].isna() | (df["calculated_density"] <= 5.0)) &
        (df["smiles_canonical"].notna())
    )
    df_filtered = df[mask].copy()
    df_filtered, _ = filter_mw_iqr(df_filtered)
    return df_filtered



def build_polymorph_table(df_filtered):
    """Return one row per inchi_single group with polymorph statistics.

    Columns returned:
      inchi_single, smiles, num_entries,
      num_polymorphs (distinct forms after deduplication),
      all_identical (True if the whole group is one crystal re-deposited),
      identifiers (pipe-separated refcodes)
    """
    df = df_filtered.dropna(subset=["inchi_single"]).copy()
    rows = []
    for inchi, group in df.groupby("inchi_single"):
        group = group.reset_index(drop=True)
        if len(group) == 1:
            clusters = [[group["identifier"].iloc[0]]]
        else:
            clusters = find_duplicate_clusters(group)
        rows.append({
            "inchi_single":        inchi,
            "smiles":              group["smiles"].iloc[0],
            "num_entries":         len(group),
            "num_distinct_forms":  len(clusters),
            "identifiers":         "|".join(group["identifier"].tolist()),
        })
    return pd.DataFrame(rows)


IQR_K = 3.0  # IQR multiplier for molecular weight upper bound


def filter_mw_iqr(df, k=IQR_K):
    """Upper-only IQR filter on molecular_weight (IQR computed on df)."""
    mw = df["molecular_weight"].dropna()
    q1, q3 = mw.quantile(0.25), mw.quantile(0.75)
    upper = q3 + k * (q3 - q1)
    mask = df["molecular_weight"].isna() | (df["molecular_weight"] <= upper)
    return df[mask].copy(), round(upper, 1)


def build_ml_dataset_naive(df_filtered):
    """Same as build_ml_dataset but counts raw entries instead of deduplicating.
    A molecule is labelled polymorphic if it has more than one entry in the CSD,
    regardless of whether those entries are genuinely distinct crystal forms.
    """
    df = df_filtered.dropna(subset=["inchi_single", "smiles_canonical"]).copy()
    counts = df.groupby("inchi_single")["identifier"].count().rename("num_entries")
    smiles_map = df.drop_duplicates("inchi_single").set_index("inchi_single")["smiles_canonical"]
    poly = pd.DataFrame({"num_entries": counts, "smiles_canonical": smiles_map})
    poly["is_polymorphic"] = (poly["num_entries"] > 1).astype(int)
    ml = (
        poly.groupby("smiles_canonical")
        .agg(is_polymorphic=("is_polymorphic", "max"))
        .reset_index()
    )
    return ml


def build_ml_dataset(df_filtered):
    """Return one row per smiles_canonical with a polymorphism label for ML.

    Groups by inchi_single to count distinct crystal forms (handles stereoisomers
    separately), then aggregates to smiles_canonical level using max:
    if any stereoisomer is polymorphic, the SMILES is labelled polymorphic.

    Columns returned:
      smiles_canonical, is_polymorphic, num_inchi_groups,
      max_distinct_forms, min_distinct_forms
    """
    poly = build_polymorph_table(df_filtered)
    poly["is_polymorphic"] = (poly["num_distinct_forms"] > 1).astype(int)

    smiles_map = (
        df_filtered[["inchi_single", "smiles_canonical"]]
        .drop_duplicates("inchi_single")
    )
    poly = poly.merge(smiles_map, on="inchi_single", how="left")

    ml = (
        poly.groupby("smiles_canonical")
        .agg(
            is_polymorphic    =("is_polymorphic",    "max"),
            num_inchi_groups  =("inchi_single",       "count"),
            max_distinct_forms=("num_distinct_forms", "max"),
            min_distinct_forms=("num_distinct_forms", "min"),
        )
        .reset_index()
    )
    return ml


def run_processing(input_path="../csd_all.csv", output_path="../csd_ml.csv"):
    df = pd.read_csv(input_path, low_memory=False)

    df_filtered = filter_entries(df)
    print(f"Filtered : {len(df):,} → {len(df_filtered):,} entries")

    ml_df = build_ml_dataset(df_filtered)
    n_poly = int(ml_df["is_polymorphic"].sum())
    print(f"ML dataset : {len(ml_df):,} SMILES — {n_poly:,} polymorphiques ({100*n_poly/len(ml_df):.1f}%)")

    ml_df.to_csv(output_path, index=False)
    print(f"Saved → {output_path}")
    return ml_df


if __name__ == "__main__":
    run_processing()
