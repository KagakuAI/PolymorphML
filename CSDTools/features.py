import os

import numpy as np
import pandas as pd
from molfeat.trans import MoleculeTransformer


def get_or_build_features(csv_file, featurizer, smiles_list=None, cache_root="feature_cache",
                           n_jobs=-1, cache_key=None, **featurizer_kwargs):
    """Load cached features for (featurizer, csv_file), computing and caching them if needed.

    cache_key defaults to featurizer, but can be set explicitly to disambiguate
    featurizer variants sharing the same molfeat name (e.g. "pharm2d" with different
    `factory` kwargs) so they don't collide in the cache directory.
    Cache is invalidated when csv_file is newer than the cached .npy files.
    Returns (x, valid_ids).
    """
    cache_key = cache_key or featurizer
    cache_dir = os.path.join(cache_root, cache_key)
    os.makedirs(cache_dir, exist_ok=True)
    csv_stem  = os.path.splitext(os.path.basename(csv_file))[0]
    cache_x   = os.path.join(cache_dir, f"{csv_stem}.npy")
    cache_ids = os.path.join(cache_dir, f"{csv_stem}_ids.npy")
    csv_mtime = os.path.getmtime(csv_file)
    cache_ok  = os.path.exists(cache_x) and os.path.getmtime(cache_x) > csv_mtime

    if cache_ok:
        x, valid_ids = np.load(cache_x), np.load(cache_ids)
        print(f"Cache chargé [{cache_key} / {csv_file}] : {x.shape}")
        return x, valid_ids

    if smiles_list is None:
        smiles_list = pd.read_csv(csv_file)["smiles_canonical"].to_list()

    desc_calc = MoleculeTransformer(featurizer=featurizer, n_jobs=n_jobs, dtype=float, **featurizer_kwargs)
    x, valid_ids = desc_calc(smiles_list, ignore_errors=True)
    np.save(cache_x, x)
    np.save(cache_ids, valid_ids)
    print(f"Featurisation [{cache_key} / {csv_file}] terminée : {x.shape}")
    return x, valid_ids


FEATURIZERS = [
    "avalon", "rdkit", "maccs", "atompair-count", "fcfp", "fcfp-count",
    "ecfp", "ecfp-count", "topological", "topological-count", "secfp", "scaffoldkeys",
    "desc2D", "estate", "erg", "cats2d", "cats",
    ("pharm2d-pmapper", "pharm2d", {}),
    ("pharm2d-gobbi", "pharm2d", {"factory": "gobbi"}),
]


def build_all_caches(csv_file="../csd_ml.csv", featurizers=FEATURIZERS):
    """Precompute and cache features for every featurizer in `featurizers`.

    Entries are either a plain molfeat name (str) or a (cache_key, featurizer, kwargs)
    tuple for variants sharing a molfeat name (e.g. pharm2d's pmapper/gobbi factories).

    Run from CSDTools/: `cd CSDTools && python features.py`.
    Each featurizer already parallelizes internally via n_jobs (joblib/loky),
    so featurizers are processed one at a time rather than in parallel.
    """
    smiles_list = pd.read_csv(csv_file)["smiles_canonical"].to_list()
    failed = []
    for entry in featurizers:
        cache_key, featurizer, kwargs = entry if isinstance(entry, tuple) else (entry, entry, {})
        try:
            get_or_build_features(csv_file, featurizer, smiles_list=smiles_list,
                                   cache_root='../feature_cache', cache_key=cache_key, **kwargs)
        except Exception as e:
            print(f"⛔ {cache_key} failed: {e}")
            failed.append(cache_key)
    if failed:
        print(f"\n{len(failed)}/{len(featurizers)} featurizers failed: {failed}")
    else:
        print(f"\nAll {len(featurizers)} featurizers cached successfully.")


if __name__ == "__main__":
    build_all_caches()
