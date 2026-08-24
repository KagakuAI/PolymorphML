import time

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

INPUT_CSV = "csd_cleaned.csv"
OUTPUT_CSV = "similarity_top10.csv"
FAILED_CSV = "similarity_failed_smiles.csv"
FP_RADIUS = 2
FP_SIZE = 2048
TOP_K = 10
CHUNK_SIZE = 500

# ── unique molecules ─────────────────────────────────────────────────────────
# Similarity depends only on the SMILES fed into RDKit, so dedup on SMILES
# (not InChI): some SMILES in this dataset map to several InChI (undefined
# stereocenters), but those rows would get bit-identical fingerprints anyway.
data = pd.read_csv(INPUT_CSV, usecols=["smiles"])
smiles = data["smiles"].dropna().unique()
print(f"{len(smiles)} unique SMILES", flush=True)

# ── fingerprints ─────────────────────────────────────────────────────────────
gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_SIZE)

valid_smiles = []
fps = []
failed = []
for s in smiles:
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        failed.append(s)
        continue
    valid_smiles.append(s)
    fps.append(gen.GetFingerprintAsNumPy(mol))

if failed:
    pd.Series(failed, name="smiles").to_csv(FAILED_CSV, index=False)
    print(f"{len(failed)} SMILES failed to parse -> {FAILED_CSV}", flush=True)

X = np.array(fps, dtype=np.float32)
onbits = X.sum(axis=1)
n = X.shape[0]
k = min(TOP_K, n - 1)
print(f"fingerprint matrix: {X.shape}, {X.nbytes / 1e9:.2f} GB", flush=True)

# ── chunked top-k Tanimoto ────────────────────────────────────────────────────
# Tanimoto(a, b) = |a & b| / (|a| + |b| - |a & b|). For 0/1 fingerprints,
# a @ b is exactly the AND-popcount, so a block matmul against the full
# matrix gives every pairwise intersection in one BLAS call per chunk.
rows_out, ranks_out, neighbors_out, scores_out = [], [], [], []

t0 = time.time()
for start in range(0, n, CHUNK_SIZE):
    end = min(start + CHUNK_SIZE, n)

    inter = X[start:end] @ X.T
    union = onbits[start:end, None] + onbits[None, :] - inter
    sim = np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)

    local_idx = np.arange(end - start)
    sim[local_idx, start + local_idx] = -1  # exclude self-match

    top_idx = np.argpartition(-sim, k, axis=1)[:, :k]
    top_scores = np.take_along_axis(sim, top_idx, axis=1)
    order = np.argsort(-top_scores, axis=1)
    top_idx = np.take_along_axis(top_idx, order, axis=1)
    top_scores = np.take_along_axis(top_scores, order, axis=1)

    for local_i in local_idx:
        global_i = start + local_i
        rows_out.extend([valid_smiles[global_i]] * k)
        ranks_out.extend(range(1, k + 1))
        neighbors_out.extend(valid_smiles[j] for j in top_idx[local_i])
        scores_out.extend(top_scores[local_i].tolist())

    elapsed = time.time() - t0
    eta = elapsed / end * (n - end)
    print(
        f"{end}/{n} ({100 * end / n:.1f}%) "
        f"elapsed={elapsed / 60:.1f}min eta={eta / 60:.1f}min",
        flush=True,
    )

out = pd.DataFrame(
    {
        "smiles": rows_out,
        "rank": ranks_out,
        "neighbor_smiles": neighbors_out,
        "tanimoto": scores_out,
    }
)
out.to_csv(OUTPUT_CSV, index=False)
print(f"saved {len(out)} rows -> {OUTPUT_CSV}", flush=True)
