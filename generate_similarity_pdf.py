"""
Generate a PDF of monomorphic candidates with their close polymorphic neighbors.

One page per group (molecules sharing the exact same polymorphic neighbor set)
or per solo candidate.

Each molecule box shows:
  - 2D structure
  - CSD identifier(s)
  - Tanimoto similarity (polymorphic neighbors only)
  - Number of crystal forms (nb_forms)
  - Deposition date(s)
  - DOI(s)
"""

import io
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# ── Config ────────────────────────────────────────────────────────────────────
TANIMOTO_THRESHOLD = 0.9
MIN_POLY_NEIGHBORS = 4
OUTPUT_PDF         = "similarity_candidates.pdf"
MOL_W, MOL_H      = 900, 630   # 2D image size in pixels
POLY_COLS          = 4          # polymorphic neighbors per row
MAX_IDS            = 3          # identifiers shown before "+N"
MAX_SMILES         = 30         # chars before truncating SMILES

COLOR_MONO       = "#dce8f5"
COLOR_POLY       = "#fde8e4"
COLOR_MONO_TITLE = "#1a4a7a"
COLOR_POLY_TITLE = "#8b1a1a"

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading data…")
data = pd.read_csv("csd_cleaned.csv")
sim  = pd.read_csv("similarity_enriched.csv")

smiles_to_inchi = data.drop_duplicates("smiles").set_index("smiles")["inchi"]

# Per-SMILES lookup: list of (identifier, deposition_date, doi)
def build_smiles_lookup(data):
    lookup = {}
    for smi, grp in data.groupby("smiles"):
        rows = []
        for _, r in grp.iterrows():
            rows.append({
                "id":   r["identifier"],
                "date": r["deposition_date"],
                "doi":  r["doi"] if pd.notna(r["doi"]) else None,
            })
        lookup[smi] = rows
    return lookup

smiles_info = build_smiles_lookup(data)

# ── Build candidate groups ────────────────────────────────────────────────────
sim["poly_neighbor"] = sim["nb_forms_2"] > 1
mono = sim[sim["nb_forms_1"] == 1]

close = mono[mono["poly_neighbor"] & (mono["tanimoto"] >= TANIMOTO_THRESHOLD)].copy()
close["same_mol"] = (
    close["smiles_1"].map(smiles_to_inchi) == close["smiles_2"].map(smiles_to_inchi)
)
close = close[~close["same_mol"]]

refined = (
    close.groupby("smiles_1")
    .agg(n_close_poly=("tanimoto", "count"), best_sim=("tanimoto", "max"))
    .sort_values(["n_close_poly", "best_sim"], ascending=False)
)
strong       = refined[refined["n_close_poly"] >= MIN_POLY_NEIGHBORS].index
strong_close = close[close["smiles_1"].isin(strong)]

neighbor_sets = strong_close.groupby("smiles_1")["smiles_2"].apply(frozenset)
bucket = defaultdict(list)
for smi, nset in neighbor_sets.items():
    bucket[nset].append(smi)

pages = []
for nset, mono_list in bucket.items():
    # Tanimoto from every mono candidate for every poly neighbor
    sub   = strong_close[strong_close["smiles_1"].isin(mono_list)]
    pivot = sub.pivot_table(index="smiles_2", columns="smiles_1", values="tanimoto")

    # Base info (nb_forms, identifiers) from first candidate
    poly_base = (
        strong_close[
            (strong_close["smiles_1"] == mono_list[0]) &
            (strong_close["smiles_2"].isin(nset))
        ]
        [["smiles_2", "nb_forms_2", "identifiers_2"]]
        .drop_duplicates("smiles_2")
        .set_index("smiles_2")
    )
    poly_base["mean_tanimoto"] = pivot.mean(axis=1)
    poly_df = poly_base.sort_values("mean_tanimoto", ascending=False).reset_index()

    pages.append((mono_list, poly_df, pivot))

pages.sort(key=lambda x: (-len(x[0]), -len(x[1])))
print(f"Pages to generate: {len(pages)}")


# ── Drawing helpers ───────────────────────────────────────────────────────────
def mol_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(MOL_W, MOL_H)
    d.drawOptions().padding = 0.15
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return np.array(Image.open(io.BytesIO(d.GetDrawingText())))

def fmt_smi(smi, maxlen=MAX_SMILES):
    return smi if len(smi) <= maxlen else smi[:maxlen] + "…"

def mol_caption(smiles, nb_forms=None, tanimoto_vals=None, max_ids=MAX_IDS):
    """Build multi-line caption string for a molecule box.

    tanimoto_vals: list of floats (one per monomorphic candidate on the page),
                   or None for monomorphic boxes.
    """
    rows = smiles_info.get(smiles, [])

    ids   = [r["id"]   for r in rows]
    dates = [r["date"] for r in rows if pd.notna(r["date"])]
    dois  = [r["doi"]  for r in rows if r["doi"]]

    id_str = ", ".join(ids[:max_ids])
    if len(ids) > max_ids:
        id_str += f"  +{len(ids) - max_ids}"

    if dates:
        dates_sorted = sorted(dates)
        date_str = dates_sorted[0][:4]
        if dates_sorted[-1][:4] != dates_sorted[0][:4]:
            date_str += f"–{dates_sorted[-1][:4]}"
    else:
        date_str = "—"

    doi_str = dois[0] if dois else "—"
    if len(doi_str) > 35:
        doi_str = doi_str[:35] + "…"

    lines = [f"ID: {id_str}", f"Date: {date_str}  ·  nb_forms: {nb_forms if nb_forms is not None else 1}"]
    if tanimoto_vals is not None:
        unique = set(round(v, 4) for v in tanimoto_vals)
        tan_str = f"{tanimoto_vals[0]:.3f}" if len(unique) == 1 else " / ".join(f"{v:.3f}" for v in tanimoto_vals)
        lines.append(f"Tanimoto: {tan_str}")
    lines.append(f"DOI: {doi_str}")
    return "\n".join(lines)

def draw_mol_box(ax, smiles, bg_color, title, caption):
    ax.set_facecolor(bg_color)
    for spine in ax.spines.values():
        spine.set_edgecolor("#aaaaaa")
        spine.set_linewidth(0.5)

    img = mol_image(smiles)
    if img is not None:
        ax.imshow(img, aspect="auto")

    title_color = COLOR_MONO_TITLE if bg_color == COLOR_MONO else COLOR_POLY_TITLE
    ax.set_title(title, fontsize=6.5, color=title_color, fontweight="bold", pad=3)
    ax.set_xlabel(caption, fontsize=5.5, color="#333333", labelpad=4, linespacing=1.6)
    ax.set_xticks([])
    ax.set_yticks([])

# ── Generate PDF ──────────────────────────────────────────────────────────────
with PdfPages(OUTPUT_PDF) as pdf:
    for page_num, (mono_list, poly_df, pivot) in enumerate(pages):
        n_mono  = len(mono_list)
        n_poly  = len(poly_df)
        n_pcols = POLY_COLS
        n_prows = (n_poly + n_pcols - 1) // n_pcols

        fig_h = 3.2 + 3.0 * n_prows
        fig   = plt.figure(figsize=(14, max(fig_h, 6)))
        fig.patch.set_facecolor("white")

        total_rows = 1 + n_prows
        n_cols = max(n_mono, n_pcols)
        gs = GridSpec(
            total_rows, n_cols, figure=fig,
            hspace=0.7, wspace=0.25,
            top=0.92, bottom=0.04, left=0.02, right=0.98,
        )

        group_label = f"{n_mono} monomorphic candidates" if n_mono > 1 else "Monomorphic candidate"
        fig.suptitle(
            f"{group_label}  ·  {n_poly} close polymorphic neighbors  (Tanimoto ≥ {TANIMOTO_THRESHOLD})",
            fontsize=10, fontweight="bold",
        )

        # Monomorphic row (centred)
        offset = (n_cols - n_mono) // 2
        for i, smi in enumerate(mono_list):
            ax = fig.add_subplot(gs[0, offset + i])
            draw_mol_box(
                ax, smi,
                bg_color=COLOR_MONO,
                title=f"MONOMORPHIC  ·  {fmt_smi(smi)}",
                caption=mol_caption(smi),
            )

        # Separator
        sep_y = gs[1, 0].get_position(fig).y1 + 0.012
        fig.add_artist(plt.Line2D(
            [0.02, 0.98], [sep_y, sep_y],
            transform=fig.transFigure,
            color="#888888", linewidth=0.8, linestyle="--",
        ))

        # Polymorphic neighbors grid
        for j, row in poly_df.iterrows():
            r = j // n_pcols
            c = j % n_pcols
            ax = fig.add_subplot(gs[1 + r, c])
            # One tanimoto per mono candidate (in display order)
            tan_vals = [pivot.loc[row["smiles_2"], smi] for smi in mono_list]
            draw_mol_box(
                ax, row["smiles_2"],
                bg_color=COLOR_POLY,
                title=fmt_smi(row["smiles_2"]),
                caption=mol_caption(
                    row["smiles_2"],
                    nb_forms=int(row["nb_forms_2"]),
                    tanimoto_vals=tan_vals,
                ),
            )

        pdf.savefig(fig, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [{page_num+1:>2}/{len(pages)}] {n_mono} mono · {n_poly} poly")

print(f"\nDone → {OUTPUT_PDF}")
