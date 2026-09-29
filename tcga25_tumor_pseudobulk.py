# ─────────────────────────────────────────────────────────────
# RUN ON: local Mac (Elements drive mounted)
#   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
#   python3 tcga25_tumor_pseudobulk.py
# Tumor-enriched pseudobulk for the 25 TCGA cetuximab patients.
# Same logic as Win4221/DeepSpotM/Code/tumor_pseudobulk.py:
#   tumor spot = top 25% by mean of KRT17, KRT14, KRT5, KRT6A, TP63, SFN
#   pseudobulk = mean of all genes over tumor spots
# Input:  /Volumes/Elements/DeepspotM/TCGA/*.h5ad.gz
#         ../Result/cetuximab_25_OS_groups.csv  (from KM_cetuximab25_OS_risk.R)
# Output: ../Result/tcga25_tumor_pseudobulk.csv
# ─────────────────────────────────────────────────────────────
import anndata as ad
import numpy as np
import pandas as pd
import os, glob, gzip, shutil, tempfile
from scipy.sparse import issparse

SLIDE_DIR = "/Volumes/Elements/DeepspotM/TCGA"
BASE = "/Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM"
GROUP_CSV = os.path.join(BASE, "Result", "cetuximab_25_OS_groups.csv")
OUTPUT_CSV = os.path.join(BASE, "Result", "tcga25_tumor_pseudobulk.csv")

TUMOR_GENES = ["KRT17", "KRT14", "KRT5", "KRT6A", "TP63", "SFN"]
QUANTILE = 0.75  # top 25%

patients = pd.read_csv(GROUP_CSV)["patient_id"].tolist()
print(f"{len(patients)} patients in {os.path.basename(GROUP_CSV)}")

results = []
for pid in patients:
    hits = sorted(glob.glob(os.path.join(SLIDE_DIR, pid + "*.h5ad.gz")))
    if not hits:
        print(f"{pid}: no slide, SKIP")
        continue
    f = hits[0]
    print(f"Processing: {pid} ...", end=" ", flush=True)

    # some ".h5ad.gz" files are really plain h5ad (HDF5), check magic bytes
    with open(f, "rb") as fh:
        is_gz = fh.read(2) == b"\x1f\x8b"
    if is_gz:
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            with gzip.open(f, "rb") as gz:
                shutil.copyfileobj(gz, tmp)
            tmp_path = tmp.name
        try:
            adata = ad.read_h5ad(tmp_path)
        finally:
            os.remove(tmp_path)
    else:
        adata = ad.read_h5ad(f)

    available = [g for g in TUMOR_GENES if g in adata.var_names]
    if len(available) == 0:
        print(f"SKIP - none of {TUMOR_GENES} found")
        continue
    print(f"({len(available)}/{len(TUMOR_GENES)} panel genes) ", end="")

    X = adata[:, available].X
    if issparse(X):
        X = X.toarray()
    score = np.mean(X, axis=1)
    tumor_mask = score >= np.quantile(score, QUANTILE)
    n_tumor, n_total = int(tumor_mask.sum()), adata.n_obs
    print(f"tumor spots: {n_tumor}/{n_total}")

    X_tumor = adata[tumor_mask].X
    if issparse(X_tumor):
        X_tumor = X_tumor.toarray()
    row = pd.Series(np.asarray(X_tumor.mean(axis=0)).flatten(),
                    index=adata.var_names, name=pid)
    row["n_tumor_spots"] = n_tumor
    row["n_total_spots"] = n_total
    results.append(row)

df = pd.DataFrame(results)
meta_cols = ["n_tumor_spots", "n_total_spots"]
df = df[meta_cols + [c for c in df.columns if c not in meta_cols]]
df.to_csv(OUTPUT_CSV)
print(f"\nSaved: {OUTPUT_CSV}")
print(f"Shape: {df.shape[0]} samples x {df.shape[1]} columns")
