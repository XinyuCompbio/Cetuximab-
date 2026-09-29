# ─────────────────────────────────────────────────────────────
# RUN ON: local Mac (Elements drive mounted)
#   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
#   python3 phoenix_resp_tumor_pseudobulk.py
# Phoenix 112um tumor-region pseudobulk for the cetuximab responder /
# non-responder patients. Same tumor logic as hnsc_slides/hpv_tumor_DE.py:
#   tumor score = mean z-score of EPCAM, COL17A1, LY6D, SOX2 (377-gene panel)
#   tumor spots = top 25% per slide; pseudobulk = mean over tumor spots
#   (values stay in Phoenix log1p space); slides averaged per patient.
# Only the labeled patients are read, so it is quick.
# Input : /Volumes/Elements/hnsc_slides/phoenix_112um/<slide>/phoenix_pred.h5
#         ../tcga_hnsc_cetuximab_patients.csv
# Output: ../Result/phoenix_resp_tumor_pseudobulk.csv  (patients x genes + group)
# ─────────────────────────────────────────────────────────────
import glob, os
import h5py
import numpy as np
import pandas as pd

RESULTS_DIR = "/Volumes/Elements/hnsc_slides/phoenix_112um"
BASE = "/Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM"
RESP_CSV = os.path.join(BASE, "tcga_hnsc_cetuximab_patients.csv")
OUT_CSV = os.path.join(BASE, "Result", "phoenix_resp_tumor_pseudobulk.csv")

TUMOR_MARKERS = ["EPCAM", "COL17A1", "LY6D", "SOX2"]
TUMOR_QUANTILE = 0.75

resp = pd.read_csv(RESP_CSV)
resp["group"] = np.where(resp.cetuximab_response.str.contains("Complete Response|Partial Response"), "R",
                np.where(resp.cetuximab_response.str.contains("Progressive Disease|Stable Disease"), "NR", ""))
resp = resp[resp.group != ""].set_index("patient_id")


def tumor_pseudobulk(h5path):
    with h5py.File(h5path, "r", locking=False) as f:
        if "X" not in f or "genes" not in f:
            return None, None
        genes = [g.decode() if isinstance(g, bytes) else str(g) for g in f["genes"][:]]
        ndone = int(f.attrs["n_done"]) if "n_done" in f.attrs else f["X"].shape[0]
        ndone = max(1, min(ndone, f["X"].shape[0]))
        X = f["X"][:ndone]
    if X.shape[1] != len(genes):
        X = X.T
    gidx = {g: i for i, g in enumerate(genes)}
    mk = [gidx[g] for g in TUMOR_MARKERS if g in gidx]
    if len(mk) < 2 or X.shape[0] < 50:
        return None, genes
    Z = X[:, mk]
    Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-8)
    score = Z.mean(1)
    tum = X[score >= np.quantile(score, TUMOR_QUANTILE)]
    if tum.shape[0] < 20:
        return None, genes
    return tum.mean(0), genes


rows, bcs, genes_ref = [], [], None
for pid in resp.index:
    for h5 in sorted(glob.glob(os.path.join(RESULTS_DIR, pid + "*", "phoenix_pred.h5"))):
        try:
            v, genes = tumor_pseudobulk(h5)
        except Exception as e:
            print("  skip", h5, e); continue
        if v is None:
            print("  skip", os.path.basename(os.path.dirname(h5))); continue
        rows.append(v); bcs.append(pid); genes_ref = genes
        print(f"{pid} ({resp.loc[pid, 'group']}) ok")

df = pd.DataFrame(np.vstack(rows), columns=genes_ref, index=bcs).groupby(level=0).mean()
df.insert(0, "group", resp.loc[df.index, "group"])
df.to_csv(OUT_CSV)
print(f"\nNR={int((df.group == 'NR').sum())}, R={int((df.group == 'R').sum())}, genes={df.shape[1] - 1}")
print("Saved:", OUT_CSV)
