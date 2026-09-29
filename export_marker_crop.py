#!/usr/bin/env python
# ============================================================
# RUN ON: local Mac   (needs: pip install openslide-python openslide-bin;
#                      external drive /Volumes/Elements mounted)
#   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/Win4221/ROSIE/Code
#   python export_marker_crop.py 8-1-30810-2      (Post RT)
#   python export_marker_crop.py 8-2-30810-1      (Post RT+Nivo)
# ============================================================
"""
Step 1 of the representative BCL2 / FoxP3 / LAG3 figure (step 2 = fig_marker_rep.R).
Same crop and H&E as export_mycaf_crop.py (template), markers instead of myCAF.
Exports to ../Result/:
  rep_<slide>_HE.png        H&E from NDPI at ds8 (skipped if already exported)
  rep_<slide>_<marker>.png  ds32: R = marker intensity in tissue, G = tissue, B = DAPI
  rep_<slide>_info.json     crop box, um/px, marker positive fractions
Assumes ROSIE map pixel (r, c) = level-0 pixel (r*32, c*32).
"""
import os, sys, glob, json
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

SLIDE = sys.argv[1] if len(sys.argv) > 1 else "8-1-30810-2"
MARKERS = ["BCL2", "FoxP3", "LAG3"]
NDPI_DIR = "/Volumes/Elements/Win4221"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
MAPS = os.path.join(BASE, "Data", "ds32", SLIDE, "stitched_maps")
OUT = os.path.join(BASE, "Result")
os.makedirs(OUT, exist_ok=True)

DS = 32
HE_UP = 4
BG, POS = 10, 128
PAD = 8
MIN_PX = 5

def load(n):
    f = "coverage_mask.png" if n == "coverage" else f"stitched_{n}.png"
    return np.array(Image.open(os.path.join(MAPS, f)).convert("L"))

asma, panck, dapi = (load(m) for m in ["aSMA", "PanCK", "DAPI"])
tissue = (load("coverage") > 0) & ((dapi > BG) | (panck > BG) | (asma > BG))   # same as myCAF crop

H, W = tissue.shape
rows = np.flatnonzero(tissue.sum(1) >= MIN_PX)
cols = np.flatnonzero(tissue.sum(0) >= MIN_PX)
r1, r2 = max(0, rows[0] - PAD), min(H, rows[-1] + PAD + 1)
c1, c2 = max(0, cols[0] - PAD), min(W, cols[-1] + PAD + 1)

pos_frac = {}
for m in MARKERS:
    mk = load(m)
    sig = np.where(tissue, mk, 0).astype(np.uint8)
    rgb = np.zeros((r2 - r1, c2 - c1, 3), np.uint8)
    rgb[..., 0] = sig[r1:r2, c1:c2]
    rgb[..., 1] = tissue[r1:r2, c1:c2] * 255
    rgb[..., 2] = dapi[r1:r2, c1:c2]
    Image.fromarray(rgb).save(os.path.join(OUT, f"rep_{SLIDE}_{m}.png"))
    pos_frac[m] = float(((mk >= POS) & tissue).sum() / tissue.sum())

he_path = os.path.join(OUT, f"rep_{SLIDE}_HE.png")
ndpi = glob.glob(os.path.join(NDPI_DIR, SLIDE + " - *.ndpi"))[0]
import openslide
sl = openslide.OpenSlide(ndpi)
mpp = float(sl.properties["openslide.mpp-x"])
out_ds = DS / HE_UP
if not os.path.exists(he_path):
    lvl = sl.get_best_level_for_downsample(out_ds)
    lds = sl.level_downsamples[lvl]
    out_w, out_h = (c2 - c1) * HE_UP, (r2 - r1) * HE_UP
    he = Image.new("RGB", (out_w, out_h), "white")
    T = 1024
    for oy in range(0, out_h, T):
        for ox in range(0, out_w, T):
            tw, th = min(T, out_w - ox), min(T, out_h - oy)
            x0 = int(c1 * DS + ox * out_ds); y0 = int(r1 * DS + oy * out_ds)
            rw, rh = int(np.ceil(tw * out_ds / lds)), int(np.ceil(th * out_ds / lds))
            tile = sl.read_region((x0, y0), lvl, (rw, rh))
            bg = Image.new("RGB", tile.size, "white"); bg.paste(tile, mask=tile.split()[3])
            he.paste(bg.resize((tw, th), Image.LANCZOS), (ox, oy))
    he.save(he_path)

info_path = os.path.join(OUT, f"rep_{SLIDE}_info.json")
info = json.load(open(info_path)) if os.path.exists(info_path) else {}
info.update(dict(slide=SLIDE, r1=int(r1), r2=int(r2), c1=int(c1), c2=int(c2),
                 um_per_px_ds32=mpp * DS, he_um_per_px=mpp * out_ds,
                 **{f"slide_{m}_pos_frac": v for m, v in pos_frac.items()}))
json.dump(info, open(info_path, "w"), indent=1)
print(f"{SLIDE}: crop {r2-r1} x {c2-c1} ds32 px; " +
      ", ".join(f"{m} {100*v:.1f}%" for m, v in pos_frac.items()))
