#!/usr/bin/env python3
"""Second standard figure -- change above detection limits -- for the NGV-corrected DoD,
and an honest look at whether the detection threshold improved.

The `stable` mask is not saved by difference_dem, so it is REBUILT here with the pipeline's
own lines (pipeline.py 515-526) from the saved z_after.npy. Same code, not an approximation.

Reported, because it is the real question: the correction carries its own uncertainty. Its
coefficient has a block-bootstrap SE of 44.4 mm per unit NGV, and that scales with NGV, so it
is largest exactly where the correction is largest. A corrected LoD must include it.

    ./lidar-icp/bin/python analysis/detect_change_ngv.py
"""
import argparse, json, os
import numpy as np
from scipy.ndimage import uniform_filter, gaussian_filter, distance_transform_edt as edt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lidar_diff_icp import coreg
from lidar_diff_icp.detect import detect_change_standard
from lidar_diff_icp.viz import hillshade

B_SE_MM_PER_NGV = 44.4      # block bootstrap, analysis/ngv.py --marks
Z95 = 1.96


def stable_mask(Z21, res):
    """Verbatim from pipeline.py: terrain masks from the reference ground."""
    Zf = Z21.copy(); nanm = np.isnan(Zf)
    if nanm.any():
        Zf = Zf[tuple(edt(nanm, return_distances=False, return_indices=True))]
    tpi = Z21 - uniform_filter(Zf, size=int(2 * 300 / res), mode="nearest")
    sdeg = np.degrees(coreg.slope_aspect(gaussian_filter(Zf, 2.0), res)[0])
    Zsm = gaussian_filter(Zf, 50 / res / 2)
    lap = (np.gradient(np.gradient(Zsm, res, axis=0), res, axis=0)
           + np.gradient(np.gradient(Zsm, res, axis=1), res, axis=1))
    convex = (sdeg > 5) & (sdeg < 35) & (tpi > -2) & (lap < 0)
    return ((sdeg < 3) & (tpi > -2)) | convex


def nmad(v):
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if v.size else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile", default="data/derived/elba_fulldensity")
    ap.add_argument("--out", default="figures/elba_change_ngv.png")
    ap.add_argument("--v", type=float, default=0.3)
    a = ap.parse_args()
    T = a.tile
    cfg = json.load(open(os.path.join(T, "corrections.json")))
    X0, Y0, X1, Y1 = cfg["bounds"]; res = float(cfg["res_m"])

    dod = np.load(f"{T}/dod.npy")
    corr = np.load(f"{T}/dod_ngv.npy")
    ngv = np.load(f"{T}/ngv.npy")
    lod = np.load(f"{T}/lod.npy")
    Z21 = np.load(f"{T}/z_after.npy")
    ny, nx = dod.shape
    st = stable_mask(Z21, res) & np.isfinite(dod)
    print(f"stable mask rebuilt from z_after.npy: {int(st.sum()):,} cells "
          f"({100 * st.sum() / np.isfinite(dod).sum():.1f}% of cells with a DoD)")

    print(f"\nDID THE CORRECTION IMPROVE THE THRESHOLD? -- scatter on STABLE ground,")
    print(f"which is what sets the LoD. If the correction removes a real bias, this falls.")
    print(f"  {'':22s} {'median':>9} {'NMAD':>9}   (mm)")
    for nm, d in (("DoD before", dod), ("DoD after NGV", corr)):
        print(f"  {nm:22s} {1000 * np.median(d[st]):+9.1f} {1000 * nmad(d[st]):9.1f}")
    print(f"  NGV on stable ground   median {np.median(ngv[st & np.isfinite(ngv)]):.3f}"
          f"   p90 {np.percentile(ngv[st & np.isfinite(ngv)], 90):.3f}")

    # corrected LoD: add the correction's own uncertainty in quadrature
    sig = lod / Z95
    sig_c = np.sqrt(sig ** 2 + (B_SE_MM_PER_NGV * np.nan_to_num(ngv) / 1000.0) ** 2)
    lod_c = Z95 * sig_c
    fin = np.isfinite(dod) & np.isfinite(lod)
    print(f"\nLoD, WITH the correction's uncertainty propagated (mm)")
    print(f"  {'':22s} {'median':>9} {'p90':>9} {'max':>9}")
    print(f"  {'lod.npy as built':22s} {1000*np.median(lod[fin]):9.1f} "
          f"{1000*np.percentile(lod[fin],90):9.1f} {1000*np.nanmax(lod[fin]):9.1f}")
    print(f"  {'+ NGV coefficient SE':22s} {1000*np.median(lod_c[fin]):9.1f} "
          f"{1000*np.percentile(lod_c[fin],90):9.1f} {1000*np.nanmax(lod_c[fin]):9.1f}")
    np.save(f"{T}/lod_ngv.npy", lod_c)

    print(f"\nDETECTION, Wheaton coherence + tau_sys floor (detect_change_standard)")
    out = {}
    for nm, d, L in (("before", dod, lod), ("after, old LoD", corr, lod),
                     ("after, LoD+SE", corr, lod_c)):
        det = detect_change_standard(d, L, st, res)
        ch = det["change"]
        vol = sum(r["volume_m3"] for r in det["regions"])
        print(f"  {nm:15s} {len(det['regions']):4d} regions  {int(ch.sum()):7,d} cells "
              f"({100*ch.sum()/fin.sum():5.2f}%)  net {vol:+11,.0f} m3  "
              f"tau_sys {1000*det['tau_sys_m']:5.1f} mm  L {det['corr_length_m']:5.0f} m")
        out[nm] = (det, ch)

    det, change = out["after, LoD+SE"]
    np.save(f"{T}/change_ngv.npy", change)
    hs = hillshade(Z21, res, X0, Y0, fill_gaps=True)
    ext = (X0, X0 + nx * res, Y0, Y0 + ny * res)
    over = np.where(change, corr, np.nan)
    net = sum(r["volume_m3"] for r in det["regions"])
    fig, ax = plt.subplots(figsize=(10, 11))
    ax.imshow(hs, extent=ext, origin="lower", cmap="gray")
    im = ax.imshow(over, extent=ext, origin="lower", cmap="RdBu", vmin=-a.v, vmax=a.v, alpha=0.7)
    ax.set_title("elba: topographic change above detection limits\n"
                 "gen2 corrected for leaf-on vegetation ($-325.2\\,$mm $\\times$ NGV); "
                 f"LoD includes the coefficient's own SE\n"
                 f"{len(det['regions'])} regions, net {net:+,.0f} m$^3$   "
                 "red = erosion, blue = deposition")
    ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
    fig.colorbar(im, ax=ax, shrink=0.6, extend="both", label="detected DoD (m)")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"\nwrote {a.out}")
    print(f"wrote {T}/change_ngv.npy and {T}/lod_ngv.npy")


if __name__ == "__main__":
    main()
