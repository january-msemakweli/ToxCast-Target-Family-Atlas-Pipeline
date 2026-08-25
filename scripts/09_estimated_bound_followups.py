"""
Follow-up analyses restricted to chemicals with an estimated cytotoxicity bound.

1) Family rank stability (actives vs cytotoxicity-separated actives).
2) Family-level chemical label-set shift (Jaccard) for the display families.

Outputs T16, T17 and derived/ESTIMATED_BOUND_FOLLOWUPS.{txt,json}.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
TOXCAST = ROOT.parent / "ToxCast" / "derived" / "summary_extract"
MC56 = TOXCAST / "mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv"
ANN = TOXCAST / "assay_annotations_invitrodb_v4_3_AUG2024.xlsx"
CYTO = TOXCAST / "cytotox_invitrodb_v4_3_AUG2024.xlsx"
TAB = ROOT / "tables"
DER = ROOT / "derived"

HIT_THRESH = 0.9
DEFAULT_BOUND_UM = 1000.0
DEFAULT_TOL = 0.1
MARGINS = {"1x": 1.0, "3x": 3.0, "10x": 10.0}
CONTROL_FAMILIES = {"background measurement", "background control"}
FAMILY_RENAME = {
    "channel 1": "ion channel (panel 1)",
    "channel 2": "ion channel (panel 2)",
}
DISPLAY_FAMS = [
    "cell cycle",
    "nuclear receptor",
    "cyp",
    "cytokine",
    "gpcr",
    "kinase",
    "protease",
    "neurodevelopment",
    "steroid hormone",
    "apoptosis",
]


def main():
    ann = pd.read_excel(ANN, usecols=["aeid", "intended_target_family"]).dropna(subset=["aeid"])
    ann["aeid"] = ann["aeid"].astype(np.int64)
    fam = ann["intended_target_family"].fillna("unannotated").astype(str).str.strip().str.lower()
    ann["family_display"] = fam.map(lambda x: FAMILY_RENAME.get(x, x))
    ann["is_ctrl"] = fam.isin(CONTROL_FAMILIES)
    fam_by_aeid = ann.set_index("aeid")["family_display"].to_dict()
    ctrl_by_aeid = ann.set_index("aeid")["is_ctrl"].to_dict()

    cy = pd.read_excel(CYTO, usecols=["chid", "cytotox_lower_bound_um"]).dropna(subset=["chid"])
    cy["chid"] = cy["chid"].astype(np.int64)
    bound = cy.set_index("chid")["cytotox_lower_bound_um"].to_dict()
    estimated = {
        int(c)
        for c, b in bound.items()
        if pd.notna(b) and b < DEFAULT_BOUND_UM - DEFAULT_TOL
    }

    fam_chem_active: dict[str, set[int]] = defaultdict(set)
    fam_chem_sep: dict[str, dict[str, set[int]]] = {m: defaultdict(set) for m in MARGINS}

    n_rows = 0
    for ci, chunk in enumerate(
        pd.read_csv(
            MC56,
            usecols=["chid", "aeid", "hitc", "ac50"],
            chunksize=500_000,
            low_memory=False,
        )
    ):
        chunk = chunk.dropna(subset=["chid", "aeid", "hitc"])
        if chunk.empty:
            continue
        n_rows += len(chunk)
        hit = chunk["hitc"].to_numpy(dtype=float) >= HIT_THRESH
        if not hit.any():
            continue
        h = chunk.loc[hit]
        chids = h["chid"].to_numpy(dtype=np.int64)
        aeids = h["aeid"].to_numpy(dtype=np.int64)
        ac = h["ac50"].to_numpy(dtype=float)
        keep = np.array([int(c) in estimated for c in chids])
        if not keep.any():
            continue
        chids, aeids, ac = chids[keep], aeids[keep], ac[keep]
        fams = np.array([str(fam_by_aeid.get(int(a), "unannotated")) for a in aeids])
        ctrl = np.array([bool(ctrl_by_aeid.get(int(a), False)) for a in aeids])
        bnd = np.array([bound.get(int(c), np.nan) for c in chids], dtype=float)
        okac = np.isfinite(ac) & (ac > 0)

        for f, c, k in zip(fams[~ctrl], chids[~ctrl], np.ones(int((~ctrl).sum()), bool)):
            fam_chem_active[f].add(int(c))
        for m, div in MARGINS.items():
            sel = okac & np.isfinite(bnd) & (ac < bnd / div) & ~ctrl
            for f, c in zip(fams[sel], chids[sel]):
                fam_chem_sep[m][f].add(int(c))
        if ci % 4 == 0:
            print(f"  streamed {n_rows:,} rows ...", flush=True)

    # ---- T16: rank stability within estimated-bound stratum ----
    t13 = pd.read_csv(TAB / "T13_family_by_bound_provenance.csv")
    bio = t13.loc[t13["family"].ne("unannotated")].copy()
    a = bio["n_active_estimated_bound"]
    s = bio["n_separated_1x_estimated_bound"]
    rho, p = stats.spearmanr(a.rank(ascending=False), s.rank(ascending=False))

    def jac(n):
        x = set(bio.nlargest(n, "n_active_estimated_bound")["family"])
        y = set(bio.nlargest(n, "n_separated_1x_estimated_bound")["family"])
        return len(x & y) / len(x | y)

    stability = {
        "n_families": int(len(bio)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "jaccard_top10": float(jac(10)),
        "jaccard_top20": float(jac(20)),
    }
    pd.DataFrame([stability]).to_csv(
        TAB / "T16_rank_stability_estimated_bound.csv", index=False
    )

    # ---- T17: label-set shift ----
    rows = []
    for f in DISPLAY_FAMS:
        act = fam_chem_active.get(f, set())
        row = {"family": f, "n_chemicals_active": len(act)}
        for m in MARGINS:
            sep = fam_chem_sep[m].get(f, set())
            union = act | sep
            row[f"n_chemicals_separated_{m}"] = len(sep)
            row[f"jaccard_{m}"] = len(act & sep) / len(union) if union else np.nan
            row[f"frac_retained_{m}"] = len(act & sep) / len(act) if act else np.nan
        rows.append(row)
    t17 = pd.DataFrame(rows)
    t17.to_csv(TAB / "T17_label_shift_estimated_bound.csv", index=False)

    out = {"rank_stability": stability, "label_shift": t17.to_dict(orient="records")}
    (DER / "ESTIMATED_BOUND_FOLLOWUPS.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )

    lines = [
        "Follow-ups restricted to estimated-bound chemicals (invitrodb v4.3)",
        "",
        "Family rank stability (actives vs cytotoxicity-separated, 1x):",
        f"  families: {stability['n_families']}",
        f"  Spearman rho: {stability['spearman_rho']:.3f} (p={stability['spearman_p']:.3g})",
        f"  Jaccard top10/top20: {stability['jaccard_top10']:.3f}/{stability['jaccard_top20']:.3f}",
        "",
        "Family chemical label-set shift (Jaccard of active vs separated chemical sets):",
    ]
    for r in rows:
        lines.append(
            f"  {r['family']}: n_active_chem={r['n_chemicals_active']}, "
            f"J1x={r['jaccard_1x']:.3f}, J3x={r['jaccard_3x']:.3f}, J10x={r['jaccard_10x']:.3f}, "
            f"retained1x={100*r['frac_retained_1x']:.1f}%"
        )
    (DER / "ESTIMATED_BOUND_FOLLOWUPS.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("\n".join(lines))


if __name__ == "__main__":
    main()
