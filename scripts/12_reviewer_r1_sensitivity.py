"""
Reviewer-driven sensitivity analyses for the Computational Toxicology R1.

1. Burst-hit support (nhit) of estimated vs default bounds, and 1x AC50
   separation after restricting to better-supported estimates.
2. Alternative potency metrics (ACC, AC20, BMD) vs the published bound.
3. How often AC50 itself lies outside the tested concentration range.
4. Full-library percentiles of the highest tested concentration.
5. Coverage-normalised promiscuity when a family counts as tested only if
   the chemical has >=k usable curves in that family.

Primary AC50 / 1x numbers are recomputed as a sanity check against T12.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

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
USECOLS = [
    "chid",
    "aeid",
    "hitc",
    "ac50",
    "acc",
    "ac20",
    "bmd",
    "conc_min",
    "conc_max",
]
FAMILY_RENAME = {
    "channel 1": "ion channel (panel 1)",
    "channel 2": "ion channel (panel 2)",
}
CONTROL_FAMILIES = {"background measurement", "background control"}
NHIT_BINS = [
    ("2-4", 2, 4),
    ("5-9", 5, 9),
    ("10-19", 10, 19),
    ("20+", 20, 10_000),
]
NHIT_MINIMA = [2, 5, 10, 20]
DEPTHS = [1, 2, 3]
METRICS = [("ac50", "AC50"), ("acc", "ACC"), ("ac20", "AC20"), ("bmd", "BMD")]


def load_ann() -> pd.DataFrame:
    ann = pd.read_excel(
        ANN,
        usecols=["aeid", "intended_target_family"],
    ).dropna(subset=["aeid"])
    ann["aeid"] = ann["aeid"].astype(np.int64)
    fam = ann["intended_target_family"].fillna("unannotated").astype(str).str.strip().str.lower()
    ann["family_display"] = fam.map(lambda x: FAMILY_RENAME.get(x, x))
    ann["is_control"] = fam.isin(CONTROL_FAMILIES)
    return ann


def load_cyto() -> pd.DataFrame:
    cy = pd.read_excel(
        CYTO,
        usecols=[
            "chid",
            "cytotox_median_um",
            "cytotox_lower_bound_um",
            "ntested",
            "nhit",
        ],
    ).dropna(subset=["chid"])
    cy["chid"] = cy["chid"].astype(np.int64)
    cy["is_default"] = cy["cytotox_lower_bound_um"] >= (DEFAULT_BOUND_UM - DEFAULT_TOL)
    return cy


def pctiles(arr: np.ndarray, qs=(5, 10, 25, 50, 75, 90, 95, 99)) -> dict:
    out = {}
    for q in qs:
        out[f"p{q}"] = float(np.percentile(arr, q))
    out["min"] = float(arr.min())
    out["max"] = float(arr.max())
    out["n"] = int(arr.size)
    return out


def stream(ann: pd.DataFrame, cy: pd.DataFrame) -> dict:
    fam_by_aeid = ann.set_index("aeid")["family_display"].to_dict()
    ctrl_by_aeid = ann.set_index("aeid")["is_control"].to_dict()
    bound = cy.set_index("chid")["cytotox_lower_bound_um"].to_dict()
    is_default = cy.set_index("chid")["is_default"].to_dict()
    nhit = cy.set_index("chid")["nhit"].to_dict()

    n_active = Counter()
    n_sep_metric = {m: Counter() for m, _ in METRICS}  # metric -> stratum -> count
    n_finite_metric = {m: Counter() for m, _ in METRICS}
    n_ac50_outside = Counter()
    n_ac50_finite = Counter()

    conc_max_vals = []
    chem_n_active = Counter()
    chem_n_sep_ac50 = Counter()

    chem_fam_n: dict[int, Counter] = defaultdict(Counter)
    chem_fam_active: dict[int, set[str]] = defaultdict(set)
    chem_fam_sep: dict[int, set[str]] = defaultdict(set)

    n_rows = 0
    for ci, chunk in enumerate(
        pd.read_csv(MC56, usecols=USECOLS, chunksize=500_000, low_memory=False)
    ):
        chunk = chunk.dropna(subset=["chid", "aeid", "hitc"]).copy()
        if chunk.empty:
            continue
        chunk["chid"] = chunk["chid"].astype(np.int64)
        chunk["aeid"] = chunk["aeid"].astype(np.int64)
        n_rows += len(chunk)

        families = chunk["aeid"].map(fam_by_aeid).fillna("unannotated").astype(str)
        is_ctrl = chunk["aeid"].map(ctrl_by_aeid).fillna(False).astype(bool)
        bio = ~is_ctrl.to_numpy()
        for chid, fam in zip(chunk["chid"].to_numpy()[bio], families.to_numpy()[bio]):
            chem_fam_n[int(chid)][str(fam)] += 1

        hit = chunk["hitc"].to_numpy(dtype=float) >= HIT_THRESH
        if not hit.any():
            continue
        h = chunk.loc[hit]
        hf = families.loc[hit].to_numpy()
        hc = is_ctrl.loc[hit].to_numpy()
        chids = h["chid"].to_numpy(dtype=np.int64)
        ac50 = h["ac50"].to_numpy(dtype=float)
        acc = h["acc"].to_numpy(dtype=float)
        ac20 = h["ac20"].to_numpy(dtype=float)
        bmd = h["bmd"].to_numpy(dtype=float)
        cmin = h["conc_min"].to_numpy(dtype=float)
        cmax = h["conc_max"].to_numpy(dtype=float)
        bnd = np.array([bound.get(int(c), np.nan) for c in chids], dtype=float)
        dflt = np.array([bool(is_default.get(int(c), False)) for c in chids])
        strata = np.where(dflt, "default", "estimated")

        for s, cnt in pd.Series(strata).value_counts().items():
            n_active[str(s)] += int(cnt)

        have_cmax = np.isfinite(cmax) & (cmax > 0)
        if have_cmax.any():
            conc_max_vals.append(cmax[have_cmax].astype(np.float64, copy=False))

        valid_ac50 = np.isfinite(ac50) & (ac50 > 0)
        have_range = np.isfinite(cmin) & np.isfinite(cmax)
        outside = valid_ac50 & have_range & ((ac50 < cmin) | (ac50 > cmax))
        for s, cnt in pd.Series(strata[valid_ac50]).value_counts().items():
            n_ac50_finite[str(s)] += int(cnt)
        for s, cnt in pd.Series(strata[outside]).value_counts().items():
            n_ac50_outside[str(s)] += int(cnt)

        metric_arr = {"ac50": ac50, "acc": acc, "ac20": ac20, "bmd": bmd}
        for key, arr in metric_arr.items():
            finite = np.isfinite(arr) & (arr > 0)
            sep = finite & np.isfinite(bnd) & (arr < bnd)
            for s, cnt in pd.Series(strata[finite]).value_counts().items():
                n_finite_metric[key][str(s)] += int(cnt)
            for s, cnt in pd.Series(strata[sep]).value_counts().items():
                n_sep_metric[key][str(s)] += int(cnt)

        sep_ac50 = valid_ac50 & np.isfinite(bnd) & (ac50 < bnd)
        for cid, s, s1 in zip(chids, strata, sep_ac50):
            chem_n_active[int(cid)] += 1
            if s1:
                chem_n_sep_ac50[int(cid)] += 1

        for cid, fam, isbio, s1 in zip(chids, hf, ~hc, sep_ac50):
            if not isbio:
                continue
            chem_fam_active[int(cid)].add(str(fam))
            if s1:
                chem_fam_sep[int(cid)].add(str(fam))

        if ci % 4 == 0:
            print(f"  streamed {n_rows:,} rows ...", flush=True)

    conc = np.concatenate(conc_max_vals) if conc_max_vals else np.array([], dtype=float)
    return {
        "n_rows": n_rows,
        "n_active": n_active,
        "n_sep_metric": n_sep_metric,
        "n_finite_metric": n_finite_metric,
        "n_ac50_outside": n_ac50_outside,
        "n_ac50_finite": n_ac50_finite,
        "conc_max": conc,
        "chem_n_active": chem_n_active,
        "chem_n_sep_ac50": chem_n_sep_ac50,
        "chem_fam_n": chem_fam_n,
        "chem_fam_active": chem_fam_active,
        "chem_fam_sep": chem_fam_sep,
        "nhit": nhit,
        "is_default": is_default,
    }


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    DER.mkdir(parents=True, exist_ok=True)

    print("Loading annotations / cytotox ...", flush=True)
    ann = load_ann()
    cy = load_cyto()
    est_cy = cy.loc[~cy["is_default"]].copy()
    dfl_cy = cy.loc[cy["is_default"]].copy()
    print(
        f"  cytotox n={len(cy):,} estimated={len(est_cy):,} default={len(dfl_cy):,}",
        flush=True,
    )

    # ----- T20: nhit distribution (no streaming needed) -----
    t20_rows = []
    for label, sub in (("estimated", est_cy), ("default", dfl_cy)):
        nh = sub["nhit"].to_numpy(dtype=float)
        row = {
            "bound_provenance": label,
            "n_chemicals": int(len(sub)),
            "nhit_min": int(np.nanmin(nh)),
            "nhit_max": int(np.nanmax(nh)),
            "nhit_median": float(np.nanmedian(nh)),
            "nhit_p25": float(np.nanpercentile(nh, 25)),
            "nhit_p75": float(np.nanpercentile(nh, 75)),
            "n_nhit_0": int((nh == 0).sum()),
            "n_nhit_1": int((nh == 1).sum()),
            "n_nhit_2_to_4": int(((nh >= 2) & (nh <= 4)).sum()),
            "n_nhit_5_to_9": int(((nh >= 5) & (nh <= 9)).sum()),
            "n_nhit_10_to_19": int(((nh >= 10) & (nh <= 19)).sum()),
            "n_nhit_20plus": int((nh >= 20).sum()),
            "ntested_median": float(sub["ntested"].median()),
        }
        t20_rows.append(row)
    t20 = pd.DataFrame(t20_rows)
    t20.to_csv(TAB / "T20_nhit_distribution.csv", index=False)

    print("Streaming mc5-6 ...", flush=True)
    agg = stream(ann, cy)

    # sanity vs T12
    act_est = agg["n_active"]["estimated"]
    sep_est = agg["n_sep_metric"]["ac50"]["estimated"]
    print(
        f"  estimated actives={act_est:,}  AC50 1x separated="
        f"{sep_est:,} ({100 * sep_est / act_est:.1f}%)",
        flush=True,
    )

    # ----- T21: nhit sensitivity of 1x AC50 separation -----
    t21_rows = []
    for lo, hi_label, hi in [(2, "4", 4), (5, "9", 9), (10, "19", 19), (20, "+", 10_000)]:
        label = f"{lo}-{hi_label}" if hi_label != "+" else "20+"
        chems = set(est_cy.loc[est_cy["nhit"].between(lo, hi), "chid"].astype(int))
        n_act = sum(agg["chem_n_active"][c] for c in chems)
        n_sep = sum(agg["chem_n_sep_ac50"][c] for c in chems)
        t21_rows.append(
            {
                "subset": f"estimated-bound, nhit {label}",
                "n_chemicals": len(chems),
                "n_active": n_act,
                "n_separated_1x_ac50": n_sep,
                "pct_separated_1x_ac50": 100 * n_sep / n_act if n_act else np.nan,
            }
        )
    for m in NHIT_MINIMA:
        chems = set(est_cy.loc[est_cy["nhit"] >= m, "chid"].astype(int))
        n_act = sum(agg["chem_n_active"][c] for c in chems)
        n_sep = sum(agg["chem_n_sep_ac50"][c] for c in chems)
        t21_rows.append(
            {
                "subset": f"estimated-bound, nhit >= {m}",
                "n_chemicals": len(chems),
                "n_active": n_act,
                "n_separated_1x_ac50": n_sep,
                "pct_separated_1x_ac50": 100 * n_sep / n_act if n_act else np.nan,
            }
        )
    t21 = pd.DataFrame(t21_rows)
    t21.to_csv(TAB / "T21_nhit_sensitivity.csv", index=False)

    # ----- T22: potency metric sensitivity -----
    t22_rows = []
    for key, label in METRICS:
        for stratum in ("estimated", "default"):
            n_act = agg["n_active"][stratum]
            n_fin = agg["n_finite_metric"][key][stratum]
            n_sep = agg["n_sep_metric"][key][stratum]
            t22_rows.append(
                {
                    "potency_metric": label,
                    "bound_provenance": stratum,
                    "n_active": n_act,
                    "n_finite_metric": n_fin,
                    "pct_actives_with_finite_metric": 100 * n_fin / n_act if n_act else np.nan,
                    "n_separated_1x": n_sep,
                    "pct_separated_of_actives": 100 * n_sep / n_act if n_act else np.nan,
                    "pct_separated_of_finite": 100 * n_sep / n_fin if n_fin else np.nan,
                }
            )
    t22 = pd.DataFrame(t22_rows)
    t22.to_csv(TAB / "T22_potency_metric_sensitivity.csv", index=False)

    # ----- T23: conc_max percentiles -----
    conc = agg["conc_max"]
    p = pctiles(conc)
    t23 = pd.DataFrame(
        [
            {
                "quantity": "highest tested concentration (uM), active curves",
                "n": p["n"],
                "min": p["min"],
                "p5": p["p5"],
                "p10": p["p10"],
                "p25": p["p25"],
                "p50": p["p50"],
                "p75": p["p75"],
                "p90": p["p90"],
                "p95": p["p95"],
                "p99": p["p99"],
                "max": p["max"],
                "pct_tested_to_1000uM": 100 * float((conc >= 999.9).mean()),
            }
        ]
    )
    t23.to_csv(TAB / "T23_conc_max_percentiles.csv", index=False)

    # AC50 outside tested range
    ac50_out = {
        s: {
            "n_finite": agg["n_ac50_finite"][s],
            "n_outside": agg["n_ac50_outside"][s],
            "pct_outside_of_finite": 100 * agg["n_ac50_outside"][s] / agg["n_ac50_finite"][s]
            if agg["n_ac50_finite"][s]
            else np.nan,
        }
        for s in ("estimated", "default")
    }

    # ----- T24: promiscuity depth -----
    est_ids = set(est_cy["chid"].astype(int))
    dfl_ids = set(dfl_cy["chid"].astype(int))
    t24_rows = []
    assays_per_tested_family = []
    for stratum, ids in (("estimated", est_ids), ("default", dfl_ids)):
        # chemicals with at least one biological active family
        chems = [c for c in ids if agg["chem_fam_active"][c]]
        for k in DEPTHS:
            fr_act = []
            fr_sep = []
            n_used = 0
            for c in chems:
                fam_n = agg["chem_fam_n"][c]
                tested = {f for f, n in fam_n.items() if n >= k}
                if not tested:
                    continue
                n_used += 1
                act = agg["chem_fam_active"][c] & tested
                sep = agg["chem_fam_sep"][c] & tested
                fr_act.append(len(act) / len(tested))
                fr_sep.append(len(sep) / len(tested))
                if k == 1 and stratum == "estimated":
                    assays_per_tested_family.extend(fam_n[f] for f in tested)
            t24_rows.append(
                {
                    "bound_provenance": stratum,
                    "min_assays_to_count_family_tested": k,
                    "n_chemicals": n_used,
                    "median_frac_tested_active": float(np.median(fr_act)) if fr_act else np.nan,
                    "median_frac_tested_separated": float(np.median(fr_sep)) if fr_sep else np.nan,
                }
            )
    t24 = pd.DataFrame(t24_rows)
    t24.to_csv(TAB / "T24_promiscuity_depth.csv", index=False)

    apf = np.array(assays_per_tested_family, dtype=float)
    depth_summary = {
        "estimated_median_assays_per_tested_family": float(np.median(apf)) if apf.size else np.nan,
        "estimated_p25_assays_per_tested_family": float(np.percentile(apf, 25)) if apf.size else np.nan,
        "estimated_p75_assays_per_tested_family": float(np.percentile(apf, 75)) if apf.size else np.nan,
        "estimated_pct_tested_families_with_1_assay": float(100 * (apf == 1).mean()) if apf.size else np.nan,
    }

    summary = {
        "n_rows": agg["n_rows"],
        "n_active": dict(agg["n_active"]),
        "ac50_1x_pct_estimated": 100 * sep_est / act_est if act_est else None,
        "t20": t20.to_dict(orient="records"),
        "t21": t21.to_dict(orient="records"),
        "t22": t22.to_dict(orient="records"),
        "t23": t23.to_dict(orient="records"),
        "t24": t24.to_dict(orient="records"),
        "ac50_outside_tested_range": ac50_out,
        "assays_per_tested_family": depth_summary,
        "nhit_estimated": {
            "min": int(est_cy["nhit"].min()),
            "median": float(est_cy["nhit"].median()),
            "p25": float(est_cy["nhit"].quantile(0.25)),
            "p75": float(est_cy["nhit"].quantile(0.75)),
            "max": int(est_cy["nhit"].max()),
            "n_nhit_eq_2": int((est_cy["nhit"] == 2).sum()),
            "n_nhit_eq_3": int((est_cy["nhit"] == 3).sum()),
            "n_nhit_eq_4": int((est_cy["nhit"] == 4).sum()),
        },
        "nhit_default": {
            "min": int(dfl_cy["nhit"].min()),
            "median": float(dfl_cy["nhit"].median()),
            "max": int(dfl_cy["nhit"].max()),
            "n_nhit_0": int((dfl_cy["nhit"] == 0).sum()),
            "n_nhit_1": int((dfl_cy["nhit"] == 1).sum()),
            "n_nhit_ge_2": int((dfl_cy["nhit"] >= 2).sum()),
        },
    }
    (DER / "R1_SENSITIVITY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n=== T20 nhit distribution ===")
    print(t20.to_string(index=False))
    print("\n=== T21 nhit sensitivity ===")
    print(t21.to_string(index=False))
    print("\n=== T22 potency metrics ===")
    print(t22.to_string(index=False))
    print("\n=== T23 conc_max percentiles ===")
    print(t23.to_string(index=False))
    print("\n=== T24 promiscuity depth ===")
    print(t24.to_string(index=False))
    print("\n=== AC50 outside tested range ===")
    print(json.dumps(ac50_out, indent=2))
    print("\n=== assays per tested family (estimated, k=1) ===")
    print(json.dumps(depth_summary, indent=2))
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
