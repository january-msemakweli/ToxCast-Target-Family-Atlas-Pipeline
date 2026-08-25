"""
Assay-format and mc6-flag stratification recomputed within the estimated-bound
stratum, so that these secondary analyses are on the same footing as the
primary ones.

Outputs T18, T19 and derived/FORMAT_FLAGS_ESTIMATED.{txt,json}.
"""
from __future__ import annotations

import json
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
MARGINS = {"1x": 1.0, "3x": 3.0, "10x": 10.0}
SEVERE = {5, 6, 7, 10, 11, 15, 17, 18, 19}


def classify_format(row) -> str:
    fmt = str(row.get("assay_format_type") or "").strip().lower()
    cell = str(row.get("cell_format") or "").strip().lower()
    if fmt in {"biochemical", "cell-free"} or cell in {"cell-free", "tissue-based cell-free"}:
        return "cell_free_biochemical"
    if fmt in {"cell-based", "organism"} or cell in {
        "cell line",
        "primary cell",
        "primary cell co-culture",
        "secondary cell",
        "cell-based",
        "whole embryo",
    }:
        return "cell_based"
    return "other_unclassified"


def parse_severe(series: pd.Series) -> np.ndarray:
    s = series.fillna("").astype(str)
    s = s.str.replace("|", ",", regex=False).str.replace(";", ",", regex=False)
    severe = np.zeros(len(s), dtype=bool)
    for fid in SEVERE:
        severe |= s.str.contains(rf"(?:^|,)\s*{fid}\s*(?:,|$)", regex=True, na=False).to_numpy()
    return severe


def main():
    ann = pd.read_excel(
        ANN, usecols=["aeid", "assay_format_type", "cell_format"]
    ).dropna(subset=["aeid"])
    ann["aeid"] = ann["aeid"].astype(np.int64)
    ann["fmt_class"] = ann.apply(classify_format, axis=1)
    fmt_by_aeid = ann.set_index("aeid")["fmt_class"].to_dict()

    cy = pd.read_excel(CYTO, usecols=["chid", "cytotox_lower_bound_um"]).dropna(subset=["chid"])
    cy["chid"] = cy["chid"].astype(np.int64)
    bound = cy.set_index("chid")["cytotox_lower_bound_um"].to_dict()
    estimated = {
        int(c) for c, b in bound.items()
        if pd.notna(b) and b < DEFAULT_BOUND_UM - DEFAULT_TOL
    }

    fmt_active = {}
    fmt_sep = {m: {} for m in MARGINS}
    n_act = 0
    n_act_nosevere = 0
    n_sep = {m: 0 for m in MARGINS}
    n_sep_nosevere = {m: 0 for m in MARGINS}

    n_rows = 0
    for ci, chunk in enumerate(
        pd.read_csv(
            MC56,
            usecols=["chid", "aeid", "hitc", "ac50", "mc6_flags"],
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
        keep = np.array([int(c) in estimated for c in chids])
        if not keep.any():
            continue
        h = h.loc[keep]
        chids = chids[keep]
        ac = h["ac50"].to_numpy(dtype=float)
        bnd = np.array([bound.get(int(c), np.nan) for c in chids], dtype=float)
        fmts = np.array(
            [str(fmt_by_aeid.get(int(a), "other_unclassified")) for a in h["aeid"].to_numpy()]
        )
        severe = parse_severe(h["mc6_flags"])
        okac = np.isfinite(ac) & (ac > 0)

        n_act += len(h)
        n_act_nosevere += int((~severe).sum())
        for f, c in pd.Series(fmts).value_counts().items():
            fmt_active[f] = fmt_active.get(f, 0) + int(c)
        for m, div in MARGINS.items():
            sel = okac & np.isfinite(bnd) & (ac < bnd / div)
            n_sep[m] += int(sel.sum())
            n_sep_nosevere[m] += int((sel & ~severe).sum())
            for f, c in pd.Series(fmts[sel]).value_counts().items():
                fmt_sep[m][f] = fmt_sep[m].get(f, 0) + int(c)
        if ci % 4 == 0:
            print(f"  streamed {n_rows:,} rows ...", flush=True)

    rows = []
    for f in sorted(fmt_active):
        row = {"format_class": f, "n_active": fmt_active[f]}
        for m in MARGINS:
            s = fmt_sep[m].get(f, 0)
            row[f"n_separated_{m}"] = s
            row[f"pct_separated_{m}"] = 100 * s / fmt_active[f] if fmt_active[f] else np.nan
        rows.append(row)
    t18 = pd.DataFrame(rows)
    t18.to_csv(TAB / "T18_format_stratification_estimated_bound.csv", index=False)

    frows = []
    for m in MARGINS:
        frows.append(
            {
                "margin": m,
                "n_active_all": n_act,
                "n_separated_all": n_sep[m],
                "pct_separated_all": 100 * n_sep[m] / n_act if n_act else np.nan,
                "n_active_no_severe": n_act_nosevere,
                "n_separated_no_severe": n_sep_nosevere[m],
                "pct_separated_no_severe": (
                    100 * n_sep_nosevere[m] / n_act_nosevere if n_act_nosevere else np.nan
                ),
            }
        )
    t19 = pd.DataFrame(frows)
    t19.to_csv(TAB / "T19_flag_sensitivity_estimated_bound.csv", index=False)

    out = {"format": rows, "flags": frows}
    (DER / "FORMAT_FLAGS_ESTIMATED.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    lines = ["Assay format and mc6 flags, estimated-bound chemicals only", "", "Format:"]
    for r in rows:
        lines.append(
            f"  {r['format_class']}: n_active={r['n_active']:,}, "
            f"sep1x={r['pct_separated_1x']:.1f}%, sep3x={r['pct_separated_3x']:.1f}%, "
            f"sep10x={r['pct_separated_10x']:.1f}%"
        )
    lines += ["", "mc6 severe-flag sensitivity:"]
    for r in frows:
        lines.append(
            f"  {r['margin']}: all {r['pct_separated_all']:.1f}% (n={r['n_active_all']:,}), "
            f"no-severe {r['pct_separated_no_severe']:.1f}% (n={r['n_active_no_severe']:,})"
        )
    (DER / "FORMAT_FLAGS_ESTIMATED.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
