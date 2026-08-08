"""
Sensitivity and stratification analyses for the target-family atlas:
1) Selectivity-margin sensitivity (1x, 3x, 10x windows)
2) Cell-free / biochemical vs cell-based stratification
3) Severe mc6-flag secondary filter
4) Family chemical-set composition shift (QSAR-relevant set change)

Requires the same invitrodb v4.3 inputs as 01_target_family_atlas.py.
Set INVITRODB_SUMMARY_DIR or INVITRODB_MC56 / INVITRODB_ANNOTATIONS /
INVITRODB_CYTOTOX.
"""
from __future__ import annotations

import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "tables"
FIG = ROOT / "figures"
DER = ROOT / "derived"


def resolve_summary_dir() -> Path:
    env_dir = os.environ.get("INVITRODB_SUMMARY_DIR")
    if env_dir:
        return Path(env_dir)
    return ROOT.parent / "ToxCast" / "derived" / "summary_extract"


def resolve_input(env_key: str, filename: str) -> Path:
    env_file = os.environ.get(env_key)
    if env_file:
        return Path(env_file)
    return resolve_summary_dir() / filename


MC56 = resolve_input("INVITRODB_MC56", "mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv")
ANN = resolve_input("INVITRODB_ANNOTATIONS", "assay_annotations_invitrodb_v4_3_AUG2024.xlsx")
CYTO = resolve_input("INVITRODB_CYTOTOX", "cytotox_invitrodb_v4_3_AUG2024.xlsx")
sys.path.insert(0, str(ROOT))
from figstyle import PALETTE, apply_style, panel_label, save_fig, style_panel  # noqa: E402

HIT_THRESH = 0.9
SEVERE = {5, 6, 7, 10, 11, 15, 17, 18, 19}
USECOLS = ["chid", "casn", "chnm", "aeid", "hitc", "ac50", "mc6_flags"]
CONTROL_FAMILIES = {"background measurement", "background control"}
FAMILY_RENAME = {
    "channel 1": "ion channel (panel 1)",
    "channel 2": "ion channel (panel 2)",
}
MARGINS = {
    "1x": 1.0,
    "3x": 3.0,
    "10x": 10.0,
}
FOCUS_FAMS = [
    "cell cycle",
    "nuclear receptor",
    "cyp",
    "neurodevelopment",
    "steroid hormone",
    "apoptosis",
    "cytokine",
    "gpcr",
    "kinase",
    "protease",
]


def parse_severe_flags(series: pd.Series) -> np.ndarray:
    s = series.fillna("").astype(str)
    s = s.str.replace("|", ",", regex=False).str.replace(";", ",", regex=False)
    severe = np.zeros(len(s), dtype=bool)
    for fid in SEVERE:
        pat = rf"(?:^|,)\s*{fid}\s*(?:,|$)"
        severe |= s.str.contains(pat, regex=True, na=False).to_numpy()
    return severe


def classify_format(row) -> str:
    """Map annotation fields to cell-free/biochemical vs cell-based vs other."""
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


def load_ann() -> pd.DataFrame:
    ann = pd.read_excel(
        ANN,
        usecols=[
            "aeid",
            "intended_target_family",
            "assay_format_type",
            "cell_format",
            "assay_source_name",
        ],
    )
    ann = ann.dropna(subset=["aeid"]).copy()
    ann["aeid"] = ann["aeid"].astype(np.int64)
    ann["intended_target_family"] = (
        ann["intended_target_family"].fillna("unannotated").astype(str).str.strip()
    )
    ann["family"] = ann["intended_target_family"].str.lower().map(
        lambda x: FAMILY_RENAME.get(x, x)
    )
    ann["is_control"] = ann["intended_target_family"].str.lower().isin(CONTROL_FAMILIES)
    ann["format_class"] = ann.apply(classify_format, axis=1)
    return ann


def main():
    print("Loading annotations/cytotox ...", flush=True)
    ann = load_ann()
    cy = pd.read_excel(CYTO, usecols=["chid", "cytotox_lower_bound_um"])
    cy = cy.dropna(subset=["chid"]).copy()
    cy["chid"] = cy["chid"].astype(np.int64)
    bound = cy.set_index("chid")["cytotox_lower_bound_um"].to_dict()
    fam = ann.set_index("aeid")["family"].to_dict()
    ctrl = ann.set_index("aeid")["is_control"].to_dict()
    fmt = ann.set_index("aeid")["format_class"].to_dict()

    # endpoint inventory
    ep_fmt = ann["format_class"].value_counts().to_dict()
    print("Endpoint format inventory:", ep_fmt)

    # Counters
    n_active = 0
    n_active_no_severe = 0
    margin_counts = {m: 0 for m in MARGINS}
    margin_counts_nsev = {m: 0 for m in MARGINS}
    # format x selectivity (1x)
    fmt_active = Counter()
    fmt_sel = {m: Counter() for m in MARGINS}
    # family x margin selective counts (biological families)
    fam_active = Counter()
    fam_sel = {m: Counter() for m in MARGINS}
    # chemical sets for focus families
    chem_active_sets = {f: set() for f in FOCUS_FAMS}
    chem_sel_sets = {m: {f: set() for f in FOCUS_FAMS} for m in MARGINS}
    # severe-flag subset family selective (1x)
    fam_active_nsev = Counter()
    fam_sel_nsev = Counter()
    # promiscuity under margins
    chem_biofam_active = defaultdict(set)
    chem_biofam_sel = {m: defaultdict(set) for m in MARGINS}

    print("Streaming", MC56, flush=True)
    for ci, chunk in enumerate(
        pd.read_csv(MC56, usecols=USECOLS, chunksize=500_000, low_memory=False)
    ):
        chunk = chunk.dropna(subset=["chid", "aeid", "hitc"]).copy()
        if chunk.empty:
            continue
        chunk["chid"] = chunk["chid"].astype(np.int64)
        chunk["aeid"] = chunk["aeid"].astype(np.int64)
        is_hit = chunk["hitc"].to_numpy(dtype=float) >= HIT_THRESH
        if not is_hit.any():
            continue
        hits = chunk.loc[is_hit].copy()
        n_active += len(hits)
        severe = parse_severe_flags(hits["mc6_flags"])
        no_sev = ~severe
        n_active_no_severe += int(no_sev.sum())

        ac = hits["ac50"].to_numpy(dtype=float)
        chids = hits["chid"].to_numpy(dtype=np.int64)
        bounds = np.array([bound.get(int(c), np.nan) for c in chids], dtype=float)
        families = hits["aeid"].map(fam).fillna("unannotated").astype(str).to_numpy()
        formats = hits["aeid"].map(fmt).fillna("other_unclassified").astype(str).to_numpy()
        is_ctrl = hits["aeid"].map(ctrl).fillna(False).to_numpy(dtype=bool)

        ok = np.isfinite(ac) & (ac > 0) & np.isfinite(bounds) & (bounds > 0)
        for mname, factor in MARGINS.items():
            thr = bounds / factor
            sel = ok & (ac < thr)
            margin_counts[mname] += int(sel.sum())
            margin_counts_nsev[mname] += int((sel & no_sev).sum())
            for fcls, cnt in pd.Series(formats[sel]).value_counts().items():
                fmt_sel[mname][str(fcls)] += int(cnt)
            for fam_s, cnt in pd.Series(families[sel]).value_counts().items():
                fam_sel[mname][str(fam_s)] += int(cnt)
            for chid, fam_s, ctrl_i, is_s in zip(chids, families, is_ctrl, sel):
                if is_s and (not ctrl_i):
                    chem_biofam_sel[mname][int(chid)].add(str(fam_s))
            for chid, fam_s, is_s in zip(chids, families, sel):
                if is_s and fam_s in chem_sel_sets[mname]:
                    chem_sel_sets[mname][fam_s].add(int(chid))

        for fcls, cnt in pd.Series(formats).value_counts().items():
            fmt_active[str(fcls)] += int(cnt)
        for fam_s, cnt in pd.Series(families).value_counts().items():
            fam_active[str(fam_s)] += int(cnt)
        for chid, fam_s, ctrl_i in zip(chids, families, is_ctrl):
            if not ctrl_i:
                chem_biofam_active[int(chid)].add(str(fam_s))
        for chid, fam_s in zip(chids, families):
            if fam_s in chem_active_sets:
                chem_active_sets[fam_s].add(int(chid))

        # no-severe subset, 1x margin
        if no_sev.any():
            fams_ns = families[no_sev]
            ac_ns = ac[no_sev]
            b_ns = bounds[no_sev]
            ok_ns = np.isfinite(ac_ns) & (ac_ns > 0) & np.isfinite(b_ns) & (b_ns > 0)
            sel_ns = ok_ns & (ac_ns < b_ns)
            for fam_s, cnt in pd.Series(fams_ns).value_counts().items():
                fam_active_nsev[str(fam_s)] += int(cnt)
            for fam_s, cnt in pd.Series(fams_ns[sel_ns]).value_counts().items():
                fam_sel_nsev[str(fam_s)] += int(cnt)

        if ci % 4 == 0:
            print(f"  streamed actives so far {n_active:,} ...", flush=True)

    # ---- tables ----
    sens_rows = []
    for mname, factor in MARGINS.items():
        nsel = margin_counts[mname]
        sens_rows.append(
            {
                "rule": f"AC50 < cytotox_lower_bound / {factor:g}",
                "margin_label": mname,
                "n_selective": nsel,
                "pct_of_actives": 100 * nsel / max(n_active, 1),
                "n_selective_no_severe": margin_counts_nsev[mname],
                "pct_of_no_severe_actives": 100
                * margin_counts_nsev[mname]
                / max(n_active_no_severe, 1),
            }
        )
    sens = pd.DataFrame(sens_rows)
    sens.to_csv(TAB / "T7_selectivity_margin_sensitivity.csv", index=False)

    fmt_rows = []
    for fcls in sorted(set(fmt_active) | set(fmt_sel["1x"])):
        na = fmt_active[fcls]
        row = {
            "format_class": fcls,
            "n_endpoints": int((ann["format_class"] == fcls).sum()),
            "n_active": na,
            "pct_of_library_actives": 100 * na / max(n_active, 1),
        }
        for mname in MARGINS:
            ns = fmt_sel[mname][fcls]
            row[f"n_selective_{mname}"] = ns
            row[f"pct_selective_of_format_actives_{mname}"] = 100 * ns / max(na, 1)
        fmt_rows.append(row)
    fmt_tab = pd.DataFrame(fmt_rows).sort_values("n_active", ascending=False)
    fmt_tab.to_csv(TAB / "T8_format_stratification.csv", index=False)

    # family retention across margins for focus + all bio families with >=100 actives
    fam_rows = []
    for fam_s, na in fam_active.items():
        if fam_s in CONTROL_FAMILIES:
            continue
        row = {"family": fam_s, "n_active": na}
        for mname in MARGINS:
            ns = fam_sel[mname][fam_s]
            row[f"n_selective_{mname}"] = ns
            row[f"pct_selective_{mname}"] = 100 * ns / max(na, 1)
        # no severe
        na2 = fam_active_nsev[fam_s]
        ns2 = fam_sel_nsev[fam_s]
        row["n_active_no_severe"] = na2
        row["n_selective_1x_no_severe"] = ns2
        row["pct_selective_1x_no_severe"] = 100 * ns2 / max(na2, 1)
        fam_rows.append(row)
    fam_tab = pd.DataFrame(fam_rows).sort_values("n_active", ascending=False)
    fam_tab.to_csv(TAB / "T9_family_margin_and_flag_sensitivity.csv", index=False)

    # chemical set Jaccard for focus families
    j_rows = []
    for fam_s in FOCUS_FAMS:
        a = chem_active_sets[fam_s]
        for mname in MARGINS:
            b = chem_sel_sets[mname][fam_s]
            inter = len(a & b)
            union = len(a | b)
            j_rows.append(
                {
                    "family": fam_s,
                    "margin": mname,
                    "n_chemicals_active": len(a),
                    "n_chemicals_selective": len(b),
                    "jaccard_active_vs_selective_chemset": inter / union if union else np.nan,
                    "pct_active_chems_retained": 100 * len(b) / max(len(a), 1),
                }
            )
    jtab = pd.DataFrame(j_rows)
    jtab.to_csv(TAB / "T10_family_chemical_set_shift.csv", index=False)

    # promiscuity under margins
    prom_rows = []
    active_chems = [c for c, s in chem_biofam_active.items() if s]
    for mname in MARGINS:
        losses = []
        ge50 = 0
        for c in active_chems:
            na = len(chem_biofam_active[c])
            ns = len(chem_biofam_sel[mname][c])
            losses.append(na - ns)
            if na > 0 and ns / na <= 0.5:
                ge50 += 1
        prom_rows.append(
            {
                "margin": mname,
                "n_chemicals_with_active": len(active_chems),
                "median_biofamilies_selective": float(
                    np.median([len(chem_biofam_sel[mname][c]) for c in active_chems])
                ),
                "median_family_loss": float(np.median(losses)),
                "pct_chemicals_losing_ge50pct_families": 100 * ge50 / max(len(active_chems), 1),
            }
        )
    pd.DataFrame(prom_rows).to_csv(TAB / "T11_promiscuity_by_margin.csv", index=False)

    # summary json/text
    summary = {
        "n_active": n_active,
        "n_active_no_severe": n_active_no_severe,
        "pct_actives_without_severe_flag": 100 * n_active_no_severe / max(n_active, 1),
        "margin_counts": margin_counts,
        "margin_pct_of_actives": {
            k: 100 * v / max(n_active, 1) for k, v in margin_counts.items()
        },
        "margin_counts_no_severe": margin_counts_nsev,
        "format_active": dict(fmt_active),
        "format_selective_1x": dict(fmt_sel["1x"]),
        "endpoint_format_inventory": ep_fmt,
    }
    (DER / "REVISION_RESULTS.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = [
        "Revision analyses",
        f"Active curves: {n_active:,}",
        f"Actives without severe mc6 flags: {n_active_no_severe:,} "
        f"({summary['pct_actives_without_severe_flag']:.1f}%)",
        "",
        "Selectivity margins (% of all actives):",
    ]
    for mname in MARGINS:
        lines.append(
            f"  {mname}: {margin_counts[mname]:,} "
            f"({summary['margin_pct_of_actives'][mname]:.1f}%)"
        )
    lines.append("")
    lines.append("Format stratification (actives / 1x selective / % selective):")
    for _, r in fmt_tab.iterrows():
        lines.append(
            f"  {r['format_class']}: {int(r['n_active']):,} / "
            f"{int(r['n_selective_1x']):,} / {r['pct_selective_of_format_actives_1x']:.1f}%"
        )
    lines.append("")
    lines.append("Focus family % selective at 1x / 3x / 10x:")
    for fam_s in FOCUS_FAMS:
        sub = fam_tab.loc[fam_tab["family"] == fam_s]
        if sub.empty:
            continue
        r = sub.iloc[0]
        lines.append(
            f"  {fam_s}: {r['pct_selective_1x']:.1f}% / "
            f"{r['pct_selective_3x']:.1f}% / {r['pct_selective_10x']:.1f}%"
        )
    lines.append("")
    lines.append("No-severe subset, library % selective at 1x:")
    lines.append(
        f"  {margin_counts_nsev['1x']:,} / {n_active_no_severe:,} = "
        f"{100 * margin_counts_nsev['1x'] / max(n_active_no_severe, 1):.1f}%"
    )
    (DER / "REVISION_RESULTS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    # ---- Figure S1 / F4 ----
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.2))

    ax = axes[0]
    style_panel(ax, grid_axis="y")
    labels = ["1x\n(AC50 < bound)", "3x\n(AC50 < bound/3)", "10x\n(AC50 < bound/10)"]
    vals = [margin_counts[m] / n_active * 100 for m in ("1x", "3x", "10x")]
    vals_ns = [
        margin_counts_nsev[m] / max(n_active_no_severe, 1) * 100 for m in ("1x", "3x", "10x")
    ]
    x = np.arange(3)
    ax.bar(x - 0.18, vals, width=0.36, color=PALETTE["navy"], label="All actives", zorder=3)
    ax.bar(
        x + 0.18,
        vals_ns,
        width=0.36,
        color=PALETTE["green"],
        label="No severe mc6 flags",
        zorder=3,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("% of active curves that are selective")
    ax.set_ylim(0, max(vals + vals_ns) * 1.25)
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "A", "Selectivity margin")

    ax = axes[1]
    style_panel(ax, grid_axis="y")
    # Only formats with actives (no empty "Other" bar)
    order = [k for k in ("cell_based", "cell_free_biochemical") if fmt_active[k] > 0]
    nice = {
        "cell_based": "Cell-based",
        "cell_free_biochemical": "Cell-free /\nbiochemical",
    }
    pcts = [100 * fmt_sel["1x"][k] / max(fmt_active[k], 1) for k in order]
    colors = [PALETTE["navy"], PALETTE["green"]][: len(order)]
    ax.bar(range(len(order)), pcts, color=colors, zorder=3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([nice[k] for k in order], fontsize=8)
    ax.set_ylabel("% selective among format actives (1x)")
    for i, k in enumerate(order):
        n = int(fmt_active[k])
        if n >= 1000:
            nlab = f"{n/1e3:.0f}k"
        else:
            nlab = str(n)
        ax.text(
            i,
            pcts[i] + 1.5,
            f"{nlab}\nactives",
            ha="center",
            va="bottom",
            fontsize=7.5,
        )
    panel_label(ax, "B", "Assay format")

    ax = axes[2]
    style_panel(ax, grid_axis="y")
    show = [
        "neurodevelopment",
        "steroid hormone",
        "cyp",
        "nuclear receptor",
        "cytokine",
        "cell cycle",
        "apoptosis",
    ]
    x = np.arange(len(show))
    w = 0.25
    for i, mname in enumerate(("1x", "3x", "10x")):
        ys = []
        for fam_s in show:
            sub = fam_tab.loc[fam_tab["family"] == fam_s]
            ys.append(float(sub.iloc[0][f"pct_selective_{mname}"]) if len(sub) else 0)
        ax.bar(
            x + (i - 1) * w,
            ys,
            width=w,
            label=mname,
            color=[PALETTE["green"], PALETTE["navy"], PALETTE["crimson"]][i],
            zorder=3,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(show, rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel("% of family actives selective")
    ax.legend(frameon=False, fontsize=8, title="Margin")
    panel_label(ax, "C", "Family retention")

    fig.tight_layout()
    save_fig(fig, "F4_sensitivity", FIG)
    print("Wrote tables T7-T11 and figures/F4_sensitivity.*")


if __name__ == "__main__":
    main()
