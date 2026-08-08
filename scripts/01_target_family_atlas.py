"""
ToxCast invitrodb v4.3 target-family bioactivity atlas with cytotoxicity filtering.

Claim: after removing cytotox-linked actives, most apparent multi-family
promiscuity collapses and only some intended target families retain
selective chemical hits.

Requires the public invitrodb v4.3 extracts. Set INVITRODB_SUMMARY_DIR to the
folder containing the winning-model summary CSV, assay annotations XLSX, and
cytotoxicity XLSX, or set INVITRODB_MC56 / INVITRODB_ANNOTATIONS /
INVITRODB_CYTOTOX to individual file paths.
"""
from __future__ import annotations

import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

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
import matplotlib.pyplot as plt  # noqa: E402

TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
DER.mkdir(parents=True, exist_ok=True)

HIT_THRESH = 0.9
USECOLS = ["chid", "casn", "chnm", "aeid", "hitc", "ac50"]

# Annotation families that are controls / non-target biology for primary ranking
CONTROL_FAMILIES = {
    "background measurement",
    "background control",
}

# Display cleanup for a few awkward EPA labels
FAMILY_RENAME = {
    "channel 1": "ion channel (panel 1)",
    "channel 2": "ion channel (panel 2)",
}


def load_annotations() -> pd.DataFrame:
    ann = pd.read_excel(ANN, usecols=[
        "aeid",
        "assay_component_endpoint_name",
        "intended_target_family",
        "intended_target_family_sub",
        "cell_viability_assay",
        "burst_assay",
        "organism",
    ])
    ann = ann.dropna(subset=["aeid"]).copy()
    ann["aeid"] = ann["aeid"].astype(np.int64)
    ann["intended_target_family"] = (
        ann["intended_target_family"].fillna("unannotated").astype(str).str.strip()
    )
    ann["family_display"] = ann["intended_target_family"].str.lower().map(
        lambda x: FAMILY_RENAME.get(x, x)
    )
    ann["is_control_family"] = ann["intended_target_family"].str.lower().isin(CONTROL_FAMILIES)
    ann["is_viability"] = ann["cell_viability_assay"].fillna(0).astype(int).eq(1)
    ann["is_burst"] = ann["burst_assay"].fillna(0).astype(int).eq(1)
    return ann


def load_cytotox() -> pd.DataFrame:
    cy = pd.read_excel(
        CYTO,
        usecols=[
            "chid",
            "casn",
            "chnm",
            "cytotox_median_um",
            "cytotox_lower_bound_um",
            "ntested",
            "nhit",
        ],
    )
    cy = cy.dropna(subset=["chid"]).copy()
    cy["chid"] = cy["chid"].astype(np.int64)
    return cy


def stream_aggregate(path: Path, ann: pd.DataFrame, cy: pd.DataFrame):
    fam_by_aeid = ann.set_index("aeid")["family_display"].to_dict()
    ctrl_by_aeid = ann.set_index("aeid")["is_control_family"].to_dict()
    viab_by_aeid = ann.set_index("aeid")["is_viability"].to_dict()

    cy_bound = cy.set_index("chid")["cytotox_lower_bound_um"].to_dict()
    chem_name = cy.set_index("chid")["chnm"].to_dict()
    chem_casn = cy.set_index("chid")["casn"].astype(str).to_dict()

    n_rows = 0
    n_active = 0
    n_selective = 0
    n_nonspecific = 0
    n_active_no_bound = 0

    # family-level: tested curves, actives, selective actives
    fam_tested = Counter()
    fam_active = Counter()
    fam_selective = Counter()
    fam_nonspecific = Counter()

    # chemical x family active / selective (for promiscuity)
    chem_fam_active: dict[int, set[str]] = defaultdict(set)
    chem_fam_selective: dict[int, set[str]] = defaultdict(set)
    chem_tested = Counter()
    chem_active_n = Counter()
    chem_selective_n = Counter()

    # also track non-control biological families only
    chem_biofam_active: dict[int, set[str]] = defaultdict(set)
    chem_biofam_selective: dict[int, set[str]] = defaultdict(set)

    unknown_aeid = 0

    for ci, chunk in enumerate(
        pd.read_csv(path, usecols=USECOLS, chunksize=500_000, low_memory=False)
    ):
        chunk = chunk.dropna(subset=["chid", "aeid", "hitc"]).copy()
        if chunk.empty:
            continue
        chunk["chid"] = chunk["chid"].astype(np.int64)
        chunk["aeid"] = chunk["aeid"].astype(np.int64)
        n_rows += len(chunk)

        families = chunk["aeid"].map(fam_by_aeid)
        unknown_aeid += int(families.isna().sum())
        families = families.fillna("unannotated")
        is_ctrl = chunk["aeid"].map(ctrl_by_aeid).fillna(False).astype(bool)
        # viability endpoints still counted in family totals but flagged elsewhere

        # tested
        for fam, cnt in families.value_counts().items():
            fam_tested[str(fam)] += int(cnt)
        vc = chunk["chid"].value_counts()
        for chid, cnt in vc.items():
            chem_tested[int(chid)] += int(cnt)
        for chid, casn, chnm in chunk[["chid", "casn", "chnm"]].drop_duplicates("chid").itertuples(
            index=False
        ):
            cid = int(chid)
            if cid not in chem_name:
                chem_name[cid] = "" if pd.isna(chnm) else str(chnm)
                chem_casn[cid] = "" if pd.isna(casn) else str(casn)

        is_hit = chunk["hitc"].to_numpy(dtype=float) >= HIT_THRESH
        n_active += int(is_hit.sum())
        if not is_hit.any():
            if ci % 4 == 0:
                print(f"  streamed {n_rows:,} rows ...", flush=True)
            continue

        hits = chunk.loc[is_hit].copy()
        hits["family"] = families.loc[is_hit].astype(str).to_numpy()
        hits["is_ctrl"] = is_ctrl.loc[is_hit].to_numpy()
        ac = hits["ac50"].to_numpy(dtype=float)
        chids = hits["chid"].to_numpy(dtype=np.int64)
        bounds = np.array([cy_bound.get(int(c), np.nan) for c in chids], dtype=float)
        selective = np.isfinite(ac) & (ac > 0) & np.isfinite(bounds) & (ac < bounds)
        nonspecific = ~selective
        n_selective += int(selective.sum())
        n_nonspecific += int(nonspecific.sum())
        n_active_no_bound += int((~np.isfinite(bounds)).sum())

        for fam, cnt in pd.Series(hits["family"]).value_counts().items():
            fam_active[str(fam)] += int(cnt)
        for fam, cnt in pd.Series(hits.loc[selective, "family"]).value_counts().items():
            fam_selective[str(fam)] += int(cnt)
        for fam, cnt in pd.Series(hits.loc[nonspecific, "family"]).value_counts().items():
            fam_nonspecific[str(fam)] += int(cnt)

        # chemical aggregates (vectorized unique pairs where possible)
        fams = hits["family"].to_numpy()
        ctrls = hits["is_ctrl"].to_numpy()
        for chid, fam, sel, ctrl in zip(chids, fams, selective, ctrls):
            cid = int(chid)
            fams_s = str(fam)
            chem_active_n[cid] += 1
            chem_fam_active[cid].add(fams_s)
            if not ctrl:
                chem_biofam_active[cid].add(fams_s)
            if sel:
                chem_selective_n[cid] += 1
                chem_fam_selective[cid].add(fams_s)
                if not ctrl:
                    chem_biofam_selective[cid].add(fams_s)

        if ci % 4 == 0:
            print(f"  streamed {n_rows:,} rows ...", flush=True)

    return {
        "n_rows": n_rows,
        "n_active": n_active,
        "n_selective": n_selective,
        "n_nonspecific": n_nonspecific,
        "n_active_no_bound": n_active_no_bound,
        "unknown_aeid": unknown_aeid,
        "fam_tested": fam_tested,
        "fam_active": fam_active,
        "fam_selective": fam_selective,
        "fam_nonspecific": fam_nonspecific,
        "chem_tested": chem_tested,
        "chem_active_n": chem_active_n,
        "chem_selective_n": chem_selective_n,
        "chem_fam_active": chem_fam_active,
        "chem_fam_selective": chem_fam_selective,
        "chem_biofam_active": chem_biofam_active,
        "chem_biofam_selective": chem_biofam_selective,
        "chem_name": chem_name,
        "chem_casn": chem_casn,
        "ann": ann,
    }


def family_table(agg) -> pd.DataFrame:
    ann = agg["ann"]
    ctrl_map = (
        ann.groupby(ann["family_display"].str.lower())["is_control_family"].max().to_dict()
    )
    families = sorted(set(agg["fam_tested"]) | set(agg["fam_active"]))
    rows = []
    for fam in families:
        n_test = agg["fam_tested"][fam]
        n_act = agg["fam_active"][fam]
        n_sel = agg["fam_selective"][fam]
        n_ns = agg["fam_nonspecific"][fam]
        rows.append(
            {
                "family": fam,
                "is_control_family": bool(ctrl_map.get(fam.lower(), False)),
                "n_tested_curves": n_test,
                "n_active": n_act,
                "n_selective": n_sel,
                "n_nonspecific": n_ns,
                "pct_active_of_tested": 100 * n_act / n_test if n_test else np.nan,
                "pct_selective_of_actives": 100 * n_sel / n_act if n_act else np.nan,
                "pct_nonspecific_of_actives": 100 * n_ns / n_act if n_act else np.nan,
            }
        )
    df = pd.DataFrame(rows)
    # ranks among biological (non-control) families
    bio = df.loc[~df["is_control_family"]].copy()
    bio["rank_active"] = bio["n_active"].rank(ascending=False, method="average")
    bio["rank_selective"] = bio["n_selective"].rank(ascending=False, method="average")
    df = df.merge(
        bio[["family", "rank_active", "rank_selective"]], on="family", how="left"
    )
    return df.sort_values(["is_control_family", "n_active"], ascending=[True, False])


def chemical_table(agg) -> pd.DataFrame:
    rows = []
    for chid, n_test in agg["chem_tested"].items():
        n_act = agg["chem_active_n"][chid]
        n_sel = agg["chem_selective_n"][chid]
        fam_a = agg["chem_biofam_active"][chid]
        fam_s = agg["chem_biofam_selective"][chid]
        rows.append(
            {
                "chid": chid,
                "casn": agg["chem_casn"].get(chid, ""),
                "chnm": agg["chem_name"].get(chid, ""),
                "n_tested": n_test,
                "n_active": n_act,
                "n_selective": n_sel,
                "n_biofamilies_active": len(fam_a),
                "n_biofamilies_selective": len(fam_s),
                "family_loss": len(fam_a) - len(fam_s),
            }
        )
    return pd.DataFrame(rows)


def family_stability(fam: pd.DataFrame) -> dict:
    bio = fam.loc[~fam["is_control_family"]].copy()
    mask = bio["rank_active"].notna() & bio["rank_selective"].notna()
    rho, p = stats.spearmanr(bio.loc[mask, "rank_active"], bio.loc[mask, "rank_selective"])
    # Jaccard top 10 / top 20 by counts
    def jaccard(n):
        a = set(bio.nlargest(n, "n_active")["family"])
        b = set(bio.nlargest(n, "n_selective")["family"])
        return len(a & b) / len(a | b) if (a or b) else np.nan

    return {
        "n_biological_families": int(len(bio)),
        "spearman_family_ranks": float(rho),
        "spearman_pvalue": float(p),
        "jaccard_top10": float(jaccard(10)),
        "jaccard_top20": float(jaccard(20)),
    }


def make_figures(fam: pd.DataFrame, chem: pd.DataFrame, totals: dict, stab: dict):
    apply_style()
    bio = fam.loc[~fam["is_control_family"]].copy()
    bio = bio.sort_values("n_active", ascending=False)

    # -------- Figure 1: three-panel summary --------
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.0))

    ax = axes[0]
    style_panel(ax, grid_axis="y")
    vals = [totals["n_active"], totals["n_selective"], totals["n_nonspecific"]]
    labels = ["All actives", "Selective\n(AC50 < cytotox\nlower bound)", "Cytotox-linked\nor unbound"]
    colors = [PALETTE["navy"], PALETTE["green"], PALETTE["crimson"]]
    ax.bar(range(3), vals, color=colors, width=0.72, zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("Active curves (hitc ≥ 0.9)")
    ymax = max(vals) * 1.18
    ax.set_ylim(0, ymax)
    for i, v in enumerate(vals):
        ax.text(i, v + ymax * 0.02, f"{v/1e3:.0f}k\n({100*v/max(totals['n_active'],1):.0f}%)",
                ha="center", va="bottom", fontsize=8)
    panel_label(ax, "A", "Library-wide retention")

    ax = axes[1]
    style_panel(ax, grid_axis="y")
    # distribution of selective fraction across biological families with >=50 actives
    sub = bio.loc[bio["n_active"] >= 50, "pct_selective_of_actives"].dropna()
    ax.hist(sub, bins=15, color=PALETTE["navy"], edgecolor="white", zorder=3)
    ax.axvline(sub.median(), color=PALETTE["crimson"], lw=1.5, ls="--",
               label=f"Median {sub.median():.0f}%")
    ax.set_xlabel("% of family actives that are selective")
    ax.set_ylabel(f"Biological families (n={len(sub)}, ≥50 actives)")
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "B", "Family selectivity")

    ax = axes[2]
    style_panel(ax)
    x = chem["n_biofamilies_active"].to_numpy()
    y = chem["n_biofamilies_selective"].to_numpy()
    # density-friendly scatter subsample
    if len(chem) > 8000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(chem), size=8000, replace=False)
        xs, ys = x[idx], y[idx]
    else:
        xs, ys = x, y
    ax.scatter(xs, ys, s=6, alpha=0.22, c=PALETTE["navy"], linewidths=0, zorder=3)
    lim = max(float(np.nanmax(x)), float(np.nanmax(y))) * 1.02
    ax.plot([0, lim], [0, lim], color="#666666", lw=0.9, ls="--", zorder=4)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Biological families with ≥1 active")
    ax.set_ylabel("Biological families with ≥1 selective active")
    # annotate median collapse
    med_loss = float(np.nanmedian(chem.loc[chem["n_biofamilies_active"] > 0, "family_loss"]))
    ax.text(0.98, 0.05, f"Median families lost: {med_loss:.0f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5)
    panel_label(ax, "C", "Chemical promiscuity")

    fig.tight_layout()
    save_fig(fig, "F1_main", FIG)

    # -------- Figure 2: family atlas (top 25 by actives) --------
    apply_style()
    top = bio.nlargest(25, "n_active").iloc[::-1]  # bottom-to-top for barh
    fig, ax = plt.subplots(figsize=(9.5, 8.2))
    style_panel(ax, grid_axis="x")
    y = np.arange(len(top))
    ax.barh(y, top["n_nonspecific"], color=PALETTE["crimson"], height=0.75,
            label="Cytotox-linked / unbound", zorder=3)
    ax.barh(y, top["n_selective"], left=top["n_nonspecific"], color=PALETTE["green"],
            height=0.75, label="Selective", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(top["family"], fontsize=9)
    ax.set_xlabel("Active curves (hitc ≥ 0.9)")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    save_fig(fig, "F2_family_atlas", FIG)

    # -------- Figure 3: promiscuity distributions --------
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    active_chems = chem.loc[chem["n_biofamilies_active"] > 0]
    ax = axes[0]
    style_panel(ax, grid_axis="y")
    bins = np.arange(0, min(40, int(active_chems["n_biofamilies_active"].max()) + 2))
    ax.hist(active_chems["n_biofamilies_active"], bins=bins, color=PALETTE["navy"],
            alpha=0.85, edgecolor="white", label="Any active", zorder=3)
    ax.hist(active_chems["n_biofamilies_selective"], bins=bins, color=PALETTE["green"],
            alpha=0.65, edgecolor="white", label="Selective only", zorder=4)
    ax.set_xlabel("Number of biological target families")
    ax.set_ylabel("Chemicals")
    ax.legend(frameon=False, fontsize=8.5)
    panel_label(ax, "A", "Promiscuity before vs after filter")

    ax = axes[1]
    style_panel(ax, grid_axis="y")
    # top 15 chemicals by family loss among those with many families active
    movers = active_chems.nlargest(15, "family_loss").iloc[::-1]
    labels = [((n[:28] + "…") if len(str(n)) > 29 else str(n)) for n in movers["chnm"]]
    y = np.arange(len(movers))
    ax.barh(y, movers["n_biofamilies_active"], color=PALETTE["navy"], height=0.7,
            label="Before", zorder=3)
    ax.barh(y, movers["n_biofamilies_selective"], color=PALETTE["green"], height=0.45,
            label="After selective filter", zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Biological families with ≥1 hit")
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "B", "Largest family collapses")

    fig.tight_layout()
    save_fig(fig, "F3_promiscuity", FIG)


def main():
    print("Loading annotations and cytotox ...", flush=True)
    if not MC56.exists():
        raise SystemExit(f"Missing {MC56}")
    ann = load_annotations()
    cy = load_cytotox()
    print(f"  endpoints={len(ann):,}  cytotox chemicals={len(cy):,}", flush=True)
    print("Streaming", MC56, flush=True)
    agg = stream_aggregate(MC56, ann, cy)
    print(
        f"Done: rows={agg['n_rows']:,} active={agg['n_active']:,} "
        f"selective={agg['n_selective']:,} nonspecific={agg['n_nonspecific']:,}",
        flush=True,
    )

    fam = family_table(agg)
    fam.to_csv(TAB / "T1_family_summary.csv", index=False)

    chem = chemical_table(agg)
    chem.to_csv(TAB / "T0_chemical_promiscuity.csv", index=False)

    stab = family_stability(fam)
    pd.DataFrame([stab]).to_csv(TAB / "T2_family_rank_stability.csv", index=False)

    # chemical summary stats
    active_chems = chem.loc[chem["n_biofamilies_active"] > 0]
    chem_summary = {
        "n_chemicals_tested": int(len(chem)),
        "n_chemicals_with_active": int((chem["n_active"] > 0).sum()),
        "n_chemicals_with_selective": int((chem["n_selective"] > 0).sum()),
        "median_biofamilies_active": float(active_chems["n_biofamilies_active"].median()),
        "median_biofamilies_selective": float(active_chems["n_biofamilies_selective"].median()),
        "median_family_loss": float(active_chems["family_loss"].median()),
        "mean_biofamilies_active": float(active_chems["n_biofamilies_active"].mean()),
        "mean_biofamilies_selective": float(active_chems["n_biofamilies_selective"].mean()),
        "pct_chemicals_losing_ge50pct_families": float(
            100
            * (
                (active_chems["n_biofamilies_selective"] / active_chems["n_biofamilies_active"])
                <= 0.5
            ).mean()
        ),
    }
    pd.DataFrame([chem_summary]).to_csv(TAB / "T3_chemical_promiscuity_summary.csv", index=False)

    movers = active_chems.nlargest(50, "family_loss")[
        [
            "chid",
            "casn",
            "chnm",
            "n_active",
            "n_selective",
            "n_biofamilies_active",
            "n_biofamilies_selective",
            "family_loss",
        ]
    ]
    movers.to_csv(TAB / "T4_top_family_collapses.csv", index=False)

    totals = {
        "n_rows": agg["n_rows"],
        "n_active": agg["n_active"],
        "n_selective": agg["n_selective"],
        "n_nonspecific": agg["n_nonspecific"],
        "n_active_no_bound": agg["n_active_no_bound"],
        "pct_selective": 100 * agg["n_selective"] / max(agg["n_active"], 1),
        "pct_nonspecific": 100 * agg["n_nonspecific"] / max(agg["n_active"], 1),
    }
    pd.DataFrame([totals]).to_csv(TAB / "T0b_library_totals.csv", index=False)

    # families with high activity but low selectivity
    bio = fam.loc[~fam["is_control_family"]].copy()
    interesting = bio.loc[bio["n_active"] >= 100].sort_values("pct_selective_of_actives")
    interesting.head(20).to_csv(TAB / "T5_least_selective_families.csv", index=False)
    bio.loc[bio["n_active"] >= 100].sort_values(
        "pct_selective_of_actives", ascending=False
    ).head(20).to_csv(TAB / "T6_most_selective_families.csv", index=False)

    print("Making figures ...", flush=True)
    make_figures(fam, chem, totals, stab)

    # RESULTS text
    n_bio = int((~fam["is_control_family"]).sum())
    n_bio_ge50 = int(((~fam["is_control_family"]) & (fam["n_active"] >= 50)).sum())
    med_sel = float(
        bio.loc[bio["n_active"] >= 50, "pct_selective_of_actives"].median()
    )
    n_fam_sel_ge50pct = int(
        ((~fam["is_control_family"]) & (fam["n_active"] >= 50) & (fam["pct_selective_of_actives"] >= 50)).sum()
    )
    lines = [
        "ToxCast target-family atlas with cytotoxicity filtering (invitrodb v4.3)",
        f"Curve rows: {agg['n_rows']:,}",
        f"Active curves (hitc>=0.9): {agg['n_active']:,}",
        f"Selective (AC50 < cytotox lower bound): {agg['n_selective']:,} ({totals['pct_selective']:.1f}%)",
        f"Cytotox-linked or unbound: {agg['n_nonspecific']:,} ({totals['pct_nonspecific']:.1f}%)",
        f"Biological families: {n_bio}",
        f"Biological families with >=50 actives: {n_bio_ge50}",
        f"Median % selective among those families: {med_sel:.1f}%",
        f"Families (>=50 actives) with >=50% selective: {n_fam_sel_ge50pct}",
        "",
        "Chemical promiscuity (biological families):",
        f"  chemicals with >=1 active: {chem_summary['n_chemicals_with_active']:,}",
        f"  median families active -> selective: "
        f"{chem_summary['median_biofamilies_active']:.0f} -> "
        f"{chem_summary['median_biofamilies_selective']:.0f}",
        f"  median families lost: {chem_summary['median_family_loss']:.0f}",
        f"  % chemicals losing >=50% of active families: "
        f"{chem_summary['pct_chemicals_losing_ge50pct_families']:.1f}%",
        "",
        "Family rank stability (biological families):",
        f"  Spearman active vs selective ranks: {stab['spearman_family_ranks']:.3f}",
        f"  Jaccard top10/top20: {stab['jaccard_top10']:.3f}/{stab['jaccard_top20']:.3f}",
        "",
        "Least selective families (>=100 actives):",
    ]
    for _, r in interesting.head(8).iterrows():
        lines.append(
            f"  {r['family']}: {int(r['n_active']):,} actives, "
            f"{r['pct_selective_of_actives']:.1f}% selective"
        )
    lines.append("")
    lines.append("Most selective families (>=100 actives):")
    for _, r in bio.loc[bio["n_active"] >= 100].sort_values(
        "pct_selective_of_actives", ascending=False
    ).head(8).iterrows():
        lines.append(
            f"  {r['family']}: {int(r['n_active']):,} actives, "
            f"{r['pct_selective_of_actives']:.1f}% selective"
        )

    (DER / "RESULTS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print("\nWrote tables/, figures/, derived/RESULTS.txt")


if __name__ == "__main__":
    main()
