"""
Figures for the revised manuscript.

F1  bound provenance, tested-concentration scale, testability of the comparison
F2  family-level retention within the estimated-bound stratum, burst assays marked
F3  coverage-normalised chemical promiscuity, estimated vs default bound
F4  sensitivity: potency window, curve-fit flags, assay format, label-set shift

All panels are built from the derived tables written by scripts 08-10, except
the tested-concentration distribution, which is sampled once from the
multi-concentration summary and cached.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TOXCAST = ROOT.parent / "ToxCast" / "derived" / "summary_extract"
MC56 = TOXCAST / "mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv"
TAB = ROOT / "tables"
FIG = ROOT / "figures"
DER = ROOT / "derived"
sys.path.insert(0, str(ROOT))
from figstyle import PALETTE, apply_style, panel_label, save_fig, style_panel  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

HIT_THRESH = 0.9
MARGINS = ["1x", "3x", "10x"]
MARGIN_LABEL = {"1x": "1×", "3x": "3×", "10x": "10×"}
CONC_CACHE = DER / "active_conc_max_sample.csv"


def load_conc_sample(n: int = 200_000) -> np.ndarray:
    if CONC_CACHE.exists():
        return pd.read_csv(CONC_CACHE)["conc_max"].to_numpy(dtype=float)
    print("Sampling tested concentrations from mc5-6 ...", flush=True)
    keep = []
    rng = np.random.default_rng(7)
    for chunk in pd.read_csv(
        MC56, usecols=["hitc", "conc_max"], chunksize=500_000, low_memory=False
    ):
        chunk = chunk.dropna(subset=["hitc", "conc_max"])
        act = chunk.loc[chunk["hitc"].astype(float) >= HIT_THRESH, "conc_max"]
        keep.append(act.to_numpy(dtype=float))
    vals = np.concatenate(keep)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if len(vals) > n:
        vals = rng.choice(vals, size=n, replace=False)
    pd.DataFrame({"conc_max": vals}).to_csv(CONC_CACHE, index=False)
    return vals


def fig1(t12: pd.DataFrame, t14: pd.DataFrame, conc: np.ndarray):
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))

    est = t12.loc[t12["bound_provenance"] == "estimated"].iloc[0]
    dfl = t12.loc[t12["bound_provenance"] == "default"].iloc[0]
    pooled_n_act = est["n_active"] + dfl["n_active"]

    # ---- A: retention by provenance ----
    ax = axes[0]
    style_panel(ax, grid_axis="y")
    x = np.arange(len(MARGINS))
    w = 0.34
    e = [est[f"pct_separated_{m}"] for m in MARGINS]
    d = [dfl[f"pct_separated_{m}"] for m in MARGINS]
    pooled = [
        100 * (est[f"n_separated_{m}"] + dfl[f"n_separated_{m}"]) / pooled_n_act
        for m in MARGINS
    ]
    ax.bar(x - w / 2, e, w, color=PALETTE["navy"], label="Estimated bound", zorder=3)
    ax.bar(x + w / 2, d, w, color=PALETTE["crimson"], label="Default 1000 µM bound", zorder=3)
    ax.plot(x, pooled, "D--", color=PALETTE["ink"], ms=6, lw=1.3,
            label="Pooled (as usually reported)", zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{MARGIN_LABEL[m]} window" for m in MARGINS])
    ax.set_xlabel("Required potency window below bound")
    ax.set_ylabel("% of actives cytotoxicity-separated")
    ax.set_ylim(0, 112)
    ax.legend(frameon=False, fontsize=8, loc="center left")
    for xi, v in zip(x - w / 2, e):
        ax.text(xi, v + 2.5, f"{v:.0f}", ha="center", fontsize=8)
    for xi, v in zip(x + w / 2, d):
        ax.text(xi, v + 2.5, f"{v:.0f}", ha="center", fontsize=8)
    panel_label(ax, "A", "Retention depends on provenance")

    # ---- B: tested concentration scale ----
    ax = axes[1]
    style_panel(ax, grid_axis="y")
    bins = np.logspace(np.log10(max(conc.min(), 1e-3)), np.log10(2000), 45)
    ax.hist(conc, bins=bins, color=PALETTE["navy"], edgecolor="white", zorder=3)
    ax.set_xscale("log")
    med = float(np.median(conc))
    ax.axvline(med, color=PALETTE["green"], lw=1.6, ls="-",
               label=f"Median top conc. {med:.0f} µM", zorder=5)
    ax.axvline(1000, color=PALETTE["crimson"], lw=1.8, ls="--",
               label="Default bound 1000 µM", zorder=5)
    ax.set_xlabel("Highest concentration tested (µM, log scale)")
    ax.set_ylabel("Active curves (sampled)")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    panel_label(ax, "B", "The default bound is off-scale")

    # ---- C: testability ----
    ax = axes[2]
    style_panel(ax, grid_axis="y")
    tot = float(
        t14.loc[t14["quantity"].str.startswith("actives with finite"), "count"].iloc[0]
    )

    def pct(sub):
        return 100 * float(t14.loc[t14["quantity"].str.contains(sub), "count"].iloc[0]) / tot

    labels = [
        "Bound above\nhighest tested\nconcentration",
        "3× window below\nlowest tested\nconcentration",
        "10× window below\nlowest tested\nconcentration",
    ]
    vals = [pct("exceeds max tested"), pct("3x window"), pct("10x window")]
    ax.bar(range(3), vals, color=PALETTE["crimson"], width=0.62, zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("% of active curves")
    ymax = max(vals) * 1.3
    ax.set_ylim(0, ymax)
    for i, v in enumerate(vals):
        ax.text(i, v + ymax * 0.025, f"{v:.1f}%", ha="center", fontsize=8.5)
    panel_label(ax, "C", "Comparisons that cannot be answered")

    fig.tight_layout()
    save_fig(fig, "F1_provenance", FIG)


def fig2(t13: pd.DataFrame):
    apply_style()
    bio = t13.loc[t13["family"].ne("unannotated")].copy()
    top = bio.nlargest(22, "n_active_estimated_bound").iloc[::-1]

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 7.6),
                             gridspec_kw={"width_ratios": [1.25, 1.0]})
    y = np.arange(len(top))

    ax = axes[0]
    style_panel(ax, grid_axis="x")
    sep = top["n_separated_1x_estimated_bound"].to_numpy(dtype=float)
    prox = top["n_active_estimated_bound"].to_numpy(dtype=float) - sep
    ax.barh(y, prox, color=PALETTE["crimson"], height=0.74,
            label="Cytotoxicity-proximal", zorder=3)
    ax.barh(y, sep, left=prox, color=PALETTE["green"], height=0.74,
            label="Cytotoxicity-separated", zorder=3)
    ax.set_yticks(y)
    labels = [
        f"{f} *" if b else f"{f}"
        for f, b in zip(top["family"], top["is_burst_related"])
    ]
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Active curves, estimated-bound chemicals")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    panel_label(ax, "A", "Filter impact on family tallies")

    ax = axes[1]
    style_panel(ax, grid_axis="x")
    pct = top["pct_separated_1x_estimated_bound"].to_numpy(dtype=float)
    colors = [
        PALETTE["crimson"] if b else PALETTE["navy"]
        for b in top["is_burst_related"]
    ]
    ax.barh(y, pct, color=colors, height=0.74, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    med = float(
        bio.loc[bio["n_active_estimated_bound"] >= 50,
                "pct_separated_1x_estimated_bound"].median()
    )
    ax.axvline(med, color=PALETTE["ink"], lw=1.3, ls="--",
               label=f"Median across families {med:.0f}%", zorder=5)
    ax.set_xlabel("% of actives cytotoxicity-separated (1× window)")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    for yi, v in zip(y, pct):
        ax.text(v + 1.2, yi, f"{v:.0f}", va="center", fontsize=8)
    ax.set_xlim(0, max(pct) * 1.2)
    panel_label(ax, "B", "Retention, burst-annotated families in red")

    fig.tight_layout()
    save_fig(fig, "F2_family_retention", FIG)


def fig3(t15: pd.DataFrame):
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))
    est = t15.loc[t15["bound_provenance"] == "estimated"]
    dfl = t15.loc[t15["bound_provenance"] == "default"]
    bins = np.linspace(0, 1, 21)

    for ax, sub, name in (
        (axes[0], est, "estimated bound"),
        (axes[1], dfl, "default 1000 µM bound"),
    ):
        style_panel(ax, grid_axis="y")
        ax.hist(sub["frac_tested_active"], bins=bins, color=PALETTE["navy"],
                edgecolor="white", label="Any active", zorder=3)
        ax.hist(sub["frac_tested_separated"], bins=bins, histtype="step",
                color=PALETTE["green"], lw=2.0,
                label="Cytotoxicity-separated", zorder=4)
        ax.set_xlabel("Fraction of tested families hit")
        ax.set_ylabel(f"Chemicals (n={len(sub):,})")
        ax.legend(frameon=False, fontsize=8.5)
        ma = float(sub["frac_tested_active"].median())
        ms = float(sub["frac_tested_separated"].median())
        ax.text(0.97, 0.62, f"median\n{ma:.2f} → {ms:.2f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=9)
    panel_label(axes[0], "A", "Estimated bound")
    panel_label(axes[1], "B", "Default bound: filter is inert")

    ax = axes[2]
    style_panel(ax)
    x = est["n_biofam_active"].to_numpy(dtype=float)
    ycount = est["n_biofam_separated"].to_numpy(dtype=float)
    ax.scatter(x, ycount, s=9, alpha=0.28, c=PALETTE["navy"], linewidths=0, zorder=3)
    lim = max(x.max(), ycount.max()) * 1.04
    ax.plot([0, lim], [0, lim], color="#666666", lw=0.9, ls="--", zorder=4)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Families with ≥1 active")
    ax.set_ylabel("Families with ≥1 separated active")
    frac = 100 * float(((est["n_biofam_separated"] / est["n_biofam_active"]) <= 0.5).mean())
    ax.text(0.97, 0.06, f"{frac:.0f}% lose ≥half\nof active families",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5)
    panel_label(ax, "C", "Raw family counts, estimated bound")

    fig.tight_layout()
    save_fig(fig, "F3_promiscuity", FIG)


def fig4(t19: pd.DataFrame, t18: pd.DataFrame, t17: pd.DataFrame):
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3))
    x = np.arange(len(MARGINS))

    # ---- A: margin x severe-flag ----
    ax = axes[0]
    style_panel(ax, grid_axis="y")
    t19i = t19.set_index("margin")
    allv = [t19i.loc[m, "pct_separated_all"] for m in MARGINS]
    nsv = [t19i.loc[m, "pct_separated_no_severe"] for m in MARGINS]
    w = 0.34
    ax.bar(x - w / 2, allv, w, color=PALETTE["navy"], label="All actives", zorder=3)
    ax.bar(x + w / 2, nsv, w, color=PALETTE["sky"],
           label="No severe curve-fit flag", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([MARGIN_LABEL[m] for m in MARGINS])
    ax.set_xlabel("Potency window below bound")
    ax.set_ylabel("% cytotoxicity-separated")
    ax.set_ylim(0, max(allv + nsv) * 1.32)
    ax.legend(frameon=False, fontsize=8.5)
    for xi, v in zip(x - w / 2, allv):
        ax.text(xi, v + 0.7, f"{v:.1f}", ha="center", fontsize=8)
    for xi, v in zip(x + w / 2, nsv):
        ax.text(xi, v + 0.7, f"{v:.1f}", ha="center", fontsize=8)
    panel_label(ax, "A", "Window and curve-fit flags")

    # ---- B: assay format ----
    ax = axes[1]
    style_panel(ax, grid_axis="y")
    t18i = t18.set_index("format_class")
    cb = [t18i.loc["cell_based", f"pct_separated_{m}"] for m in MARGINS]
    cf = [t18i.loc["cell_free_biochemical", f"pct_separated_{m}"] for m in MARGINS]
    ax.bar(x - w / 2, cb, w, color=PALETTE["navy"], label="Cell-based", zorder=3)
    ax.bar(x + w / 2, cf, w, color=PALETTE["purple"],
           label="Cell-free / biochemical", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([MARGIN_LABEL[m] for m in MARGINS])
    ax.set_xlabel("Potency window below bound")
    ax.set_ylabel("% cytotoxicity-separated")
    ax.set_ylim(0, max(cb + cf) * 1.32)
    ax.legend(frameon=False, fontsize=8.5)
    for xi, v in zip(x - w / 2, cb):
        ax.text(xi, v + 0.9, f"{v:.1f}", ha="center", fontsize=8)
    for xi, v in zip(x + w / 2, cf):
        ax.text(xi, v + 0.9, f"{v:.1f}", ha="center", fontsize=8)
    panel_label(ax, "B", "Assay format")

    # ---- C: label-set shift ----
    ax = axes[2]
    style_panel(ax, grid_axis="y")
    t = t17.sort_values("jaccard_1x", ascending=False)
    xi = np.arange(len(t))
    for m, col, mk in zip(
        MARGINS,
        [PALETTE["navy"], PALETTE["sky"], PALETTE["crimson"]],
        ["o", "s", "^"],
    ):
        ax.plot(xi, t[f"jaccard_{m}"], mk + "-", color=col, ms=5, lw=1.4,
                label=f"{MARGIN_LABEL[m]} window", zorder=4)
    ax.set_xticks(xi)
    ax.set_xticklabels(t["family"], rotation=42, ha="right", fontsize=8)
    ax.set_ylabel("Jaccard overlap of chemical sets")
    ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, fontsize=8.5)
    panel_label(ax, "C", "Chemical label sets before vs after filtering")

    fig.tight_layout()
    save_fig(fig, "F4_sensitivity", FIG)


def main():
    t12 = pd.read_csv(TAB / "T12_bound_provenance_stratification.csv")
    t13 = pd.read_csv(TAB / "T13_family_by_bound_provenance.csv")
    t14 = pd.read_csv(TAB / "T14_window_testability.csv")
    t15 = pd.read_csv(TAB / "T15_promiscuity_coverage_normalised.csv")
    t17 = pd.read_csv(TAB / "T17_label_shift_estimated_bound.csv")
    t18 = pd.read_csv(TAB / "T18_format_stratification_estimated_bound.csv")
    t19 = pd.read_csv(TAB / "T19_flag_sensitivity_estimated_bound.csv")

    conc = load_conc_sample()
    print(f"conc sample n={len(conc):,} median={np.median(conc):.1f} µM", flush=True)

    fig1(t12, t14, conc)
    print("  F1_provenance", flush=True)
    fig2(t13)
    print("  F2_family_retention", flush=True)
    fig3(t15)
    print("  F3_promiscuity", flush=True)
    fig4(t19, t18, t17)
    print("  F4_sensitivity", flush=True)
    print("Done.")


if __name__ == "__main__":
    main()
