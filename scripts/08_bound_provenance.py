"""
Reviewer-driven analyses for the revised manuscript.

Addresses three technical objections to using the invitrodb cytotoxicity
lower bound as a uniform library-wide filter:

1. Bound provenance. In v4.3 most chemicals carry the 1000 uM ceiling rather
   than an estimated cytotoxicity point. Actives from those chemicals pass any
   "AC50 below bound" test trivially, because ToxCast screening rarely reaches
   1000 uM. We stratify every result by whether the chemical has an estimated
   or a default bound.

2. Testability of the comparison. Using conc_min / conc_max of each fitted
   curve we quantify how often the bound sits above the tested range (the
   comparison is vacuous) and how often a 3x or 10x window falls below the
   lowest tested concentration (the window is unobservable).

3. Coverage-normalised promiscuity. Chemicals are not screened in identical
   panels, so raw family counts confound promiscuity with panel coverage. We
   express promiscuity as the fraction of the biological families a chemical
   was actually tested in.

Outputs T12-T15, derived/BOUND_PROVENANCE.{txt,json} and figures/F5_bound_provenance.
"""
from __future__ import annotations

import json
import sys
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
FIG = ROOT / "figures"
DER = ROOT / "derived"
sys.path.insert(0, str(ROOT))
from figstyle import PALETTE, apply_style, panel_label, save_fig, style_panel  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

for d in (TAB, FIG, DER):
    d.mkdir(parents=True, exist_ok=True)

HIT_THRESH = 0.9
DEFAULT_BOUND_UM = 1000.0
DEFAULT_TOL = 0.1  # bound >= 999.9 counts as the ceiling
MARGINS = {"1x": 1.0, "3x": 3.0, "10x": 10.0}
USECOLS = ["chid", "casn", "chnm", "aeid", "hitc", "ac50", "conc_min", "conc_max"]
CONTROL_FAMILIES = {"background measurement", "background control"}
FAMILY_RENAME = {
    "channel 1": "ion channel (panel 1)",
    "channel 2": "ion channel (panel 2)",
}


def load_ann() -> pd.DataFrame:
    ann = pd.read_excel(
        ANN,
        usecols=["aeid", "intended_target_family", "burst_assay", "cell_viability_assay"],
    )
    ann = ann.dropna(subset=["aeid"]).copy()
    ann["aeid"] = ann["aeid"].astype(np.int64)
    ann["intended_target_family"] = (
        ann["intended_target_family"].fillna("unannotated").astype(str).str.strip()
    )
    fam = ann["intended_target_family"].str.lower()
    ann["family_display"] = fam.map(lambda x: FAMILY_RENAME.get(x, x))
    ann["is_control_family"] = fam.isin(CONTROL_FAMILIES)
    ann["is_burst"] = ann["burst_assay"].fillna(0).astype(int).eq(1)
    ann["is_viability"] = ann["cell_viability_assay"].fillna(0).astype(int).eq(1)
    return ann


def load_cyto() -> pd.DataFrame:
    cy = pd.read_excel(
        CYTO,
        usecols=[
            "chid",
            "chnm",
            "cytotox_median_um",
            "cytotox_lower_bound_um",
            "ntested",
            "nhit",
        ],
    )
    cy = cy.dropna(subset=["chid"]).copy()
    cy["chid"] = cy["chid"].astype(np.int64)
    cy["is_default_bound"] = cy["cytotox_lower_bound_um"] >= (DEFAULT_BOUND_UM - DEFAULT_TOL)
    cy["no_burst_hit"] = cy["nhit"].fillna(0).astype(float).eq(0)
    return cy


def stream(ann: pd.DataFrame, cy: pd.DataFrame) -> dict:
    fam_by_aeid = ann.set_index("aeid")["family_display"].to_dict()
    ctrl_by_aeid = ann.set_index("aeid")["is_control_family"].to_dict()
    burst_by_aeid = ann.set_index("aeid")["is_burst"].to_dict()

    bound = cy.set_index("chid")["cytotox_lower_bound_um"].to_dict()
    is_default = cy.set_index("chid")["is_default_bound"].to_dict()
    chem_name = cy.set_index("chid")["chnm"].to_dict()

    # counters keyed by stratum: "estimated" | "default"
    n_active = Counter()
    n_sep = {m: Counter() for m in MARGINS}
    n_chem_active: dict[str, set[int]] = defaultdict(set)

    # testability
    n_bound_above_max = 0            # bound > conc_max -> comparison vacuous
    n_bound_above_max_est = 0
    n_window_below_min = {m: 0 for m in MARGINS}
    n_active_with_conc = 0

    # family x stratum
    fam_active = {"estimated": Counter(), "default": Counter()}
    fam_sep = {m: {"estimated": Counter(), "default": Counter()} for m in MARGINS}
    fam_burst_flag: dict[str, bool] = {}

    # coverage-normalised promiscuity (biological families only)
    chem_fam_tested: dict[int, set[str]] = defaultdict(set)
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
        is_burst = chunk["aeid"].map(burst_by_aeid).fillna(False).astype(bool)

        # tested coverage (biological families only)
        bio = ~is_ctrl.to_numpy()
        for chid, fam in zip(
            chunk["chid"].to_numpy()[bio], families.to_numpy()[bio]
        ):
            chem_fam_tested[int(chid)].add(str(fam))

        for fam, b in zip(families.to_numpy(), is_burst.to_numpy()):
            if str(fam) not in fam_burst_flag:
                fam_burst_flag[str(fam)] = bool(b)
            elif b:
                fam_burst_flag[str(fam)] = True

        hit = chunk["hitc"].to_numpy(dtype=float) >= HIT_THRESH
        if not hit.any():
            continue
        h = chunk.loc[hit]
        hf = families.loc[hit].to_numpy()
        hc = is_ctrl.loc[hit].to_numpy()
        chids = h["chid"].to_numpy(dtype=np.int64)
        ac = h["ac50"].to_numpy(dtype=float)
        cmin = h["conc_min"].to_numpy(dtype=float)
        cmax = h["conc_max"].to_numpy(dtype=float)
        bnd = np.array([bound.get(int(c), np.nan) for c in chids], dtype=float)
        dflt = np.array([bool(is_default.get(int(c), False)) for c in chids])

        strata = np.where(dflt, "default", "estimated")
        for s, cnt in pd.Series(strata).value_counts().items():
            n_active[str(s)] += int(cnt)
        for s, cid in zip(strata, chids):
            n_chem_active[str(s)].add(int(cid))

        valid_ac = np.isfinite(ac) & (ac > 0)
        have_conc = np.isfinite(cmax) & np.isfinite(cmin)
        n_active_with_conc += int(have_conc.sum())
        vac = have_conc & np.isfinite(bnd) & (bnd > cmax)
        n_bound_above_max += int(vac.sum())
        n_bound_above_max_est += int((vac & ~dflt).sum())

        for m, div in MARGINS.items():
            thresh = bnd / div
            sep = valid_ac & np.isfinite(thresh) & (ac < thresh)
            for s, cnt in pd.Series(np.where(dflt, "default", "estimated")[sep]).value_counts().items():
                n_sep[m][str(s)] += int(cnt)
            n_window_below_min[m] += int(
                (have_conc & np.isfinite(thresh) & (thresh < cmin)).sum()
            )
            for fam, s in zip(hf[sep & ~hc], strata[sep & ~hc]):
                fam_sep[m][str(s)][str(fam)] += 1

        for fam, s in zip(hf[~hc], strata[~hc]):
            fam_active[str(s)][str(fam)] += 1

        sep1 = valid_ac & np.isfinite(bnd) & (ac < bnd)
        for chid, fam, isbio, s1 in zip(chids, hf, ~hc, sep1):
            if not isbio:
                continue
            cid = int(chid)
            chem_fam_active[cid].add(str(fam))
            if s1:
                chem_fam_sep[cid].add(str(fam))

        if ci % 4 == 0:
            print(f"  streamed {n_rows:,} rows ...", flush=True)

    return {
        "n_rows": n_rows,
        "n_active": n_active,
        "n_sep": n_sep,
        "n_chem_active": {k: len(v) for k, v in n_chem_active.items()},
        "n_bound_above_max": n_bound_above_max,
        "n_bound_above_max_est": n_bound_above_max_est,
        "n_window_below_min": n_window_below_min,
        "n_active_with_conc": n_active_with_conc,
        "fam_active": fam_active,
        "fam_sep": fam_sep,
        "fam_burst_flag": fam_burst_flag,
        "chem_fam_tested": chem_fam_tested,
        "chem_fam_active": chem_fam_active,
        "chem_fam_sep": chem_fam_sep,
        "chem_name": chem_name,
    }


def main():
    print("Loading annotations / cytotox ...", flush=True)
    ann = load_ann()
    cy = load_cyto()

    n_chem = len(cy)
    n_default = int(cy["is_default_bound"].sum())
    n_noburst = int(cy["no_burst_hit"].sum())
    print(
        f"  chemicals={n_chem:,} default-bound={n_default:,} "
        f"({100*n_default/n_chem:.1f}%) no-burst-hit={n_noburst:,}",
        flush=True,
    )

    agg = stream(ann, cy)

    # ---------- T12: stratified retention ----------
    rows = []
    for stratum in ("estimated", "default"):
        act = agg["n_active"][stratum]
        row = {
            "bound_provenance": stratum,
            "n_chemicals_with_active": agg["n_chem_active"].get(stratum, 0),
            "n_active": act,
        }
        for m in MARGINS:
            s = agg["n_sep"][m][stratum]
            row[f"n_separated_{m}"] = s
            row[f"pct_separated_{m}"] = 100 * s / act if act else np.nan
        rows.append(row)
    t12 = pd.DataFrame(rows)
    t12.to_csv(TAB / "T12_bound_provenance_stratification.csv", index=False)

    # ---------- T14: testability ----------
    t14 = pd.DataFrame(
        [
            {
                "quantity": "actives with finite conc_min and conc_max",
                "count": agg["n_active_with_conc"],
            },
            {
                "quantity": "actives whose cytotox bound exceeds max tested conc",
                "count": agg["n_bound_above_max"],
            },
            {
                "quantity": "same, estimated-bound chemicals only",
                "count": agg["n_bound_above_max_est"],
            },
        ]
        + [
            {
                "quantity": f"actives whose {m} window falls below min tested conc",
                "count": agg["n_window_below_min"][m],
            }
            for m in MARGINS
        ]
    )
    t14.to_csv(TAB / "T14_window_testability.csv", index=False)

    # ---------- T13: family table, estimated-bound chemicals only ----------
    fams = sorted(set(agg["fam_active"]["estimated"]) | set(agg["fam_active"]["default"]))
    frows = []
    for fam in fams:
        a_est = agg["fam_active"]["estimated"][fam]
        a_def = agg["fam_active"]["default"][fam]
        row = {
            "family": fam,
            "is_burst_related": bool(agg["fam_burst_flag"].get(fam, False)),
            "n_active_estimated_bound": a_est,
            "n_active_default_bound": a_def,
            "pct_actives_from_default_bound": (
                100 * a_def / (a_est + a_def) if (a_est + a_def) else np.nan
            ),
        }
        for m in MARGINS:
            s_est = agg["fam_sep"][m]["estimated"][fam]
            row[f"n_separated_{m}_estimated_bound"] = s_est
            row[f"pct_separated_{m}_estimated_bound"] = (
                100 * s_est / a_est if a_est else np.nan
            )
        frows.append(row)
    t13 = pd.DataFrame(frows).sort_values(
        "n_active_estimated_bound", ascending=False
    )
    t13.to_csv(TAB / "T13_family_by_bound_provenance.csv", index=False)

    # ---------- T15: coverage-normalised promiscuity ----------
    default_set = set(cy.loc[cy["is_default_bound"], "chid"].astype(int))
    prows = []
    for cid, tested in agg["chem_fam_tested"].items():
        n_t = len(tested)
        n_a = len(agg["chem_fam_active"].get(cid, ()))
        n_s = len(agg["chem_fam_sep"].get(cid, ()))
        if n_a == 0:
            continue
        prows.append(
            {
                "chid": cid,
                "chnm": agg["chem_name"].get(cid, ""),
                "bound_provenance": "default" if cid in default_set else "estimated",
                "n_biofam_tested": n_t,
                "n_biofam_active": n_a,
                "n_biofam_separated": n_s,
                "frac_tested_active": n_a / n_t if n_t else np.nan,
                "frac_tested_separated": n_s / n_t if n_t else np.nan,
                "family_loss": n_a - n_s,
            }
        )
    t15 = pd.DataFrame(prows)
    t15.to_csv(TAB / "T15_promiscuity_coverage_normalised.csv", index=False)

    est = t15.loc[t15["bound_provenance"] == "estimated"]
    dfl = t15.loc[t15["bound_provenance"] == "default"]

    summary = {
        "n_chemicals_cytotox_table": n_chem,
        "n_default_bound_chemicals": n_default,
        "pct_default_bound_chemicals": 100 * n_default / n_chem,
        "n_no_burst_hit_chemicals": n_noburst,
        "n_rows_streamed": agg["n_rows"],
        "actives_estimated_bound": agg["n_active"]["estimated"],
        "actives_default_bound": agg["n_active"]["default"],
        "pct_actives_from_default_bound": 100
        * agg["n_active"]["default"]
        / max(agg["n_active"]["estimated"] + agg["n_active"]["default"], 1),
        "n_active_with_conc": agg["n_active_with_conc"],
        "n_bound_above_max_conc": agg["n_bound_above_max"],
        "pct_bound_above_max_conc": 100
        * agg["n_bound_above_max"]
        / max(agg["n_active_with_conc"], 1),
        "n_bound_above_max_conc_estimated": agg["n_bound_above_max_est"],
    }
    for m in MARGINS:
        e = agg["n_active"]["estimated"]
        d = agg["n_active"]["default"]
        summary[f"pct_separated_{m}_estimated_bound"] = (
            100 * agg["n_sep"][m]["estimated"] / e if e else np.nan
        )
        summary[f"pct_separated_{m}_default_bound"] = (
            100 * agg["n_sep"][m]["default"] / d if d else np.nan
        )
        summary[f"pct_separated_{m}_pooled"] = (
            100 * (agg["n_sep"][m]["estimated"] + agg["n_sep"][m]["default"]) / max(e + d, 1)
        )
        summary[f"n_window_below_min_conc_{m}"] = agg["n_window_below_min"][m]
        summary[f"pct_window_below_min_conc_{m}"] = (
            100 * agg["n_window_below_min"][m] / max(agg["n_active_with_conc"], 1)
        )
    summary.update(
        {
            "median_frac_tested_active_estimated": float(est["frac_tested_active"].median()),
            "median_frac_tested_separated_estimated": float(
                est["frac_tested_separated"].median()
            ),
            "median_frac_tested_active_default": float(dfl["frac_tested_active"].median()),
            "median_frac_tested_separated_default": float(
                dfl["frac_tested_separated"].median()
            ),
            "median_biofam_tested_estimated": float(est["n_biofam_tested"].median()),
            "median_biofam_active_estimated": float(est["n_biofam_active"].median()),
            "median_biofam_separated_estimated": float(est["n_biofam_separated"].median()),
            "n_active_chemicals_estimated": int(len(est)),
            "n_active_chemicals_default": int(len(dfl)),
            "pct_est_chem_losing_half_families": float(
                100 * ((est["n_biofam_separated"] / est["n_biofam_active"]) <= 0.5).mean()
            ),
        }
    )

    (DER / "BOUND_PROVENANCE.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    lines = [
        "Bound provenance and testability analyses (invitrodb v4.3)",
        "",
        f"Chemicals in cytotox table: {n_chem:,}",
        f"  with default 1000 uM lower bound: {n_default:,} ({100*n_default/n_chem:.1f}%)",
        f"  with no active burst assay (nhit=0): {n_noburst:,}",
        "",
        f"Active curves from estimated-bound chemicals: {agg['n_active']['estimated']:,}",
        f"Active curves from default-bound chemicals:   {agg['n_active']['default']:,}"
        f" ({summary['pct_actives_from_default_bound']:.1f}% of actives)",
        "",
        "Cytotoxicity-separated fraction by bound provenance:",
    ]
    for m in MARGINS:
        lines.append(
            f"  {m}: estimated {summary[f'pct_separated_{m}_estimated_bound']:.1f}%  "
            f"default {summary[f'pct_separated_{m}_default_bound']:.1f}%  "
            f"pooled {summary[f'pct_separated_{m}_pooled']:.1f}%"
        )
    lines += [
        "",
        "Testability of the comparison:",
        f"  actives with tested-range info: {agg['n_active_with_conc']:,}",
        f"  bound above max tested conc: {agg['n_bound_above_max']:,} "
        f"({summary['pct_bound_above_max_conc']:.1f}%)",
        f"    of which estimated-bound: {agg['n_bound_above_max_est']:,}",
    ]
    for m in MARGINS:
        lines.append(
            f"  {m} window below min tested conc: {agg['n_window_below_min'][m]:,} "
            f"({summary[f'pct_window_below_min_conc_{m}']:.1f}%)"
        )
    lines += [
        "",
        "Coverage-normalised promiscuity (estimated-bound chemicals):",
        f"  n chemicals with >=1 active: {len(est):,}",
        f"  median biological families tested: {est['n_biofam_tested'].median():.0f}",
        f"  median families active -> separated: "
        f"{est['n_biofam_active'].median():.0f} -> {est['n_biofam_separated'].median():.0f}",
        f"  median fraction of tested families active: "
        f"{est['frac_tested_active'].median():.3f}",
        f"  median fraction of tested families separated: "
        f"{est['frac_tested_separated'].median():.3f}",
        f"  % losing >=half of active families: "
        f"{summary['pct_est_chem_losing_half_families']:.1f}%",
        "",
        "Same, default-bound chemicals (for contrast):",
        f"  n chemicals with >=1 active: {len(dfl):,}",
        f"  median fraction of tested families active: "
        f"{dfl['frac_tested_active'].median():.3f}",
        f"  median fraction of tested families separated: "
        f"{dfl['frac_tested_separated'].median():.3f}",
    ]
    (DER / "BOUND_PROVENANCE.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    make_figure(t12, t13, t15, summary)
    print("\nWrote T12-T15, derived/BOUND_PROVENANCE.*, figures/F5_bound_provenance")


def make_figure(t12: pd.DataFrame, t13: pd.DataFrame, t15: pd.DataFrame, summary: dict):
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.1))

    # Panel A: separated fraction by provenance and margin
    ax = axes[0]
    style_panel(ax, grid_axis="y")
    margins = list(MARGINS)
    est = [summary[f"pct_separated_{m}_estimated_bound"] for m in margins]
    dfl = [summary[f"pct_separated_{m}_default_bound"] for m in margins]
    x = np.arange(len(margins))
    w = 0.36
    ax.bar(x - w / 2, est, w, color=PALETTE["navy"], label="Estimated bound", zorder=3)
    ax.bar(x + w / 2, dfl, w, color=PALETTE["crimson"], label="Default 1000 µM", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{m} window" for m in margins])
    ax.set_xlabel("Required potency window below bound")
    ax.set_ylabel("% of actives cytotoxicity-separated")
    ax.set_ylim(0, 105)
    ax.legend(frameon=False, fontsize=8)
    for xi, v in zip(x - w / 2, est):
        ax.text(xi, v + 2, f"{v:.0f}", ha="center", fontsize=8)
    for xi, v in zip(x + w / 2, dfl):
        ax.text(xi, v + 2, f"{v:.0f}", ha="center", fontsize=8)
    panel_label(ax, "A", "Bound provenance drives retention")

    # Panel B: testability
    ax = axes[1]
    style_panel(ax, grid_axis="y")
    labels = ["Bound above\nmax tested\nconc.", "3× window\nbelow min\ntested conc.",
              "10× window\nbelow min\ntested conc."]
    vals = [
        summary["pct_bound_above_max_conc"],
        summary["pct_window_below_min_conc_3x"],
        summary["pct_window_below_min_conc_10x"],
    ]
    ax.bar(range(3), vals, color=PALETTE["crimson"], width=0.65, zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("% of actives with tested-range info")
    ymax = max(vals) * 1.25 if max(vals) > 0 else 1
    ax.set_ylim(0, ymax)
    for i, v in enumerate(vals):
        ax.text(i, v + ymax * 0.02, f"{v:.1f}%", ha="center", fontsize=8)
    panel_label(ax, "B", "When the comparison is untestable")

    # Panel C: coverage-normalised promiscuity, estimated-bound chemicals
    ax = axes[2]
    style_panel(ax, grid_axis="y")
    e = t15.loc[t15["bound_provenance"] == "estimated"]
    bins = np.linspace(0, 1, 21)
    ax.hist(e["frac_tested_active"], bins=bins, color=PALETTE["navy"], alpha=0.85,
            edgecolor="white", label="Any active", zorder=3)
    ax.hist(e["frac_tested_separated"], bins=bins, color=PALETTE["green"], alpha=0.65,
            edgecolor="white", label="Cytotox-separated", zorder=4)
    ax.set_xlabel("Fraction of tested biological families hit")
    ax.set_ylabel(f"Chemicals (n={len(e):,}, estimated bound)")
    ax.legend(frameon=False, fontsize=8)
    panel_label(ax, "C", "Coverage-normalised promiscuity")

    fig.tight_layout()
    save_fig(fig, "F5_bound_provenance", FIG)


if __name__ == "__main__":
    main()
