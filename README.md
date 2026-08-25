# Auditing the ToxCast cytotoxicity filter (invitrodb v4.3)

Full research article package for invitrodb v4.3. The study asks what the published chemical-level cytotoxicity lower bound actually does when applied uniformly across a release: how often the bound is an estimate rather than a default ceiling, how often the comparison is answerable at the concentrations screened, and how much of a chemical's apparent multi-family breadth survives the filter once assay coverage is accounted for.

Public pipeline archive: https://doi.org/10.5281/zenodo.22100042 (all versions: https://doi.org/10.5281/zenodo.21855656)

## Findings

- 7,593 of 10,487 chemicals (72.4%) in the v4.3 cytotoxicity table carry the 1000 uM ceiling rather than an estimated bound; 6,203 have no active burst assay.
- Actives from those chemicals clear the bound at every window (99.9% at 1x, 97.2% at 10x), because the median top tested concentration is 76.6 uM.
- Among actives from chemicals with an estimated bound, only 19.8% separate from cytotoxicity at 1x (9.0% at 3x, 4.4% at 10x). The pooled figures of 32.1%, 22.9% and 18.6% are mixtures of the two regimes.
- The bound exceeds the highest tested concentration for 15.9% of actives; a 10x window falls below the lowest tested concentration for 8.9%.
- All 90 burst-annotated endpoints sit in two target families (cell cycle, cell morphology), and those two retain the least. Family retention tracks distance from the burst set, not mechanistic specificity, so it is not used here to rank families.
- Coverage-normalised promiscuity is the filter's defensible use: among estimated-bound chemicals the median fraction of tested families hit falls from 0.53 to 0.21, and 68.3% lose at least half of their active families. Among default-bound chemicals it does not move.

## Reproduce

```bash
python scripts/01_target_family_atlas.py        # pooled library, family and chemical tallies
python scripts/07_revision_analyses.py          # potency-window, assay-format and mc6-flag analyses (pooled)
python scripts/08_bound_provenance.py           # bound provenance, testability, coverage-normalised promiscuity
python scripts/09_estimated_bound_followups.py  # rank stability and label-set shift, estimated-bound stratum
python scripts/10_format_and_flags_estimated.py # format and flag analyses recomputed within stratum
python scripts/11_figures_revised.py            # manuscript figures F1-F4
```

Requires local invitrodb v4.3 extracts under `../ToxCast/derived/summary_extract/`:

- `mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv`
- `assay_annotations_invitrodb_v4_3_AUG2024.xlsx`
- `cytotox_invitrodb_v4_3_AUG2024.xlsx`

Scripts 08-11 stream the multi-concentration summary in chunks, so a typical workstation is sufficient. Each streaming pass takes roughly one minute.

## Outputs

- `tables/T1`-`T11` pooled family, chemical and sensitivity summaries
- `tables/T12`-`T15` bound provenance, family tallies by provenance, window testability, coverage-normalised promiscuity
- `tables/T16`-`T19` rank stability, label-set shift, format and flag analyses within the estimated-bound stratum
- `figures/F1_provenance`, `F2_family_retention`, `F3_promiscuity`, `F4_sensitivity` (png/pdf/eps/tif)
- `derived/BOUND_PROVENANCE.{txt,json}`, `derived/ESTIMATED_BOUND_FOLLOWUPS.{txt,json}`, `derived/FORMAT_FLAGS_ESTIMATED.{txt,json}`
- `manuscript/manuscript.tex` (elsarticle; *Computational Toxicology*)

## Definitions

- Active: `hitc >= 0.9`
- Cytotoxicity-separated: active with a finite, positive AC50 strictly below the chemical's cytotoxicity lower bound (1x window); 3x and 10x windows divide the bound by 3 and 10
- Bound provenance: `estimated` when the released lower bound is below the 1000 uM ceiling, `default` when it equals it
- Untestable comparison: bound above the curve's highest tested concentration, or required threshold below its lowest tested concentration
- Family: EPA `intended_target_family` (control and background families excluded from biological tallies)
- Coverage-normalised promiscuity: families hit divided by families the chemical was actually tested in

The cytotoxicity bound is the median AC50 across the release's burst-annotated endpoints, offset by three times a global median absolute deviation estimated once across a broadly screened core chemical set. Released values are used as distributed and are not recomputed here.
