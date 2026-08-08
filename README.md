# ToxCast target-family bioactivity atlas (invitrodb v4.3)

Pipeline to build a cytotoxicity-filtered atlas of ToxCast intended target
families: which families retain selective chemical activity after actives at or
above the chemical cytotoxicity lower bound are removed, and how chemical
multi-family promiscuity changes under that filter.

Archive: https://doi.org/10.5281/zenodo.21855657

Data source (not redistributed here):

- US EPA (2025). ToxCast invitrodb version 4.3.
  https://doi.org/10.23645/epacomptox.6062623.v14

## Contents

```
scripts/01_target_family_atlas.py   Primary atlas analysis
scripts/07_revision_analyses.py     Selectivity margins, assay format, mc6 flags
figstyle.py                         Figure styling helpers
tables/                             Aggregate result tables
figures/                            Figures F1-F4 (png/pdf/eps/tif)
derived/RESULTS.txt                 Key numeric summary
derived/REVISION_RESULTS.*          Sensitivity-analysis summary
requirements.txt                    Python dependencies
```

## Requirements

- Python 3.10+
- Public invitrodb v4.3 extracts:
  - `mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv`
  - `assay_annotations_invitrodb_v4_3_AUG2024.xlsx`
  - `cytotox_invitrodb_v4_3_AUG2024.xlsx`

## Reproduce

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
pip install -r requirements.txt

# Point to the EPA summary extract folder:
set INVITRODB_SUMMARY_DIR=C:\path\to\summary_extract

python scripts/01_target_family_atlas.py
python scripts/07_revision_analyses.py
```

On Unix shells use `export` instead of `set`. Individual file paths can also be
set with `INVITRODB_MC56`, `INVITRODB_ANNOTATIONS`, and `INVITRODB_CYTOTOX`.

## Definitions

| Term | Rule |
|------|------|
| Active | Continuous hit call `hitc >= 0.9` |
| Selective (1x) | Active and AC50 strictly below the chemical cytotoxicity lower bound |
| Selective (3x / 10x) | AC50 below bound/3 or bound/10 |
| Cytotoxicity-linked | Active but not selective under the chosen window |
| Family | EPA intended target family (control/background families excluded from biological rankings) |

Outputs refresh `tables/`, `figures/F1_main.*` through `F4_sensitivity.*`, and
`derived/RESULTS.txt` / `REVISION_RESULTS.*`.
