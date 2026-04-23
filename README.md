# TB Screening Cost-Effectiveness Analysis — UK Migrants

Cost-effectiveness analysis of 43 tuberculosis screening strategies for migrants entering the UK. Coupled decision tree–Markov cohort model in R.

## Requirements

- R ≥ 4.4
- Packages: `dplyr`, `ggplot2`, `patchwork`, `readxl`, `scales`, `tibble`, `parallel`

Install all at once:
```r
install.packages(c("dplyr", "ggplot2", "patchwork", "readxl", "scales", "tibble"))
```

## Files

| File | Purpose |
|---|---|
| `MasterTBModel.R` | Main orchestrator — runs all model components |
| `scripts/model_functions.R` | Config loading and all function definitions |
| `scripts/run_sa.R` | Deterministic and probabilistic sensitivity analyses |
| `scripts/run_transmission.R` | Societal perspective and transmission prevention |
| `scripts/run_plots_and_tables.R` | All figures and Excel output tables |
| `input/config.csv` | All 124 model parameters (values, distributions, sources) |
| `input/Screening_Decision_Tree.xlsx` | Decision tree initial conditions (Zenner et al. 2025 ERJ) |
| `output/` | Generated figures (TIFF) and results tables (CSV, Excel) |

## How to run

```bash
Rscript MasterTBModel.R
```

Or open in RStudio and source `MasterTBModel.R`. Runtime approximately 10–15 minutes (PSA: 1,000 simulations, parallelised).

All outputs are written to `output/body/` (figures and main tables) and `output/csv/` (results CSVs).

## Model structure

- 16 health states, 43 screening strategies (34 sequential + 8 parallel IGRA + passive case finding)
- 100,000 cohort, 660 monthly cycles (55-year lifetime horizon, entry age ~25)
- Primary perspective: societal (NHS costs + secondary transmission prevention savings)
- NHS perspective analysed as robustness check
- 3.5% discount rate for costs and QALYs (NICE PMG36 reference case)
- PSA: 1,000 Monte Carlo simulations (beta/gamma distributions)
- DSA: tornado diagrams for three frontier strategies
- Background mortality: age-varying ONS national life tables (25–79 years, 5-year bands)

## Efficient frontier (base case, societal perspective)

| Strategy | Cost per person (£) | Sequential ICER |
|---|---|---|
| Passive case finding | £203 | Reference |
| Cough+CXR (TB sx) | £209 | £2,205/QALY |
| Symptom screen+CXR | £220 | £7,963/QALY |
| Parallel Sx+QFT (Ultra) | £345 (societal) | £12,738/QALY |

All frontier strategies cost-effective at NICE £25,000–£35,000/QALY threshold.

## Citation

Kuzminova K, Allel Henriquez K, Zenner D, Haghparast-Bidgoli H. Cost-effectiveness of tuberculosis screening strategies for UK migrants: a coupled decision tree–Markov model. *The Lancet Regional Health – Europe*. 2026 [in preparation].
