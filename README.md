# Fused-Lasso
Change points detection for All Of Us

# Age-varying biomarker effects

Detect age-dependent change points in biomarker effects on CVD risk.
A per-age sliding-window logistic GLM gives a coefficient curve
$\beta_k(\text{age})$ per biomarker; two lasso variants then recover
piecewise-constant structure:

- **Step 1** — per-biomarker fused lasso (`flam::flamCV`)
- **Step 2** — shared change points across biomarkers (`gglasso::cv.gglasso`)

Each step refines age into bins at the detected change points, fits a Cox
model, and reports the held-out C-index.

## Layout

```
R/
├── 00_prepare_data.R     # merge PRS, 60/20/20 split
├── step1_fused_lasso.R   # per-biomarker fused lasso
└── step2_group_lasso.R   # shared group lasso
```

## Run

```r
install.packages(c("dplyr","purrr","tibble","broom","tidyr",
                   "Matrix","survival","flam","gglasso"))

source("R/00_prepare_data.R")   # once
source("R/step1_fused_lasso.R")
source("R/step2_group_lasso.R")
```

## Inputs

- `data/relax_prevent_cohort_time_prs.csv` — cohort (`time2`, `cvd`, `age`, biomarkers `sbp/bmi/hdl/tc/scre`, covariates)
- `data/prs_top100.tsv` — PRS (`IID`, `PRS_STD`)

## Outputs

`result/` — change point tables and fused coefficient curves for each step.
