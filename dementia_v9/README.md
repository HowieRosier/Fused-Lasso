# Dementia change points on All of Us v9

Scripts for the dementia analysis on the All of Us Controlled Tier v9 data (C2025Q4R6): the cohort, the
change-point method, the two baselines, and the checks. They run in R inside the Researcher Workbench.
This folder holds code only. Participant-level data stay in the Workbench.

## What the analysis does

1. **Per-age effects.** Participants are grouped by baseline age (one-year groups; neighbouring ages are merged
   until a group has at least 40 dementia cases). In each group a logistic model relates dementia to the
   biomarkers jointly, each per SD, with adjustment. This gives an effect and a standard error per biomarker
   and age group.
2. **Change points.** Three lines search these effects for ages at which an effect jumps to a new level:
   - line 1: group lasso (the biomarkers change at shared ages), the change-point code of `../dementia/step2_group_lasso.R`, unchanged;
   - line 2: fused lasso (each biomarker on its own), the change-point code of `../dementia/step1_fused_lasso.R`, unchanged;
   - line 3: sparse group lasso, fine-tuned. The penalty is re-weighted step by step. At each stage the mixing
     weight and the penalty strength are chosen along the path by a BIC-type criterion,
     misfit + log(n) x (1.5 per change age + 1 per jump), scaled by a multiplier.
3. **Penalty multiplier.** Set by rule: the smallest value at which data simulated with no change produce a
   change point in at most 5% of data sets.
4. **Cox evaluation.** Each line's change points enter a Cox model in which a biomarker has its own effect
   between its change points; the model is compared with one whose effects are constant in age.

## Scripts, in run order

| Script | Purpose |
|---|---|
| `dementia_cp_02_v9_cohort.R` | Builds the cohort from the OMOP tables and writes two files: the wide file (ages 30-79) and the analysis set (ages 45-79) |
| `dementia_cp_11_extra_labs.R` | Adds 19 further biomarkers (heart rate, glucose, eGFR, blood count, liver panel and others) to the cohort files; the cohort itself is not rebuilt |
| `dementia_cp_07_tune.R` | Per-age models and line 3; tests of whether an effect varies with age; step against straight line; figure of the per-age estimates |
| `dementia_cp_08_compare.R` | Lines 1, 2 and 3 on the same per-age estimates (uses the functions of `dementia_cp_07_tune.R`; the code of lines 1 and 2 is copied into the script) |
| `dementia_cp_10_cox.R` | Cox evaluation of the change points of each line (run after `dementia_cp_08_compare.R` in the same R session) |
| `dementia_cp_12_panel.R` | Screen of the further biomarkers: each added to the five original ones on its own, and joint panels |
| `dementia_cp_13_sensitivity.R` | Line 3 under other settings: age grouping, penalty multiplier, fine-tuning |

Each script starts with a header that describes its rules, options and output. Scripts print model-level results
only (ages, coefficients, rates).

## How to run

Copy the scripts to the home folder of an R app in the Workbench. The cohort scripts find the dataset and the
billing project from the workspace.

```r
source("~/dementia_cp_02_v9_cohort.R")       # cohort files in ~/dementia_v9/
source("~/dementia_cp_11_extra_labs.R")      # adds the further biomarkers

# five original biomarkers
ADJUST <- c("sex_at_birth", "anti", "imp_bmi", "imp_ldlc", "imp_hdl", "imp_trig")
source("~/dementia_cp_08_compare.R"); source("~/dementia_cp_10_cox.R")

# six biomarkers (five + heart rate)
PB <- c("sbp", "bmi", "ldlc", "hdl", "trig", "hr")
BIOMARKERS <- PB; ADJUST <- c("sex_at_birth", "anti", paste0("imp_", PB[-1])); BIN_MIN_EVENTS <- 40; BIC_STRICT <- 2.5; RUN_ORIGINAL <- FALSE
source("~/dementia_cp_08_compare.R")
BIOMARKERS <- PB; source("~/dementia_cp_10_cox.R")
BIOMARKERS <- PB; ADJUST <- c("sex_at_birth", "anti", paste0("imp_", PB[-1])); BIN_MIN_EVENTS <- 40; SEARCH <- "ft"; BIC_C <- 2.5; NULL_SIMS <- 200; N_BOOT <- 0; PRINT_PERAGE <- TRUE
source("~/dementia_cp_07_tune.R")

source("~/dementia_cp_12_panel.R")           # screen and joint panels (about 2 h on 4 cores)
source("~/dementia_cp_13_sensitivity.R")     # settings grid (about 1 h on 4 cores)
```

R packages: bigrquery, jsonlite, lubridate, dplyr, tidyr, purrr, tibble, broom, sparsegl, gglasso, flam, survival.

## Notes

- Settings of the main analysis (age grouping, penalty multiplier, fine-tuning) and the rules of the screen were
  fixed before the corresponding results were seen; each script states its rules in the header.
- Follow the All of Us Data and Statistics Dissemination Policy when taking results out of the Workbench
  (no participant-level data; no counts of 1 to 20).
