## step1_fused_lasso.R
## Step 1: Per-biomarker Fused Lasso for age-related change point detection
## Pipeline: adaptive Gaussian kernel (matched with Step 3)
##           + flam::flamCV (per-biomarker, no shared structure)
##
## Role: independent verification of biomarker-specific transitions
##       that Step 3's shared-CP framework may miss.
##
## Cohort: All of Us dementia relax cohort
## Outcome: incident dementia
## Biomarkers: SBP, BMI, LDL-C, HDL-C, Triglycerides

library(dplyr)
library(tidyverse)
library(broom)
library(survival)
library(flam)

# ---- 1. Data ----------------------------------------------------------------

system("mkdir -p ~/data")
system(paste0("gsutil cp ", Sys.getenv("WORKSPACE_BUCKET"),
              "/dementia_relax_cohort_time_all_final1.csv ~/data/"))

data <- read.csv("~/data/dementia_relax_cohort_time_all_final1.csv")
data <- data %>% filter(is.finite(time2), time2 > 0, is.finite(age))

set.seed(42)
train_idx  <- sample(seq_len(nrow(data)), size = floor(0.60 * nrow(data)))
train_data <- data[train_idx, ]
test_data  <- data[-train_idx, ]

# ---- 2. Per-age GLM (adaptive Gaussian kernel, same as Steps 2 & 3) -------

biomarkers <- c("sbp", "bmi", "ldlc", "hdl", "trig")
age_grid   <- 45:79   # SAP main analysis (age >= 45)

fit_one_adaptive_kernel <- function(df, a_g, targets, min_events = 100) {
  df <- df %>% filter(is.finite(age))
  for (h in seq(1, 10, by = 0.5)) {
    w    <- exp(-((df$age - a_g)^2) / (2 * h^2))
    keep <- w > 0.01 * max(w)
    if (sum(df$Dementia[keep]) >= min_events) break
  }
  dwin <- df[keep, ]; wk <- w[keep]
  if (nrow(dwin) < 50)
    return(tibble(age_grid = a_g, term = targets,
                  estimate = NA_real_, se = NA_real_))
  form <- as.formula(paste("Dementia ~", paste(targets, collapse = " + ")))
  mod  <- try(glm(form, data = dwin, family = binomial(), weights = wk),
              silent = TRUE)
  if (inherits(mod, "try-error"))
    return(tibble(age_grid = a_g, term = targets,
                  estimate = NA_real_, se = NA_real_))
  co <- broom::tidy(mod) %>%
    filter(term %in% targets) %>%
    select(term, estimate, se = std.error)
  tibble(age_grid = a_g, term = targets) %>%
    left_join(co, by = "term")
}

cat("Building perage (adaptive kernel, min_events=100)...\n")
perage <- map_dfr(age_grid,
  ~ fit_one_adaptive_kernel(train_data, .x, biomarkers, min_events = 100)
) %>% arrange(term, age_grid)
cat("Done.\n\n")

# ---- 3. Per-biomarker fused lasso (flamCV) ---------------------------------
# Note: flam package does not support adaptive lasso weights.
# Step 1 uses standard CV-selected lambda only.

fuse_one_predictor <- function(pred, perage, age_grid_all,
                               n.fold = 10, seed = 123,
                               within1SE = FALSE, tol_jump = 1e-6) {
  dat <- perage %>% filter(term == pred) %>% arrange(age_grid)
  y  <- dat$estimate
  se <- dat$se
  keep <- is.finite(y) & is.finite(se) & se > 0

  B <- length(age_grid_all)
  if (sum(keep) < max(10, ceiling(B / 3))) {
    return(tibble(term = pred, age_grid = age_grid_all,
                  beta_raw = NA_real_, beta_fused = NA_real_,
                  cp = FALSE, lambda = NA_real_))
  }

  xk <- matrix(dat$age_grid[keep], ncol = 1)
  yk <- y[keep]

  cvfit <- flam::flamCV(
    x = xk, y = yk, family = "gaussian",
    alpha = 1, n.fold = n.fold,
    lambda.min.ratio = 0.1, seed = seed, within1SE = within1SE
  )
  lam_star <- cvfit$lambda.cv
  yhat <- as.numeric(predict(cvfit$flam.out, new.x = xk,
                             lambda = lam_star, alpha = cvfit$alpha))

  out <- tibble(term = pred, age_grid = dat$age_grid,
                beta_raw = y, beta_fused = NA_real_,
                cp = FALSE, lambda = NA_real_)
  out$beta_fused[keep] <- yhat

  jumps <- which(abs(diff(yhat)) > tol_jump)
  if (length(jumps) > 0) {
    kept_ages <- out$age_grid[keep]
    cp_age <- kept_ages[jumps + 1]
    out$cp[out$age_grid %in% cp_age]     <- TRUE
    out$lambda[out$age_grid %in% cp_age] <- lam_star
  }
  out
}

cat("Running per-biomarker fused lasso (flamCV)...\n")
fused_all <- map_dfr(biomarkers,
                     ~ fuse_one_predictor(.x, perage, age_grid))
cat("Done.\n\n")

# ---- 4. Results -------------------------------------------------------------

cp_table <- fused_all %>%
  filter(cp) %>%
  select(term, age_grid, beta_fused, lambda) %>%
  arrange(term, age_grid)

cat("===== Per-biomarker CPs =====\n")
for (bm in biomarkers) {
  cps <- cp_table %>% filter(term == bm) %>% pull(age_grid)
  cat(sprintf("  %-5s (%2d): %s\n", toupper(bm), length(cps),
              if (length(cps) == 0) "(none)" else paste(cps, collapse = ", ")))
}

union_cuts <- fused_all %>% filter(cp) %>% pull(age_grid) %>% unique() %>% sort()
cat(sprintf("\nUnion of all CPs (%d): %s\n", length(union_cuts),
            paste(union_cuts, collapse = ", ")))

# ---- 5. Cox model (using union of CPs) -------------------------------------

predictors_cat <- c("diabete", "gender", "race", "smoker",
                    "anti", "hearing_loss", "ethnicity",
                    "brain", "depression")

if (length(union_cuts) > 0) {
  data$age_regime <- cut(data$age, breaks = c(-Inf, union_cuts, Inf),
                         right = FALSE, include.lowest = TRUE)
} else {
  data$age_regime <- 1L
}

form_cox <- as.formula(paste(
  "Surv(time2, Dementia) ~ factor(age_regime) + sbp + bmi + ldlc + hdl + trig +",
  paste0("factor(", predictors_cat, ")", collapse = " + ")
))
fit_cox <- suppressWarnings(coxph(form_cox, data = data, x = TRUE, y = TRUE))

set.seed(42)
tr_idx <- sample(seq_len(nrow(data)), floor(0.60 * nrow(data)))
fit_tr <- suppressWarnings(coxph(form_cox, data = data[tr_idx, ]))
lp_te  <- predict(fit_tr, newdata = data[-tr_idx, ], type = "lp")
c_te   <- concordance(Surv(time2, Dementia) ~ lp_te,
                      data = data[-tr_idx, ], reverse = TRUE)$concordance

cat(sprintf("\nCox BIC: %.1f\n", BIC(fit_cox)))
cat(sprintf("Test C-index: %.4f\n", c_te))
