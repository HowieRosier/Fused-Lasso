## step1_fused_lasso.R
## Step 1: Per-biomarker Fused Lasso for age-related change point detection
## Pipeline: adaptive Gaussian kernel (matched with Step 3)
##           + flam::flamCV (per-biomarker, no shared structure)
##
## Role: independent verification of biomarker-specific transitions
##       that Step 3's shared-CP framework may miss.

library(dplyr)
library(tidyverse)
library(broom)
library(survival)
library(flam)

# ---- 1. Data ----------------------------------------------------------------

system("mkdir -p ~/data")
system(paste0("gsutil cp ", Sys.getenv("WORKSPACE_BUCKET"),
              "/relax_prevent_cohort_time_prs.csv ~/data/"))
system(paste0("gsutil cp ", Sys.getenv("WORKSPACE_BUCKET"),
              "/prs_top100.tsv ~/"))

data <- read.csv("~/data/relax_prevent_cohort_time_prs.csv")
prs  <- read.delim("~/prs_top100.tsv")
prs$person_id <- prs$IID
data <- data %>%
  inner_join(prs %>% mutate(prs = PRS_STD), by = "person_id") %>%
  select(-PRS_STD)

data$ethnicity[data$ethnicity %in% c(
  "PMI: Skip", "What Race Ethnicity: Race Ethnicity None Of These",
  "PMI: Prefer Not To Answer")] <- "No information"

data$race[data$race %in% c(
  "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander",
  "I prefer not to answer", "None Indicated", "None of these",
  "PMI: Skip", "Middle Eastern or North African",
  "More than one population")] <- "Others"

data$gender_clean <- case_when(
  data$gender == "Female" ~ "Female",
  data$gender == "Male"   ~ "Male",
  TRUE                    ~ "Other"
)

set.seed(42)
train_idx  <- sample(seq_len(nrow(data)), size = floor(0.60 * nrow(data)))
train_data <- data[train_idx, ]
test_data  <- data[-train_idx, ]

# ---- 2. Per-age GLM (adaptive Gaussian kernel, same as Steps 2 & 3) -------

biomarkers <- c("sbp", "bmi", "hdl", "tc", "scre")
age_grid   <- 30:79

fit_one_adaptive_kernel <- function(df, a_g, targets, min_events = 100) {
  df <- df %>% filter(is.finite(age))
  for (h in seq(1, 10, by = 0.5)) {
    w    <- exp(-((df$age - a_g)^2) / (2 * h^2))
    keep <- w > 0.01 * max(w)
    if (sum(df$cvd[keep]) >= min_events) break
  }
  dwin <- df[keep, ]; wk <- w[keep]
  if (nrow(dwin) < 50)
    return(tibble(age_grid = a_g, term = targets,
                  estimate = NA_real_, se = NA_real_))
  form <- as.formula(paste("cvd ~", paste(targets, collapse = " + ")))
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

# ---- 5. Verify golden standard transitions --------------------------------

golden <- c(37, 44, 51, 61, 68, 78)
cat("\n===== Golden standard verification =====\n")
for (g in golden) {
  hits <- cp_table %>%
    filter(abs(age_grid - g) <= 1) %>%
    select(term, age_grid)
  if (nrow(hits) > 0) {
    bms <- paste(sprintf("%s@%d", toupper(hits$term), hits$age_grid), collapse = ", ")
    cat(sprintf("  Golden %d: FOUND by %s\n", g, bms))
  } else {
    cat(sprintf("  Golden %d: not found\n", g))
  }
}

# ---- 6. Cox model (using union of CPs) -------------------------------------

predictors_cat <- c("diabete", "gender_clean", "race", "smoker",
                    "anti", "statin", "ethnicity")

if (length(union_cuts) > 0) {
  data$age_regime <- cut(data$age, breaks = c(-Inf, union_cuts, Inf),
                         right = FALSE, include.lowest = TRUE)
} else {
  data$age_regime <- 1L
}

form_cox <- as.formula(paste(
  "Surv(time2, cvd) ~ factor(age_regime) + sbp + bmi + hdl + tc + scre +",
  paste0("factor(", predictors_cat, ")", collapse = " + ")
))
fit_cox <- suppressWarnings(coxph(form_cox, data = data, x = TRUE, y = TRUE))

set.seed(42)
tr_idx <- sample(seq_len(nrow(data)), floor(0.60 * nrow(data)))
fit_tr <- suppressWarnings(coxph(form_cox, data = data[tr_idx, ]))
lp_te  <- predict(fit_tr, newdata = data[-tr_idx, ], type = "lp")
c_te   <- concordance(Surv(time2, cvd) ~ lp_te,
                      data = data[-tr_idx, ], reverse = TRUE)$concordance

cat(sprintf("\nCox BIC: %.1f\n", BIC(fit_cox)))
cat(sprintf("Test C-index: %.4f\n", c_te))
