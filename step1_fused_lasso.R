# Step 1: Fused Lasso (per-biomarker change points)
#
# For each biomarker, fit a per-age-window logistic GLM to get a
# coefficient curve beta(a), then run a fused lasso (via flamCV) along
# the age grid to detect change points. The union of change points is
# used to bin age for a downstream Cox model.

library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(broom)
library(survival)
library(flam)

train_data <- read.csv("data/relax_prevent_cohort_time_prs_train.csv")
test_data  <- read.csv("data/relax_prevent_cohort_time_prs_test.csv")

OUTCOME    <- "cvd"
age_var    <- "age"
biomarkers <- c("sbp", "bmi", "hdl", "tc", "scre")

train_data[[age_var]] <- as.numeric(train_data[[age_var]])

rhs_full <- paste(biomarkers, collapse = " + ")

# ---------- helpers ----------

make_age_grid <- function(df, age_var = "age", grid_step = 1) {
  rng <- range(df[[age_var]], na.rm = TRUE)
  seq(floor(rng[1]), ceiling(rng[2]), by = grid_step)
}

# fit one GLM on the window [a_g - win, a_g + win] and return estimates
# + SEs for the target biomarkers
fit_one_window <- function(df, a_g, win, rhs_full, targets,
                           age_var = "age", outcome = "cvd") {
  dwin <- df %>%
    filter(is.finite(.data[[age_var]]),
           abs(.data[[age_var]] - a_g) <= win)

  if (nrow(dwin) < 50) {
    return(tibble(age_grid = a_g, term = targets,
                  estimate = NA_real_, se = NA_real_))
  }

  form <- as.formula(paste(outcome, "~", rhs_full))
  mod  <- try(glm(form, data = dwin, family = binomial()), silent = TRUE)
  if (inherits(mod, "try-error")) {
    return(tibble(age_grid = a_g, term = targets,
                  estimate = NA_real_, se = NA_real_))
  }

  co <- broom::tidy(mod) %>%
    filter(term %in% targets) %>%
    select(term, estimate, se = std.error)

  tibble(age_grid = a_g, term = targets) %>%
    left_join(co, by = "term")
}

# ---------- per-age GLMs ----------

age_grid <- make_age_grid(train_data, age_var, grid_step = 1)
win      <- 2  # +/- 2 years

perage <- map_dfr(
  age_grid,
  ~ fit_one_window(train_data, a_g = .x, win = win,
                   rhs_full = rhs_full, targets = biomarkers,
                   age_var = age_var, outcome = OUTCOME)
) %>%
  arrange(term, age_grid)

# ---------- fused lasso per biomarker (flam-CV picks lambda) ----------

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
                  cp = FALSE, cp_after = NA_character_, lambda = NA_real_))
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
                cp = FALSE, cp_after = NA_character_, lambda = NA_real_)
  out$beta_fused[keep] <- yhat

  jumps <- which(abs(diff(yhat)) > tol_jump)
  if (length(jumps) > 0) {
    kept_ages <- out$age_grid[keep]
    cp_age <- kept_ages[jumps + 1]
    out$cp[out$age_grid %in% cp_age]     <- TRUE
    out$lambda[out$age_grid %in% cp_age] <- lam_star
    for (a in cp_age) {
      prev_a <- kept_ages[which(kept_ages == a) - 1]
      out$cp_after[out$age_grid == a] <- paste0(prev_a, " -> ", a)
    }
  }
  out
}

fused_all_cont <- map_dfr(biomarkers,
                          ~ fuse_one_predictor(.x, perage, age_grid))

cp_table_cont <- fused_all_cont %>%
  filter(cp) %>%
  select(term, age_grid, cp_after, beta_fused, lambda) %>%
  arrange(term, age_grid)

print(cp_table_cont)

dir.create("result", showWarnings = FALSE)
write.csv(fused_all_cont, "result/fused_all_cont_flamCV.csv", row.names = FALSE)
write.csv(cp_table_cont,  "result/change_point_age_flamCV.csv", row.names = FALSE)

# ---------- union cuts -> Cox on age_refined ----------

union_cuts <- fused_all_cont %>%
  filter(cp) %>% pull(age_grid) %>% unique() %>% sort()
print(union_cuts)

make_age_refined <- function(df, cuts) {
  if (length(cuts) == 0) return(rep(1L, nrow(df)))
  as.integer(cut(df[[age_var]],
                 breaks = c(-Inf, cuts, Inf),
                 right = FALSE, include.lowest = TRUE))
}

train_data$age_refined <- make_age_refined(train_data, union_cuts)
test_data$age_refined  <- make_age_refined(test_data,  union_cuts)

fit_cox <- coxph(
  Surv(time2, cvd) ~ factor(age_refined) +
    sbp + bmi + hdl + tc + scre +
    factor(gender) + factor(race) + factor(smoker) +
    factor(anti) + factor(statin) + factor(diabete) + factor(ethnicity),
  data = train_data, x = TRUE, y = TRUE
)

# held-out C-index
lp_test <- predict(fit_cox, newdata = test_data, type = "lp")
c_test  <- survConcordance(Surv(time2, cvd) ~ lp_test, data = test_data)$concordance
cat("FLS individual C-index:", round(c_test, 6), "\n")

cat("logLik:", logLik(fit_cox), "\n")
cat("AIC:   ", AIC(fit_cox),    "\n")
cat("BIC:   ", BIC(fit_cox),    "\n")
