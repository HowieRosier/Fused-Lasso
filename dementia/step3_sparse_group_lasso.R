## step3_sparse_group_lasso.R
## Step 3: Sparse Group Lasso for shared change point detection
## Pipeline: adaptive Gaussian kernel + column standardization
##           + sparsegl with BIC-selected (alpha, lambda) under df_cap
##           + iterated adaptive lasso refit (until convergence)
##
## Cohort: All of Us dementia relax cohort
## Outcome: incident dementia
## Biomarkers: SBP, BMI, LDL-C, HDL-C, Triglycerides

library(dplyr)
library(tidyverse)
library(broom)
library(survival)
library(sparsegl)

# ---- 1. Data ----------------------------------------------------------------

system("mkdir -p ~/data")
system(paste0("gsutil cp ", Sys.getenv("WORKSPACE_BUCKET"),
              "/dementia_relax_cohort_time_all_final1.csv ~/data/"))

data <- read.csv("~/data/dementia_relax_cohort_time_all_final1.csv")
data <- data %>% filter(is.finite(time2), time2 > 0, is.finite(age))

# Filter to SAP main analysis age range
data <- data %>% filter(age >= 45, age < 80)

# Standardize biomarkers (mean=0, SD=1)
biomarkers <- c("sbp", "bmi", "ldlc", "hdl", "trig")
for (bm in biomarkers) {
  data[[bm]] <- as.numeric(scale(data[[bm]]))
}

# ---- 2. Per-age GLM (adaptive Gaussian kernel) -----------------------------

age_grid <- 45:79

fit_one_adaptive_kernel <- function(df, a_g, targets, min_events = 100) {
  df <- df %>% filter(is.finite(age))
  for (h in seq(1, 10, by = 0.5)) {
    w    <- exp(-((df$age - a_g)^2) / (2 * h^2))
    keep <- w > 0.01 * max(w)
    if (sum(df$Dementia[keep]) >= min_events) break
  }
  dwin <- df[keep, ]; wk <- w[keep]
  if (nrow(dwin) < 50)
    return(tibble(age_grid = a_g, term = targets, bandwidth = h,
                  estimate = NA_real_, se = NA_real_))
  form <- as.formula(paste("Dementia ~", paste(targets, collapse = " + ")))
  mod  <- try(glm(form, data = dwin, family = binomial(), weights = wk),
              silent = TRUE)
  if (inherits(mod, "try-error"))
    return(tibble(age_grid = a_g, term = targets, bandwidth = h,
                  estimate = NA_real_, se = NA_real_))
  co <- broom::tidy(mod) %>%
    filter(term %in% targets) %>%
    select(term, estimate, std.error) %>%
    rename(se = std.error)
  tibble(age_grid = a_g, term = targets, bandwidth = h) %>%
    left_join(co, by = "term")
}

perage <- map_dfr(age_grid,
  ~ fit_one_adaptive_kernel(data, .x, biomarkers, min_events = 100)
) %>% arrange(term, age_grid)

# ---- 3. Design matrix + column standardization ----------------------------

wide_est <- perage %>%
  select(age_grid, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  arrange(age_grid)
wide_se <- perage %>%
  select(age_grid, term, se) %>%
  pivot_wider(names_from = term, values_from = se) %>%
  arrange(age_grid)

Y_all  <- as.matrix(wide_est[, biomarkers])
SE_all <- as.matrix(wide_se[, biomarkers])
keep_row <- apply(Y_all, 1, function(r) all(is.finite(r))) &
            apply(SE_all, 1, function(r) all(is.finite(r) & r > 0))

age_keep <- wide_est$age_grid[keep_row]
Y  <- Y_all[keep_row, ]
SE <- SE_all[keep_row, ]
b  <- nrow(Y)
K  <- ncol(Y)

# Column standardization: balance biomarker scales
col_sd <- apply(Y, 2, sd)
col_sd[col_sd < 1e-8] <- 1
Y_std  <- sweep(Y,  2, col_sd, "/")
SE_std <- sweep(SE, 2, col_sd, "/")

# Reparameterize: beta_t = alpha + sum of delta_1..delta_{t-1}
y_vec   <- as.numeric(t(Y_std))
X_alpha <- kronecker(rep(1, b), diag(K))
J <- matrix(0, b, b - 1)
for (t in 2:b) J[t, 1:(t - 1)] <- 1
X_delta <- kronecker(J, diag(K))
X       <- cbind(X_alpha, X_delta)
group   <- c(1:K, rep((K + 1):(K + b - 1), each = K))

# SE-based group penalty
SE_by_age <- apply(SE_std, 1, mean)
n_groups  <- max(group)
pf_g <- numeric(n_groups); pf_g[1:K] <- 1
for (m in 1:(b - 1)) pf_g[K + m] <- 1 / SE_by_age[m]
pf_g[(K+1):n_groups] <- pf_g[(K+1):n_groups] /
  mean(pf_g[(K+1):n_groups])

# ---- 4. BIC selection over (alpha, lambda) under df_cap ------------------

df_cap     <- b - 1   # prevent saturated overfitting at n=p edge
alpha_grid <- seq(0.05, 0.95, by = 0.05)

pick_best_bic <- function(X, y_vec, group, pf_group,
                          pf_sparse = NULL, df_cap, alpha_grid) {
  best <- list(bic = Inf, alpha = NA, lam = NA, df = NA, coef = NULL)
  for (a in alpha_grid) {
    fit <- if (is.null(pf_sparse))
      sparsegl(X, y_vec, group = group, asparse = a,
               pf_group = pf_group, intercept = FALSE)
    else
      sparsegl(X, y_vec, group = group, asparse = a,
               pf_group = pf_group, pf_sparse = pf_sparse, intercept = FALSE)
    er <- estimate_risk(fit, X, type = "BIC", approx_df = FALSE)
    valid <- er$df <= df_cap
    if (!any(valid)) next
    j <- which.min(er$BIC[valid])
    if (er$BIC[valid][j] < best$bic) {
      best$bic   <- er$BIC[valid][j]
      best$alpha <- a
      best$lam   <- er$lambda[valid][j]
      best$df    <- er$df[valid][j]
      cv <- as.numeric(coef(fit, s = best$lam))
      if (length(cv) == ncol(X) + 1) cv <- cv[-1]
      best$coef <- cv
    }
  }
  best
}

# ---- 5. Stage 0: Initial SGL fit ------------------------------------------

best <- pick_best_bic(X, y_vec, group, pf_group = pf_g,
                      df_cap = df_cap, alpha_grid = alpha_grid)
cat(sprintf("Stage 0: alpha=%.2f lambda=%.4g df=%.1f BIC=%.3f\n",
            best$alpha, best$lam, best$df, best$bic))

# ---- 6. Iterated adaptive lasso refit (until convergence) ----------------
# Adaptive weights: w_g = pf_g[g] / ||beta_g||^gamma  (group)
#                   w_j = 1 / |beta_j|^gamma           (individual)

max_stages <- 5
gamma      <- 1
prev_cps   <- NULL
prev_coef  <- best$coef

for (stage in 1:max_stages) {
  apf_g <- numeric(n_groups)
  for (g in 1:n_groups) {
    idx <- which(group == g)
    ng  <- sqrt(sum(best$coef[idx]^2))
    apf_g[g] <- if (ng > 1e-8) pf_g[g] / ng^gamma else 1e6
  }
  apf_s <- ifelse(abs(best$coef) > 1e-8, 1 / abs(best$coef)^gamma, 1e6)

  best_new <- pick_best_bic(X, y_vec, group,
                            pf_group = apf_g, pf_sparse = apf_s,
                            df_cap = df_cap, alpha_grid = alpha_grid)
  delta_tmp <- matrix(best_new$coef[(K + 1):length(best_new$coef)],
                      nrow = b - 1, byrow = TRUE)
  jn_tmp <- apply(delta_tmp, 1, function(v) sqrt(sum(v^2)))
  cps_new <- age_keep[which(jn_tmp > 1e-6) + 1]

  cat(sprintf("Stage %d: alpha=%.2f lambda=%.4g df=%.1f BIC=%.3f | CPs (%d): %s\n",
              stage, best_new$alpha, best_new$lam, best_new$df, best_new$bic,
              length(cps_new), paste(cps_new, collapse = ",")))

  coef_change <- sqrt(sum((best_new$coef - best$coef)^2)) /
                 max(sqrt(sum(best$coef^2)), 1e-8)
  same_cps <- identical(sort(cps_new), sort(prev_cps))
  if (same_cps && coef_change < 1e-3) {
    cat(sprintf("  Converged (relative coef change %.4f).\n", coef_change))
    best <- best_new
    break
  }
  prev_cps <- cps_new
  best     <- best_new
}

# ---- 7. Recover fused curves & extract change points ----------------------

coef_final  <- best$coef
delta_final <- matrix(coef_final[(K + 1):length(coef_final)],
                      nrow = b - 1, byrow = TRUE)
colnames(delta_final) <- biomarkers

# Un-standardize coefficients
alpha_hat_orig <- coef_final[1:K] * col_sd
delta_orig     <- sweep(delta_final, 2, col_sd, "*")

B_fused <- matrix(NA_real_, nrow = b, ncol = K,
                  dimnames = list(age_keep, biomarkers))
B_fused[1, ] <- alpha_hat_orig
for (t in 2:b) {
  B_fused[t, ] <- alpha_hat_orig + colSums(delta_orig[1:(t - 1), , drop = FALSE])
}

jump_norms    <- apply(delta_final, 1, function(v) sqrt(sum(v^2)))
cp_pos        <- which(jump_norms > 1e-6)
change_points <- if (length(cp_pos)) age_keep[cp_pos + 1] else integer(0)

cat(sprintf("\nFinal alpha:  %.2f\n", best$alpha))
cat(sprintf("Final lambda: %g\n", best$lam))
cat(sprintf("CPs (%d): %s\n", length(change_points),
            paste(change_points, collapse = ", ")))

# Per-CP biomarker decomposition (standardized scale)
if (length(cp_pos) > 0) {
  cat("\nPer-CP biomarker contributions (standardized jumps):\n")
  for (i in seq_along(cp_pos)) {
    dv  <- delta_final[cp_pos[i], ]
    pct <- 100 * dv^2 / sum(dv^2)
    cat(sprintf("  age %d (||delta||=%.3f): %s\n",
                change_points[i], sqrt(sum(dv^2)),
                paste(sprintf("%s=%+.2f (%.0f%%)",
                              biomarkers, dv, pct), collapse = ", ")))
  }
}

# ---- 8. Cox model ----------------------------------------------------------

predictors_cat <- c("diabete", "gender", "race", "smoker",
                    "anti", "hearing_loss", "ethnicity",
                    "brain", "depression")

if (length(change_points) > 0) {
  data$age_regime <- cut(data$age, breaks = c(-Inf, change_points, Inf),
                         right = FALSE, include.lowest = TRUE)
} else {
  data$age_regime <- 1L
}

print(data %>% group_by(age_regime) %>%
      summarise(n = n(), n_dem = sum(Dementia), rate = round(mean(Dementia)*100, 2)))

form_cox <- as.formula(paste(
  "Surv(time2, Dementia) ~ factor(age_regime) + sbp + bmi + ldlc + hdl + trig +",
  paste0("factor(", predictors_cat, ")", collapse = " + ")
))
fit_cox <- suppressWarnings(coxph(form_cox, data = data, x = TRUE, y = TRUE))

# Test C-index
set.seed(42)
tr_idx <- sample(seq_len(nrow(data)), floor(0.60 * nrow(data)))
fit_tr <- suppressWarnings(coxph(form_cox, data = data[tr_idx, ]))
lp_te  <- predict(fit_tr, newdata = data[-tr_idx, ], type = "lp")
c_te   <- concordance(Surv(time2, Dementia) ~ lp_te,
                      data = data[-tr_idx, ], reverse = TRUE)$concordance

cat(sprintf("\nCox BIC: %.1f\n", BIC(fit_cox)))
cat(sprintf("Test C-index: %.4f\n", c_te))
