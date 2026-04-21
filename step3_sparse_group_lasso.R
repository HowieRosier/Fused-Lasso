## step3_sparse_group_lasso_v2.R
## Step 3: Sparse Group Lasso with adaptive kernel + column standardization
##
## Changes from v1:
##   1. Adaptive Gaussian kernel (instead of fixed ±2yr window)
##   2. Column standardization of biomarker coefficients before SGL
##   3. Alpha selected by BIC grid search (not manually fixed)
##   4. Weighted CP merge by jump norm (not simple mean)

library(dplyr)
library(tidyverse)
library(broom)
library(sparsegl)
library(survival)

# 1. Data

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

# 2. Per-age GLM (adaptive kernel)

biomarkers <- c("sbp", "bmi", "hdl", "tc", "scre")
age_grid   <- 30:79

# Gaussian kernel GLM at each age. Bandwidth h expands until
# the window has >= min_events CVD events.
fit_one_adaptive_kernel <- function(df, a_g, targets, min_events = 100) {
  df <- df %>% filter(is.finite(age))
  for (h in seq(1, 10, by = 0.5)) {
    w    <- exp(-((df$age - a_g)^2) / (2 * h^2))
    keep <- w > 0.01 * max(w)
    if (sum(df$cvd[keep]) >= min_events) break
  }
  dwin <- df[keep, ]; wk <- w[keep]
  if (nrow(dwin) < 50)
    return(tibble(age_grid = a_g, term = targets, bandwidth = h,
                  estimate = NA_real_, se = NA_real_))

  form <- as.formula(paste("cvd ~", paste(targets, collapse = " + ")))
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
  ~ fit_one_adaptive_kernel(train_data, .x, biomarkers, min_events = 100)
) %>% arrange(term, age_grid)

# 3. Design matrix

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

# [NEW] Column standardization: prevent SCRE from dominating group norms.
# SCRE coefficients are ~17x larger than other biomarkers; without this,
# SGL detects only SCRE-specific CPs rather than consensus CPs.
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

# 4. Sparse group lasso

# SE-based group penalty (on standardized scale)
SE_by_age   <- apply(SE_std, 1, mean)
n_groups    <- max(group)
pf_group_se <- numeric(n_groups)
pf_group_se[1:K] <- 1
for (m in 1:(b - 1)) pf_group_se[K + m] <- 1 / SE_by_age[m]
pf_group_se[(K+1):n_groups] <- pf_group_se[(K+1):n_groups] /
  mean(pf_group_se[(K+1):n_groups])

# [NEW] BIC grid search over alpha (data-driven, not manually fixed)
best_bic <- Inf; best_alpha <- NA; best_lam <- NA
for (a in seq(0.05, 0.95, by = 0.05)) {
  fit <- sparsegl(x = X, y = y_vec, group = group,
                  asparse = a, pf_group = pf_group_se, intercept = FALSE)
  er <- estimate_risk(fit, X, type = "BIC", approx_df = TRUE)
  if (min(er$BIC) < best_bic) {
    best_bic   <- min(er$BIC)
    best_alpha <- a
    best_lam   <- er$lambda[which.min(er$BIC)]
  }
}

# Adaptive lasso: initial fit -> derive weights -> refit
fit_init  <- sparsegl(x = X, y = y_vec, group = group,
                      asparse = best_alpha, pf_group = pf_group_se,
                      intercept = FALSE)
coef_init <- as.numeric(coef(fit_init, s = best_lam))
if (length(coef_init) == ncol(X) + 1) coef_init <- coef_init[-1]

apf_group <- numeric(n_groups)
for (g in 1:n_groups) {
  idx <- which(group == g)
  ng  <- sqrt(sum(coef_init[idx]^2))
  apf_group[g] <- if (ng > 1e-8) pf_group_se[g] / ng else 1e6
}
apf_sparse <- ifelse(abs(coef_init) > 1e-8, 1 / abs(coef_init), 1e6)

fit_final <- sparsegl(x = X, y = y_vec, group = group,
                      asparse = best_alpha,
                      pf_group = apf_group, pf_sparse = apf_sparse,
                      intercept = FALSE)
er_final  <- estimate_risk(fit_final, X, type = "BIC", approx_df = TRUE)
lam_final <- er_final$lambda[which.min(er_final$BIC)]

# Extract change points
coef_final  <- as.numeric(coef(fit_final, s = lam_final))
if (length(coef_final) == ncol(X) + 1) coef_final <- coef_final[-1]
delta_final <- matrix(coef_final[(K+1):length(coef_final)], nrow = b-1, byrow = TRUE)
jump_norms  <- apply(delta_final, 1, function(v) sqrt(sum(v^2)))
cp_pos      <- which(jump_norms > 1e-6)
change_points <- age_keep[cp_pos + 1]
cp_jn         <- jump_norms[cp_pos]

# [NEW] Reconstruct coefficients in original scale
alpha_hat_std  <- coef_final[1:K]
alpha_hat_orig <- alpha_hat_std * col_sd
delta_orig     <- sweep(delta_final, 2, col_sd, "*")

B_fused <- matrix(NA, nrow = b, ncol = K,
                  dimnames = list(age_keep, biomarkers))
B_fused[1, ] <- alpha_hat_orig
for (t in 2:b) {
  B_fused[t, ] <- alpha_hat_orig + colSums(delta_orig[1:(t-1), , drop = FALSE])
}

cat("Alpha:", best_alpha, "| Lambda:", lam_final, "\n")
cat("Raw CPs:", paste(change_points, collapse = ", "), "\n")
cat("Jump norms:", paste(round(cp_jn, 4), collapse = ", "), "\n")

# 5. Weighted merge

# [NEW] Merge adjacent CPs weighted by jump norm (not simple mean).
# Stronger jumps pull the representative age toward them.
merge_weighted <- function(ages, jump_norms, min_gap = 3) {
  if (length(ages) <= 1) return(ages)
  ord <- order(ages); ages <- ages[ord]; jn <- jump_norms[ord]
  groups <- list(); g_ages <- c(ages[1]); g_jn <- c(jn[1])
  for (i in 2:length(ages)) {
    if (ages[i] - ages[i-1] < min_gap) {
      g_ages <- c(g_ages, ages[i]); g_jn <- c(g_jn, jn[i])
    } else {
      groups <- c(groups, list(list(ages = g_ages, jn = g_jn)))
      g_ages <- c(ages[i]); g_jn <- c(jn[i])
    }
  }
  groups <- c(groups, list(list(ages = g_ages, jn = g_jn)))
  sapply(groups, function(g) round(weighted.mean(g$ages, g$jn)))
}

cp_merged <- merge_weighted(change_points, cp_jn, min_gap = 3)
cat("Merged CPs:", paste(cp_merged, collapse = ", "), "\n")

# 6. Cox model

predictors_cat <- c("diabete", "gender_clean", "race", "smoker",
                    "anti", "statin", "ethnicity")

data$age_regime <- cut(data$age, breaks = c(-Inf, cp_merged, Inf),
                       right = FALSE, include.lowest = TRUE)

form_cox <- as.formula(paste(
  "Surv(time2, cvd) ~ age_regime + sbp + bmi + hdl + tc + scre +",
  paste0("factor(", predictors_cat, ")", collapse = " + ")
))

# Full data: BIC
fit_cox <- coxph(form_cox, data = data, x = TRUE, y = TRUE)

# Test data: C-index
set.seed(42)
train_idx  <- sample(seq_len(nrow(data)), size = floor(0.60 * nrow(data)))
train_data <- data[train_idx, ]
test_data  <- data[-train_idx, ]
fit_train  <- coxph(form_cox, data = train_data)
lp_test    <- predict(fit_train, newdata = test_data, type = "lp")
c_test     <- concordance(Surv(time2, cvd) ~ lp_test,
                          data = test_data, reverse = TRUE)$concordance

cat("\n Final Results \n")
cat(sprintf("Alpha (BIC-selected): %.2f\n", best_alpha))
cat(sprintf("Lambda: %g\n", lam_final))
cat(sprintf("Raw CPs (%d): %s\n", length(change_points),
            paste(change_points, collapse = ", ")))
cat(sprintf("Merged CPs (%d): %s\n", length(cp_merged),
            paste(cp_merged, collapse = ", ")))
cat(sprintf("Cox BIC: %.1f\n", BIC(fit_cox)))
cat(sprintf("Test C-index: %.4f\n", c_test))

cat("\nAge regimes:\n")
print(data %>% group_by(age_regime) %>%
      summarise(n = n(), n_cvd = sum(cvd),
                rate = round(mean(cvd) * 100, 2)))
