## step2_group_lasso.R
## Step 2: Group Lasso for shared change point detection
## Pipeline: adaptive Gaussian kernel + column standardization
##           + BIC-selected lambda (matched with Step 3 SGL framework)
##           + adaptive lasso refit (gamma=1)
##
## Variable-controlled comparison with Step 3 (only difference: no L1).
##
## Cohort: All of Us dementia relax cohort
## Outcome: incident dementia
## Biomarkers: SBP, BMI, LDL-C, HDL-C, Triglycerides

library(dplyr)
library(tidyverse)
library(broom)
library(survival)
library(gglasso)

# ---- 1. Data ----------------------------------------------------------------

system("mkdir -p ~/data")
system(paste0("gsutil cp ", Sys.getenv("WORKSPACE_BUCKET"),
              "/dementia_relax_cohort_time_all_final1.csv ~/data/"))

data <- read.csv("~/data/dementia_relax_cohort_time_all_final1.csv")
data <- data %>% filter(is.finite(time2), time2 > 0, is.finite(age))

# Filter to SAP main analysis age range
data <- data %>% filter(age >= 45, age < 80)

set.seed(42)
train_idx  <- sample(seq_len(nrow(data)), size = floor(0.60 * nrow(data)))
train_data <- data[train_idx, ]
test_data  <- data[-train_idx, ]

# ---- 2. Per-age GLM (adaptive Gaussian kernel, same as Step 3) -------------

biomarkers <- c("sbp", "bmi", "ldlc", "hdl", "trig")
age_grid   <- 45:79

# Standardize biomarkers (mean=0, SD=1) — cohort uses raw clinical-scale values
for (bm in biomarkers) {
  train_data[[bm]] <- as.numeric(scale(train_data[[bm]]))
  data[[bm]]       <- as.numeric(scale(data[[bm]]))
}

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

perage <- map_dfr(age_grid,
  ~ fit_one_adaptive_kernel(train_data, .x, biomarkers, min_events = 100)
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

# Column standardization: prevent any biomarker from dominating jump norms
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

# SE-based penalty factor
SE_by_age <- apply(SE_std, 1, mean)
n_groups  <- max(group)
pf_g <- numeric(n_groups); pf_g[1:K] <- 1
for (m in 1:(b - 1)) pf_g[K + m] <- 1 / SE_by_age[m]
pf_g[(K+1):n_groups] <- pf_g[(K+1):n_groups] /
  mean(pf_g[(K+1):n_groups])

# ---- 4. BIC computation with exact df (matches Step 3 sparsegl) ----------
# Exact df for group lasso: total parameters in active groups

compute_bic_gglasso <- function(fit, X, y, group) {
  n  <- length(y)
  betas <- as.matrix(fit$beta)
  preds <- predict(fit, newx = X, type = "link")
  if (is.null(dim(preds))) preds <- matrix(preds, ncol = 1)
  rss <- colSums((y - preds)^2)
  df_exact <- sapply(seq_along(fit$lambda), function(j) {
    beta_j <- betas[, j]
    active_groups <- unique(group[abs(beta_j) > 1e-8])
    if (length(active_groups) == 0) return(0)
    sum(group %in% active_groups & abs(beta_j) > 1e-8)
  })
  list(bic = log(rss / n) + log(n) * df_exact / n,
       df  = df_exact, lambda = fit$lambda, rss = rss)
}

# ---- 5. Stage 1: Initial group lasso (BIC + df_cap) ----------------------

df_cap <- b - 1   # prevent saturated overfitting at n=p edge

fit_init <- gglasso(x = X, y = y_vec, group = group,
                    pf = pf_g, loss = "ls",
                    intercept = FALSE,
                    nlambda = 200,
                    lambda.factor = 1e-5)
er_init <- compute_bic_gglasso(fit_init, X, y_vec, group)
valid_init <- er_init$df <= df_cap
j_init  <- which.min(er_init$bic[valid_init])
lam_init <- er_init$lambda[valid_init][j_init]
coef_init <- as.numeric(coef(fit_init, s = lam_init))
if (length(coef_init) == ncol(X) + 1) coef_init <- coef_init[-1]

# ---- 6. Stage 2: Adaptive lasso refit (gamma = 1) ------------------------
# Adaptive weights: w_g = pf_g[g] / ||beta_init_g||^gamma
# Strong initial signals get small weights (preserved),
# weak signals get large weights (suppressed).

gamma <- 1
apf_g <- numeric(n_groups)
for (g in 1:n_groups) {
  idx <- which(group == g)
  ng  <- sqrt(sum(coef_init[idx]^2))
  apf_g[g] <- if (ng > 1e-8) pf_g[g] / ng^gamma else 1e6
}

fit_final <- gglasso(x = X, y = y_vec, group = group,
                     pf = apf_g, loss = "ls",
                     intercept = FALSE,
                     nlambda = 200,
                     lambda.factor = 1e-5)
er_final <- compute_bic_gglasso(fit_final, X, y_vec, group)
valid_final <- er_final$df <= df_cap
j_final  <- which.min(er_final$bic[valid_final])
lam_final <- er_final$lambda[valid_final][j_final]
coef_final <- as.numeric(coef(fit_final, s = lam_final))
if (length(coef_final) == ncol(X) + 1) coef_final <- coef_final[-1]

# ---- 7. Recover fused curves & extract change points ----------------------

alpha_hat <- coef_final[1:K]
delta_hat <- matrix(coef_final[(K + 1):length(coef_final)],
                    nrow = b - 1, byrow = TRUE)
colnames(delta_hat) <- biomarkers

# Un-standardize coefficients
alpha_hat_orig <- alpha_hat * col_sd
delta_orig     <- sweep(delta_hat, 2, col_sd, "*")

B_fused <- matrix(NA_real_, nrow = b, ncol = K,
                  dimnames = list(age_keep, biomarkers))
B_fused[1, ] <- alpha_hat_orig
for (t in 2:b) {
  B_fused[t, ] <- alpha_hat_orig + colSums(delta_orig[1:(t - 1), , drop = FALSE])
}

jump_norms <- apply(delta_hat, 1, function(v) sqrt(sum(v^2)))
cp_pos     <- which(jump_norms > 1e-6)
change_points <- if (length(cp_pos)) age_keep[cp_pos + 1] else integer(0)

cat(sprintf("Stage 1 lambda: %g\n", lam_init))
cat(sprintf("Stage 2 lambda: %g (gamma=%d)\n", lam_final, gamma))
cat(sprintf("CPs (%d): %s\n", length(change_points),
            paste(change_points, collapse = ", ")))

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
