## step2_group_lasso.R
## Step 2: Group Lasso for shared change point detection
## Pipeline: adaptive Gaussian kernel + column standardization
##           + cv.gglasso (initial) + adaptive lasso refit (gamma=3)
##
## Variable-controlled comparison with Step 3 (only difference: no L1).

library(dplyr)
library(tidyverse)
library(broom)
library(survival)
library(gglasso)

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

# ---- 2. Per-age GLM (adaptive Gaussian kernel, same as Step 3) -------------

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

# Column standardization: prevent SCRE from dominating jump norms
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

# ---- 4. Stage 1: Initial group lasso (CV lambda.min) ----------------------

set.seed(123)
cvfit_init <- cv.gglasso(x = X, y = y_vec, group = group,
                         pf = pf_g,
                         pred.loss = "L2", nfolds = 10,
                         loss = "ls", intercept = FALSE)
lam_init <- cvfit_init$lambda.min

fit_init <- gglasso(x = X, y = y_vec, group = group,
                    pf = pf_g, loss = "ls",
                    lambda = lam_init, intercept = FALSE)
coef_init <- as.numeric(coef(fit_init))
if (length(coef_init) == ncol(X) + 1) coef_init <- coef_init[-1]

# ---- 5. Stage 2: Adaptive lasso refit (gamma = 3) -------------------------
# Adaptive weights: w_g = pf_g[g] / ||beta_init_g||^gamma
# Strong initial signals get small weights (preserved),
# weak signals get large weights (suppressed).

gamma <- 3
apf_g <- numeric(n_groups)
for (g in 1:n_groups) {
  idx <- which(group == g)
  ng  <- sqrt(sum(coef_init[idx]^2))
  apf_g[g] <- if (ng > 1e-8) pf_g[g] / ng^gamma else 1e6
}

set.seed(123)
cvfit_final <- cv.gglasso(x = X, y = y_vec, group = group,
                          pf = apf_g,
                          pred.loss = "L2", nfolds = 10,
                          loss = "ls", intercept = FALSE)
lam_final <- cvfit_final$lambda.min

fit_final <- gglasso(x = X, y = y_vec, group = group,
                     pf = apf_g, loss = "ls",
                     lambda = lam_final, intercept = FALSE)
coef_final <- as.numeric(coef(fit_final))
if (length(coef_final) == ncol(X) + 1) coef_final <- coef_final[-1]

# ---- 6. Recover fused curves & extract change points ----------------------

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

# ---- 7. Cox model ----------------------------------------------------------

predictors_cat <- c("diabete", "gender_clean", "race", "smoker",
                    "anti", "statin", "ethnicity")

if (length(change_points) > 0) {
  data$age_regime <- cut(data$age, breaks = c(-Inf, change_points, Inf),
                         right = FALSE, include.lowest = TRUE)
} else {
  data$age_regime <- 1L
}

print(data %>% group_by(age_regime) %>%
      summarise(n = n(), n_cvd = sum(cvd), rate = round(mean(cvd)*100, 2)))

form_cox <- as.formula(paste(
  "Surv(time2, cvd) ~ factor(age_regime) + sbp + bmi + hdl + tc + scre +",
  paste0("factor(", predictors_cat, ")", collapse = " + ")
))
fit_cox <- suppressWarnings(coxph(form_cox, data = data, x = TRUE, y = TRUE))

# Test C-index
set.seed(42)
tr_idx <- sample(seq_len(nrow(data)), floor(0.60 * nrow(data)))
fit_tr <- suppressWarnings(coxph(form_cox, data = data[tr_idx, ]))
lp_te  <- predict(fit_tr, newdata = data[-tr_idx, ], type = "lp")
c_te   <- concordance(Surv(time2, cvd) ~ lp_te,
                      data = data[-tr_idx, ], reverse = TRUE)$concordance

cat(sprintf("\nCox BIC: %.1f\n", BIC(fit_cox)))
cat(sprintf("Test C-index: %.4f\n", c_te))
