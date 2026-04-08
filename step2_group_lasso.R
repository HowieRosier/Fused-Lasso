# Step 2: Group Lasso (shared change points across biomarkers)
#
# Reparameterize the per-age coefficient curves as
#     beta_t = alpha + sum_{m < t} delta_m
# where each delta_m is a K-vector (one entry per biomarker). A
# shared change point at boundary m corresponds to delta_m != 0.
# We solve it via cv.gglasso with groups = (baseline cols, each delta_m).

library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(broom)
library(Matrix)
library(survival)
library(gglasso)

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
win      <- 2

perage <- map_dfr(
  age_grid,
  ~ fit_one_window(train_data, a_g = .x, win = win,
                   rhs_full = rhs_full, targets = biomarkers,
                   age_var = age_var, outcome = OUTCOME)
) %>%
  arrange(term, age_grid)

# ---------- wide coefficient matrix (b age points x K biomarkers) ----------

wide_est <- perage %>%
  select(age_grid, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  arrange(age_grid)

wide_se <- perage %>%
  select(age_grid, term, se) %>%
  pivot_wider(names_from = term, values_from = se) %>%
  arrange(age_grid)

Y_all  <- as.matrix(wide_est[, biomarkers, drop = FALSE])
SE_all <- as.matrix(wide_se[, biomarkers, drop = FALSE])

keep_row <- apply(Y_all,  1, function(r) all(is.finite(r))) &
            apply(SE_all, 1, function(r) all(is.finite(r) & r > 0))
if (sum(keep_row) < 15) stop("Too few complete age grid points.")

age_keep <- wide_est$age_grid[keep_row]
Y <- Y_all[keep_row, , drop = FALSE]
b <- nrow(Y)
K <- ncol(Y)

# ---------- design matrix: baseline + jump groups ----------
# y_vec = vec(Y^T) stacks (Y[1, ], Y[2, ], ..., Y[b, ])
# X_alpha picks the baseline alpha (K cols)
# X_delta picks the cumulative sum of jumps delta_1..delta_{b-1}

y_vec   <- as.numeric(t(Y))
X_alpha <- kronecker(rep(1, b), diag(K))

J <- matrix(0, nrow = b, ncol = b - 1)
for (t in 2:b) J[t, 1:(t - 1)] <- 1
X_delta <- kronecker(J, diag(K))

X <- cbind(X_alpha, X_delta)

# groups: baseline cols are singletons, each delta_m is one group of K
group_alpha <- 1:K
group_delta <- rep((K + 1):(K + b - 1), each = K)
group <- c(group_alpha, group_delta)

# ---------- CV + final fit ----------

set.seed(123)
cvfit <- cv.gglasso(
  x = X, y = y_vec, group = group,
  pred.loss = "L2", nfolds = 10,
  loss = "ls", intercept = FALSE
)
lam_star <- cvfit$lambda.1se
print(lam_star)

fit_final <- gglasso(
  x = X, y = y_vec, group = group,
  loss = "ls", lambda = lam_star, intercept = FALSE
)

coef_hat <- as.numeric(coef(fit_final))
if (length(coef_hat) == ncol(X) + 1) coef_hat <- coef_hat[-1]
stopifnot(length(coef_hat) == ncol(X))

# ---------- recover fused curves & shared change points ----------

alpha_hat <- coef_hat[1:K]
delta_hat <- matrix(coef_hat[(K + 1):length(coef_hat)],
                    nrow = b - 1, byrow = TRUE)
colnames(delta_hat) <- biomarkers
rownames(delta_hat) <- paste0("jump_", age_keep[-length(age_keep)])

B_fused <- matrix(NA_real_, nrow = b, ncol = K,
                  dimnames = list(age_keep, biomarkers))
B_fused[1, ] <- alpha_hat
for (t in 2:b) {
  B_fused[t, ] <- alpha_hat + colSums(delta_hat[1:(t - 1), , drop = FALSE])
}

jump_norm <- apply(delta_hat, 1, function(v) sqrt(sum(v^2)))
jump_pos  <- which(jump_norm > 1e-6)
cp_age     <- if (length(jump_pos)) age_keep[jump_pos + 1] else numeric(0)
union_cuts <- sort(unique(cp_age))

cp_table_shared <- tibble(
  age_grid  = cp_age,
  cp_after  = if (length(jump_pos))
                paste0(age_keep[jump_pos], " -> ", age_keep[jump_pos + 1])
              else character(0),
  jump_norm = jump_norm[jump_pos],
  lambda    = lam_star
)
print(cp_table_shared)

beta_fused_long <- as.data.frame(B_fused) %>%
  tibble::rownames_to_column("age_grid") %>%
  mutate(age_grid = as.numeric(age_grid)) %>%
  pivot_longer(-age_grid, names_to = "term", values_to = "beta_fused")

dir.create("result", showWarnings = FALSE)
write.csv(cp_table_shared,
          "result/change_point_age_cvgglasso_shared_1se.csv", row.names = FALSE)
write.csv(beta_fused_long,
          "result/beta_fused_cvgglasso_shared_1se.csv", row.names = FALSE)
print(union_cuts)

# ---------- Cox on age_refined ----------

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

lp_test <- predict(fit_cox, newdata = test_data, type = "lp")
c_test  <- survConcordance(Surv(time2, cvd) ~ lp_test, data = test_data)$concordance
cat("FLS group C-index:", round(c_test, 6), "\n")

cat("logLik:", logLik(fit_cox), "\n")
cat("AIC:   ", AIC(fit_cox),    "\n")
cat("BIC:   ", BIC(fit_cox),    "\n")
