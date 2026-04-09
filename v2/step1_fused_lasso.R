## step1_fused_lasso.R
## Step 1: Per-variable Fused Lasso for age-related change point detection
## Uses modified BIC (mult=3.4) to select lambda + stability analysis

library(dplyr)
library(tidyverse)
library(broom)
library(genlasso)
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

no_info_vals <- c(
  "PMI: Skip", "What Race Ethnicity: Race Ethnicity None Of These",
  "PMI: Prefer Not To Answer")
data$ethnicity[data$ethnicity %in% no_info_vals] <- "No information"

no_info_vals2 <- c(
  "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander",
  "I prefer not to answer", "None Indicated", "None of these",
  "PMI: Skip", "Middle Eastern or North African",
  "More than one population")
data$race[data$race %in% no_info_vals2] <- "Others"

data$gender_clean <- case_when(
  data$gender == "Female" ~ "Female",
  data$gender == "Male"   ~ "Male",
  TRUE                    ~ "Other"
)

set.seed(42)
train_idx  <- sample(seq_len(nrow(data)), size = floor(0.60 * nrow(data)))
train_data <- data[train_idx, ]
test_data  <- data[-train_idx, ]

# 2. Per-age GLM (fixed ±2 window)

biomarkers <- c("sbp", "bmi", "hdl", "tc", "scre")
age_grid   <- 30:79
win        <- 2

fit_one_window <- function(df, a_g, win, targets) {
  dwin <- df %>% filter(is.finite(age), abs(age - a_g) <= win)
  if (nrow(dwin) < 50)
    return(tibble(age_grid = a_g, term = targets, estimate = NA_real_, se = NA_real_))
  form <- as.formula(paste("cvd ~", paste(targets, collapse = " + ")))
  mod  <- try(glm(form, data = dwin, family = binomial()), silent = TRUE)
  if (inherits(mod, "try-error"))
    return(tibble(age_grid = a_g, term = targets, estimate = NA_real_, se = NA_real_))
  co <- broom::tidy(mod) %>%
    filter(term %in% targets) %>%
    select(term, estimate, std.error) %>%
    rename(se = std.error)
  tibble(age_grid = a_g, term = targets) %>% left_join(co, by = "term")
}

perage <- map_dfr(age_grid,
  ~ fit_one_window(train_data, .x, win, biomarkers)
) %>% arrange(term, age_grid)

# 3. Per-variable fused lasso with modified BIC

# For each biomarker separately: fit fused lasso on the per-age coefficient
# curve, select lambda by modified BIC with adjustable multiplier.
# Standard BIC (mult=1) is too weak because inverse-variance weights are large
# (thousands of patients per window -> small SE -> large 1/SE^2).

fuse_one_predictor_mBIC <- function(pred, perage, mult = 3.4) {
  dat  <- perage %>% filter(term == pred) %>% arrange(age_grid)
  y    <- dat$estimate; se <- dat$se
  keep <- is.finite(y) & is.finite(se) & se > 0
  if (sum(keep) < 10)
    return(tibble(term = pred, age_grid = dat$age_grid,
                  beta_raw = y, beta_fused = NA_real_, cp = FALSE))

  yk <- y[keep]; sek <- se[keep]; wk <- 1 / (sek^2)
  X <- diag(length(yk)); D <- diff(diag(length(yk)))
  Xw <- diag(sqrt(wk)) %*% X; yw <- as.numeric(sqrt(wk) * yk)

  fit <- fusedlasso(y = yw, X = Xw, D = D)

  # Modified BIC: n*log(WRSS/n) + mult*df*log(n)
  lam   <- fit$lambda
  coefs <- coef(fit, lambda = lam)$beta
  n     <- length(yk)
  wrss  <- apply(coefs, 2, function(b) sum(wk * (yk - b)^2))
  jumps_count <- apply(apply(coefs, 2, diff), 2, function(d) sum(abs(d) > 1e-6))
  df       <- jumps_count + 1
  IC_vals  <- n * log(wrss / n) + mult * df * log(n)
  lam_star <- lam[which.min(IC_vals)]
  beta_hat <- as.numeric(coef(fit, lambda = lam_star)$beta)

  out <- tibble(term = pred, age_grid = dat$age_grid,
                beta_raw = y, beta_fused = NA_real_, cp = FALSE, lambda = lam_star)
  out$beta_fused[keep] <- beta_hat
  jumps <- which(abs(diff(beta_hat)) > 1e-6)
  if (length(jumps) > 0) {
    kept_ages <- out$age_grid[keep]
    cp_age    <- kept_ages[jumps + 1]
    out$cp[out$age_grid %in% cp_age] <- TRUE
  }
  out
}

# 4. Run + sensitivity analysis

# Main result at mult=3.4
res_final <- map_dfr(biomarkers, ~ fuse_one_predictor_mBIC(.x, perage, mult = 3.4))
cp_summary <- res_final %>% filter(cp) %>%
  group_by(term) %>%
  summarise(n_cp = n(), ages = paste(age_grid, collapse = ", "))
cat("=== Step 1 Result (mult=3.4) ===\n")
print(cp_summary)
cat("Total CPs:", sum(cp_summary$n_cp), "\n")

# Sensitivity: sweep mult from 3.0 to 4.0
cat("\n=== Sensitivity Analysis ===\n")
for (m in c(3.0, 3.2, 3.4, 3.6, 3.8, 4.0)) {
  res <- map_dfr(biomarkers, ~ fuse_one_predictor_mBIC(.x, perage, mult = m))
  cp_sum <- res %>% filter(cp) %>%
    group_by(term) %>%
    summarise(n_cp = n(), ages = paste(age_grid, collapse = ", "))
  cat("mult =", m, "| Total CPs:", sum(cp_sum$n_cp), "\n")
}

# 5. Plot

library(ggplot2)
ggplot(res_final, aes(x = age_grid)) +
  geom_point(aes(y = beta_raw), alpha = 0.4, size = 1, color = "blue") +
  geom_step(aes(y = beta_fused), linewidth = 0.8, color = "red") +
  geom_vline(data = res_final %>% filter(cp),
             aes(xintercept = age_grid),
             linetype = "dashed", color = "red", alpha = 0.5) +
  facet_wrap(~term, scales = "free_y", ncol = 2) +
  labs(title = "Fused Lasso (mBIC, mult=3.4)", x = "Age", y = "Coefficient") +
  theme_minimal()
