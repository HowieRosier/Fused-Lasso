## step2_group_lasso.R
## Step 2: Group Lasso for shared change point detection
## Uses CV lambda.min + jump norm filtering + adjacent merging

library(dplyr)
library(tidyverse)
library(broom)
library(gglasso)
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

# 3. Design matrix

wide_est <- perage %>%
  select(age_grid, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  arrange(age_grid)

Y_all <- as.matrix(wide_est[, biomarkers])
keep_row <- apply(Y_all, 1, function(r) all(is.finite(r)))

age_keep <- wide_est$age_grid[keep_row]
Y <- Y_all[keep_row, ]
b <- nrow(Y); K <- ncol(Y)

# Reparameterize: beta_t = alpha + sum of delta_1..delta_{t-1}
y_vec   <- as.numeric(t(Y))
X_alpha <- kronecker(rep(1, b), diag(K))
J <- matrix(0, b, b - 1)
for (t in 2:b) J[t, 1:(t - 1)] <- 1
X_delta <- kronecker(J, diag(K))
X       <- cbind(X_alpha, X_delta)
group   <- c(1:K, rep((K + 1):(K + b - 1), each = K))

# 4. Group lasso

# CV to get lambda path
set.seed(123)
cvfit <- cv.gglasso(x = X, y = y_vec, group = group,
                    pred.loss = "L2", nfolds = 10, loss = "ls",
                    intercept = FALSE)

# Use lambda.min, then filter by jump norm
fit_cv  <- gglasso(x = X, y = y_vec, group = group,
                   loss = "ls", lambda = cvfit$lambda.min, intercept = FALSE)
coef_cv <- as.numeric(coef(fit_cv))
if (length(coef_cv) == ncol(X) + 1) coef_cv <- coef_cv[-1]

delta_cv  <- matrix(coef_cv[(K+1):length(coef_cv)], nrow = b-1, byrow = TRUE)
jump_norm <- apply(delta_cv, 1, function(v) sqrt(sum(v^2)))

# Keep top 25% by jump norm magnitude
threshold  <- quantile(jump_norm[jump_norm > 1e-6], 0.75)
cp_raw     <- age_keep[which(jump_norm >= threshold) + 1]

# Merge adjacent change points (within 3 years)
merge_adjacent <- function(ages, min_gap = 3) {
  ages <- sort(ages)
  groups <- list(); g <- c(ages[1])
  for (i in 2:length(ages)) {
    if (ages[i] - ages[i-1] < min_gap) g <- c(g, ages[i])
    else { groups <- c(groups, list(g)); g <- c(ages[i]) }
  }
  groups <- c(groups, list(g))
  sapply(groups, mean)
}

cp_merged <- merge_adjacent(cp_raw, min_gap = 3)

cat("Raw CPs (top 25% jump norm):", paste(cp_raw, collapse = ", "), "\n")
cat("Merged CPs:", round(cp_merged, 1), "\n")

# 5. Cox model

predictors_cat <- c("diabete", "gender_clean", "race", "smoker",
                     "anti", "statin", "ethnicity")

data$age_regime <- cut(data$age,
                       breaks = c(-Inf, round(cp_merged, 1), Inf),
                       right = FALSE, include.lowest = TRUE)
print(table(data$age_regime))

form_cox <- as.formula(paste(
  "Surv(time2, cvd) ~ age_regime + sbp + bmi + hdl + tc + scre +",
  paste0("factor(", predictors_cat, ")", collapse = " + ")
))
fit_cox <- coxph(form_cox, data = data, x = TRUE, y = TRUE)

cat("Cox BIC:", BIC(fit_cox), "\n")
