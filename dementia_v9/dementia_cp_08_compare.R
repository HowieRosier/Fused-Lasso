## =============================================================================
## dementia_cp_08_compare.R
## All change-point methods on the same cohort file, in one run, with one comparison table at the end.
##
##   A. The paper's three lines as written in github.com/HowieRosier/Fused-Lasso (dementia/):
##        line 1  step2_group_lasso.R         group lasso, all biomarkers share the change points   (package gglasso)
##        line 2  step1_fused_lasso.R         fused lasso, each biomarker on its own                (package flam)
##        line 3  step3_sparse_group_lasso.R  sparse group lasso                                    (package sparsegl)
##      Sections 2-7 of step2 and 2-3 of step1 are copied unchanged below (kernel per-age models on the 60% training
##      sample with seed 42, as in those scripts); line 3 is run through dementia_cp_07_tune.R with METHOD = "paper",
##      which reproduces dementia_cp_03_sgl_v9.R. The Cox models of the three scripts are not run here.
##
##   B. THE COMPARISON FOR THE PAPER: the three lines on ONE table of per-age estimates (the one line 3 uses: separate one-year
##      age groups, adjusted as set by ADJUST), so that only the method differs:
##        line 1  standard group lasso   change-point code of step2_group_lasso.R (sections 3-7), unchanged
##        line 2  standard fused lasso   change-point code of step1_fused_lasso.R (section 3), unchanged
##        line 3  combined, fine-tuned   sparse group lasso re-weighted step by step (see below)
##      SUPPLEMENT (printed after the comparison): the same three lines and the exact search on that table (separate one-year age groups,
##      every estimate weighted by its standard error, starting level free), each judged the same way:
##      a set of jumps is refitted without shrinkage and scored by misfit + BIC x log(n) per jump.
##        line 1  group lasso          every biomarker jumps at the same ages
##        line 2  fused lasso          each biomarker on its own, with its own penalty
##        line 3  sparse group lasso   shared and individual jumps mixed
##      The three lines are also run FINE-TUNED IN THE SAME WAY, so that they differ only in the penalty (the fair comparison):
##        line 1 fine-tuned            group penalty only (mixing weight 0): a change age moves every biomarker;
##                                     criterion misfit + log(n) x 1.5 x BIC_STRICT per change age
##        line 2 fine-tuned            each biomarker on its own (single-jump penalty, mixing weight 1);
##                                     criterion misfit + log(n) x BIC_STRICT per jump
##        line 3 fine-tuned            the sparse group lasso, re-weighted step by step: after each stage every jump gets a new
##                                     penalty weight 1 / (|jump| + eps) and every age 1 / (size of its jumps + eps), eps halved
##                                     every stage; each stage picks the mixing weight and the penalty strength by the combined
##                                     criterion below. Every jump column is rescaled so that a jump is judged by its evidence.
##        line 3 rebuilt (combined)    the sparse group lasso, with its mixing weight running from 0 (= line 1) to 1 (= line 2),
##                                     proposes the candidate ages. A change age then costs 1.5 x log(n) once, plus 1 x log(n) for
##                                     each biomarker that jumps there. A biomarker that jumps on its own pays 2.5 x log(n), as in
##                                     line 2 and the exact search; biomarkers that jump together share the price of the age.
##        exact search                 every possible set of jump ages, each biomarker on its own (script 07's main method)
##      The estimates come from dementia_cp_07_tune.R, so ADJUST and MEASURED_ONLY work here exactly as there.
##
## HOW TO RUN (dementia_cp_07_tune.R, version 7.1 or later, must be in the home folder):
##     source("~/dementia_cp_08_compare.R")
##     ADJUST <- c("sex_at_birth", "anti", "imp_bmi", "imp_ldlc", "imp_hdl", "imp_trig"); source("~/dementia_cp_08_compare.R")
## Wider panel (after dementia_cp_11_extra_labs.R): BIOMARKERS and BIN_MIN_EVENTS go to part B, e.g. with the panel and multiplier that
## dementia_cp_12_panel.R chose:  BIOMARKERS <- c("sbp", "bmi", "ldlc", "hdl", "trig", "dbp", ...); BIN_MIN_EVENTS <- 100; BIC_STRICT <- 2.5;
##     ADJUST <- c(... , "imp_dbp", ...); RUN_ORIGINAL <- FALSE; source("~/dementia_cp_08_compare.R")   (part A always uses the paper's five)
## Other settings: IN_CSV (cohort file; default = newest ~/dementia_v9/dementia_v9_analysis_*.csv), BIC_STRICT (2.5),
## RUN_ORIGINAL (TRUE), RUN_SAME (TRUE), ORIGINAL_LINES (1:3; e.g. c(1, 2) skips the slow line 3 of part A), TUNE_FILE ("~/dementia_cp_07_tune.R").
## Output: change-point ages only. No participant-level data and no counts of people are printed or saved.
## The table is also saved as ~/dementia_v9/compare_<date-time>.txt.
## =============================================================================

.c8 <- list(IN_CSV = if (exists("IN_CSV")) IN_CSV else "", ADJUST = if (exists("ADJUST")) as.character(ADJUST) else character(0),
            MEASURED_ONLY = if (exists("MEASURED_ONLY")) as.character(MEASURED_ONLY) else character(0),
            BIOMARKERS = if (exists("BIOMARKERS")) as.character(BIOMARKERS) else character(0),
            BIN_MIN_EVENTS = if (exists("BIN_MIN_EVENTS")) as.numeric(BIN_MIN_EVENTS) else numeric(0),
            BIC_STRICT = if (exists("BIC_STRICT")) as.numeric(BIC_STRICT) else 2.5, RUN_ORIGINAL = if (exists("RUN_ORIGINAL")) isTRUE(RUN_ORIGINAL) else TRUE,
            RUN_SAME = if (exists("RUN_SAME")) isTRUE(RUN_SAME) else TRUE,
            ORIGINAL_LINES = if (exists("ORIGINAL_LINES")) as.integer(ORIGINAL_LINES) else 1:3, TUNE_FILE = path.expand(if (exists("TUNE_FILE")) TUNE_FILE else "~/dementia_cp_07_tune.R"))
suppressWarnings(rm(list = intersect(c("ADJUST", "MEASURED_ONLY", "BIC_STRICT", "RUN_ORIGINAL", "RUN_SAME", "ORIGINAL_LINES", "BIOMARKERS", "BIN_MIN_EVENTS"),
                                     ls(envir = globalenv())), envir = globalenv()))
.c8$t0 <- Sys.time(); options(width = 170)
if (!file.exists(.c8$TUNE_FILE)) stop("dementia_cp_07_tune.R (version 7.1 or later) was not found at ", .c8$TUNE_FILE, call. = FALSE)
if (!any(grepl('SCRIPT_VERSION <- "(7\\.[1-9]|[89]|[1-9][0-9])', readLines(.c8$TUNE_FILE, warn = FALSE))))
  stop("dementia_cp_07_tune.R is older than version 7.1: upload the 7.1 file first (its METHOD = \"paper\" branch is needed here).", call. = FALSE)
for (p in c("dplyr", "tidyr", "purrr", "tibble", "broom", "sparsegl")) if (!requireNamespace(p, quietly = TRUE)) utils::install.packages(p, quiet = TRUE)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(tibble); library(broom); library(sparsegl) })
hdr8 <- function(x) cat("\n==================== ", x, " ====================\n", sep = "")
fmt_bm <- function(lst, nm) { has <- vapply(lst, length, integer(1)) > 0
  if (!any(has)) "none" else paste(sprintf("%s %s", nm[has], vapply(lst[has], function(a) paste(sort(a), collapse = ", "), "")), collapse = " | ") }
N_HIDDEN_WARN <- 0L
hide_weight_warning <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (grepl("non-integer #successes", conditionMessage(w), fixed = TRUE)) { N_HIDDEN_WARN <<- N_HIDDEN_WARN + 1L; invokeRestart("muffleWarning") } })
CMP_RESULT <- list(original = list(), same = list())
## Install a package if it is missing. Tried in this order: CRAN; the package's last source file in the CRAN archive (for packages
## that CRAN has archived); the read-only GitHub mirror of CRAN (github.com/cran/<package>). Only public package code is fetched.
get_pkg <- function(pkg, archive_url = NULL, github = NULL, needs = character(0)) {
  ok <- function() requireNamespace(pkg, quietly = TRUE)
  if (ok()) return(TRUE)
  quiet <- function(expr) try(suppressWarnings(suppressMessages(utils::capture.output(expr, type = "message"))), silent = TRUE)
  for (q in needs) if (!requireNamespace(q, quietly = TRUE)) quiet(utils::install.packages(q, quiet = TRUE))
  quiet(utils::install.packages(pkg, quiet = TRUE)); if (ok()) return(TRUE)
  if (!is.null(archive_url)) { quiet(utils::install.packages(archive_url, repos = NULL, type = "source", quiet = TRUE)); if (ok()) return(TRUE) }
  if (!is.null(github)) {
    if (nzchar(Sys.which("git"))) { d <- file.path(tempdir(), paste0(pkg, "_src")); unlink(d, recursive = TRUE)
      st <- suppressWarnings(system2("git", c("clone", "--depth", "1", "-q", paste0("https://github.com/", github), shQuote(d)), stdout = FALSE, stderr = FALSE))
      if (identical(as.integer(st), 0L)) quiet(utils::install.packages(d, repos = NULL, type = "source", quiet = TRUE))
      if (ok()) return(TRUE) }
    if (!requireNamespace("remotes", quietly = TRUE)) quiet(utils::install.packages("remotes", quiet = TRUE))
    if (requireNamespace("remotes", quietly = TRUE)) quiet(remotes::install_github(github, quiet = TRUE, upgrade = "never"))
    if (ok()) return(TRUE) }
  FALSE
}

## ---- cohort file (the same file and the same row rules as scripts 03 and 07) ---------------------------------------
.c8$loc <- if (nzchar(.c8$IN_CSV)) path.expand(.c8$IN_CSV) else sort(Sys.glob(path.expand("~/dementia_v9/dementia_v9_analysis_*.csv")), decreasing = TRUE)[1]
if (is.na(.c8$loc) || !file.exists(.c8$loc)) stop("Cohort file not found. Run dementia_cp_02_v9_cohort.R first, or set IN_CSV.", call. = FALSE)
hdr8("Settings")
cat(sprintf("  script: dementia_cp_08_compare.R, version 1.6 (2026-11-20) | cohort file: %s\n", basename(.c8$loc)))
cat(sprintf("  part A (the paper's scripts as written): %s | part B (same estimates for every method): %s\n", if (.c8$RUN_ORIGINAL) "yes" else "no", if (.c8$RUN_SAME) "yes" else "no"))
if (length(.c8$ADJUST) || length(.c8$MEASURED_ONLY)) cat("  ADJUST / MEASURED_ONLY apply to part B only; the paper's scripts have no such option and are run as written.\n")

## =====================================================================================================================
## PART A: the paper's scripts as written
## =====================================================================================================================
if (.c8$RUN_ORIGINAL) {
  hdr8("Part A: the paper's three lines as written")
  .c8$raw <- utils::read.csv(.c8$loc, stringsAsFactors = FALSE)
  .c8$raw <- .c8$raw %>% filter(is.finite(time2), time2 > 0, is.finite(age), age >= 45, age < 80)
  if (all(c("ldlc", "trig") %in% names(.c8$raw))) .c8$raw <- .c8$raw %>% filter(ldlc < 400, trig < 700)
  cat(sprintf("  rows used: N %s, dementia cases %s (lines 1 and 2 estimate on a random 60%% of them, seed 42, as in the scripts)\n",
              format(nrow(.c8$raw), big.mark = ","), format(sum(.c8$raw$Dementia), big.mark = ",")))

  ## ---- line 1: step2_group_lasso.R ------------------------------------------------------------------------------
  CMP_RESULT$original$line1 <- if (!1L %in% .c8$ORIGINAL_LINES) list(error = "skipped (ORIGINAL_LINES)") else tryCatch({
    if (!get_pkg("gglasso", github = "cran/gglasso")) stop("package gglasso could not be installed (tried CRAN and github.com/cran/gglasso)")
    suppressPackageStartupMessages(library(gglasso))
    local({
      data <- .c8$raw
      data <- data %>% filter(is.finite(time2), time2 > 0, is.finite(age))
      data <- data %>% filter(age >= 45, age < 80)
      set.seed(42)
      train_idx  <- sample(seq_len(nrow(data)), size = floor(0.60 * nrow(data)))
      train_data <- data[train_idx, ]
      invisible(utils::capture.output({
## >>>>>>>>>> copied unchanged from step2_group_lasso.R, sections 2-7 >>>>>>>>>>
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

perage <- hide_weight_warning(map_dfr(age_grid,
  ~ fit_one_adaptive_kernel(train_data, .x, biomarkers, min_events = 100)
) %>% arrange(term, age_grid))

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
## <<<<<<<<<< end of the copied sections <<<<<<<<<<
      }))
      level_only <- all(abs(alpha_hat) < 1e-8)                       # the jump at the second age then only sets the starting level
      cps <- if (level_only) setdiff(change_points, age_keep[2]) else change_points
      list(cps = cps, level_jump = if (level_only && age_keep[2] %in% change_points) age_keep[2] else NA, lambda = lam_final)
    })
  }, error = function(e) list(error = conditionMessage(e)))

  ## ---- line 2: step1_fused_lasso.R ------------------------------------------------------------------------------
  CMP_RESULT$original$line2 <- if (!2L %in% .c8$ORIGINAL_LINES) list(error = "skipped (ORIGINAL_LINES)") else tryCatch({
    if (!get_pkg("flam", archive_url = "https://cran.r-project.org/src/contrib/Archive/flam/flam_3.2.tar.gz", github = "cran/flam", needs = c("Rcpp", "MASS")))
      stop("package flam could not be installed (tried CRAN, the CRAN archive and github.com/cran/flam)")
    suppressPackageStartupMessages(library(flam))
    local({
      data <- .c8$raw
      data <- data %>% filter(is.finite(time2), time2 > 0, is.finite(age))
      set.seed(42)
      train_idx  <- sample(seq_len(nrow(data)), size = floor(0.60 * nrow(data)))
      train_data <- data[train_idx, ]
      invisible(utils::capture.output({
## >>>>>>>>>> copied unchanged from step1_fused_lasso.R, sections 2-3 >>>>>>>>>>
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
perage <- hide_weight_warning(map_dfr(age_grid,
  ~ fit_one_adaptive_kernel(train_data, .x, biomarkers, min_events = 100)
) %>% arrange(term, age_grid))
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
## <<<<<<<<<< end of the copied sections <<<<<<<<<<
      }))
      per <- lapply(biomarkers, function(bm) fused_all$age_grid[fused_all$term == bm & fused_all$cp]); names(per) <- biomarkers
      list(per = per, union = sort(unique(unlist(per))))
    })
  }, error = function(e) list(error = conditionMessage(e)))

  ## ---- line 3: step3_sparse_group_lasso.R, through script 07 with METHOD = "paper" ----------------------------------
  CMP_RESULT$original$line3 <- if (!3L %in% .c8$ORIGINAL_LINES) list(error = "skipped (ORIGINAL_LINES)") else tryCatch({
    IN_CSV <- .c8$loc; METHOD <- "paper"; LABEL <- "compare: paper line 3"; N_BOOT <- 0
    invisible(utils::capture.output(source(.c8$TUNE_FILE, local = FALSE)))
    list(cps = CP, level_jump = if (isTRUE(TUNE_RESULT$level_only) && length(TUNE_RESULT$cps)) TUNE_RESULT$age_keep[2] else NA)
  }, error = function(e) list(error = conditionMessage(e)))
}

## =====================================================================================================================
## PART B: every method on the same table of per-age estimates
## =====================================================================================================================
cmp_setup <- function(Y, SE) {
  b <- nrow(Y); K <- ncol(Y); J <- matrix(0, b, b - 1); for (t in 2:b) J[t, 1:(t - 1)] <- 1
  rw <- as.numeric(t(1 / SE)); XA <- kronecker(rep(1, b), diag(K)) * rw; XD <- kronecker(J, diag(K)) * rw
  P <- XA %*% solve(crossprod(XA), t(XA)); X <- XD - P %*% XD; y <- as.numeric(t(Y / SE)); y <- as.numeric(y - P %*% y)
  list(X = X, y = y, b = b, K = K, n = b * K, group = rep(seq_len(b - 1), each = K), jm = rep(seq_len(b - 1), each = K), jk = rep(seq_len(K), b - 1))
}
cmp_ridge <- function(X, y) {       # ridge start for the adaptive weights (penalty by generalized cross-validation), as in script 07
  s <- svd(X); d <- s$d; uy <- as.numeric(crossprod(s$u, y)); lam <- exp(seq(log(max(d)^2 * 1e-6), log(max(d)^2 * 10), length.out = 60))
  gcv <- vapply(lam, function(l) { f <- d^2 / (d^2 + l); sum((y - s$u %*% (f * uy))^2) / (length(y) - sum(f))^2 }, numeric(1))
  l <- lam[which.min(gcv)]; as.numeric(s$v %*% (d / (d^2 + l) * uy))
}
## candidate sets of jumps from sparse-group-lasso paths (adaptive weights from the ridge start); each distinct set once
cmp_sets <- function(X, y, group, alphas) {
  ini <- cmp_ridge(X, y); ng <- max(group); norm2 <- function(v) sqrt(sum(v^2))
  pfg <- vapply(seq_len(ng), function(g) { v <- norm2(ini[group == g]); if (v > 1e-8) 1 / v else 1e6 }, numeric(1))
  pfs <- ifelse(abs(ini) > 1e-8, 1 / abs(ini), 1e6); sets <- list(); seen <- character(0); note <- character(0)
  for (a in alphas) {
    fit <- NULL
    for (aa in unique(c(a, min(max(a, 0.001), 0.999)))) {          # some sparsegl versions refuse exactly 0 or 1
      fit <- try(suppressWarnings(sparsegl::sparsegl(X, y, group = group, asparse = aa, pf_group = pfg, pf_sparse = pfs, intercept = FALSE, standardize = FALSE)), silent = TRUE)
      if (!inherits(fit, "try-error")) { if (aa != a) note <- c(note, sprintf("asparse %g used in place of %g", aa, a)); break }
    }
    if (inherits(fit, "try-error")) next
    Bm <- as.matrix(fit$beta)
    for (j in seq_len(ncol(Bm))) { S <- which(abs(Bm[, j]) > 1e-8); key <- paste(S, collapse = ","); if (length(S) && !key %in% seen) { seen <- c(seen, key); sets[[length(sets) + 1]] <- S } }
  }
  attr(sets, "note") <- unique(note); sets
}
cmp_best <- function(X, y, sets, pen, df_cap) {                      # the candidate set with the smallest misfit + pen per jump (refitted)
  best <- list(crit = sum(y^2), S = integer(0))
  for (S in sets) { if (length(S) > df_cap) next
    r <- sum(stats::lm.fit(X[, S, drop = FALSE], y)$residuals^2) + pen * length(S); if (r < best$crit - 1e-9) best <- list(crit = r, S = S) }
  best$S
}
cmp_seg <- function(y, w, pen, allowed = NULL, kmax = 6) {           # exact search for one biomarker; allowed = positions where a new level may start
  Tn <- length(y); allowed <- if (is.null(allowed)) 2:Tn else sort(intersect(allowed, 2:Tn)); kmax <- min(kmax, length(allowed))
  cw <- c(0, cumsum(w)); cwy <- c(0, cumsum(w * y)); cwyy <- c(0, cumsum(w * y * y))
  cost <- function(i, j) { sw <- cw[j + 1] - cw[i]; swy <- cwy[j + 1] - cwy[i]; pmax(0, (cwyy[j + 1] - cwyy[i]) - swy^2 / sw) }
  Fm <- matrix(Inf, kmax + 1, Tn); Bk <- matrix(0L, kmax + 1, Tn); Fm[1, ] <- cost(1, seq_len(Tn))
  for (k in seq_len(kmax)) for (j in 2:Tn) { i <- allowed[allowed <= j & allowed > k]; if (!length(i)) next
    v <- Fm[k, i - 1] + cost(i, j); m <- which.min(v); if (is.finite(v[m])) { Fm[k + 1, j] <- v[m]; Bk[k + 1, j] <- i[m] } }
  k <- which.min(Fm[, Tn] + pen * (0:kmax)) - 1L; cuts <- integer(0); j <- Tn
  for (kk in rev(seq_len(k))) { i <- Bk[kk + 1, j]; cuts <- c(i, cuts); j <- i - 1L }
  cuts
}
## ---- line 3 fine-tuned: sparse group lasso, re-weighted step by step ----------------------------------------------------
## Y, SE: per-age estimates and standard errors (age groups x biomarkers). Every estimate is divided by its own SE and each
## biomarker's starting level is left free; every jump column is then rescaled to length 1, so a jump is judged by its evidence
## (like a z statistic) and one near either end of the age range is not at a disadvantage.
## Stage 0: the plain sparse group lasso. Stage s = 1, 2, ...: every jump gets a new penalty weight from the stage before,
##     weight of one jump = 1 / (|jump| + eps_s),  weight of one age (group) = 1 / (size of all jumps at that age + eps_s),
##     eps_s = eps0 x decay^(s - 1)   (halves every stage: the weights are tightened step by step),
## so jumps that keep being supported are penalised less and less and the others more and more.
## At every stage the mixing weight alpha (0 = all biomarkers jump together, 1 = each on its own) and the penalty strength are
## chosen by the combined criterion on refits without shrinkage:
##     misfit + log(n) x (age_cost x number of change ages + jump_cost x number of jumps).
## Stops when two stages in a row select the same jumps (at stage 2 at the earliest), or after max_stage stages.
std_perage <- function(Y, SE, ages, biomarkers)
  tibble(age_grid = rep(ages, times = length(biomarkers)), term = rep(biomarkers, each = length(ages)), estimate = as.numeric(Y), se = as.numeric(SE))
## line 2 standard: fused lasso of each biomarker's per-age estimates, penalty by 10-fold cross-validation (flamCV)
std_line2 <- function(Y, SE, ages, biomarkers) {
  perage <- std_perage(Y, SE, ages, biomarkers); age_grid <- ages
  invisible(utils::capture.output({
## >>>>>>>>>> copied unchanged from step1_fused_lasso.R, section 3 >>>>>>>>>>
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
## <<<<<<<<<< end of the copied section <<<<<<<<<<
  }))
  lapply(biomarkers, function(bm) fused_all$age_grid[fused_all$term == bm & fused_all$cp])
}
## line 1 standard: group lasso, a change age moves every biomarker; BIC with a df cap, one adaptive refit (gglasso)
std_line1 <- function(Y, SE, ages, biomarkers) {
  perage <- std_perage(Y, SE, ages, biomarkers)
  invisible(utils::capture.output({
## >>>>>>>>>> copied unchanged from step2_group_lasso.R, sections 3-7 >>>>>>>>>>
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
## <<<<<<<<<< end of the copied sections <<<<<<<<<<
  }))
  level_only <- all(abs(alpha_hat) < 1e-8)                 # then the jump at the second age only sets the starting level
  cps <- if (level_only) setdiff(change_points, age_keep[2]) else change_points
  replicate(length(biomarkers), cps, simplify = FALSE)     # a group-lasso change point moves every biomarker
}

c8_ft <- function(Y, SE, age_cost = 1.5, jump_cost = 1, eps0 = 1, decay = 0.5, max_stage = 6,
                     alphas = c(0, 0.25, 0.5, 0.75, 1), nlambda = 60, logn = NULL) {   # logn: set it when one biomarker is run on its own
  b <- nrow(Y); K <- ncol(Y); n <- b * K; if (is.null(logn)) logn <- log(n); W <- 1 / SE^2; norm2 <- function(v) sqrt(sum(v^2))
  J <- matrix(0, b, b - 1); for (t in 2:b) J[t, 1:(t - 1)] <- 1
  rw <- as.numeric(t(1 / SE)); XA <- kronecker(rep(1, b), diag(K)) * rw; XD <- kronecker(J, diag(K)) * rw
  P <- XA %*% solve(crossprod(XA), t(XA)); X <- XD - P %*% XD; y <- as.numeric(t(Y / SE)); y <- as.numeric(y - P %*% y)
  X <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  group <- rep(seq_len(b - 1), each = K); jm <- group; jk <- rep(seq_len(K), b - 1); ng <- b - 1; p <- ncol(X)
  rss_seg <- function(yy, w, cuts) { st <- c(1L, cuts); en <- c(cuts - 1L, length(yy)); s <- 0
    for (i in seq_along(st)) { r <- st[i]:en[i]; m <- sum(w[r] * yy[r]) / sum(w[r]); s <- s + sum(w[r] * (yy[r] - m)^2) }; s }
  to_cuts <- function(S) lapply(seq_len(K), function(k) sort(jm[S[jk[S] == k]] + 1L))
  pen_of <- function(S) logn * (age_cost * length(unique(jm[S])) + jump_cost * length(S))
  cache <- new.env()
  rss_refit <- function(S) { key <- if (length(S)) paste(S, collapse = ",") else "none"; v <- cache[[key]]
    if (is.null(v)) { cu <- to_cuts(S); v <- sum(vapply(seq_len(K), function(k) rss_seg(Y[, k], W[, k], cu[[k]]), numeric(1))); cache[[key]] <- v }; v }
  run_stage <- function(pfg, pfs) {
    best <- list(S = integer(0), v = rss_refit(integer(0)), alpha = NA, coef = numeric(p))
    for (a in alphas) for (aa in unique(c(a, min(max(a, 0.001), 0.999)))) {      # some sparsegl versions refuse exactly 0 or 1
      fit <- try(suppressWarnings(sparsegl::sparsegl(X, y, group = group, asparse = aa, pf_group = pfg, pf_sparse = pfs,
                                                     intercept = FALSE, standardize = FALSE, nlambda = nlambda)), silent = TRUE)
      if (inherits(fit, "try-error")) next
      Bm <- as.matrix(fit$beta)
      for (j in seq_len(ncol(Bm))) { S <- which(abs(Bm[, j]) > 1e-8); if (length(S) > 4 * K) next
        v <- rss_refit(S) + pen_of(S); if (v < best$v - 1e-9) best <- list(S = S, v = v, alpha = a, coef = Bm[, j]) }
      break }
    best }
  s0 <- run_stage(rep(1, ng), rep(1, p)); d <- s0$coef
  log <- list(list(stage = 0, eps = NA, alpha = s0$alpha, crit = s0$v, cuts = to_cuts(s0$S)))
  prev <- NULL; same <- 0; cur <- s0
  for (s in seq_len(max_stage)) {
    eps <- eps0 * decay^(s - 1); gn <- vapply(seq_len(ng), function(g) norm2(d[group == g]), numeric(1))
    cur <- run_stage(1 / (gn + eps), 1 / (abs(d) + eps)); d <- cur$coef
    log[[length(log) + 1]] <- list(stage = s, eps = eps, alpha = cur$alpha, crit = cur$v, cuts = to_cuts(cur$S))
    if (!is.null(prev) && identical(prev, cur$S)) same <- same + 1 else same <- 0
    prev <- cur$S; if (same >= 1 && s >= 2) break
  }
  list(cuts = to_cuts(cur$S), crit = cur$v, misfit = rss_refit(cur$S), n_ages = length(unique(jm[cur$S])), n_jumps = length(cur$S), log = log)
}

## ---- line 3 rebuilt: shared and separate jumps in one criterion ---------------------------------------------------------
## cost of a set of staircases = misfit + log(n) x sum over change ages of sched[m], m = number of biomarkers jumping at that age
rss_seg <- function(y, w, cuts) { st <- c(1L, cuts); en <- c(cuts - 1L, length(y)); s <- 0
  for (i in seq_along(st)) { r <- st[i]:en[i]; m <- sum(w[r] * y[r]) / sum(w[r]); s <- s + sum(w[r] * (y[r] - m)^2) }; s }
ga2_search <- function(Y, W, sched, logn, cand = 2:nrow(Y), max_iter = 200) {
  K <- ncol(Y); cuts <- replicate(K, integer(0), simplify = FALSE); rss <- vapply(seq_len(K), function(k) rss_seg(Y[, k], W[, k], integer(0)), numeric(1))
  price <- function(m) if (m <= 0) 0 else sched[m] * logn
  for (it in seq_len(max_iter)) {
    best <- list(gain = 1e-9)
    for (a in cand) {
      mem <- which(vapply(cuts, function(cc) a %in% cc, logical(1))); non <- setdiff(seq_len(K), mem); m <- length(mem)
      if (length(non)) {                                             # add the biomarkers that gain most from a jump at this age
        g <- vapply(non, function(k) rss[k] - rss_seg(Y[, k], W[, k], sort(c(cuts[[k]], a))), numeric(1)); o <- order(g, decreasing = TRUE)
        imp <- cumsum(g[o]) - (vapply(m + seq_along(o), price, numeric(1)) - price(m)); j <- which.max(imp)
        if (imp[j] > best$gain) best <- list(gain = imp[j], type = "add", a = a, ks = non[o[seq_len(j)]]) }
      if (m) {                                                       # drop the biomarkers that lose least when their jump at this age goes
        l <- vapply(mem, function(k) rss_seg(Y[, k], W[, k], setdiff(cuts[[k]], a)) - rss[k], numeric(1)); o <- order(l)
        imp <- (price(m) - vapply(m - seq_along(o), price, numeric(1))) - cumsum(l[o]); j <- which.max(imp)
        if (imp[j] > best$gain) best <- list(gain = imp[j], type = "drop", a = a, ks = mem[o[seq_len(j)]]) }
    }
    if (is.null(best$type)) break
    for (k in best$ks) { cuts[[k]] <- if (best$type == "add") sort(c(cuts[[k]], best$a)) else setdiff(cuts[[k]], best$a); rss[k] <- rss_seg(Y[, k], W[, k], cuts[[k]]) }
  }
  cuts
}
## candidate ages from sparse-group-lasso paths: for each mixing weight (0 = group lasso, 1 = lasso on single jumps), every age that
## is in a solution with at most max_ages change ages; plain paths and paths with adaptive weights
sgl_propose <- function(Y, SE, max_ages = 6, alphas = c(0, 0.25, 0.5, 0.75, 1)) {
  st <- cmp_setup(Y, SE); ini <- cmp_ridge(st$X, st$y); norm2 <- function(v) sqrt(sum(v^2)); ng <- max(st$group)
  pfg <- vapply(seq_len(ng), function(g) { v <- norm2(ini[st$group == g]); if (v > 1e-8) 1 / v else 1e6 }, numeric(1)); pfs <- ifelse(abs(ini) > 1e-8, 1 / abs(ini), 1e6)
  pos <- integer(0)
  for (a in alphas) for (w in 1:2) for (aa in unique(c(a, min(max(a, 0.001), 0.999)))) {      # some sparsegl versions refuse exactly 0 or 1
    args <- list(st$X, st$y, group = st$group, asparse = aa, intercept = FALSE, standardize = FALSE, pf_group = if (w == 1) rep(1, ng) else pfg); if (w == 2) args$pf_sparse <- pfs
    fit <- try(suppressWarnings(do.call(sparsegl::sparsegl, args)), silent = TRUE); if (inherits(fit, "try-error")) next
    Bm <- as.matrix(fit$beta)
    for (j in seq_len(ncol(Bm))) { g <- unique(st$jm[abs(Bm[, j]) > 1e-8]); if (length(g) && length(g) <= max_ages) pos <- union(pos, g + 1L) }
    break
  }
  sort(pos)
}
## all methods on one table of estimates; returns, per method, the list of change-point ages for every biomarker
cmp_all <- function(Y, SE, ages, bic_std = 1, bic_strict = 2.5, which = c("line1", "line2", "line3", "line3check", "exact")) {
  st <- cmp_setup(Y, SE); b <- st$b; K <- st$K; logn <- log(st$n); out <- list(); notes <- character(0)
  by_bm <- function(S) lapply(seq_len(K), function(k) ages[st$jm[S[st$jk[S] == k]] + 1])
  safe <- function(expr) tryCatch(expr, error = function(e) { notes <<- c(notes, conditionMessage(e)); NULL })
  if ("line1" %in% which) { s1 <- safe(cmp_sets(st$X, st$y, st$group, 0)); notes <- c(notes, attr(s1, "note"))
    if (!is.null(s1)) { s1 <- lapply(s1, function(S) which(st$jm %in% unique(st$jm[S])))          # all biomarkers jump at a chosen age
      for (bc in c(bic_std, bic_strict)) out[[sprintf("line1|%g", bc)]] <- by_bm(cmp_best(st$X, st$y, s1, bc * logn, b - 1)) } }
  if ("line2" %in% which) for (bc in c(bic_std, bic_strict)) out[[sprintf("line2|%g", bc)]] <- safe(lapply(seq_len(K), function(k) {
      cols <- which(st$jk == k); sk <- cmp_sets(st$X[, cols, drop = FALSE], st$y, seq_along(cols), 0.5)   # groups of size 1 = plain lasso on this biomarker's jumps
      ages[st$jm[cols[cmp_best(st$X[, cols, drop = FALSE], st$y, sk, bc * logn, b - 1)]] + 1] }))
  if (any(c("line3", "line3check") %in% which)) { s3 <- safe(cmp_sets(st$X, st$y, st$group, seq(0.05, 0.95, by = 0.05)))
    if (!is.null(s3)) { for (bc in c(bic_std, bic_strict)) out[[sprintf("line3|%g", bc)]] <- by_bm(cmp_best(st$X, st$y, s3, bc * logn, b - 1))
      if ("line3check" %in% which) { prop <- sort(unique(unlist(out[[sprintf("line3|%g", bic_std)]]))); pos <- match(prop, ages)
        out[[sprintf("line3check|%g", bic_strict)]] <- lapply(seq_len(K), function(k) ages[cmp_seg(Y[, k], 1 / SE[, k]^2, bic_strict * logn, allowed = pos)])
        attr(out, "proposed") <- prop } } }
  if ("exact" %in% which) for (bc in c(bic_std, bic_strict)) out[[sprintf("exact|%g", bc)]] <- lapply(seq_len(K), function(k) ages[cmp_seg(Y[, k], 1 / SE[, k]^2, bc * logn)])
  attr(out, "notes") <- unique(notes); out
}

if (.c8$RUN_SAME) {
  hdr8("Part B: every method on the same per-age estimates")
  .c8$ok <- tryCatch({
    IN_CSV <- .c8$loc; LABEL <- "compare: estimates"; N_BOOT <- 0; NULL_SIMS <- 0; BIC_C <- .c8$BIC_STRICT
    if (length(.c8$ADJUST)) ADJUST <- .c8$ADJUST
    if (length(.c8$MEASURED_ONLY)) MEASURED_ONLY <- .c8$MEASURED_ONLY
    if (length(.c8$BIOMARKERS)) BIOMARKERS <- .c8$BIOMARKERS
    if (length(.c8$BIN_MIN_EVENTS)) BIN_MIN_EVENTS <- .c8$BIN_MIN_EVENTS
    .c8$out07 <- utils::capture.output(source(.c8$TUNE_FILE, local = FALSE)); TRUE
  }, error = function(e) { cat("  the per-age estimates could not be made:", conditionMessage(e), "\n"); FALSE })
  if (.c8$ok) {
    cat(paste(grep("^  (input:|METHOD|only people|every per-age model|per-age estimates:)", .c8$out07, value = TRUE), collapse = "\n"), "\n")
    .c8$bm <- colnames(TUNE_RESULT$Y)
    CMP_RESULT$same <- cmp_all(TUNE_RESULT$Y, TUNE_RESULT$SE, TUNE_RESULT$age_keep, 1, .c8$BIC_STRICT)
    CMP_RESULT$same_exact07 <- lapply(seq_along(.c8$bm), function(k) TUNE_RESULT$cps[abs(TUNE_RESULT$jumps[, k]) > 1e-8])
    CMP_RESULT$line3new <- tryCatch({ Yb <- TUNE_RESULT$Y; Sb <- TUNE_RESULT$SE; ag <- TUNE_RESULT$age_keep; Kb <- ncol(Yb); sched <- (.c8$BIC_STRICT - 1) + seq_len(Kb)
      prop <- sgl_propose(Yb, Sb)
      list(per = lapply(ga2_search(Yb, 1 / Sb^2, sched, log(length(Yb)), cand = prop), function(q) ag[q]), proposed = ag[prop],
           every_age = lapply(ga2_search(Yb, 1 / Sb^2, sched, log(length(Yb))), function(q) ag[q])) }, error = function(e) list(error = conditionMessage(e)))
    CMP_RESULT$line3ft <- tryCatch({ r <- c8_ft(TUNE_RESULT$Y, TUNE_RESULT$SE, age_cost = .c8$BIC_STRICT - 1, jump_cost = 1)
      r$per <- lapply(r$cuts, function(q) TUNE_RESULT$age_keep[q]); r$ages <- TUNE_RESULT$age_keep; r }, error = function(e) list(error = conditionMessage(e)))
    ## lines 1 and 2 with the same fine-tuning; log(n) is that of the whole table in every line, so the three criteria are on one scale
    ## the paper's standard baselines on the same per-age estimates (their change-point code unchanged)
    CMP_RESULT$std1 <- tryCatch({ if (!get_pkg("gglasso", github = "cran/gglasso")) stop("package gglasso could not be installed (tried CRAN and github.com/cran/gglasso)")
      suppressPackageStartupMessages(library(gglasso)); list(per = std_line1(TUNE_RESULT$Y, TUNE_RESULT$SE, TUNE_RESULT$age_keep, .c8$bm)) },
      error = function(e) list(error = conditionMessage(e)))
    CMP_RESULT$std2 <- tryCatch({ if (!get_pkg("flam", archive_url = "https://cran.r-project.org/src/contrib/Archive/flam/flam_3.2.tar.gz", github = "cran/flam", needs = c("Rcpp", "MASS")))
        stop("package flam could not be installed (tried CRAN, the CRAN archive and github.com/cran/flam)")
      list(per = std_line2(TUNE_RESULT$Y, TUNE_RESULT$SE, TUNE_RESULT$age_keep, .c8$bm)) }, error = function(e) list(error = conditionMessage(e)))
    CMP_RESULT$line1ft <- tryCatch({ r <- c8_ft(TUNE_RESULT$Y, TUNE_RESULT$SE, age_cost = 1.5 * .c8$BIC_STRICT, jump_cost = 0, alphas = 0)
      list(per = lapply(r$cuts, function(q) TUNE_RESULT$age_keep[q])) }, error = function(e) list(error = conditionMessage(e)))
    CMP_RESULT$line2ft <- tryCatch({ Yb <- TUNE_RESULT$Y; Sb <- TUNE_RESULT$SE
      list(per = lapply(seq_len(ncol(Yb)), function(k) TUNE_RESULT$age_keep[c8_ft(Yb[, k, drop = FALSE], Sb[, k, drop = FALSE], age_cost = 0, jump_cost = .c8$BIC_STRICT,
                                                                                      alphas = 1, logn = log(length(Yb)))$cuts[[1]]])) }, error = function(e) list(error = conditionMessage(e)))
    CMP_RESULT$biomarkers <- .c8$bm
  }
}
IN_CSV <- .c8$loc

## =====================================================================================================================
## the comparison
## =====================================================================================================================
.c8$txt <- utils::capture.output({
  hdr8("Comparison: change points found by each method")
  cat(sprintf("  cohort file: %s\n", basename(.c8$loc)))
  if (.c8$RUN_ORIGINAL) {
    cat("\n  A. The paper's three scripts as written, each with its own per-age estimates (kernel windows, 60% sample for lines 1 and 2, no adjustment)\n")
    o <- CMP_RESULT$original; show <- function(x, f) if (!is.null(x$error)) { if (startsWith(x$error, "skipped")) "skipped" else paste("could not be run:", x$error) } else f(x)
    lv <- function(x) if (!is.na(x$level_jump)) sprintf("   (a jump at %d only sets the starting level and is not counted)", x$level_jump) else ""
    cat(sprintf("     %-38s %s\n", "line 1  group lasso (shared)", show(o$line1, function(x) paste0(if (length(x$cps)) paste(x$cps, collapse = ", ") else "none", lv(x)))))
    if (!is.null(o$line2$error)) cat(sprintf("     %-38s %s\n", "line 2  fused lasso (per biomarker)", show(o$line2, identity))) else {
      cat(sprintf("     %-38s %d different ages in all\n", "line 2  fused lasso (per biomarker)", length(o$line2$union)))
      for (v in names(o$line2$per)) cat(sprintf("     %-38s %s\n", paste0("          ", v, sprintf(" (%d)", length(o$line2$per[[v]]))),
                                                if (length(o$line2$per[[v]])) paste(o$line2$per[[v]], collapse = ", ") else "none")) }
    cat(sprintf("     %-38s %s\n", "line 3  sparse group lasso", show(o$line3, function(x) paste0(if (length(x$cps)) paste(x$cps, collapse = ", ") else "none", lv(x)))))
  }
  if (.c8$RUN_SAME && length(CMP_RESULT$same)) {
    r <- CMP_RESULT$same; bm <- CMP_RESULT$biomarkers; g <- function(key) if (is.null(r[[key]])) "could not be run" else fmt_bm(r[[key]], bm)
    rows <- list(c("line1", "line 1  group lasso (shared)"), c("line2", "line 2  fused lasso (per biomarker)"), c("line3", "line 3  sparse group lasso"))
    cat(sprintf("\n  B. THE COMPARISON FOR THE PAPER: the three lines on the same per-age estimates%s\n",
                if (length(.c8$ADJUST)) paste0(" (adjusted for ", paste(.c8$ADJUST, collapse = ", "), ")") else " (no adjustment)"))
    n3 <- CMP_RESULT$line3new; f3 <- CMP_RESULT$line3ft; f1 <- CMP_RESULT$line1ft; f2 <- CMP_RESULT$line2ft; s1 <- CMP_RESULT$std1; s2 <- CMP_RESULT$std2
    shw <- function(x, all_same = FALSE) { if (!is.null(x$error)) return(paste("could not be run:", x$error))
      if (all_same && length(x$per) && all(vapply(x$per, function(a) identical(sort(a), sort(x$per[[1]])), logical(1))) && length(x$per[[1]]))
        return(sprintf("%s (all %d biomarkers jump there)", paste(sort(x$per[[1]]), collapse = ", "), length(x$per)))
      fmt_bm(x$per, bm) }
    cat(sprintf("     %-56s %s\n", "line 1  overall: group lasso, standard (step2)", shw(s1, TRUE)))
    cat(sprintf("     %-56s %s\n", "line 2  each on its own: fused lasso, standard (step1)", shw(s2)))
    cat(sprintf("     %-56s %s\n", "line 3  combined: sparse group lasso, fine-tuned", shw(f3)))
    cat("     lines 1 and 2 run the change-point code of step2 / step1 unchanged; only their per-age estimates are the ones line 3 uses.\n")
    cat("     line 1 penalizes the starting level, so it tends to put a change near the youngest ages even when nothing changes (in simulations\n")
    cat("     with no change at all, mostly at 48-50); line 2 picks its penalty by cross-validation and reports about 10 jumps when nothing changes.\n")
    if (is.null(f3$error)) {
      cat("\n     line 3 stage by stage (stage 0 = the plain sparse group lasso; later stages re-weighted, eps halved each time):\n")
      for (l in f3$log) cat(sprintf("       stage %d %-11s alpha %-5s criterion %7.1f | %s\n", l$stage, if (l$stage == 0) "" else sprintf("(eps %g)", l$eps),
                                    if (is.na(l$alpha)) "-" else format(l$alpha), l$crit, fmt_bm(lapply(l$cuts, function(q) f3$ages[q]), bm)))
      cat(sprintf("       criterion = misfit + log(n) x (%g per change age + 1 per jump): one biomarker on its own pays %g x log(n); biomarkers that jump\n", .c8$BIC_STRICT - 1, .c8$BIC_STRICT))
      cat("       at the same age share the price of that age. alpha 0 = all biomarkers together (as line 1), 1 = each on its own (as line 2).\n")
    }
    cat("\n  SUPPLEMENT (not the main comparison)\n")
    cat("     the three lines fine-tuned in the same way, so that they differ only in the penalty:\n")
    cat(sprintf("     %-56s %s\n", "  line 1  overall (group lasso)", shw(f1, TRUE)))
    cat(sprintf("     %-56s %s\n", "  line 2  each on its own (fused lasso)", shw(f2)))
    cat(sprintf("     %-56s %s\n", "  line 3  combined (sparse group lasso)", shw(f3)))
    cat(sprintf("     %-56s %s\n", "exact search, each biomarker on its own", g(sprintf("exact|%g", .c8$BIC_STRICT))))
    cat(sprintf("     %-56s %s\n", "line 3 rebuilt (greedy search, same criterion)", if (!is.null(n3$error)) paste("could not be run:", n3$error) else fmt_bm(n3$per, bm)))
    cat("     the three lines as plain lasso paths with the refit criterion, without the fine-tuning:\n")
    for (m in rows) cat(sprintf("     %-56s %s\n", paste0("  ", m[2]), g(sprintf("%s|%g", m[1], .c8$BIC_STRICT))))
    if (length(attr(r, "notes"))) cat("     notes:", paste(attr(r, "notes"), collapse = "; "), "\n")
  }
})
cat(paste(.c8$txt, collapse = "\n"), "\n")
.c8$out_dir <- path.expand("~/dementia_v9"); dir.create(.c8$out_dir, showWarnings = FALSE)
.c8$f <- file.path(.c8$out_dir, sprintf("compare_%s.txt", format(Sys.time(), "%Y%m%d-%H%M")))
writeLines(c(sprintf("dementia_cp_08_compare.R, run %s", format(Sys.time(), "%Y-%m-%d %H:%M")), .c8$txt), .c8$f)
cat(sprintf("\n  saved: %s  (change-point ages only)\n", sub(path.expand("~"), "~", .c8$f, fixed = TRUE)))
if (N_HIDDEN_WARN > 0) cat(sprintf("  (%d expected 'non-integer #successes' warnings of the kernel-weighted models were hidden)\n", N_HIDDEN_WARN))
cat(sprintf("\nruntime %.1f min\n", as.numeric(difftime(Sys.time(), .c8$t0, units = "mins"))))
