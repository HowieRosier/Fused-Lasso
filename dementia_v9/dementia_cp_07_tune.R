## =============================================================================
## dementia_cp_07_tune.R
## Change-point search with switchable settings, for trying variants of the method on one cohort file.
##
## METHOD <- "fixed" (default): the paper's two-step idea with its weak points repaired
##     - per-age estimates from one-year age groups that do not overlap (so they are independent);
##     - the fit is weighted by each estimate's standard error;
##     - the starting level of each curve is not penalized (only jumps are);
##     - sparsegl does not rescale the design (a jump costs the same at every age);
##     - a staircase is judged by BIC with the known error variance: misfit + BIC_C x log(n) for each jump;
##     - SEARCH <- "exact" (default): for each biomarker, the staircase with the smallest criterion is found exactly
##       (dynamic programming over every possible set of jump ages). Nothing is left to a search heuristic.
##       SEARCH <- "ft": line 3 fine-tuned. The sparse group lasso, re-weighted step by step: after each stage every jump gets a
##       new penalty weight 1 / (|jump| + eps) and every age 1 / (size of its jumps + eps), with eps halved every stage, so the
##       jumps that keep being supported are penalised less and less. Each stage picks the mixing weight and the penalty strength
##       by the combined criterion: misfit + log(n) x ((BIC_C - 1) per change age + 1 per jump), so one biomarker jumping on its
##       own pays BIC_C x log(n), as in the exact search, and biomarkers jumping at the same age share the price of that age.
##       Settings: FT_EPS0 (1), FT_DECAY (0.5), FT_STAGES (6).
##       SEARCH <- "sgl": the sparse-group-lasso path proposes sets of jumps, each set is refitted without shrinkage
##       (versions 1-5 of this script). That search only proposes jumps in which several biomarkers move at the same
##       age, so it can miss a staircase in which one biomarker jumps on its own;
##     - the penalty multiplier BIC_C was set from simulations with no true change (see FIXED_BIC_C_EXACT below);
##     - stability is checked on random halves of the cohort (N_BOOT resamples), not on bootstrap samples.
## METHOD <- "paper": exactly the method of dementia_cp_03_sgl_v9.R (repo dementia/step3).
## Two robustness checks (both use columns already in the cohort file; nothing is re-extracted):
##     ADJUST <- c("sex_at_birth", "anti", "imp_hdl")   extra columns added to every per-age model (not searched for change points).
##         imp_<lab> = 1 marks a lab value that was missing and filled in; "anti" = on blood-pressure drugs at baseline.
##     MEASURED_ONLY <- c("hdl")                        keep only people whose value of these biomarkers was really measured.
## PRINT_PERAGE <- TRUE also prints the estimate of every age group and saves a figure of them (points with 95% intervals,
##     the fitted steps and the straight line) in ~/dementia_v9/. Estimates only: no counts of people are printed or drawn.
## Each switch can also be set by hand (see Settings). Set what you want, then source() the file:
##     N_BOOT <- 50; source("~/dementia_cp_07_tune.R")
##     METHOD <- "paper"; LABEL <- "paper method"; source("~/dementia_cp_07_tune.R")
##     BIOMARKERS <- c("sbp", "ldlc", "hdl"); BIC_C <- 0.75; source("~/dementia_cp_07_tune.R")
## After a run, the settings you set are cleared (except IN_CSV), so the next run starts from the defaults.
##
## Every run is appended to ~/dementia_v9/tune_log.csv: settings, change points, bootstrap stability.
## That file is the record of what was tried. Output is model-level only (coefficients, change points).
## =============================================================================

.opt <- function(name, default) if (exists(name, envir = globalenv(), inherits = FALSE)) get(name, envir = globalenv()) else default
.TUNE_NAMES <- c("METHOD", "LABEL", "BIOMARKERS", "N_BOOT", "N_CORES", "AGE_EST", "WEIGHT", "BASE_PEN", "STANDARDIZE", "CRIT",
                 "ADAPT_INIT", "MAX_STAGES", "DF_CAP", "MIN_EVENTS", "H_MAX", "BIN_MIN_EVENTS", "BIC_C", "GAMMA", "ALPHA_STEP", "RESAMPLE", "REFIT", "SHAPE",
                 "PRINT_PERAGE", "KEEP_SETTINGS", "SEARCH", "MAX_JUMPS", "NULL_SIMS", "ADJUST", "MEASURED_ONLY",
                 "FT_EPS0", "FT_DECAY", "FT_STAGES")

## ---- Settings (shown with their defaults) ----------------------------------------------
IN_FILE <- .opt("IN_CSV", "")            # IN_CSV: "" = newest ~/dementia_v9/dementia_v9_analysis_*.csv (from dementia_cp_02_v9_cohort.R)
cfg <- list(METHOD = match.arg(.opt("METHOD", "fixed"), c("fixed", "paper")))
## Default penalty multiplier of the "fixed" method. 2 x log(n) is close to 2 x log(number of possible jumps), the usual
## penalty when many candidate positions are searched. In simulated cohorts of this size with no true change it gave
## a false change point in 0 of 20 cohorts (BIC_C = 1.5: 2 of 20; BIC_C = 1: 3 of 33). With BIC_C = 1 the method also reports
## change points when the effects only drift gradually with age (25 of 26 simulated cohorts): see the note under the stability table.
FIXED_BIC_C <- 2
FIXED_BIC_C_SLOPE <- 2   # same for changes of slope (checked on simulated cohorts whose effects drift in a straight line)
.shape <- as.character(.opt("SHAPE", "step"))
if (!identical(.shape, "step")) suppressWarnings(rm(list = intersect("SHAPE", ls(envir = globalenv())), envir = globalenv()))
if (!identical(.shape, "step")) stop("SHAPE = \"slope\" is switched off: in simulated cohorts of this size it found 0 of 10 true changes of slope ",
                                     "and placed a very strong one 8 years off. Use SHAPE = \"step\".", call. = FALSE)
## Default penalty multiplier of the exact search. The exact search looks at every possible jump age, so it needs a larger
## penalty than the lasso-path search. With the standard errors of the v9 run and no true change, it reported at least one
## change point in 41% of simulated data sets at BIC_C = 1.5, 11% at 2, 3% at 2.5 and 1% at 3. Every run re-computes these
## rates for its own data (table 'Change points at other penalty multipliers').
FIXED_BIC_C_EXACT <- 2.5
.search <- match.arg(as.character(.opt("SEARCH", if (cfg$METHOD == "paper") "sgl" else "exact")), c("exact", "sgl", "ft"))
if (cfg$METHOD == "paper" && .search %in% c("exact", "ft")) { suppressWarnings(rm(list = intersect("SEARCH", ls(envir = globalenv())), envir = globalenv()))
  stop("METHOD = \"paper\" has its own search; SEARCH = \"exact\" and \"ft\" belong to METHOD = \"fixed\".", call. = FALSE) }
.pre <- if (cfg$METHOD == "paper")
  list(AGE_EST = "kernel", WEIGHT = "sd", BASE_PEN = TRUE,  STANDARDIZE = TRUE,  CRIT = "logmse", ADAPT_INIT = "sgl",   MAX_STAGES = 5, RESAMPLE = "boot", REFIT = FALSE, BIC_C = 1) else
  list(AGE_EST = "bins",   WEIGHT = "se", BASE_PEN = FALSE, STANDARDIZE = FALSE, CRIT = "known",  ADAPT_INIT = "ridge", MAX_STAGES = 1, RESAMPLE = "half", REFIT = TRUE,
       BIC_C = if (.shape == "slope") FIXED_BIC_C_SLOPE else if (.search %in% c("exact", "ft")) FIXED_BIC_C_EXACT else FIXED_BIC_C)
cfg$LABEL          <- as.character(.opt("LABEL", ""))                   # free text for the log
cfg$SHAPE          <- .shape                                            # a change point = the effect jumps to a new level
cfg$SEARCH         <- .search                                           # "exact": best staircase found exactly / "ft": line 3 fine-tuned / "sgl": lasso path (versions 1-5)
cfg$FT_EPS0        <- as.numeric(.opt("FT_EPS0", 1))                    # line 3 fine-tuned: eps of the first re-weighting stage
cfg$FT_DECAY       <- as.numeric(.opt("FT_DECAY", 0.5))                 # line 3 fine-tuned: eps is multiplied by this at every stage
cfg$FT_STAGES      <- as.integer(.opt("FT_STAGES", 6))                  # line 3 fine-tuned: most re-weighting stages
cfg$ADJUST         <- as.character(.opt("ADJUST", character(0)))        # extra columns of the cohort file added to every per-age model
cfg$MEASURED_ONLY  <- as.character(.opt("MEASURED_ONLY", character(0))) # biomarkers that must be really measured (imp_<lab> = 0)
cfg$MAX_JUMPS      <- as.integer(.opt("MAX_JUMPS", 6))                  # exact search: most jumps allowed for one biomarker
cfg$NULL_SIMS      <- as.integer(.opt("NULL_SIMS", 1000))               # exact search: simulations behind the 'false alarm' and 'step or gradual' numbers (0 = skip)
cfg$PRINT_PERAGE   <- isTRUE(.opt("PRINT_PERAGE", FALSE))               # TRUE: also print the estimates of the separate age groups
cfg$BIOMARKERS     <- .opt("BIOMARKERS", c("sbp", "bmi", "ldlc", "hdl", "trig"))   # columns of the cohort file
cfg$N_BOOT         <- as.integer(.opt("N_BOOT", 0))                     # resamples for the stability check (0 = none)
cfg$RESAMPLE       <- match.arg(.opt("RESAMPLE", .pre$RESAMPLE), c("half", "boot"))  # each resample: a random half of the people, or a bootstrap sample
cfg$N_CORES        <- as.integer(.opt("N_CORES", max(1L, parallel::detectCores())))
cfg$AGE_EST        <- match.arg(.opt("AGE_EST", .pre$AGE_EST), c("bins", "kernel"))     # per-age estimates: one-year groups / Gaussian kernel
cfg$WEIGHT         <- match.arg(.opt("WEIGHT", .pre$WEIGHT), c("se", "sd"))            # fit weighted by SE / each curve scaled by its SD over ages
cfg$BASE_PEN       <- isTRUE(.opt("BASE_PEN", .pre$BASE_PEN))                          # penalize the starting level (paper) or not
cfg$STANDARDIZE    <- isTRUE(.opt("STANDARDIZE", .pre$STANDARDIZE))                    # let sparsegl rescale the design columns
cfg$CRIT           <- match.arg(.opt("CRIT", .pre$CRIT), c("known", "logmse"))         # BIC with known variance (needs WEIGHT = "se") / log(MSE) BIC
cfg$ADAPT_INIT     <- match.arg(.opt("ADAPT_INIT", .pre$ADAPT_INIT), c("ridge", "sgl")) # where the first adaptive weights come from
cfg$MAX_STAGES     <- as.integer(.opt("MAX_STAGES", .pre$MAX_STAGES))                  # adaptive refits
cfg$DF_CAP         <- as.numeric(.opt("DF_CAP", 1))                                    # largest model: DF_CAP x (number of age points - 1)
cfg$MIN_EVENTS     <- as.numeric(.opt("MIN_EVENTS", 100))                              # kernel: cases needed in the window
cfg$H_MAX          <- as.numeric(.opt("H_MAX", 10))                                    # kernel: largest bandwidth
cfg$BIN_MIN_EVENTS <- as.numeric(.opt("BIN_MIN_EVENTS", 40))                           # bins: neighbouring ages are merged until this many cases
cfg$REFIT          <- isTRUE(.opt("REFIT", .pre$REFIT))                                # judge each candidate set of jumps after refitting it without shrinkage
cfg$BIC_C          <- as.numeric(.opt("BIC_C", .pre$BIC_C))                            # multiplier on the BIC penalty (smaller = more change points)
cfg$GAMMA          <- as.numeric(.opt("GAMMA", 1))                                     # adaptive-weight exponent
cfg$ALPHA_STEP     <- as.numeric(.opt("ALPHA_STEP", 0.05))                             # grid step for the lasso / group-lasso mix
KEEP <- isTRUE(.opt("KEEP_SETTINGS", FALSE))
if (cfg$CRIT == "known" && cfg$WEIGHT != "se") stop("CRIT = \"known\" needs WEIGHT = \"se\".", call. = FALSE)
if (cfg$REFIT && (cfg$CRIT != "known" || cfg$BASE_PEN)) stop("REFIT = TRUE needs CRIT = \"known\" and BASE_PEN = FALSE.", call. = FALSE)
if (cfg$SHAPE == "slope" && cfg$BASE_PEN) stop("SHAPE = \"slope\" needs BASE_PEN = FALSE (use METHOD = \"fixed\").", call. = FALSE)
if (cfg$ADAPT_INIT == "ridge") cfg$MAX_STAGES <- max(1L, cfg$MAX_STAGES)
if (cfg$SEARCH %in% c("exact", "ft") && (cfg$AGE_EST != "bins" || cfg$WEIGHT != "se" || cfg$BASE_PEN || cfg$CRIT != "known"))
  stop("SEARCH = \"exact\" or \"ft\" needs AGE_EST = \"bins\", WEIGHT = \"se\", BASE_PEN = FALSE and CRIT = \"known\" (the METHOD = \"fixed\" defaults).", call. = FALSE)

SCRIPT_VERSION <- "8 (2026-11-20)"
t0 <- Sys.time()
options(width = 170)
for (p in c("dplyr", "tidyr", "purrr", "tibble", "broom", "sparsegl")) {
  if (!requireNamespace(p, quietly = TRUE)) utils::install.packages(p, quiet = TRUE)
}
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(tibble); library(sparsegl) })
if (cfg$SEARCH == "sgl" && !cfg$STANDARDIZE && !"standardize" %in% names(formals(sparsegl::sparsegl)))
  stop("This sparsegl version has no `standardize` switch; set STANDARDIZE <- TRUE or update sparsegl.", call. = FALSE)
hdr <- function(s) cat("\n==================== ", s, " ====================\n", sep = "")
biomarkers <- cfg$BIOMARKERS; K <- length(biomarkers)

## ---- data ---------------------------------------------------------------------------------
loc <- if (nzchar(IN_FILE)) path.expand(IN_FILE) else {
  f <- sort(Sys.glob(path.expand("~/dementia_v9/dementia_v9_analysis_*.csv")), decreasing = TRUE); if (length(f)) f[1] else ""
}
if (!nzchar(loc) || !file.exists(loc)) stop("Cohort file not found: run dementia_cp_02_v9_cohort.R first, or set IN_CSV.", call. = FALSE)
raw <- utils::read.csv(loc)
if (any(cfg$ADJUST %in% c(biomarkers, "Dementia", "age"))) stop("ADJUST must not repeat a biomarker, Dementia or age.", call. = FALSE)
if (length(setdiff(cfg$MEASURED_ONLY, biomarkers))) stop("MEASURED_ONLY must name biomarkers of this run.", call. = FALSE)
flag_cols <- if (length(cfg$MEASURED_ONLY)) paste0("imp_", cfg$MEASURED_ONLY) else character(0)
miss <- setdiff(c("Dementia", "age", "time2", biomarkers, cfg$ADJUST, flag_cols), names(raw))
if (length(miss)) stop("The cohort file has no column(s): ", paste(miss, collapse = ", "), call. = FALSE)
raw <- raw %>% filter(is.finite(time2), time2 > 0, is.finite(age), age >= 45, age < 80)
if (all(c("ldlc", "trig") %in% names(raw))) raw <- raw %>% filter(ldlc < 400, trig < 700)
raw <- raw[stats::complete.cases(raw[, biomarkers, drop = FALSE]), unique(c("Dementia", "age", biomarkers, cfg$ADJUST, flag_cols))]
N_MAIN <- nrow(raw); EV_MAIN <- sum(raw$Dementia)
for (fc in flag_cols) raw <- raw[!is.na(raw[[fc]]) & raw[[fc]] == 0, , drop = FALSE]
ADJ_NOTE <- character(0)
for (v in cfg$ADJUST) {             # text columns: missing = its own group, groups under 200 people are pooled; number columns: rows with a missing value are left out
  if (is.numeric(raw[[v]]) || is.logical(raw[[v]])) { raw[[v]] <- as.numeric(raw[[v]]); raw <- raw[is.finite(raw[[v]]), , drop = FALSE]
    ADJ_NOTE <- c(ADJ_NOTE, sprintf("%s (number)", v))
  } else { x <- as.character(raw[[v]]); x[is.na(x) | !nzchar(x)] <- "unknown"; tb <- table(x); x[x %in% names(tb)[tb < 200]] <- "other (pooled)"
    raw[[v]] <- factor(x); ADJ_NOTE <- c(ADJ_NOTE, sprintf("%s (%d groups)", v, nlevels(raw[[v]]))) }
}
raw <- raw[, unique(c("Dementia", "age", biomarkers, cfg$ADJUST)), drop = FALSE]
N_ALL <- nrow(raw); EV_ALL <- sum(raw$Dementia)

## ---- step 1: per-age estimates ------------------------------------------------------------
form_txt <- paste("Dementia ~", paste(biomarkers, collapse = " + "))   # the formula is built inside each fit so glm finds the weights
form_for <- function(dwin) {        # adjusters are added only where they vary inside this age group; a group seen fewer than 10 times here is pooled
  if (!length(cfg$ADJUST)) return(list(txt = form_txt, data = dwin))
  use <- character(0)
  for (v in cfg$ADJUST) {
    if (is.factor(dwin[[v]])) { x <- as.character(dwin[[v]]); tb <- table(x); x[x %in% names(tb)[tb < 10]] <- "other (pooled)"; dwin[[v]] <- factor(x) }
    if (length(unique(dwin[[v]])) > 1) use <- c(use, v)
  }
  list(txt = paste(c(form_txt, sprintf("`%s`", use)), collapse = " + "), data = dwin)
}
fit_kernel <- function(df, a_g) {   # as in the paper's code
  for (h in seq(1, cfg$H_MAX, by = 0.5)) {
    w <- exp(-((df$age - a_g)^2) / (2 * h^2)); keep <- w > 0.01 * max(w)
    if (sum(df$Dementia[keep]) >= cfg$MIN_EVENTS) break
  }
  dwin <- df[keep, ]; wk <- w[keep]
  out <- tibble(age_grid = a_g, age_to = a_g, term = biomarkers, bandwidth = h, estimate = NA_real_, se = NA_real_)
  if (nrow(dwin) < 50) return(out)
  ff <- form_for(dwin); form <- stats::as.formula(ff$txt); dwin <- ff$data   # the formula is created here so that glm finds the weights
  mod <- try(suppressWarnings(stats::glm(form, data = dwin, family = stats::binomial(), weights = wk)), silent = TRUE)
  if (inherits(mod, "try-error")) return(out)
  co <- summary(mod)$coefficients; i <- match(biomarkers, rownames(co))
  out$estimate <- co[i, "Estimate"]; out$se <- co[i, "Std. Error"]; out
}
make_bins <- function(df) {         # one-year age groups; neighbours are merged until each has BIN_MIN_EVENTS cases
  ev <- vapply(45:79, function(a) sum(df$Dementia[floor(df$age) == a]), numeric(1))
  bins <- list(); cur <- integer(0); cum <- 0
  for (i in seq_along(ev)) {
    cur <- c(cur, 44L + i); cum <- cum + ev[i]
    if (cum >= cfg$BIN_MIN_EVENTS) { bins[[length(bins) + 1]] <- cur; cur <- integer(0); cum <- 0 }
  }
  if (length(cur)) { if (length(bins)) bins[[length(bins)]] <- c(bins[[length(bins)]], cur) else bins[[1]] <- cur }
  bins
}
fit_bin <- function(df, a_floor, bn) {
  dwin <- df[a_floor %in% bn, ]
  out <- tibble(age_grid = min(bn), age_to = max(bn), term = biomarkers, bandwidth = NA_real_, estimate = NA_real_, se = NA_real_)
  if (nrow(dwin) < 50 || sum(dwin$Dementia) < 5) return(out)
  ff <- form_for(dwin); form <- stats::as.formula(ff$txt); dwin <- ff$data   # the formula is created here so that glm finds the weights
  mod <- try(suppressWarnings(stats::glm(form, data = dwin, family = stats::binomial())), silent = TRUE)
  if (inherits(mod, "try-error") || !isTRUE(mod$converged)) return(out)
  co <- summary(mod)$coefficients; i <- match(biomarkers, rownames(co))
  out$estimate <- co[i, "Estimate"]; out$se <- co[i, "Std. Error"]; out
}
step1 <- function(d, bins) {
  if (cfg$AGE_EST == "kernel") return(map_dfr(45:79, ~ fit_kernel(d, .x)))
  a_floor <- floor(d$age); map_dfr(bins, ~ fit_bin(d, a_floor, .x))
}

## ---- step 2: piecewise-constant curves (staircases) ---------------------------------------------
norm2 <- function(v) sqrt(sum(v^2))
## Exact search for one biomarker. For every number of jumps 0..kmax: the smallest weighted misfit any staircase with that
## many jumps can reach, and where its jumps are (cuts = position of the first age group of each new level).
seg_exact <- function(y, w, kmax) {
  Tn <- length(y); kmax <- max(0L, min(kmax, Tn - 1L)); cw <- c(0, cumsum(w)); cwy <- c(0, cumsum(w * y)); cwyy <- c(0, cumsum(w * y * y))
  cost <- function(i, j) { sw <- cw[j + 1] - cw[i]; swy <- cwy[j + 1] - cwy[i]; pmax(0, (cwyy[j + 1] - cwyy[i]) - swy^2 / sw) }   # one level on points i..j
  Fm <- matrix(Inf, kmax + 1, Tn); Bk <- matrix(0L, kmax + 1, Tn); Fm[1, ] <- cost(1, seq_len(Tn))
  for (k in seq_len(kmax)) for (j in (k + 1):Tn) {
    i <- (k + 1):j; v <- Fm[k, i - 1] + cost(i, j); m <- which.min(v); Fm[k + 1, j] <- v[m]; Bk[k + 1, j] <- i[m]
  }
  lapply(0:kmax, function(k) { cuts <- integer(0); j <- Tn
    for (kk in rev(seq_len(k))) { i <- Bk[kk + 1, j]; cuts <- c(i, cuts); j <- i - 1L }
    list(rss = Fm[k + 1, Tn], cuts = cuts) })
}
seg_pick <- function(tab, pen) tab[[which.min(vapply(tab, function(x) x$rss, numeric(1)) + pen * (seq_along(tab) - 1))]]
## best single step and straight line for one biomarker: how much each lowers the misfit of 'no change'
scan_step_line <- function(y, w, mid) {
  n <- length(y); cw <- cumsum(w); cwy <- cumsum(w * y); tot <- cwy[n]; totw <- cw[n]
  gain <- cwy[-n]^2 / cw[-n] + (tot - cwy[-n])^2 / (totw - cw[-n]) - tot^2 / totw; j <- which.max(gain)
  m <- tot / totw; am <- sum(w * mid) / totw; sxx <- sum(w * (mid - am)^2); sl <- sum(w * (mid - am) * (y - m)) / sxx
  list(step_gain = gain[j], at = j + 1L, line_gain = sl^2 * sxx, level = m, slope = sl, am = am)
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
ft_line3 <- function(Y, SE, age_cost = 1.5, jump_cost = 1, eps0 = 1, decay = 0.5, max_stage = 6,
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
ridge_init <- function(X, y) {      # ridge fit, penalty chosen by generalized cross-validation
  s <- svd(X); d2 <- s$d^2; uy <- as.numeric(crossprod(s$u, y)); n <- length(y)
  taus <- exp(seq(log(1e-6), log(1e3), length.out = 61)) * mean(d2)
  gcv <- vapply(taus, function(tau) {
    sh <- d2 / (d2 + tau); rss <- sum(y^2) - sum(uy^2) + sum(((1 - sh) * uy)^2); n * rss / (n - sum(sh))^2
  }, numeric(1))
  tau <- taus[which.min(gcv)]
  as.numeric(s$v %*% (s$d / (d2 + tau) * uy))
}
step2 <- function(pa) {
  mid_all <- pa %>% distinct(age_grid, age_to) %>% arrange(age_grid) %>% mutate(mid = if (cfg$AGE_EST == "bins") (age_grid + age_to + 1) / 2 else age_grid)
  wide_e <- pa %>% select(age_grid, term, estimate) %>% pivot_wider(names_from = term, values_from = estimate) %>% arrange(age_grid)
  wide_s <- pa %>% select(age_grid, term, se) %>% pivot_wider(names_from = term, values_from = se) %>% arrange(age_grid)
  Y_all <- as.matrix(wide_e[, biomarkers, drop = FALSE]); SE_all <- as.matrix(wide_s[, biomarkers, drop = FALSE])
  keep_row <- apply(Y_all, 1, function(r) all(is.finite(r))) & apply(SE_all, 1, function(r) all(is.finite(r) & r > 0))
  age_keep <- wide_e$age_grid[keep_row]
  Y <- Y_all[keep_row, , drop = FALSE]; SE <- SE_all[keep_row, , drop = FALSE]; b <- nrow(Y)
  none <- list(cps = integer(0), level_only = FALSE, age_keep = age_keep, b = b, ok = FALSE)
  if (b < 5) return(none)

  if (cfg$WEIGHT == "sd") {         # paper: each biomarker's curve divided by its SD over ages
    unscale <- apply(Y, 2, stats::sd); unscale[unscale < 1e-8] <- 1
    y_vec <- as.numeric(t(sweep(Y, 2, unscale, "/"))); rw <- rep(1, b * K)
    se_age <- apply(sweep(SE, 2, unscale, "/"), 1, mean)
    pf_d <- 1 / se_age[seq_len(b - 1)]; pf_d <- pf_d / mean(pf_d)
  } else {                          # every estimate divided by its own SE: errors have variance 1
    unscale <- rep(1, K); y_vec <- as.numeric(t(Y / SE)); rw <- as.numeric(t(1 / SE)); pf_d <- rep(1, b - 1)
  }
  a_mid <- mid_all$mid[match(age_keep, mid_all$age_grid)]
  if (cfg$SEARCH %in% c("exact", "ft")) {   # exact: each biomarker separately, best staircase found exactly / ft: line 3 fine-tuned
    n <- b * K; pen <- cfg$BIC_C * log(n); Wm <- 1 / SE^2
    tabs <- lapply(seq_len(K), function(k) seg_exact(Y[, k], Wm[, k], cfg$MAX_JUMPS))      # also used by the reference tables
    ftr <- if (cfg$SEARCH == "ft") ft_line3(Y, SE, age_cost = cfg$BIC_C - 1, jump_cost = 1, eps0 = cfg$FT_EPS0, decay = cfg$FT_DECAY, max_stage = cfg$FT_STAGES) else NULL
    delta <- matrix(0, b - 1, K); B <- matrix(NA_real_, b, K, dimnames = list(age_keep, biomarkers)); rss <- 0; nj <- 0
    for (k in seq_len(K)) {
      cuts_k <- if (is.null(ftr)) seg_pick(tabs[[k]], pen)$cuts else ftr$cuts[[k]]; nj <- nj + length(cuts_k)
      st <- c(1L, cuts_k); en <- c(cuts_k - 1L, b)
      lev <- vapply(seq_along(st), function(i) sum(Wm[st[i]:en[i], k] * Y[st[i]:en[i], k]) / sum(Wm[st[i]:en[i], k]), numeric(1))
      B[, k] <- rep(lev, en - st + 1L); rss <- rss + sum(Wm[, k] * (Y[, k] - B[, k])^2); if (length(cuts_k)) delta[cuts_k - 1L, k] <- diff(lev)
    }
    colnames(delta) <- biomarkers; cp_pos <- which(apply(delta, 1, norm2) > 0)
    sl <- lapply(seq_len(K), function(k) scan_step_line(Y[, k], Wm[, k], a_mid))
    crit <- rss + pen * (nj + K)
    fmt <- function(cu) { has <- vapply(cu, length, integer(1)) > 0
      if (!any(has)) "none" else paste(sprintf("%s %s", biomarkers[has], vapply(cu[has], function(q) paste(age_keep[q], collapse = ","), "")), collapse = " | ") }
    trace <- if (is.null(ftr)) sprintf("Exact search: best staircase has %d jump(s), misfit %.1f, criterion %.3f | CPs (%d): %s", nj, rss, crit, length(cp_pos),
                                       paste(age_keep[cp_pos + 1], collapse = ",")) else
      c(sprintf("Line 3 fine-tuned: sparse group lasso, re-weighted step by step (combined criterion: misfit + log(n) x (%g per change age + 1 per jump))", cfg$BIC_C - 1),
        vapply(ftr$log, function(l) sprintf("  stage %d %-26s alpha %-5s criterion %7.1f | %s", l$stage,
                                            if (l$stage == 0) "(plain sparse group lasso)" else sprintf("(eps %g)", l$eps), if (is.na(l$alpha)) "-" else format(l$alpha), l$crit, fmt(l$cuts)), ""),
        sprintf("Final (last stage): %d change age(s), %d jump(s), misfit %.1f, combined criterion %.1f | CPs (%d): %s", ftr$n_ages, ftr$n_jumps, ftr$misfit, ftr$crit,
                length(cp_pos), paste(age_keep[cp_pos + 1], collapse = ",")))
    return(list(cps = age_keep[cp_pos + 1], level_only = FALSE, jumps = delta[cp_pos, , drop = FALSE], B = B, age_keep = age_keep, a_mid = a_mid, b = b,
                sel = list(crit = crit, alpha = NA, lam = NA, df = nj + K), ok = TRUE, n = n, df_free = K, cands = NULL, tabs = tabs, Y = Y, SE = SE, scan = sl,
                trace = trace, ft = ftr,
                desc = list(rss0 = sum(vapply(tabs, function(t) t[[1]]$rss, numeric(1))), chi_line = vapply(sl, function(x) x$line_gain, numeric(1)))))
  }
  if (cfg$SHAPE == "step") {        # level + jumps: column m is 1 from age point m + 1 on
    A0 <- kronecker(rep(1, b), diag(K))
    J <- matrix(0, b, b - 1); for (t in 2:b) J[t, 1:(t - 1)] <- 1
  } else {                          # level + slope + changes of slope: column j is (age - age of point j + 1), 0 before it
    A0 <- cbind(kronecker(rep(1, b), diag(K)), kronecker(a_mid - a_mid[1], diag(K)))
    J <- vapply(2:(b - 1), function(m) pmax(0, a_mid - a_mid[m]), numeric(b))
    pf_d <- pf_d[seq_len(ncol(J))]
  }
  D0 <- kronecker(J, diag(K)); n_d <- ncol(J)
  XA <- A0 * rw; XD <- D0 * rw
  n <- length(y_vec)
  if (cfg$BASE_PEN) {               # starting level penalized like any other group (paper)
    X <- cbind(XA, XD); y <- y_vec; group <- c(seq_len(K), rep(K + seq_len(n_d), each = K)); pf_g <- c(rep(1, K), pf_d); df_free <- 0
  } else {                          # starting level free: project it out, fit the jumps, recover it afterwards
    P <- XA %*% solve(crossprod(XA), t(XA))
    X <- XD - P %*% XD; y <- as.numeric(y_vec - P %*% y_vec); group <- rep(seq_len(n_d), each = K); pf_g <- pf_d; df_free <- ncol(XA)
  }
  ## the same estimates described without steps: one level per biomarker, or a level plus a straight line in age.
  ## chi_line[k] = how much the misfit drops when biomarker k gets a straight line (1 extra parameter)
  desc <- NULL
  if (cfg$WEIGHT == "se" && !cfg$BASE_PEN && cfg$SHAPE == "step") {
    L <- kronecker(a_mid, diag(K)) * rw; L <- L - P %*% L
    desc <- list(rss0 = sum(y^2), chi_line = as.numeric(crossprod(L, y))^2 / colSums(L^2))
  }
  n_groups <- max(group); df_cap <- cfg$DF_CAP * (b - 1)
  alpha_grid <- seq(cfg$ALPHA_STEP, 1 - cfg$ALPHA_STEP, by = cfg$ALPHA_STEP)

  cache <- new.env()                # refits of candidate jump sets (REFIT): each distinct set is fitted once
  refit <- function(S) {
    key <- if (length(S)) paste(S, collapse = ",") else "none"
    if (!is.null(cache[[key]])) return(cache[[key]])
    out <- if (!length(S)) list(rss = sum(y^2), coef = numeric(0), S = S) else {
      cf <- qr.coef(qr(X[, S, drop = FALSE]), y); cf[is.na(cf)] <- 0
      list(rss = sum((y - X[, S, drop = FALSE] %*% cf)^2), coef = cf, S = S)
    }
    cache[[key]] <- out; out
  }
  pick <- function(pfg, pfs = NULL) {
    best <- list(crit = Inf, alpha = NA, lam = NA, df = NA, coef = NULL)
    for (a in alpha_grid) {
      args <- list(X, y, group = group, asparse = a, pf_group = pfg, intercept = FALSE)
      if (!is.null(pfs)) args$pf_sparse <- pfs
      if (!cfg$STANDARDIZE) args$standardize <- FALSE
      fit <- try(do.call(sparsegl, args), silent = TRUE)
      if (inherits(fit, "try-error")) next
      if (cfg$REFIT) {              # the lasso path proposes sets of jumps; each set is refitted and judged by BIC
        Bm <- as.matrix(fit$beta)
        for (j in c(0L, seq_len(ncol(Bm)))) {
          S <- if (j == 0L) integer(0) else which(abs(Bm[, j]) > 1e-8)
          df <- length(S) + df_free
          if (df > df_cap) next
          r <- refit(S); crit <- r$rss + cfg$BIC_C * log(n) * df
          if (crit < best$crit - 1e-9) {
            cv <- numeric(ncol(X)); cv[S] <- r$coef
            best <- list(crit = crit, alpha = a, lam = if (j == 0L) NA else fit$lambda[j], df = df, coef = cv)
          }
        }
        next
      }
      er <- estimate_risk(fit, X, type = "BIC", approx_df = FALSE)
      df <- er$df + df_free
      crit <- if (cfg$CRIT == "logmse") log(fit$mse) + cfg$BIC_C * log(n) * df / n
              else colSums((y - as.matrix(X %*% fit$beta))^2) + cfg$BIC_C * log(n) * df
      valid <- which(df <= df_cap & is.finite(crit))
      if (!length(valid)) next
      j <- valid[which.min(crit[valid])]
      if (crit[j] < best$crit) {
        cv <- as.numeric(stats::coef(fit, s = er$lambda[j])); if (length(cv) == ncol(X) + 1) cv <- cv[-1]
        best <- list(crit = crit[j], alpha = a, lam = er$lambda[j], df = df[j], coef = cv)
      }
    }
    best
  }
  jumps_of <- function(cf) matrix(if (cfg$BASE_PEN) cf[-seq_len(K)] else cf, nrow = n_d, byrow = TRUE)
  cps_of <- function(cf) age_keep[which(apply(jumps_of(cf), 1, norm2) > 1e-6) + 1]

  trace <- character(0)
  if (cfg$ADAPT_INIT == "sgl") {
    best <- pick(pf_g)
    if (is.null(best$coef)) return(none)
    trace <- c(trace, sprintf("Stage 0: alpha=%.2f lambda=%.4g df=%.1f crit=%.3f", best$alpha, best$lam, best$df, best$crit))
  } else best <- list(coef = ridge_init(X, y), crit = NA, alpha = NA, lam = NA, df = NA)
  prev_cps <- NULL
  for (stage in seq_len(cfg$MAX_STAGES)) {
    apf_g <- vapply(seq_len(n_groups), function(g) {
      ng <- norm2(best$coef[group == g]); if (ng > 1e-8) pf_g[g] / ng^cfg$GAMMA else 1e6 }, numeric(1))
    apf_s <- ifelse(abs(best$coef) > 1e-8, 1 / abs(best$coef)^cfg$GAMMA, 1e6)
    best_new <- pick(apf_g, apf_s)
    if (is.null(best_new$coef)) break
    cps_new <- cps_of(best_new$coef)
    trace <- c(trace, sprintf("Stage %d: alpha=%.2f lambda=%.4g df=%.1f crit=%.3f | CPs (%d): %s", stage, best_new$alpha, best_new$lam,
                              best_new$df, best_new$crit, length(cps_new), paste(cps_new, collapse = ",")))
    change <- norm2(best_new$coef - best$coef) / max(norm2(best$coef), 1e-8)
    same <- identical(sort(cps_new), sort(prev_cps))
    best <- best_new
    if (same && change < 1e-3) break
    prev_cps <- cps_new
  }
  delta <- jumps_of(best$coef)
  alpha_hat <- if (cfg$BASE_PEN) best$coef[seq_len(K)] else as.numeric(solve(crossprod(XA), crossprod(XA, y_vec - XD %*% as.numeric(t(delta)))))
  delta_o <- sweep(delta, 2, unscale, "*")
  B <- sweep(matrix(as.numeric(A0 %*% alpha_hat + D0 %*% as.numeric(t(delta))), nrow = b, byrow = TRUE), 2, unscale, "*")
  dimnames(B) <- list(age_keep, biomarkers)
  cp_pos <- which(apply(delta, 1, norm2) > 1e-6)
  colnames(delta_o) <- biomarkers
  list(cps = age_keep[cp_pos + 1], level_only = cfg$BASE_PEN && all(abs(alpha_hat[seq_len(K)]) < 1e-8), jumps = delta_o[cp_pos, , drop = FALSE],
       B = B, age_keep = age_keep, a_mid = a_mid, b = b, sel = best[c("crit", "alpha", "lam", "df")], trace = trace, ok = TRUE, n = n, df_free = df_free, desc = desc,
       cands = if (cfg$REFIT) lapply(mget(ls(cache), envir = cache), function(r) {
         jm <- numeric(ncol(X)); jm[r$S] <- 1
         list(rss = r$rss, n_jump = length(r$S), cps = age_keep[which(rowSums(matrix(jm, nrow = n_d, byrow = TRUE)) > 0) + 1]) }) else NULL)
}
run_once <- function(d, bins) {
  for (bm in biomarkers) d[[bm]] <- as.numeric(scale(d[[bm]]))
  pa <- step1(d, bins); r <- step2(pa); r$perage <- pa; r
}
## a jump at the second age point with the level shrunk to zero only sets the level (paper method): not a change point
real_cps <- function(r) { v <- r$cps; if (isTRUE(r$level_only) && length(v) && r$b >= 2) v <- setdiff(v, r$age_keep[2]); v }

## ---- run ----------------------------------------------------------------------------------
hdr("Settings")
BINS <- if (cfg$AGE_EST == "bins") make_bins(raw) else NULL
merged <- if (is.null(BINS)) character(0) else vapply(Filter(function(x) length(x) > 1, BINS), function(x) paste0(min(x), "-", max(x)), "")
cat(sprintf("  script: dementia_cp_07_tune.R, version %s\n", SCRIPT_VERSION))
cat(sprintf("  input: %s | N %s, dementia cases %s\n", basename(loc), format(N_ALL, big.mark = ","), format(EV_ALL, big.mark = ",")))
cat(sprintf("  METHOD %s%s | biomarkers: %s\n", cfg$METHOD, if (nzchar(cfg$LABEL)) paste0(" [", cfg$LABEL, "]") else "", paste(biomarkers, collapse = ", ")))
if (length(cfg$MEASURED_ONLY)) cat(sprintf("  only people with a measured value of: %s  (main analysis: N %s, cases %s)\n", paste(cfg$MEASURED_ONLY, collapse = ", "),
                                           format(N_MAIN, big.mark = ","), format(EV_MAIN, big.mark = ",")))
if (length(cfg$ADJUST)) cat(sprintf("  every per-age model is also adjusted for: %s\n", paste(ADJ_NOTE, collapse = ", ")))
cat(sprintf("  a change point is: %s\n", if (cfg$SHAPE == "step") "a jump of the effect to a new level (SHAPE = \"step\")"
            else "a change in the slope of the effect over age (SHAPE = \"slope\")"))
cat(sprintf("  per-age estimates: %s\n", if (cfg$AGE_EST == "kernel") sprintf("Gaussian kernel (at least %g cases, bandwidth up to %g)", cfg$MIN_EVENTS, cfg$H_MAX)
            else sprintf("%d separate age groups (at least %g cases each%s)", length(BINS), cfg$BIN_MIN_EVENTS,
                         if (length(merged)) paste0("; merged: ", paste(merged, collapse = ", ")) else "")))
if (cfg$SEARCH == "ft") cat(sprintf("  fit: weights by SE | starting level free | search: line 3 fine-tuned (sparse group lasso re-weighted step by step; eps %g, x %g per stage, at most %d stages) | criterion: misfit + log(n) x (%g per change age + 1 per jump)\n",
                                    cfg$FT_EPS0, cfg$FT_DECAY, cfg$FT_STAGES, cfg$BIC_C - 1)) else
if (cfg$SEARCH == "exact") cat(sprintf("  fit: weights by SE | starting level free | criterion BIC (known variance) x %g | search: exact, each biomarker on its own (up to %d jumps)\n",
                                       cfg$BIC_C, cfg$MAX_JUMPS)) else
cat(sprintf("  fit: weights %s | starting level %s | sparsegl rescaling %s | criterion %s x %g%s | first weights from %s | adaptive refits %d | search: lasso path\n",
            if (cfg$WEIGHT == "se") "by SE" else "none (curves scaled by SD)", if (cfg$BASE_PEN) "penalized" else "free",
            if (cfg$STANDARDIZE) "on" else "off", if (cfg$CRIT == "known") "BIC (known variance)" else "BIC (log MSE)", cfg$BIC_C,
            if (cfg$REFIT) ", judged on refits without shrinkage" else "", cfg$ADAPT_INIT, cfg$MAX_STAGES))

R0 <- run_once(raw, BINS)
if (!isTRUE(R0$ok)) stop("The fit failed (too few usable age points).", call. = FALSE)
if (cfg$AGE_EST == "bins") {
  hdr("Does each biomarker's effect vary with age at all?")
  cat("  From the separate age groups, before any change-point search. 'same at every age' = test that the effect is constant;\n")
  cat("  'steady trend' = test of a straight-line change with age. Small p = evidence that the effect varies. A biomarker\n")
  cat("  whose effect does not vary cannot have a change point.\n\n")
  cat(sprintf("    %-8s %14s %22s %18s %22s %30s\n", "", "average effect", "same at every age: p", "steady trend: p", "change per 10 years",
              "age where the effect is zero"))
  for (bm in biomarkers) {
    x <- R0$perage[R0$perage$term == bm & is.finite(R0$perage$estimate) & is.finite(R0$perage$se) & R0$perage$se > 0, ]
    w <- 1 / x$se^2; m <- sum(w * x$estimate) / sum(w); Q <- sum(w * (x$estimate - m)^2)
    mid <- (x$age_grid + x$age_to) / 2 + 0.5; am <- sum(w * mid) / sum(w)
    slope <- sum(w * (mid - am) * (x$estimate - m)) / sum(w * (mid - am)^2); se_s <- sqrt(1 / sum(w * (mid - am)^2)); se_m <- sqrt(1 / sum(w))
    p_trend <- 2 * stats::pnorm(-abs(slope / se_s))
    ## straight-line fit: effect(age) = m + slope x (age - am); it is zero at am - m / slope (delta-method 95% interval)
    a0 <- am - m / slope; se0 <- sqrt(se_m^2 / slope^2 + m^2 * se_s^2 / slope^4)
    cross <- if (p_trend < 0.01 && is.finite(a0) && a0 > 45 && a0 < 80) sprintf("%.0f (%.0f to %.0f)", a0, a0 - 1.96 * se0, a0 + 1.96 * se0) else "-"
    pv <- function(p) if (p < 0.0001) "<0.0001" else sprintf("%.4f", p)
    cat(sprintf("    %-8s %14.3f %22s %18s %22s %30s\n", bm, m, pv(stats::pchisq(Q, nrow(x) - 1, lower.tail = FALSE)), pv(p_trend),
                sprintf("%+.3f", 10 * slope), cross))
  }
  cat("\n  'age where the effect is zero': for biomarkers with a clear steady trend (p < 0.01), the age at which a straight line through the\n")
  cat("  estimates crosses zero, with a 95% interval. Before that age the effect has one sign, after it the other.\n")
}
if (cfg$PRINT_PERAGE) {
  hdr("Estimate of every age group (log-odds per SD, standard error in brackets)")
  pa <- R0$perage; ag <- sort(unique(pa$age_grid))
  cat(sprintf("    %-8s%s\n", "ages", paste(sprintf("%17s", biomarkers), collapse = "")))
  for (a in ag) {
    x <- pa[pa$age_grid == a, ]; x <- x[match(biomarkers, x$term), ]
    cat(sprintf("    %-8s%s\n", if (x$age_to[1] > a) paste0(a, "-", x$age_to[1]) else as.character(a),
                paste(ifelse(is.finite(x$estimate), sprintf("%8.3f (%5.3f)", x$estimate, x$se), sprintf("%17s", "-")), collapse = " ")))
  }
}
hdr("Result")
if (cfg$AGE_EST == "bins") { lost <- setdiff(vapply(BINS, min, numeric(1)), R0$age_keep)
  if (length(lost)) cat(sprintf("  NOTE: %d age group(s) could not be fitted and are left out (starting at age %s).\n", length(lost), paste(lost, collapse = ", "))) }
for (l in R0$trace) cat("  ", l, "\n", sep = "")
CP <- real_cps(R0)
cat(sprintf("\n  CHANGE POINTS%s (%d): %s\n", if (cfg$SHAPE == "slope") " IN SLOPE" else "", length(CP), if (length(CP)) paste(CP, collapse = ", ") else "none"))
if (isTRUE(R0$level_only) && length(R0$cps)) cat(sprintf("  (the jump at %d only sets the starting level and is not counted)\n", R0$age_keep[2]))
row_of <- function(a) match(a, R0$age_keep)
if (cfg$SHAPE == "step") {
  if (length(CP)) {
    cat("\n  Jump at each change point (log-odds per SD; 0 = that biomarker does not change there):\n")
    cat(sprintf("    %-6s%s\n", "age", paste(sprintf("%9s", biomarkers), collapse = "")))
    for (a in CP) cat(sprintf("    %-6d%s\n", a, paste(sprintf("%9.3f", R0$jumps[match(a, R0$cps), ]), collapse = "")))
  }
  seg_start <- c(R0$age_keep[1], CP); seg_end <- c(CP - 1, 79)
  val_at <- seg_start; if (isTRUE(R0$level_only) && R0$age_keep[2] %in% R0$cps) val_at[1] <- R0$age_keep[2]   # level after the level-setting jump
  cat("\n  Fitted effect in each age segment (log-odds per SD):\n")
  cat(sprintf("    %-9s%s\n", "ages", paste(sprintf("%9s", biomarkers), collapse = "")))
  for (i in seq_along(seg_start)) cat(sprintf("    %-9s%s\n", paste0(seg_start[i], "-", seg_end[i]),
                                            paste(sprintf("%9.3f", R0$B[row_of(val_at[i]), ]), collapse = "")))
} else {
  if (length(CP)) {
    cat("\n  Change of slope at each change point (log-odds per SD per 10 years of age; 0 = that biomarker's slope does not change there):\n")
    cat(sprintf("    %-6s%s\n", "age", paste(sprintf("%9s", biomarkers), collapse = "")))
    for (a in CP) cat(sprintf("    %-6d%s\n", a, paste(sprintf("%9.3f", 10 * R0$jumps[match(a, R0$cps), ]), collapse = "")))
  }
  knots <- unique(c(R0$age_keep[1], CP, utils::tail(R0$age_keep, 1)))
  cat("\n  Fitted effect at the first age group, at each change point and at the last age group (log-odds per SD):\n")
  cat(sprintf("    %-9s%s\n", "age", paste(sprintf("%9s", biomarkers), collapse = "")))
  for (a in knots) cat(sprintf("    %-9d%s\n", a, paste(sprintf("%9.3f", R0$B[row_of(a), ]), collapse = "")))
  cat("\n  Slope within each segment (change of the effect per 10 years of age):\n")
  cat(sprintf("    %-9s%s\n", "ages", paste(sprintf("%9s", biomarkers), collapse = "")))
  for (i in seq_len(length(knots) - 1)) {
    j1 <- row_of(knots[i]); j2 <- row_of(knots[i + 1])
    cat(sprintf("    %-9s%s\n", paste0(knots[i], "-", knots[i + 1]),
                paste(sprintf("%9.3f", 10 * (R0$B[j2, ] - R0$B[j1, ]) / (R0$a_mid[j2] - R0$a_mid[j1])), collapse = "")))
  }
}

## the same per-age estimates described three ways, all judged by the criterion of this run
DESC <- c(none = NA, line = NA, step = NA); line_bm <- character(0)
if (!is.null(R0$desc) && cfg$REFIT && cfg$CRIT == "known" && is.finite(R0$sel$crit)) {
  pen <- cfg$BIC_C * log(R0$n); use <- R0$desc$chi_line > pen     # a straight line is kept by the same rule as a jump: it must pay for its parameter
  line_bm <- biomarkers[use]
  rws <- list(list("no change with age", K, R0$desc$rss0),
              list(if (any(use)) paste("straight line for", paste(line_bm, collapse = ", ")) else "straight line (none pays for itself)", K + sum(use),
                   R0$desc$rss0 - sum(R0$desc$chi_line[use])))
  if (length(CP)) rws[[3]] <- list(paste("steps at", paste(CP, collapse = ", ")), R0$sel$df, R0$sel$crit - pen * R0$sel$df)
  cat("\n  The same per-age estimates described in different ways. 'misfit' = distance between the description and the estimates;\n")
  cat(sprintf("  criterion = misfit + %.1f for each parameter used (BIC x %g). The description with the smallest criterion is the best.\n", pen, cfg$BIC_C))
  cat(sprintf("    %-44s %10s %10s %11s\n", "description", "parameters", "misfit", "criterion"))
  cr <- vapply(rws, function(r) r[[3]] + pen * r[[2]], numeric(1))
  for (i in seq_along(rws)) cat(sprintf("    %-44s %10.0f %10.1f %11.1f%s\n", rws[[i]][[1]], rws[[i]][[2]], rws[[i]][[3]], cr[i], if (i == which.min(cr)) "  <- best" else ""))
  DESC[c("none", "line", "step")[seq_along(cr)]] <- cr
  cat("  'straight line' = the effect changes by the same amount every year, so there is no particular age at which it changes.\n")
  if (length(CP) && which.min(cr) != 3) cat("  Here the steps are not the best description: they are one way to summarise a change that the straight line describes with fewer parameters.\n")
}

## exact search: the same search at other penalties, with simulated false-alarm rates; and step versus straight line
PEN_TAB <- ""; STEP_LINE <- ""
if (cfg$SEARCH %in% c("exact", "ft")) {
  grid <- sort(unique(c(1.5, 2, 2.5, 3, 4, cfg$BIC_C))); logn <- log(R0$n); Wm <- 1 / R0$SE^2
  cps_at <- function(tabs, cc) lapply(tabs, function(t) R0$age_keep[seg_pick(t, cc * logn)$cuts])
  fa <- rep(NA_real_, length(grid))
  if (cfg$NULL_SIMS > 0) {          # data sets with no change with age: estimate = constant + SE x noise
    old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(20261119); hit <- numeric(length(grid))
    for (s in seq_len(cfg$NULL_SIMS)) {
      tb <- lapply(seq_len(K), function(k) seg_exact(R0$SE[, k] * stats::rnorm(R0$b), Wm[, k], cfg$MAX_JUMPS))
      hit <- hit + vapply(grid, function(cc) any(vapply(tb, function(t) length(seg_pick(t, cc * logn)$cuts) > 0, logical(1))), logical(1))
    }
    fa <- hit / cfg$NULL_SIMS
  }
  cat(sprintf("\n  Change points at other penalty multipliers (same data, %s). 'false alarm' = share of simulated data sets in which\n",
              if (cfg$SEARCH == "ft") "EXACT search, shown for reference; '<- this run' marks the multiplier of this run" else "same exact search"))
  cat(sprintf("  nothing changes with age but the search still reports a change point (%s simulations with the standard errors of this run).\n",
              if (cfg$NULL_SIMS > 0) format(cfg$NULL_SIMS, big.mark = ",") else "no"))
  cat(sprintf("    %-9s %12s   %s\n", "BIC x", "false alarm", "change points"))
  for (i in seq_along(grid)) {
    cc <- cps_at(R0$tabs, grid[i]); has <- vapply(cc, length, integer(1)) > 0
    txt <- if (any(has)) paste(sprintf("%s %s", biomarkers[has], vapply(cc[has], paste, "", collapse = ", ")), collapse = " | ") else "none"
    cat(sprintf("    %-9s %12s   %s%s\n", format(grid[i]), if (is.na(fa[i])) "-" else sprintf("%.1f%%", 100 * fa[i]), txt, if (isTRUE(all.equal(grid[i], cfg$BIC_C))) "   <- this run" else ""))
    PEN_TAB <- paste0(PEN_TAB, if (nzchar(PEN_TAB)) " ; " else "", format(grid[i]), ": ", txt)
  }
  has <- which(vapply(cps_at(R0$tabs, min(grid)), length, integer(1)) > 0)   # every biomarker with a change point at any penalty in the table above
  if (length(has)) {
    cat("\n  A jump, or a gradual change? For each biomarker with a change point in the table above: how much the best single step and a straight line\n")
    cat("  each lower the misfit of 'no change' (both use one extra parameter), and how often a step would beat the line by the observed\n")
    cat("  margin if the truth were that straight line (simulated).\n")
    cat(sprintf("    %-8s %14s %10s %14s %16s %30s\n", "", "best step at", "step", "straight line", "step better by", "chance if truth is the line"))
    nsim <- max(2000L, cfg$NULL_SIMS * 2L); dsum_obs <- 0; dsum_sim <- numeric(nsim); set.seed(20261120)
    pct <- function(p) if (is.na(p)) "-" else if (p < 1 / nsim) sprintf("under %.2f%%", 100 / nsim) else sprintf("%.1f%%", 100 * p)
    for (k in has) {
      sc <- R0$scan[[k]]; d_obs <- sc$step_gain - sc$line_gain; tr <- sc$level + sc$slope * (R0$a_mid - sc$am)
      d_sim <- if (cfg$NULL_SIMS > 0) vapply(seq_len(nsim), function(i) { z <- scan_step_line(tr + R0$SE[, k] * stats::rnorm(R0$b), Wm[, k], R0$a_mid)
        z$step_gain - z$line_gain }, numeric(1)) else rep(NA_real_, nsim)
      pv <- mean(d_sim >= d_obs); dsum_obs <- dsum_obs + d_obs; dsum_sim <- dsum_sim + d_sim
      cat(sprintf("    %-8s %14d %10.1f %14.1f %16.1f %30s\n", biomarkers[k], R0$age_keep[sc$at], sc$step_gain, sc$line_gain, d_obs, pct(pv)))
      STEP_LINE <- paste0(STEP_LINE, if (nzchar(STEP_LINE)) " ; " else "", sprintf("%s: step %.1f line %.1f p=%s", biomarkers[k], sc$step_gain, sc$line_gain, if (is.na(pv)) "-" else sprintf("%.3f", pv)))
    }
    if (length(has) > 1 && cfg$NULL_SIMS > 0) { pj <- mean(dsum_sim >= dsum_obs)
      cat(sprintf("    %-8s %14s %10s %14s %16.1f %30s\n", "together", "", "", "", dsum_obs, pct(pj)))
      STEP_LINE <- paste0(STEP_LINE, sprintf(" ; together p=%.3f", pj)) }
    cat("  reading: a chance under 5% is clear evidence for a jump. Between 5% and 20% the step fits better, but a gradual change is not ruled out.\n")
  }
  if (cfg$NULL_SIMS > 0 && !is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())
}

## figure: per-age estimates with 95% intervals, the fitted steps (solid) and the straight line (dashed)
if (cfg$PRINT_PERAGE && cfg$SHAPE == "step") {
  fig_dir <- path.expand("~/dementia_v9"); dir.create(fig_dir, showWarnings = FALSE)
  fig <- file.path(fig_dir, sprintf("tune_fig_%s_%s.png", if (nzchar(cfg$LABEL)) gsub("[^A-Za-z0-9]+", "-", cfg$LABEL) else cfg$METHOD, format(Sys.time(), "%m%d-%H%M")))
  ok <- tryCatch({
    opened <- tryCatch({ grDevices::png(fig, width = 2400, height = 1500, res = 200, type = "cairo"); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!opened) grDevices::png(fig, width = 2400, height = 1500, res = 200)
    nr <- ceiling(K / 3); graphics::par(mfrow = c(nr, min(K, 3)), mar = c(4, 4.2, 2.6, 0.8), mgp = c(2.3, 0.7, 0), cex = 0.8)
    for (k in seq_len(K)) {
      x <- R0$perage[R0$perage$term == biomarkers[k] & is.finite(R0$perage$estimate) & is.finite(R0$perage$se) & R0$perage$se > 0, ]
      mid <- if (cfg$AGE_EST == "bins") (x$age_grid + x$age_to + 1) / 2 else x$age_grid
      lo <- x$estimate - 1.96 * x$se; hi <- x$estimate + 1.96 * x$se
      yl <- range(c(stats::quantile(lo, 0.03), stats::quantile(hi, 0.97), R0$B[, k]))
      graphics::plot(mid, x$estimate, type = "n", xlim = c(45, 80), ylim = yl, xlab = "age at baseline", ylab = "effect on dementia (log-odds per SD)",
                     main = toupper(biomarkers[k]))
      graphics::abline(h = 0, col = "grey50"); if (length(CP)) graphics::abline(v = CP, col = "grey75", lty = 3)
      graphics::segments(mid, pmax(lo, yl[1] - 1), mid, pmin(hi, yl[2] + 1), col = "grey60")
      graphics::points(mid, x$estimate, pch = 16, cex = 0.7, col = "grey25")
      w <- 1 / x$se^2; lf <- stats::lm(x$estimate ~ mid, weights = w)
      graphics::lines(c(45, 80), stats::predict(lf, data.frame(mid = c(45, 80))), lty = 2, lwd = 1.6, col = "#1f77b4")
      st <- c(R0$age_keep, 80); for (i in seq_along(R0$age_keep)) graphics::segments(st[i], R0$B[i, k], st[i + 1], R0$B[i, k], lwd = 2.2, col = "#d62728")
    }
    if (K < nr * min(K, 3)) { graphics::plot.new()
      graphics::legend("center", bty = "n", cex = 1.05, lwd = c(NA, 2.2, 1.6, 1), lty = c(NA, 1, 2, 3), pch = c(16, NA, NA, NA),
                       col = c("grey25", "#d62728", "#1f77b4", "grey60"),
                       legend = c("estimate of one age group, 95% interval", "fitted steps", "straight line", "change point found")) }
    grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
  cat(if (ok) sprintf("\n  figure saved: %s\n  (open it from the Files pane; it shows estimates only, no counts of people)\n", sub(path.expand("~"), "~", fig, fixed = TRUE))
      else "\n  (the figure could not be drawn on this machine; the table of estimates above has the same numbers)\n")
}

stab_txt <- ""; zero_share <- NA
if (cfg$N_BOOT > 0) {
  hdr(sprintf("Stability: %d resamples (%s) on %d core(s)", cfg$N_BOOT,
              if (cfg$RESAMPLE == "half") "each a random half of the people" else "bootstrap", cfg$N_CORES))
  boot_one <- function(i) { set.seed(20261119 + i)
    idx <- if (cfg$RESAMPLE == "half") sample.int(nrow(raw), floor(nrow(raw) / 2)) else sample.int(nrow(raw), replace = TRUE)
    tryCatch({ r <- run_once(raw[idx, ], BINS); v <- real_cps(r)
               list(cps = v, jumps = r$jumps[match(v, r$cps), , drop = FALSE]) }, error = function(e) NULL) }
  res <- vector("list", cfg$N_BOOT)
  for (s in seq(1, cfg$N_BOOT, by = max(1L, cfg$N_CORES * 2L))) {
    idx <- s:min(cfg$N_BOOT, s + cfg$N_CORES * 2L - 1L)
    res[idx] <- if (cfg$N_CORES > 1) parallel::mclapply(idx, boot_one, mc.cores = cfg$N_CORES, mc.set.seed = FALSE) else lapply(idx, boot_one)
    cat(sprintf("  done %d/%d  (%.1f min)\n", max(idx), cfg$N_BOOT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  full <- Filter(function(r) is.list(r) && !is.null(r$cps), res); res <- lapply(full, `[[`, "cps"); nb <- length(res)
  if (nb) {
    ncp <- vapply(res, length, integer(1)); zero_share <- mean(ncp == 0)
    tab <- table(factor(pmin(ncp, 6), levels = 0:6, labels = c(0:5, "6+")))
    cat(sprintf("  resamples finished: %d | number of change points: %s\n", nb, paste(sprintf("%s: %.0f%%", names(tab), 100 * as.numeric(tab) / nb), collapse = " | ")))
    win <- function(a) mean(vapply(res, function(v) any(abs(v - a) <= 1), logical(1)))
    if (length(CP) && cfg$SEARCH %in% c("exact", "ft")) {      # per biomarker: does the same biomarker jump near the same age in a random half?
      cat("\n  each change point found above, for its own biomarker: share of resamples in which that biomarker jumps within 1 year / within 2 years of it\n")
      wk <- function(k, a, tol) mean(vapply(full, function(r) any(abs(r$cps[abs(r$jumps[, k]) > 1e-8] - a) <= tol), logical(1)))
      st <- character(0)
      for (k in seq_len(K)) for (a in CP[abs(R0$jumps[match(CP, R0$cps), k]) > 1e-8]) {
        cat(sprintf("    %-5s %d   %3.0f%% / %3.0f%%  %s\n", biomarkers[k], a, 100 * wk(k, a, 1), 100 * wk(k, a, 2), strrep("#", round(40 * wk(k, a, 1)))))
        st <- c(st, sprintf("%s %d:%.0f%%", biomarkers[k], a, 100 * wk(k, a, 1)))
      }
      stab_txt <- paste(st, collapse = " ")
    } else if (length(CP)) {
      cat("\n  each change point found above: share of resamples with a change point within 1 year of it\n")
      for (a in CP) cat(sprintf("    age %d   %3.0f%%  %s\n", a, 100 * win(a), strrep("#", round(40 * win(a)))))
      stab_txt <- paste(sprintf("%d:%.0f%%", CP, 100 * vapply(CP, win, numeric(1))), collapse = " ")
    }
    freq <- vapply(46:79, function(a) mean(vapply(res, function(v) a %in% v, logical(1))), numeric(1))
    show <- which(freq >= 0.10)
    cat("\n  ages with a change point in at least 10% of resamples (exact age):\n")
    if (length(show)) for (i in show) cat(sprintf("    age %d   %3.0f%%  %s\n", 45 + i, 100 * freq[i], strrep("#", round(40 * freq[i])))) else cat("    none\n")
    wins <- list("45-49" = 45:49, "50-54" = 50:54, "55-59" = 55:59, "60-64" = 60:64, "65-69" = 65:69, "70-74" = 70:74, "75-79" = 75:79)
    cat("\n  share of resamples in which each biomarker has a jump inside a 5-year window (this allows for a change point\n")
    cat("  landing a year or two apart from one resample to the next):\n")
    cat(sprintf("    %-8s%s\n", "", paste(sprintf("%8s", names(wins)), collapse = "")))
    for (k in seq_along(biomarkers)) cat(sprintf("    %-8s%s\n", biomarkers[k], paste(vapply(wins, function(w)
      sprintf("%7.0f%%", 100 * mean(vapply(full, function(r) any(r$cps %in% w & abs(r$jumps[, k]) > 1e-8), logical(1)))), ""), collapse = "")))
    cat(sprintf("    %-8s%s\n", "any", paste(vapply(wins, function(w) sprintf("%7.0f%%", 100 * mean(vapply(res, function(v) any(v %in% w), logical(1)))), ""), collapse = "")))
    cat("  reading: this table says how firmly the AGES are pinned down, not whether the effects change with age (that is the test table above).\n")
    if (cfg$SEARCH %in% c("exact", "ft")) {
      if (cfg$SEARCH == "ft") cat("  (the reference numbers below were simulated with the exact search)\n")
      cat("  With a strict penalty, a random half of the cohort often does not reach the bar at all, even when the jump is real. In simulations\n")
      cat("  of the v9 size (exact search, BIC x 2.5; share of half-samples with that biomarker's jump within 1 year of the one found):\n")
      cat("    true effects jump (sbp at 53, hdl at 69, sizes as fitted) ..... 0-80%, typically 20-40%\n")
      cat("    true effects change gradually (straight lines) ................. 5-70%, typically 20%\n")
      cat("  The two overlap, so these shares do not settle 'jump or gradual'; the table 'A jump, or a gradual change?' above is the direct check.\n")
    } else if (cfg$BIC_C >= 1.75) {
      cat("  With the default penalty, in simulated cohorts of this size: true jumps of about 0.2 came back in about half or more of the\n")
      cat("  half-samples, false ones in about 10% or fewer.\n")
    } else {
      cat("  With BIC x 1 and the five default biomarkers, in simulated cohorts of this size (share of half-samples with a change point\n")
      cat("  within 1 year of each change point found in the full cohort):\n")
      cat("    true effects jump (the staircase fitted on v9 with BIC x 1) .. 15-80%, typically 40%   (9 cohorts)\n")
      cat("    true effects change gradually, no jump anywhere .............. mostly 10-40%, typically 22%   (16 cohorts; change points were\n")
      cat("                                                                   still reported in 15 of them)\n")
      cat("    true effects do not change with age .......................... no change point reported in 30 of 33 cohorts\n")
      cat("  So with this penalty, shares of 20-30% do not separate real jumps from a gradual change.\n")
      if (!setequal(biomarkers, c("sbp", "bmi", "ldlc", "hdl", "trig")))
        cat("  NOTE: these reference numbers are for the five default biomarkers. With fewer biomarkers each jump is penalized less, so more\n",
            "  change points and higher shares appear even when the effects only change gradually (sbp + hdl alone, gradual truth: 12-83%).\n", sep = "")
    }
  } else cat("  no resample finished\n")
}

## ---- log ------------------------------------------------------------------------------------
out_dir <- path.expand("~/dementia_v9"); dir.create(out_dir, showWarnings = FALSE)
log_file <- file.path(out_dir, "tune_log.csv")
row <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M", tz = "UTC"), label = cfg$LABEL, input = basename(loc), n = N_ALL, cases = EV_ALL,
                  method = cfg$METHOD, shape = cfg$SHAPE, biomarkers = paste(biomarkers, collapse = "+"), age_est = cfg$AGE_EST, weight = cfg$WEIGHT, base_pen = cfg$BASE_PEN,
                  standardize = cfg$STANDARDIZE, crit = cfg$CRIT, refit = cfg$REFIT, adapt_init = cfg$ADAPT_INIT, max_stages = cfg$MAX_STAGES, df_cap = cfg$DF_CAP,
                  min_events = cfg$MIN_EVENTS, h_max = cfg$H_MAX, bin_min_events = cfg$BIN_MIN_EVENTS, bic_c = cfg$BIC_C, gamma = cfg$GAMMA,
                  alpha_step = cfg$ALPHA_STEP, change_points = paste(CP, collapse = " "), n_boot = cfg$N_BOOT, resample = cfg$RESAMPLE, stability = stab_txt,
                  crit_no_change = if (is.na(DESC[["none"]])) "" else sprintf("%.1f", DESC[["none"]]), crit_line = if (is.na(DESC[["line"]])) "" else sprintf("%.1f", DESC[["line"]]),
                  crit_steps = if (is.na(DESC[["step"]])) "" else sprintf("%.1f", DESC[["step"]]), line_for = paste(line_bm, collapse = "+"), version = SCRIPT_VERSION,
                  search = cfg$SEARCH, ft_settings = if (cfg$SEARCH == "ft") sprintf("eps %g, x %g, %d stages", cfg$FT_EPS0, cfg$FT_DECAY, cfg$FT_STAGES) else "", other_penalties = PEN_TAB, step_vs_line = STEP_LINE, adjust = paste(cfg$ADJUST, collapse = "+"), measured_only = paste(cfg$MEASURED_ONLY, collapse = "+"),
                  share_no_cp = if (is.na(zero_share)) "" else sprintf("%.0f%%", 100 * zero_share), stringsAsFactors = FALSE)
old <- if (file.exists(log_file)) tryCatch(utils::read.csv(log_file, stringsAsFactors = FALSE, colClasses = "character"), error = function(e) NULL)
if (!is.null(old) && nrow(old) && !"search" %in% names(old)) old$search <- "sgl"   # runs logged before version 6 all used the lasso-path search
row[] <- lapply(row, as.character)
LOG <- if (!is.null(old) && nrow(old)) dplyr::bind_rows(old, row) else row   # columns added by a newer script version are filled in
LOG[is.na(LOG)] <- ""
utils::write.csv(LOG, log_file, row.names = FALSE)
hdr(sprintf("Log: %d run(s) so far in %s", nrow(LOG), sub(path.expand("~"), "~", log_file, fixed = TRUE)))
tl <- utils::tail(LOG, 8)
dflt <- list(fixed = c(age_est = "bins", weight = "se", base_pen = "FALSE", standardize = "FALSE", crit = "known", refit = "TRUE", adapt_init = "ridge",
                       max_stages = "1", search = "exact", bic_c = as.character(FIXED_BIC_C_EXACT)),
             paper = c(age_est = "kernel", weight = "sd", base_pen = "TRUE", standardize = "TRUE", crit = "logmse", refit = "FALSE", adapt_init = "sgl",
                       max_stages = "5", search = "sgl", bic_c = "1"))
common <- c(shape = "step", df_cap = "1", min_events = "100", h_max = "10", bin_min_events = "40", gamma = "1", alpha_step = "0.05",
            biomarkers = "sbp+bmi+ldlc+hdl+trig", adjust = "", measured_only = "", ft_settings = "eps 1, x 0.5, 6 stages")
for (i in seq_len(nrow(tl))) {
  dd <- c(dflt[[tl$method[i]]], common)
  d <- paste(unlist(lapply(names(dd), function(k) if (nzchar(tl[[k]][i]) && !identical(tl[[k]][i], dd[[k]])) paste0(k, "=", tl[[k]][i]))), collapse = ", ")
  cat(sprintf("  %2d. %-6s%s%s | CPs: %s%s\n", nrow(LOG) - nrow(tl) + i, tl$method[i], if (nzchar(tl$label[i])) paste0(" [", tl$label[i], "]") else "",
              if (nzchar(d)) paste0(" (", d, ")") else "", if (nzchar(tl$change_points[i])) tl$change_points[i] else "none",
              if (tl$n_boot[i] != "0") sprintf(" | stability: %s; no CP in %s of resamples", if (nzchar(tl$stability[i])) tl$stability[i] else "-", tl$share_no_cp[i]) else ""))
}
cat(sprintf("\nruntime %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
TUNE_RESULT <- R0
if (!KEEP) suppressWarnings(rm(list = intersect(.TUNE_NAMES, ls(envir = globalenv())), envir = globalenv()))
