## =============================================================================
## dementia_cp_13_sensitivity.R
## Line 3 (sparse group lasso, fine-tuned) on the paper's five biomarkers: how much do the change points depend on the
## settings of the method?
##
## The main analysis is NOT changed by this script. Its settings were fixed before its results were seen:
##   age groups of at least 40 dementia cases; fine-tuning eps 1, halved at every stage, at most 6 stages; penalty
##   multiplier 2.5 = the smallest with a false alarm of at most 5% in data with no change (rule 3 of
##   dementia_cp_12_panel.R; with the v9 standard errors: 36% at 1.5, 14% at 2, 6.6% at 2.25, 3.8% at 2.5; 500 data sets).
## This script repeats the analysis with other settings, as a sensitivity table. The grids below were written down before
## the run; every setting is reported, none is picked afterwards.
##   A. Age groups: at least 30, 40, 60, 80, 100, 170 cases per group (40 = main analysis; 100 and 170 are the age groups
##      of the 10- and 18-biomarker panels). For each: the multiplier by rule 3 and the change points; also the change
##      points with the multiplier held at 2.5.
##   B. Penalty multiplier (main age groups): 1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4. For each: false alarm and change points.
##   C. Fine-tuning (main age groups, main multiplier): starting eps 0.5, 1, 2 x decay 0.25, 0.5, 0.75; no fine-tuning at
##      all (the plain sparse group lasso); a finer search grid (mixing weight in steps of 0.1, 100 penalty strengths).
##      For each: change points, criterion (lower = better fit for the same penalty) and false alarm.
## Already searched inside line 3 in every run (not a setting to choose by hand): the mixing weight alpha (0, 0.25, 0.5,
## 0.75, 1) and the penalty strength lambda (60 values) - at every stage the pair with the smallest criterion is taken.
## 'false alarm' = share of simulated data sets with no change at all (estimate = standard error x noise) in which the
## method reports at least one change point.
## Per-age models and line 3 are the functions of dementia_cp_07_tune.R (read from TUNE_FILE). Every per-age model is
## adjusted as the paper's adjusted analysis (sex at birth, antihypertensive drugs, imp_bmi, imp_ldlc, imp_hdl, imp_trig).
## Output: model-level only (ages, jumps, rates). A CSV is saved in ~/dementia_v9/.
## Run:  source("~/dementia_cp_13_sensitivity.R")        (about 40 min on 4 cores; progress is printed)
## Options: NULL_SIMS (400), N_CORES, IN_CSV, TUNE_FILE ("~/dementia_cp_07_tune.R").
## =============================================================================

.o <- function(n, d) if (exists(n, envir = globalenv(), inherits = FALSE)) get(n, envir = globalenv()) else d
P <- list(IN_CSV = .o("IN_CSV", ""), TUNE_FILE = .o("TUNE_FILE", "~/dementia_cp_07_tune.R"),
          N_CORES = as.integer(.o("N_CORES", max(1L, parallel::detectCores()))), NULL_SIMS = as.integer(.o("NULL_SIMS", 400)))
suppressWarnings(rm(list = intersect("NULL_SIMS", ls(envir = globalenv())), envir = globalenv()))
## fixed before the run (see the header)
ORIG <- c("sbp", "bmi", "ldlc", "hdl", "trig"); ADJ_ORIG <- c("sex_at_birth", "anti", "imp_bmi", "imp_ldlc", "imp_hdl", "imp_trig")
BIN_GRID <- c(30, 40, 60, 80, 100, 170); BIN_MAIN <- 40; BIC_MAIN <- 2.5
BIC_GRID <- c(1.5, 2, 2.5, 3, 3.5, 4, 5, 6); FA_MAX <- 0.05                  # rule 3, as dementia_cp_12_panel.R
MULT_GRID <- c(1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4)
EPS_GRID <- c(0.5, 1, 2); DECAY_GRID <- c(0.25, 0.5, 0.75); FT_MAIN <- list(eps0 = 1, decay = 0.5, stages = 6L)
ALPHA_MAIN <- c(0, 0.25, 0.5, 0.75, 1); NLAM_MAIN <- 60L

t0 <- Sys.time(); options(width = 170)
for (p in c("dplyr", "tidyr", "purrr", "tibble", "sparsegl")) if (!requireNamespace(p, quietly = TRUE)) utils::install.packages(p, quiet = TRUE)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(tibble); library(sparsegl) })
hdr <- function(s) cat("\n==================== ", s, " ====================\n", sep = "")
mins <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

## ---- the functions of 07 ---------------------------------------------------------------------------
tf <- path.expand(P$TUNE_FILE)
if (!file.exists(tf)) stop("dementia_cp_07_tune.R not found at ", tf, " (set TUNE_FILE).", call. = FALSE)
E <- new.env(parent = globalenv())
WANT <- c("form_for", "fit_kernel", "make_bins", "fit_bin", "step1", "ft_line3")
for (e in parse(tf)) if (is.call(e) && is.name(e[[1]]) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) && as.character(e[[2]]) %in% WANT) eval(e, E)
if (length(setdiff(WANT, ls(E)))) stop("These functions are missing in ", tf, ": ", paste(setdiff(WANT, ls(E)), collapse = ", "), call. = FALSE)
V07 <- sub('^SCRIPT_VERSION <- "(.*)".*$', "\\1", grep("^SCRIPT_VERSION <-", readLines(tf), value = TRUE)[1])
set_run <- function(bin_min) {
  E$biomarkers <- ORIG; E$K <- length(ORIG); E$form_txt <- paste("Dementia ~", paste(ORIG, collapse = " + "))
  E$cfg <- list(AGE_EST = "bins", ADJUST = ADJ_ORIG, BIN_MIN_EVENTS = bin_min, MIN_EVENTS = 100, H_MAX = 10)
}

## ---- data -----------------------------------------------------------------------------------------
loc <- if (nzchar(P$IN_CSV)) path.expand(P$IN_CSV) else {
  f <- sort(Sys.glob(path.expand("~/dementia_v9/dementia_v9_analysis_*.csv")), decreasing = TRUE); if (length(f)) f[1] else "" }
if (!nzchar(loc) || !file.exists(loc)) stop("Analysis file not found: run dementia_cp_02_v9_cohort.R first (or set IN_CSV).", call. = FALSE)
A0 <- utils::read.csv(loc, stringsAsFactors = FALSE)
miss <- setdiff(c("Dementia", "age", ORIG, ADJ_ORIG), names(A0))
if (length(miss)) stop("The analysis file has no column(s): ", paste(miss, collapse = ", "), call. = FALSE)

prep <- function(bms, adjust) {     # as 07: complete biomarkers; text adjusters -> groups (missing = own group, < 200 pooled)
  d <- A0[stats::complete.cases(A0[, bms, drop = FALSE]), unique(c("Dementia", "age", bms, adjust)), drop = FALSE]
  for (v in adjust) {
    if (is.numeric(d[[v]]) || is.logical(d[[v]])) { d[[v]] <- as.numeric(d[[v]]); d <- d[is.finite(d[[v]]), , drop = FALSE]
    } else { x <- as.character(d[[v]]); x[is.na(x) | !nzchar(x)] <- "unknown"; tb <- table(x); x[x %in% names(tb)[tb < 200]] <- "other (pooled)"
      d[[v]] <- factor(x) }
  }
  d
}
perage <- function(d, bins) { for (bm in E$biomarkers) d[[bm]] <- as.numeric(scale(d[[bm]])); E$step1(d, bins) }   # as 07 run_once
mats <- function(pa) {              # as 07 step2: age groups x biomarkers; an age group with any estimate missing is left out
  bms <- E$biomarkers; ag <- sort(unique(pa$age_grid))
  grab <- function(col) vapply(bms, function(b) { x <- pa[pa$term == b, ]; x[[col]][match(ag, x$age_grid)] }, numeric(length(ag)))
  Y <- matrix(grab("estimate"), length(ag)); SE <- matrix(grab("se"), length(ag)); colnames(Y) <- colnames(SE) <- bms
  keep <- apply(Y, 1, function(r) all(is.finite(r))) & apply(SE, 1, function(r) all(is.finite(r) & r > 0))
  list(Y = Y[keep, , drop = FALSE], SE = SE[keep, , drop = FALSE], age = ag[keep], lost = ag[!keep])
}
## line 3; the defaults are those of 07 (and of the main analysis)
ft <- function(Y, SE, bic, eps0 = FT_MAIN$eps0, decay = FT_MAIN$decay, stages = FT_MAIN$stages, alphas = ALPHA_MAIN, nlambda = NLAM_MAIN)
  E$ft_line3(Y, SE, age_cost = bic - 1, jump_cost = 1, eps0 = eps0, decay = decay, max_stage = stages, alphas = alphas, nlambda = nlambda)
jumps_of <- function(M, cuts) lapply(seq_len(ncol(M$Y)), function(k) { cu <- cuts[[k]]; b <- nrow(M$Y); st <- c(1L, cu); en <- c(cu - 1L, b); w <- 1 / M$SE[, k]^2
  lev <- vapply(seq_along(st), function(i) sum(w[st[i]:en[i]] * M$Y[st[i]:en[i], k]) / sum(w[st[i]:en[i]]), numeric(1))
  if (length(cu)) data.frame(age = M$age[cu], jump = diff(lev)) else data.frame(age = numeric(0), jump = numeric(0)) })
## a jump at the second age group only moves the first group, one at the last group only the last: one age group on its own, marked *
edge_mark <- function(a, M) ifelse(a %in% c(M$age[2], M$age[length(M$age)]), "*", "")
cp_txt <- function(M, cuts) { J <- jumps_of(M, cuts); has <- vapply(J, nrow, integer(1)) > 0
  if (!any(has)) "none" else paste(vapply(which(has), function(k) sprintf("%s %s", ORIG[k], paste(sprintf("%d%s (%+.2f)", J[[k]]$age, edge_mark(J[[k]]$age, M), J[[k]]$jump),
                                                                                                    collapse = ", ")), ""), collapse = " | ") }
false_alarm <- function(SE, bic, ...) mean(unlist(parallel::mclapply(seq_len(P$NULL_SIMS), function(i) {   # data with no change: estimate = SE x noise
  set.seed(20261121 + i); r <- ft(SE * matrix(stats::rnorm(length(SE)), nrow(SE)), SE, bic, ...); any(lengths(r$cuts) > 0) },
  mc.cores = P$N_CORES, mc.set.seed = FALSE)))
choose_bic <- function(SE) {        # rule 3, as dementia_cp_12_panel.R
  tab <- data.frame(bic = numeric(0), fa = numeric(0))
  for (b in BIC_GRID) { tab <- rbind(tab, data.frame(bic = b, fa = false_alarm(SE, b))); if (tail(tab$fa, 1) <= FA_MAX) break }
  i <- match(tail(tab$bic, 1), BIC_GRID)
  if (tail(tab$fa, 1) <= FA_MAX && i > 1) { mid <- (BIC_GRID[i] + BIC_GRID[i - 1]) / 2
    tab <- rbind(tab, data.frame(bic = mid, fa = false_alarm(SE, mid))) }
  tab <- tab[order(tab$bic), ]
  pick <- if (any(tab$fa <= FA_MAX)) min(tab$bic[tab$fa <= FA_MAX]) else max(BIC_GRID)
  list(bic = pick, fa = tab$fa[tab$bic == pick], tab = tab)
}
pc <- function(x) sprintf("%.1f%%", 100 * x)
OUT <- list(); add <- function(...) OUT[[length(OUT) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

hdr("Settings")
set_run(BIN_MAIN); D <- prep(ORIG, ADJ_ORIG)
cat(sprintf("  script: dementia_cp_13_sensitivity.R, version 1 (2026-11-21; grids fixed before the run) | input: %s | N %s, dementia cases %s\n",
            basename(loc), format(nrow(D), big.mark = ","), format(sum(D$Dementia), big.mark = ",")))
cat(sprintf("  functions of dementia_cp_07_tune.R version %s | biomarkers: %s | adjusted for: %s | %d null data sets per false-alarm rate | %d core(s)\n",
            V07, paste(ORIG, collapse = ", "), paste(ADJ_ORIG, collapse = ", "), P$NULL_SIMS, P$N_CORES))
cat("  '<- main' marks the setting of the main analysis. Jumps in log-odds per SD.\n")

## ---- A. age groups ----------------------------------------------------------------------------------
hdr("A. Age groups (how many cases each group must have)")
cat(sprintf("    %-10s %-11s %-26s %s\n", "min cases", "age groups", "multiplier (false alarm)", "change points   [in brackets: with the multiplier held at 2.5, where rule 3 gives another one]"))
MA <- list()
for (bm in BIN_GRID) {
  set_run(bm); bins <- E$make_bins(D); M <- mats(perage(D, bins)); cb <- choose_bic(M$SE)
  r <- ft(M$Y, M$SE, cb$bic); r25 <- if (isTRUE(all.equal(cb$bic, 2.5))) r else ft(M$Y, M$SE, 2.5)
  MA[[as.character(bm)]] <- list(M = M, cb = cb, r = r)
  same <- isTRUE(all.equal(cb$bic, 2.5))
  cat(sprintf("    %-10d %-11d %-26s %s%s%s  (%s)\n", bm, nrow(M$Y), sprintf("%g (%s)", cb$bic, pc(cb$fa)), cp_txt(M, r$cuts),
              if (same) "" else sprintf("   [at 2.5: %s]", cp_txt(M, r25$cuts)), if (bm == BIN_MAIN) "   <- main" else "", mins()))
  add(part = "A age groups", setting = sprintf("at least %d cases", bm), main = bm == BIN_MAIN, age_groups = nrow(M$Y), multiplier = cb$bic, false_alarm = cb$fa,
      change_points = cp_txt(M, r$cuts), criterion = r$crit, stages = length(r$log) - 1L, change_points_at_2.5 = cp_txt(M, r25$cuts))
}
cat("  coarser age groups place a change less exactly (a jump can only sit between two groups) and average over more ages.\n")
cat(sprintf("  the main analysis uses the multiplier 2.5 (set earlier from 500 simulated data sets); with %d data sets per rate, rule 3 can land\n", P$NULL_SIMS))
cat("  one step beside it, and the brackets then show the result at 2.5.\n")

## ---- B. penalty multiplier ----------------------------------------------------------------------------
hdr(sprintf("B. Penalty multiplier (age groups of at least %d cases)", BIN_MAIN))
M0 <- MA[[as.character(BIN_MAIN)]]$M; CB0 <- MA[[as.character(BIN_MAIN)]]$cb; BIC0 <- BIC_MAIN; FA_MAIN <- NA_real_
cat(sprintf("    %-11s %-13s %s\n", "multiplier", "false alarm", "change points"))
for (b in sort(unique(c(MULT_GRID, BIC0)))) {
  hit <- which(abs(CB0$tab$bic - b) < 1e-9); fa <- if (length(hit)) CB0$tab$fa[hit[1]] else false_alarm(M0$SE, b)
  r <- ft(M0$Y, M0$SE, b)
  if (isTRUE(all.equal(b, BIC0))) FA_MAIN <- fa
  cat(sprintf("    %-11g %-13s %s%s\n", b, pc(fa), cp_txt(M0, r$cuts), if (isTRUE(all.equal(b, BIC0))) "   <- main" else ""))
  add(part = "B multiplier", setting = sprintf("multiplier %g", b), main = isTRUE(all.equal(b, BIC0)), age_groups = nrow(M0$Y), multiplier = b, false_alarm = fa,
      change_points = cp_txt(M0, r$cuts), criterion = r$crit, stages = length(r$log) - 1L, change_points_at_2.5 = NA_character_)
}
cat(sprintf("  the false-alarm rates here come from %d simulated data sets each, so they move by about 1 point from run to run.\n", P$NULL_SIMS))
cat(sprintf("  A smaller multiplier finds more change points, but more of them are false; a larger one keeps only the strongest. (%s)\n", mins()))

## ---- C. fine-tuning -----------------------------------------------------------------------------------
hdr(sprintf("C. Fine-tuning of line 3 (age groups of at least %d cases, multiplier %g)", BIN_MAIN, BIC0))
cat(sprintf("    %-44s %-7s %-10s %-12s %s\n", "setting", "stages", "criterion", "false alarm", "change points"))
rowC <- function(label, main = FALSE, ...) {
  r <- ft(M0$Y, M0$SE, BIC0, ...); fa <- if (main) FA_MAIN else false_alarm(M0$SE, BIC0, ...)
  cat(sprintf("    %-44s %-7d %-10.1f %-12s %s%s\n", label, length(r$log) - 1L, r$crit, pc(fa), cp_txt(M0, r$cuts), if (main) "   <- main" else ""))
  add(part = "C fine-tuning", setting = label, main = main, age_groups = nrow(M0$Y), multiplier = BIC0, false_alarm = fa, change_points = cp_txt(M0, r$cuts),
      criterion = r$crit, stages = length(r$log) - 1L, change_points_at_2.5 = NA_character_)
}
rowC("no fine-tuning (plain sparse group lasso)", stages = 0L)
for (e0 in EPS_GRID) for (dc in DECAY_GRID)
  rowC(sprintf("eps %g, x %g at every stage", e0, dc), main = isTRUE(all.equal(e0, FT_MAIN$eps0)) && isTRUE(all.equal(dc, FT_MAIN$decay)), eps0 = e0, decay = dc)
rowC("finer search grid (alpha by 0.1, 100 lambdas)", alphas = seq(0, 1, by = 0.1), nlambda = 100L)
cat("  criterion = misfit + penalty; with the same data and multiplier, a lower value is a better staircase. Settings that reach the\n")
cat("  same criterion and the same change points do not matter for the result.\n")

## ---- summary ------------------------------------------------------------------------------------------
hdr("Summary")
OUT <- do.call(rbind, OUT)
main_cp <- OUT$change_points[OUT$part == "B multiplier" & OUT$main]
ages_of <- function(txt, b) { m <- regmatches(txt, regexpr(paste0("(^| )", b, " [^|]*"), txt)); if (!length(m)) integer(0) else as.integer(regmatches(m, gregexpr("[0-9]+(?=\\*? \\()", m, perl = TRUE))[[1]]) }
cat(sprintf("  main analysis: %s\n", main_cp))
for (b in ORIG) { a0 <- ages_of(main_cp, b); if (!length(a0)) next
  for (a in a0) { near <- function(part) { x <- OUT$change_points[OUT$part == part]; sprintf("%d of %d", sum(vapply(x, function(t) any(abs(ages_of(t, b) - a) <= 2), logical(1))), length(x)) }
    cat(sprintf("    %s %d (within 2 years): age groups %s | multipliers %s | fine-tuning settings %s\n", b, a, near("A age groups"), near("B multiplier"), near("C fine-tuning"))) } }
other <- unique(unlist(lapply(seq_len(nrow(OUT)), function(i) unlist(lapply(ORIG, function(b) { a <- ages_of(OUT$change_points[i], b); a0 <- ages_of(main_cp, b)
  a <- a[!vapply(a, function(x) any(abs(x - a0) <= 2), logical(1))]; if (length(a)) paste(b, a) else NULL })))))
cat(sprintf("  change points that appear only with other settings: %s\n", if (length(other)) paste(sort(other), collapse = ", ") else "none"))
if (any(grepl("*", OUT$change_points, fixed = TRUE))) cat("  * = the jump moves only the first or only the last age group: one group on its own, read with caution.\n")
fcsv <- file.path(path.expand("~/dementia_v9"), sprintf("sensitivity_%s.csv", format(Sys.time(), "%Y%m%d-%H%M")))
dir.create(dirname(fcsv), showWarnings = FALSE, recursive = TRUE); utils::write.csv(OUT, fcsv, row.names = FALSE)
SENS_RESULT <- OUT
cat(sprintf("\n  saved: %s (model-level numbers only) | total %s\n", basename(fcsv), mins()))
