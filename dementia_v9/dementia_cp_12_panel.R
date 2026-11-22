## =============================================================================
## dementia_cp_12_panel.R
## Which biomarkers have an effect on dementia that changes with age, and at what age? Line 3 (fine-tuned) on a
## wider panel of biomarkers.
##
## Input: the analysis file after dementia_cp_11_extra_labs.R (first batch dbp, hr, a1c, glu, egfr, hgb, alb; second batch below).
## The rules below were written down before any change point of the new biomarkers was looked at. The script
## applies them in this order and prints each decision:
##   1. Panel. The five original biomarkers always stay. New biomarkers are taken in the order dbp, hr, a1c, glu,
##      egfr, hgb, alb; one enters the main analysis if (a) it is measured (not filled) in at least 40% of the
##      analysis set, and (b) its measured values are not strongly correlated (|r| > 0.5) with a biomarker already in
##      the panel: two strongly correlated biomarkers in one per-age model split their effect between them, and both
##      estimates lose much of their precision (expected case: diastolic with systolic blood pressure). Hemoglobin
##      A1c, or glucose if A1c falls short (never both: they measure the same thing). A biomarker left out by (a) or
##      (b) is still analysed in the one-at-a-time runs (rule 4). Total cholesterol is not used (about LDL + HDL + TG/5).
##   2. Age groups. With about twice as many biomarkers, each per-age model has about twice as many coefficients,
##      so neighbouring ages are merged until each group has at least 100 dementia cases (40 in the paper's
##      five-biomarker analysis).
##   3. Penalty. The same standard as the paper's analysis: the multiplier is the smallest (search below) at which
##      simulated data sets with no change at all (built from the standard errors of the run,
##      nothing else) show a change point in at most 5% of cases. With the five biomarkers and the v9 standard errors
##      this gives 2.5, the value used so far (false alarm 36% at 1.5, 14% at 2, 6.6% at 2.25, 3.8% at 2.5; 500 data
##      sets each). With more biomarkers the same multiplier can be stricter or looser, so it is recalculated in every
##      run. It is fixed before the search on the data. (Search: 1.5, 2, 2.5, 3, 3.5, 4, 5, 6 until the first value
##      at or under 5%, then the point halfway back to the value before it.)
##   4. Main result = line 3 on the whole panel. Robustness: (a) the five original biomarkers plus ONE new biomarker
##      at a time, set up as the paper's main analysis (40-case age groups), for every new biomarker, also those
##      under 40%; (b) the five original biomarkers alone with the main panel's age groups (shows what the coarser
##      age groups change on their own). The penalty follows rule 3 in every run.
## Every per-age model is adjusted for sex at birth, antihypertensive drugs and the 'filled' flags of its biomarkers
## (as the paper's adjusted analysis: imp_bmi, imp_ldlc, imp_hdl, imp_trig; plus imp_<new> for each new biomarker).
## Per-age models and the line 3 search are the functions of dementia_cp_07_tune.R (read from TUNE_FILE), so the
## numbers are those 07 gives with the same settings. Stability: random halves of the people, as in 07.
## Output: model-level numbers only (coefficients, p-values, ages; no counts of people under 21). A figure and a CSV
## are saved in ~/dementia_v9/.
##
## Second batch (BATCH = 2, the default since 2026-11-21; rules fixed before any of its results were seen): bun, na, ca, alt,
## ast, alp, bili, wbc, plt, rdw, mcv, tsh are considered after the first batch, in this order, with the same rule 1. The main
## panel is then the full panel; its age groups keep the cases per coefficient of the 10-biomarker panel (100 cases there), so a
## per-age model with more coefficients gets larger age groups (rule 2). One-at-a-time runs: the second batch only (the first
## batch was done with BATCH = 1, which reproduces the first run).
## Run:  source("~/dementia_cp_12_panel.R")            (BATCH = 2: about 1.5-2 h on 4 cores; progress is printed)
## Options (set before source):  N_HALF (50), N_HALF_ONE (20), NULL_SIMS (400), RUN_ONE (TRUE), RUN_BRIDGE (TRUE),
##   N_CORES, IN_CSV, TUNE_FILE ("~/dementia_cp_07_tune.R").
## =============================================================================

.o <- function(n, d) if (exists(n, envir = globalenv(), inherits = FALSE)) get(n, envir = globalenv()) else d
P <- list(IN_CSV = .o("IN_CSV", ""), TUNE_FILE = .o("TUNE_FILE", "~/dementia_cp_07_tune.R"),
          N_CORES = as.integer(.o("N_CORES", max(1L, parallel::detectCores()))), NULL_SIMS = as.integer(.o("NULL_SIMS", 400)),
          N_HALF = as.integer(.o("N_HALF", 50)), N_HALF_ONE = as.integer(.o("N_HALF_ONE", 20)),
          RUN_ONE = isTRUE(.o("RUN_ONE", TRUE)), RUN_BRIDGE = isTRUE(.o("RUN_BRIDGE", TRUE)), BATCH = as.integer(.o("BATCH", 2)))
suppressWarnings(rm(list = intersect(c("N_HALF", "N_HALF_ONE", "NULL_SIMS", "RUN_ONE", "RUN_BRIDGE", "BATCH"), ls(envir = globalenv())), envir = globalenv()))
## fixed by the rules in the header
MIN_MEASURED <- 0.40; MAX_R <- 0.5; BIN_MAIN <- 100; BIN_ONE <- 40; BIC_GRID <- c(1.5, 2, 2.5, 3, 3.5, 4, 5, 6); FA_MAX <- 0.05
ORIG <- c("sbp", "bmi", "ldlc", "hdl", "trig"); NEWB1 <- c("dbp", "hr", "a1c", "glu", "egfr", "hgb", "alb")
NEWB2 <- c("bun", "na", "ca", "alt", "ast", "alp", "bili", "wbc", "plt", "rdw", "mcv", "tsh")
NEWB <- if (P$BATCH == 1) NEWB1 else c(NEWB1, NEWB2)
ADJ_ORIG <- c("sex_at_birth", "anti", "imp_bmi", "imp_ldlc", "imp_hdl", "imp_trig")
LBL <- c(sbp = "systolic BP", bmi = "BMI", ldlc = "LDL cholesterol", hdl = "HDL cholesterol", trig = "triglycerides",
         dbp = "diastolic BP", hr = "heart rate", a1c = "hemoglobin A1c", glu = "glucose", egfr = "eGFR", hgb = "hemoglobin", alb = "albumin",
         bun = "urea nitrogen", na = "sodium", ca = "calcium", alt = "ALT", ast = "AST", alp = "alk. phosphatase", bili = "bilirubin",
         wbc = "white cells", plt = "platelets", rdw = "RDW", mcv = "MCV", tsh = "TSH")

t0 <- Sys.time(); options(width = 170)
for (p in c("dplyr", "tidyr", "purrr", "tibble", "sparsegl")) if (!requireNamespace(p, quietly = TRUE)) utils::install.packages(p, quiet = TRUE)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(tibble); library(sparsegl) })
hdr <- function(s) cat("\n==================== ", s, " ====================\n", sep = "")
pv <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.0001, "<0.0001", sprintf("%.4f", p)))
mins <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

## ---- the functions of 07 ---------------------------------------------------------------------------
tf <- path.expand(P$TUNE_FILE)
if (!file.exists(tf)) stop("dementia_cp_07_tune.R not found at ", tf, " (set TUNE_FILE).", call. = FALSE)
E <- new.env(parent = globalenv())
WANT <- c("form_for", "fit_kernel", "make_bins", "fit_bin", "step1", "ft_line3")
for (e in parse(tf)) if (is.call(e) && is.name(e[[1]]) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) && as.character(e[[2]]) %in% WANT) eval(e, E)
if (length(setdiff(WANT, ls(E)))) stop("These functions are missing in ", tf, ": ", paste(setdiff(WANT, ls(E)), collapse = ", "), call. = FALSE)
V07 <- sub('^SCRIPT_VERSION <- "(.*)".*$', "\\1", grep("^SCRIPT_VERSION <-", readLines(tf), value = TRUE)[1])
set_run <- function(bms, adjust, bin_min) {
  E$biomarkers <- bms; E$K <- length(bms); E$form_txt <- paste("Dementia ~", paste(bms, collapse = " + "))
  E$cfg <- list(AGE_EST = "bins", ADJUST = adjust, BIN_MIN_EVENTS = bin_min, MIN_EVENTS = 100, H_MAX = 10)
}

## ---- data -----------------------------------------------------------------------------------------
loc <- if (nzchar(P$IN_CSV)) path.expand(P$IN_CSV) else {
  f <- sort(Sys.glob(path.expand("~/dementia_v9/dementia_v9_analysis_*.csv")), decreasing = TRUE); if (length(f)) f[1] else "" }
if (!nzchar(loc) || !file.exists(loc)) stop("Analysis file not found: run dementia_cp_02_v9_cohort.R and dementia_cp_11_extra_labs.R first.", call. = FALSE)
A0 <- utils::read.csv(loc, stringsAsFactors = FALSE)
miss <- setdiff(c("Dementia", "age", "time2", ORIG, ADJ_ORIG), names(A0))
if (length(miss)) stop("The analysis file has no column(s): ", paste(miss, collapse = ", "), call. = FALSE)
have_new <- NEWB[NEWB %in% names(A0) & paste0("imp_", NEWB) %in% names(A0)]
if (!length(have_new)) stop("The analysis file has none of the new biomarkers: run dementia_cp_11_extra_labs.R first.", call. = FALSE)
if (P$BATCH == 2 && !any(NEWB2 %in% have_new)) stop("The analysis file has none of the second batch: run dementia_cp_11_extra_labs.R (version 1.6) first.", call. = FALSE)
if (P$BATCH == 2 && length(setdiff(NEWB2, have_new))) cat("NOTE: not in the analysis file, left out:", paste(setdiff(NEWB2, have_new), collapse = ", "), "\n")
A0 <- A0 %>% filter(is.finite(time2), time2 > 0, is.finite(age), age >= 45, age < 80, ldlc < 400, trig < 700)   # as 07
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
  to <- vapply(ag, function(a) pa$age_to[pa$age_grid == a][1], numeric(1))
  keep <- apply(Y, 1, function(r) all(is.finite(r))) & apply(SE, 1, function(r) all(is.finite(r) & r > 0))
  list(Y = Y[keep, , drop = FALSE], SE = SE[keep, , drop = FALSE], age = ag[keep], to = to[keep], mid = (ag[keep] + to[keep] + 1) / 2, lost = ag[!keep])
}
ft <- function(Y, SE, bic) E$ft_line3(Y, SE, age_cost = bic - 1, jump_cost = 1, eps0 = 1, decay = 0.5, max_stage = 6)   # 07's defaults
levels_of <- function(M, cuts) {    # refit without shrinkage: weighted mean in each segment
  b <- nrow(M$Y); B <- M$Y * NA; J <- list()
  for (k in seq_len(ncol(M$Y))) { cu <- cuts[[k]]; st <- c(1L, cu); en <- c(cu - 1L, b); w <- 1 / M$SE[, k]^2
    lev <- vapply(seq_along(st), function(i) sum(w[st[i]:en[i]] * M$Y[st[i]:en[i], k]) / sum(w[st[i]:en[i]]), numeric(1))
    B[, k] <- rep(lev, en - st + 1L); J[[k]] <- if (length(cu)) data.frame(age = M$age[cu], jump = diff(lev)) else data.frame(age = numeric(0), jump = numeric(0)) }
  list(B = B, J = J)
}
tests <- function(M) {              # as 07 'Does each biomarker's effect vary with age at all?'
  do.call(rbind, lapply(seq_len(ncol(M$Y)), function(k) {
    x <- M$Y[, k]; w <- 1 / M$SE[, k]^2; m <- sum(w * x) / sum(w); se_m <- sqrt(1 / sum(w)); Q <- sum(w * (x - m)^2)
    am <- sum(w * M$mid) / sum(w); sxx <- sum(w * (M$mid - am)^2); sl <- sum(w * (M$mid - am) * (x - m)) / sxx
    data.frame(biomarker = colnames(M$Y)[k], effect = m, se = se_m, p_effect = 2 * stats::pnorm(-abs(m / se_m)),
               p_vary = stats::pchisq(Q, length(x) - 1, lower.tail = FALSE), slope10 = 10 * sl, p_trend = 2 * stats::pnorm(-abs(sl / sqrt(1 / sxx))))
  }))
}
par_map <- function(n, f) unlist(parallel::mclapply(seq_len(n), f, mc.cores = P$N_CORES, mc.set.seed = FALSE), recursive = FALSE)
false_alarm <- function(SE, bic, nsim) mean(unlist(parallel::mclapply(seq_len(nsim), function(i) {   # data with no change: estimate = SE x noise
  set.seed(20261121 + i); r <- ft(SE * matrix(stats::rnorm(length(SE)), nrow(SE)), SE, bic); any(lengths(r$cuts) > 0) },
  mc.cores = P$N_CORES, mc.set.seed = FALSE)))
choose_bic <- function(SE, show = TRUE) {      # rule 3
  tab <- data.frame(bic = numeric(0), fa = numeric(0))
  for (b in BIC_GRID) { tab <- rbind(tab, data.frame(bic = b, fa = false_alarm(SE, b, P$NULL_SIMS))); if (tail(tab$fa, 1) <= FA_MAX) break }
  i <- match(tail(tab$bic, 1), BIC_GRID)
  if (tail(tab$fa, 1) <= FA_MAX && i > 1) { mid <- (BIC_GRID[i] + BIC_GRID[i - 1]) / 2           # halfway back to the value before
    tab <- rbind(tab, data.frame(bic = mid, fa = false_alarm(SE, mid, P$NULL_SIMS))) }
  tab <- tab[order(tab$bic), ]
  pick <- if (any(tab$fa <= FA_MAX)) min(tab$bic[tab$fa <= FA_MAX]) else max(BIC_GRID)
  if (show) cat(sprintf("    %s   -> multiplier %g%s\n", paste(sprintf("x %g: %.1f%%", tab$bic, 100 * tab$fa), collapse = " | "), pick,
                        if (all(tab$fa > FA_MAX)) "  (NOTE: none of the grid reaches 5%; the largest is used)" else ""))
  list(bic = pick, fa = tab$fa[tab$bic == pick], tab = tab)
}
halves <- function(d, bins, bic, n) {          # as 07 RESAMPLE = "half": a random half of the people, same age groups
  if (n <= 0) return(list())
  par_map(n, function(i) { set.seed(20261119 + i)
    list(tryCatch({ M <- mats(perage(d[sample.int(nrow(d), floor(nrow(d) / 2)), , drop = FALSE], bins)); r <- ft(M$Y, M$SE, bic)
                    lapply(r$cuts, function(cu) M$age[cu]) }, error = function(e) NULL)) })
}
near_share <- function(hs, k, a, starts) {     # share of halves in which biomarker k jumps in the same or a neighbouring age group
  hs <- Filter(Negate(is.null), hs); if (!length(hs)) return(NA_real_)
  ia <- match(a, starts); mean(vapply(hs, function(h) any(abs(match(h[[k]], starts) - ia) <= 1), logical(1)))
}
any_share <- function(hs, k) { hs <- Filter(Negate(is.null), hs); if (!length(hs)) NA_real_ else mean(vapply(hs, function(h) length(h[[k]]) > 0, logical(1))) }
## a jump at the second age group only moves the first group, one at the last group only the last: one age group on its own, marked *
edge_mark <- function(a, M) ifelse(a %in% c(M$age[2], M$age[length(M$age)]), "*", "")
cp_txt <- function(J, bms, M = NULL, digits = 2) { has <- vapply(J, nrow, integer(1)) > 0
  if (!any(has)) "none" else paste(vapply(which(has), function(k) sprintf("%s %s", bms[k], paste(sprintf(paste0("%d%s (%+.", digits, "f)"), J[[k]]$age,
    if (is.null(M)) "" else edge_mark(J[[k]]$age, M), J[[k]]$jump), collapse = ", ")), ""), collapse = " | ") }
EDGE_NOTE <- "  * = the jump moves only the first or only the last age group: one group on its own, read with caution.\n"

adj_for <- function(bms) { extra <- setdiff(bms, ORIG); unique(c(ADJ_ORIG, if (length(extra)) paste0("imp_", extra))) }
n_coef <- function(bms) { adjust <- adj_for(bms); d <- prep(bms, adjust)          # coefficients of a per-age model (without the intercept)
  ncol(stats::model.matrix(stats::as.formula(paste("~", paste(c(bms, adjust), collapse = " + "))), data = d)) - 1L }
## one complete run: per-age models -> penalty by rule 3 -> line 3 -> tests -> halves
run_panel <- function(bms, bin_min, n_half, label) {
  adjust <- adj_for(bms)
  set_run(bms, adjust, bin_min); d <- prep(bms, adjust); bins <- E$make_bins(d)
  M <- mats(perage(d, bins))
  cb <- choose_bic(M$SE, show = FALSE)
  r <- ft(M$Y, M$SE, cb$bic); L <- levels_of(M, r$cuts); te <- tests(M)
  hs <- halves(d, bins, cb$bic, n_half)
  cat(sprintf("  [%s] %d age groups (>= %d cases) | multiplier %g (false alarm %.1f%%) | change points: %s  (%s)\n",
              label, nrow(M$Y), bin_min, cb$bic, 100 * cb$fa, cp_txt(L$J, bms, M), mins()))
  list(bms = bms, adjust = adjust, d_n = nrow(d), d_ev = sum(d$Dementia), bins = bins, starts = vapply(bins, min, numeric(1)), M = M, cb = cb, r = r, L = L,
       te = te, hs = hs, label = label, bin_min = bin_min)
}

## ---- 1. coverage and the panel (rule 1) ---------------------------------------------------------------
hdr("1. Which biomarkers enter the panel (rule 1)")
cat(sprintf("  script: dementia_cp_12_panel.R, version 2 (2026-11-21; batch %d; rules fixed before the run) | input: %s | N %s, dementia cases %s | functions of dementia_cp_07_tune.R version %s | %d core(s)\n", P$BATCH, basename(loc),
            format(nrow(A0), big.mark = ","), format(sum(A0$Dementia), big.mark = ","), V07, P$N_CORES))
meas <- vapply(c(ORIG, have_new), function(b) mean(A0[[paste0("imp_", b)]] == 0), numeric(1))
r_meas <- function(a, b) { u <- A0[[paste0("imp_", a)]] == 0 & A0[[paste0("imp_", b)]] == 0     # measured values only
  if (sum(u) > 1000) stats::cor(A0[[a]][u], A0[[b]][u]) else NA_real_ }
PANEL <- ORIG; WHY <- setNames(character(length(have_new)), have_new); RTXT <- WHY; PANEL1 <- ORIG
for (b in have_new) {
  rr <- vapply(PANEL, function(a) r_meas(a, b), numeric(1)); top <- which.max(abs(rr))
  WHY[b] <- if (meas[[b]] < MIN_MEASURED) "no: measured in under 40%" else
    if (b == "glu" && "a1c" %in% PANEL) "no: A1c is used" else
    if (length(top) && isTRUE(abs(rr[top]) > MAX_R)) sprintf("no: r = %.2f with %s", rr[top], names(rr)[top]) else "yes"
  RTXT[b] <- if (length(top)) sprintf("%+.2f (%s)", rr[top], names(rr)[top]) else "-"
  if (WHY[[b]] == "yes") PANEL <- c(PANEL, b)
  if (b == tail(intersect(NEWB1, have_new), 1)) PANEL1 <- PANEL       # the first-batch panel (reference for rule 2)
}
rows1 <- if (P$BATCH == 2) intersect(NEWB2, have_new) else c(ORIG, have_new)
if (P$BATCH == 2) { out1 <- setdiff(intersect(NEWB1, have_new), PANEL1)
  cat(sprintf("  first batch (as in the first run): panel %s%s\n", paste(PANEL1, collapse = ", "),
              if (length(out1)) paste0("; one-at-a-time only: ", paste(sprintf("%s (%s)", out1, sub("^no: ", "", WHY[out1])), collapse = ", ")) else "")) }
cat(sprintf("    %-5s %-18s %10s %30s   %s\n", "", "", "measured", "largest r with the panel so far", "in the main panel"))
for (b in rows1) cat(sprintf("    %-5s %-18s %9.1f%% %30s   %s\n", b, LBL[[b]], 100 * meas[[b]], if (b %in% ORIG) "" else RTXT[[b]],
  if (b %in% ORIG) "yes (original five)" else if (WHY[[b]] == "yes") "yes" else paste(WHY[[b]], "-> one at a time only")))
cat(sprintf("  main panel (%d): %s\n", length(PANEL), paste(PANEL, collapse = ", ")))

## rule 2: the first-batch panel has 100 cases per age group; a larger panel keeps the same cases per coefficient
P1 <- n_coef(PANEL1); PM <- n_coef(PANEL)
BIN_USE <- if (P$BATCH == 1) BIN_MAIN else max(BIN_MAIN, ceiling(BIN_MAIN * PM / P1))
cat(sprintf("  rule 2: at least %d cases per age group (%d coefficients per age-group model; the first-batch panel has %d with 100 cases)\n",
            BIN_USE, PM, P1))

## ---- 2. main result ----------------------------------------------------------------------------------
hdr("2. Main result: line 3 on the whole panel")
MAIN <- run_panel(PANEL, BIN_USE, P$N_HALF, "main panel")
M <- MAIN$M; K <- length(PANEL)
merged <- vapply(Filter(function(x) length(x) > 1, MAIN$bins), function(x) paste0(min(x), "-", max(x)), "")
cat(sprintf("  adjusted for: %s\n", paste(MAIN$adjust, collapse = ", ")))
cat(sprintf("  age groups: %d%s%s\n", nrow(M$Y), if (length(merged)) paste0("; merged: ", paste(merged, collapse = ", ")) else "",
            if (length(M$lost)) sprintf("; NOTE: %d group(s) could not be fitted, left out (%s)", length(M$lost), paste(M$lost, collapse = ", ")) else ""))
cat(sprintf("  penalty (rule 3), false alarm in %d data sets with no change: %s -> %g\n", P$NULL_SIMS,
            paste(sprintf("x %g: %.1f%%", MAIN$cb$tab$bic, 100 * MAIN$cb$tab$fa), collapse = " | "), MAIN$cb$bic))
fmtc <- function(cu) { has <- lengths(cu) > 0; if (!any(has)) "none" else paste(sprintf("%s %s", PANEL[has], vapply(cu[has], function(q) paste(M$age[q], collapse = ","), "")), collapse = " | ") }
crit <- vapply(MAIN$r$log, function(l) l$crit, numeric(1)); ctx <- vapply(MAIN$r$log, function(l) fmtc(l$cuts), "")
cat(sprintf("  line 3 re-weighted in %d stages (criterion %.1f -> %.1f): %s\n", length(crit), crit[1], tail(crit, 1),
            if (length(unique(ctx)) == 1) "the same change points at every stage" else paste(sprintf("stage %d: %s", seq_along(ctx) - 1L, ctx), collapse = "; ")))
CPS <- sort(unique(unlist(lapply(MAIN$L$J, function(j) j$age))))
cat(sprintf("  CHANGE POINTS (%d): %s\n", length(CPS), if (length(CPS)) paste0(CPS, edge_mark(CPS, M), collapse = ", ") else "none"))
if (any(nzchar(edge_mark(CPS, M)))) cat(EDGE_NOTE)
chg <- which(vapply(MAIN$L$J, nrow, integer(1)) > 0)
if (length(chg)) {
  seg_s <- c(M$age[1], CPS); seg_e <- c(CPS - 1, 79)
  cat("  effect in each age segment (log-odds per SD; only the biomarkers that change):\n")
  cat(sprintf("    %-9s%s\n", "ages", paste(sprintf("%8s", PANEL[chg]), collapse = "")))
  for (i in seq_along(seg_s)) cat(sprintf("    %-9s%s\n", paste0(seg_s[i], "-", seg_e[i]), paste(sprintf("%8.3f", MAIN$L$B[match(seg_s[i], M$age), chg]), collapse = "")))
}
te <- MAIN$te; te$p_vary_holm <- stats::p.adjust(te$p_vary, "holm"); nh <- sum(!vapply(MAIN$hs, is.null, logical(1)))
cat(sprintf("\n  each biomarker (Holm = corrected for %d biomarkers; halves = share of %d random halves with a jump at any age / near the age found):\n", K, nh))
cat(sprintf("    %-5s %-16s %24s %8s %13s %18s  %-20s %s\n", "", "", "average (95% CI)", "p", "varies: Holm", "trend /10 y (p)", "change points", "halves"))
for (k in seq_len(K)) { j <- MAIN$L$J[[k]]
  hv <- if (nh > 0) sprintf("%.0f%% / %s", 100 * any_share(MAIN$hs, k), if (nrow(j)) paste(sprintf("%.0f%%", 100 * vapply(j$age, function(a) near_share(MAIN$hs, k, a, MAIN$starts), numeric(1))), collapse = ", ") else "-") else ""
  cat(sprintf("    %-5s %-16s %24s %8s %13s %18s  %-20s %s\n", PANEL[k], LBL[[PANEL[k]]],
              sprintf("%+.3f (%+.3f, %+.3f)", te$effect[k], te$effect[k] - 1.96 * te$se[k], te$effect[k] + 1.96 * te$se[k]), pv(te$p_effect[k]),
              pv(te$p_vary_holm[k]), sprintf("%+.3f (%s)", te$slope10[k], pv(te$p_trend[k])),
              if (nrow(j)) paste(sprintf("%d%s (%+.2f)", j$age, edge_mark(j$age, M), j$jump), collapse = ", ") else "none", hv)) }

## figure: estimates with 95% intervals, line 3 steps (solid), straight line (dashed)
fig <- file.path(path.expand("~/dementia_v9"), sprintf("panel_fig_%s.png", format(Sys.time(), "%m%d-%H%M")))
fig_ok <- tryCatch({
  dir.create(dirname(fig), showWarnings = FALSE)
  opened <- tryCatch({ grDevices::png(fig, width = 2600, height = 650 * ceiling((K + 1) / 4), res = 200, type = "cairo"); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!opened) grDevices::png(fig, width = 2600, height = 650 * ceiling((K + 1) / 4), res = 200)
  graphics::par(mfrow = c(ceiling((K + 1) / 4), 4), mar = c(4, 4.2, 2.6, 0.8), mgp = c(2.3, 0.7, 0), cex = 0.75)
  for (k in seq_len(K)) {
    y <- M$Y[, k]; s <- M$SE[, k]; lo <- y - 1.96 * s; hi <- y + 1.96 * s
    yl <- range(c(stats::quantile(lo, 0.03), stats::quantile(hi, 0.97), MAIN$L$B[, k]))
    graphics::plot(M$mid, y, type = "n", xlim = c(45, 80), ylim = yl, xlab = "age at baseline", ylab = "log-odds per SD", main = LBL[[PANEL[k]]])
    graphics::abline(h = 0, col = "grey50"); if (nrow(MAIN$L$J[[k]])) graphics::abline(v = MAIN$L$J[[k]]$age, col = "grey70", lty = 3)
    graphics::segments(M$mid, pmax(lo, yl[1] - 1), M$mid, pmin(hi, yl[2] + 1), col = "grey60"); graphics::points(M$mid, y, pch = 16, cex = 0.7, col = "grey25")
    lf <- stats::lm(y ~ M$mid, weights = 1 / s^2); graphics::lines(c(45, 80), stats::coef(lf)[1] + stats::coef(lf)[2] * c(45, 80), lty = 2, lwd = 1.5, col = "#1f77b4")
    st <- c(M$age, 80); for (i in seq_along(M$age)) graphics::segments(st[i], MAIN$L$B[i, k], st[i + 1], MAIN$L$B[i, k], lwd = 2.2, col = "#d62728")
  }
  graphics::plot.new(); graphics::legend("center", bty = "n", cex = 1, lwd = c(NA, 2.2, 1.5, 1), lty = c(NA, 1, 2, 3), pch = c(16, NA, NA, NA),
                                         col = c("grey25", "#d62728", "#1f77b4", "grey70"),
                                         legend = c("age group, 95% interval", "line 3 steps", "straight line", "change point"))
  grDevices::dev.off(); TRUE }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
cat(sprintf("  figure: %s\n", if (fig_ok) fig else "could not be drawn"))

## ---- 3. robustness ---------------------------------------------------------------------------------
ONE <- list(); BR <- NULL
if (P$RUN_BRIDGE) {
  hdr(sprintf("3a. The five original biomarkers with %d-case age groups (what the coarser age groups alone change)", BIN_USE))
  BR <- run_panel(ORIG, BIN_USE, 0, sprintf("original five, %d-case groups", BIN_USE))
}
if (P$RUN_ONE) {
  newb_one <- intersect(if (P$BATCH == 1) NEWB1 else NEWB2, have_new)
  hdr(sprintf("3b. The original five + one new biomarker at a time (40-case age groups, as the paper's main analysis; %d runs)", length(newb_one)))
  for (b in newb_one) ONE[[b]] <- run_panel(c(ORIG, b), BIN_ONE, P$N_HALF_ONE, sprintf("five + %s", b))
  note5 <- vapply(names(ONE), function(b) { rr <- vapply(ORIG, function(a) r_meas(a, b), numeric(1)); top <- which.max(abs(rr))
    if (length(top) && isTRUE(abs(rr[top]) > MAX_R)) sprintf("r = %.2f with %s", rr[top], names(rr)[top]) else "" }, "")
  cat(sprintf("\n  the new biomarker in each run (halves: %d random halves, jump at any age / near the age found):\n", P$N_HALF_ONE))
  cat(sprintf("    %-5s %-16s %9s %24s %8s %10s  %-20s %-14s %s\n", "", "", "measured", "average (95% CI)", "p", "varies: p", "change points", "halves", "note"))
  for (b in names(ONE)) { R1 <- ONE[[b]]; k <- length(R1$bms); t1 <- R1$te[k, ]; j <- R1$L$J[[k]]
    st <- if (P$N_HALF_ONE > 0) sprintf("%.0f%% / %s", 100 * any_share(R1$hs, k),
                                         if (nrow(j)) paste(sprintf("%.0f%%", 100 * vapply(j$age, function(a) near_share(R1$hs, k, a, R1$starts), numeric(1))), collapse = ", ") else "-") else ""
    cat(sprintf("    %-5s %-16s %8.1f%% %24s %8s %10s  %-20s %-14s %s\n", b, LBL[[b]], 100 * meas[[b]],
                sprintf("%+.3f (%+.3f, %+.3f)", t1$effect, t1$effect - 1.96 * t1$se, t1$effect + 1.96 * t1$se), pv(t1$p_effect), pv(t1$p_vary),
                if (nrow(j)) paste(sprintf("%d%s (%+.2f)", j$age, edge_mark(j$age, R1$M), j$jump), collapse = ", ") else "none", st, note5[[b]])) }
  if (any(nzchar(note5))) cat("  note: r > 0.5 with one of the five: the two share the effect in that run.\n")
  if (any(vapply(ONE, function(R1) any(nzchar(edge_mark(unlist(lapply(R1$L$J, function(j) j$age)), R1$M))), logical(1)))) cat(EDGE_NOTE)
}

## ---- 4. summary ------------------------------------------------------------------------------------
hdr("4. Summary")
lst <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"
nochg <- !(seq_len(K) %in% chg)
cat(sprintf("  main panel (%d biomarkers, %d age groups, multiplier %g): change points %s\n", K, nrow(M$Y), MAIN$cb$bic, cp_txt(MAIN$L$J, PANEL, M)))
cat(sprintf("    effect varies with age (Holm p < 0.05): %s\n", lst(PANEL[te$p_vary_holm < 0.05])))
cat(sprintf("    associated (p < 0.05), no change point: %s\n", lst(PANEL[te$p_effect < 0.05 & nochg])))
cat(sprintf("    not associated (p >= 0.05), no change point: %s\n", lst(PANEL[te$p_effect >= 0.05 & nochg])))
if (length(ONE)) cat(sprintf("  one at a time with the five (change points of the new one): %s\n",
                             paste(vapply(names(ONE), function(b) { j <- ONE[[b]]$L$J[[6]]; sprintf("%s %s", b, if (nrow(j)) paste(j$age, collapse = ",") else "-") }, ""), collapse = " | ")))
if (!is.null(BR)) cat(sprintf("  original five, %d-case groups: %s (multiplier %g)\n", BIN_USE, cp_txt(BR$L$J, ORIG, BR$M), BR$cb$bic))
cat("  a change point is credible when most random halves find it near the same age and it survives the one-at-a-time runs.\n")

## CSV (model-level only)
row_of <- function(R, k, run) { j <- R$L$J[[k]]; t1 <- R$te[k, ]
  data.frame(run = run, biomarker = R$bms[k], measured = meas[[R$bms[k]]], n_biomarkers = length(R$bms), age_groups = nrow(R$M$Y), min_cases = R$bin_min,
             multiplier = R$cb$bic, false_alarm = R$cb$fa, effect = t1$effect, se = t1$se, p_effect = t1$p_effect, p_vary = t1$p_vary, slope10 = t1$slope10,
             p_trend = t1$p_trend, change_ages = paste(j$age, collapse = " "), jumps = paste(sprintf("%.4f", j$jump), collapse = " "),
             halves = length(Filter(Negate(is.null), R$hs)), any_jump_halves = any_share(R$hs, k),
             near_halves = paste(sprintf("%.2f", vapply(j$age, function(a) near_share(R$hs, k, a, R$starts), numeric(1))), collapse = " ")) }
OUT <- do.call(rbind, c(lapply(seq_len(K), function(k) row_of(MAIN, k, "main panel")),
                        if (!is.null(BR)) lapply(1:5, function(k) row_of(BR, k, sprintf("original five, %d-case groups", BIN_USE))),
                        unlist(lapply(names(ONE), function(b) lapply(seq_along(ONE[[b]]$bms), function(k) row_of(ONE[[b]], k, paste("original five +", b)))), recursive = FALSE)))
OUT$p_vary_holm_main <- ifelse(OUT$run == "main panel", stats::p.adjust(ifelse(OUT$run == "main panel", OUT$p_vary, NA), "holm"), NA)
fcsv <- file.path(path.expand("~/dementia_v9"), sprintf("panel_results_%s.csv", format(Sys.time(), "%Y%m%d-%H%M")))
utils::write.csv(OUT, fcsv, row.names = FALSE)
pa_out <- data.frame(age_from = rep(M$age, K), age_to = rep(M$to, K), biomarker = rep(PANEL, each = nrow(M$Y)), estimate = as.numeric(M$Y),
                     se = as.numeric(M$SE), fitted = as.numeric(MAIN$L$B))
utils::write.csv(pa_out, sub("panel_results_", "panel_perage_", fcsv), row.names = FALSE)
PANEL_RESULT <- list(main = MAIN, bridge = BR, one = ONE, table = OUT)
cat(sprintf("\n  saved: %s and %s (model-level numbers only) | total %s\n", basename(fcsv), basename(sub("panel_results_", "panel_perage_", fcsv)), mins()))
