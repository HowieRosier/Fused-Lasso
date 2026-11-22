## =============================================================================
## dementia_cp_10_cox.R
## The last step of the paper: the change points of each line used in Cox models of time to dementia.
##   1. Age regimes, as in step 8 of the paper's scripts:
##        Surv(time2, Dementia) ~ age regimes + 5 biomarkers (per SD) + 9 categorical covariates
##      Cox BIC on all people, and C-index on a held-out 40% (model fitted on a random 60%, seed 42, as in the scripts).
##   2. Biomarker effects before and after the change points: the same Cox model with age adjusted by a smooth curve
##      (natural spline, 4 df), and each biomarker's effect allowed to differ between the age segments its own change points
##      define. Hazard ratios per SD in each segment, and a likelihood-ratio test against effects that stay the same at every age.
## Change points: those of the comparison for the paper (dementia_cp_08_compare.R, part B). If that script was run in this R
## session they are taken from it; otherwise the values of the v9 run of 2026-11-20 (adjusted per-age models) are used.
## Set CP_SETS to use others, e.g. CP_SETS <- list("my set" = list(sbp = 53, hdl = 74)).
## Wider panel (after dementia_cp_11_extra_labs.R): set BIOMARKERS. The change points then come from CP_SETS, or from
## dementia_cp_12_panel.R (main panel) or dementia_cp_08_compare.R run with the same biomarkers in this R session. The
## 'filled' flags of the new biomarkers (imp_dbp, ...) are added to the covariates.
##
## HOW TO RUN:  source("~/dementia_cp_10_cox.R")          (a few minutes)
## Output: model results only (hazard ratios, BIC, C-index). Counts of people are shown only when above 20.
## =============================================================================

.c10 <- list(IN_CSV = if (exists("IN_CSV")) IN_CSV else "", CP_SETS = if (exists("CP_SETS")) CP_SETS else NULL,
             SHOW_ALL_HR = if (exists("SHOW_ALL_HR")) isTRUE(SHOW_ALL_HR) else FALSE,
             BIOMARKERS = if (exists("BIOMARKERS")) as.character(BIOMARKERS) else c("sbp", "bmi", "ldlc", "hdl", "trig"))
suppressWarnings(rm(list = intersect(c("CP_SETS", "SHOW_ALL_HR", "BIOMARKERS"), ls(envir = globalenv())), envir = globalenv()))
.c10$five <- identical(.c10$BIOMARKERS, c("sbp", "bmi", "ldlc", "hdl", "trig"))
.c10$t0 <- Sys.time(); options(width = 170)
for (p in c("dplyr", "survival", "splines")) if (!requireNamespace(p, quietly = TRUE)) utils::install.packages(p, quiet = TRUE)
suppressPackageStartupMessages({ library(dplyr); library(survival); library(splines) })
hdr10 <- function(x) cat("\n==================== ", x, " ====================\n", sep = "")

## small counts (1-20) are never shown, also not when they could be worked out from other cells
SMALL <- function(n) !is.na(n) & n >= 1 & n <= 20
msk <- function(n) vapply(suppressWarnings(as.numeric(n)), function(v)
  if (is.na(v)) "NA" else if (v >= 1 && v <= 20) "<=20" else format(round(v), big.mark = ",", trim = TRUE), character(1))
part_mask <- function(cells) {
  cells <- as.numeric(cells); cells[is.na(cells)] <- 0; hide <- SMALL(cells)
  while (any(hide) && (sum(hide) == 1 || SMALL(sum(cells[hide]))) && any(!hide & cells > 0)) {
    cand <- which(!hide & cells > 0); hide[cand[which.min(cells[cand])]] <- TRUE
  }
  ifelse(hide, ifelse(SMALL(cells), "<=20", "hidden"), msk(cells))
}
exact_shown <- function(x) !grepl("<=|hidden", x)

## ---- data (the same rows as in dementia_cp_07_tune.R, plus complete covariates) ------------------------------------------------
.c10$loc <- if (nzchar(.c10$IN_CSV)) path.expand(.c10$IN_CSV) else sort(Sys.glob(path.expand("~/dementia_v9/dementia_v9_analysis_*.csv")), decreasing = TRUE)[1]
if (is.na(.c10$loc) || !file.exists(.c10$loc)) stop("Cohort file not found. Run dementia_cp_02_v9_cohort.R first, or set IN_CSV.", call. = FALSE)
biomarkers <- .c10$BIOMARKERS
covs <- c("diabete", "gender", "race", "smoker", "anti", "hearing_loss", "ethnicity", "brain", "depression")   # as in the paper's scripts
imps <- setdiff(biomarkers, c("sbp", "bmi", "ldlc", "hdl", "trig")); imps <- if (length(imps)) paste0("imp_", imps) else character(0)   # 'filled' flags of new biomarkers
dat <- utils::read.csv(.c10$loc, stringsAsFactors = FALSE)
miss <- setdiff(c("Dementia", "age", "time2", biomarkers, covs, imps), names(dat))
if (length(miss)) stop("The cohort file has no column(s): ", paste(miss, collapse = ", "), call. = FALSE)
dat <- dat %>% filter(is.finite(time2), time2 > 0, is.finite(age), age >= 45, age < 80, ldlc < 400, trig < 700)
dat <- dat[stats::complete.cases(dat[, c(biomarkers, covs, imps)]), c("Dementia", "age", "time2", biomarkers, covs, imps)]
for (v in biomarkers) dat[[v]] <- as.numeric(scale(dat[[v]]))                                   # per SD, as in the scripts
for (v in covs) dat[[v]] <- factor(dat[[v]])
set.seed(42); .c10$tr <- sample(seq_len(nrow(dat)), floor(0.60 * nrow(dat)))                    # the same split for every model

## ---- change points --------------------------------------------------------------------------------------------------------------
.c10$src <- "values of the v9 run of 2026-11-20 (dementia_cp_08_compare.R part B, adjusted per-age models)"
sets <- list("line 1  group lasso, standard"           = setNames(rep(list(48), 5), c("sbp", "bmi", "ldlc", "hdl", "trig")),
             "line 2  fused lasso, standard"           = list(sbp = c(53, 57, 67, 77), hdl = c(70, 74)),
             "line 3  sparse group lasso, fine-tuned"  = list(sbp = 53, hdl = 74))
if (!is.null(.c10$CP_SETS)) { sets <- .c10$CP_SETS; .c10$src <- "CP_SETS" } else
if (exists("CMP_RESULT") && all(c("std1", "std2", "line3ft") %in% names(CMP_RESULT)) && is.null(CMP_RESULT$std1$error) &&
    is.null(CMP_RESULT$std2$error) && is.null(CMP_RESULT$line3ft$error) && identical(CMP_RESULT$biomarkers, biomarkers)) {
  nm <- CMP_RESULT$biomarkers; as_set <- function(per) { x <- setNames(per, nm); x[vapply(x, length, 1L) > 0] }
  sets <- list("line 1  group lasso, standard" = as_set(CMP_RESULT$std1$per), "line 2  fused lasso, standard" = as_set(CMP_RESULT$std2$per),
               "line 3  sparse group lasso, fine-tuned" = as_set(CMP_RESULT$line3ft$per))
  .c10$src <- "the dementia_cp_08_compare.R run in this R session (part B)"
} else if (!.c10$five && exists("PANEL_RESULT") && identical(PANEL_RESULT$main$bms, biomarkers)) {
  J <- PANEL_RESULT$main$L$J; x <- setNames(lapply(J, function(j) j$age), biomarkers)
  sets <- list("line 3 fine-tuned, wider panel" = x[vapply(x, length, 1L) > 0])
  .c10$src <- "the dementia_cp_12_panel.R run in this R session (main panel)"
} else if (!.c10$five) stop("BIOMARKERS is not the paper's five: set CP_SETS, or run dementia_cp_12_panel.R (or 08 with the same BIOMARKERS) first.", call. = FALSE)
sets <- c(list("no change points" = list()), sets)

hdr10("Settings")
cat(sprintf("  script: dementia_cp_10_cox.R, version 1.2 (2026-11-21) | cohort file: %s | biomarkers: %s%s\n", basename(.c10$loc), paste(biomarkers, collapse = ", "), if (length(imps)) paste0(" | extra covariates: ", paste(imps, collapse = ", ")) else ""))
cat(sprintf("  people with complete biomarkers and covariates: N %s, dementia cases %s | held-out 40%% for the C-index (seed 42)\n",
            msk(nrow(dat)), msk(sum(dat$Dementia))))
cat(sprintf("  change points from: %s\n", .c10$src))
for (s in names(sets)[-1]) cat(sprintf("    %-40s %s\n", s, if (!length(sets[[s]])) "none" else
  paste(sprintf("%s %s", names(sets[[s]]), vapply(sets[[s]], function(a) paste(sort(a), collapse = ","), "")), collapse = " | ")))

cindex <- function(form) { fit <- suppressWarnings(coxph(form, data = dat[.c10$tr, ]))
  lp <- predict(fit, newdata = dat[-.c10$tr, ], type = "lp")
  concordance(Surv(time2, Dementia) ~ lp, data = dat[-.c10$tr, ], reverse = TRUE)$concordance }
cov_txt <- paste(c(paste0("factor(", covs, ")"), imps), collapse = " + ")

## ---- 1. age regimes ------------------------------------------------------------------------------------------------------------
hdr10("1. Age regimes (as in the paper's scripts): Cox BIC and C-index")
cat(sprintf("  Surv(time2, Dementia) ~ age regimes + %s + 9 covariates%s; regimes are cut at every change point of the line.\n",
            paste(biomarkers, collapse = " + "), if (length(imps)) paste0(" + ", length(imps), " filled flags") else ""))
reg_out <- list()
for (s in names(sets)) {
  cps <- sort(unique(unlist(sets[[s]])))
  if (length(cps)) { edges <- c(45, cps, 80)
    dat$age_regime <- cut(dat$age, breaks = c(-Inf, cps, Inf), right = FALSE, labels = sprintf("ages %g-%g", head(edges, -1), tail(edges, -1) - 1))
    rt <- "factor(age_regime) + " } else rt <- ""
  form <- stats::as.formula(paste("Surv(time2, Dementia) ~", rt, paste(biomarkers, collapse = " + "), "+", cov_txt))
  fit <- suppressWarnings(coxph(form, data = dat))
  reg_out[[s]] <- data.frame(line = s, regimes = length(cps) + 1, cox_bic = BIC(fit), c_index_test = cindex(form))
  if (length(cps)) {
    rt_tab <- dat %>% group_by(age_regime) %>% summarise(n = n(), ev = sum(Dementia), .groups = "drop")
    pn <- part_mask(rt_tab$n); pe <- part_mask(rt_tab$ev)
    reg_out[[s]]$regime_table <- paste(sprintf("%s N %s, cases %s%s", as.character(rt_tab$age_regime), pn, pe,
                                               ifelse(exact_shown(pn) & exact_shown(pe), sprintf(" (%.1f%%)", 100 * rt_tab$ev / rt_tab$n), "")), collapse = " ; ")
  } else reg_out[[s]]$regime_table <- ""
}
reg_out <- do.call(rbind, reg_out)
cat(sprintf("\n  %-40s %8s %12s %22s\n", "", "regimes", "Cox BIC", "C-index (held-out 40%)"))
for (i in seq_len(nrow(reg_out))) cat(sprintf("  %-40s %8d %12.1f %22.4f\n", reg_out$line[i], reg_out$regimes[i], reg_out$cox_bic[i], reg_out$c_index_test[i]))
cat("\n  Regimes (age at index; dementia cases and share):\n")
for (i in which(nzchar(reg_out$regime_table))) cat(sprintf("    %s:\n      %s\n", reg_out$line[i], gsub(" ; ", "\n      ", reg_out$regime_table[i])))
cat("  reading: these models use the change points only to cut age into regimes, so they reward cutting age into more pieces (dementia\n")
cat("  risk rises steeply with age), whether or not the biomarker effects really change there. In a simulated cohort whose true change\n")
cat("  points were sbp 53 and hdl 74, line 2's seven regimes (four of its six change points false) had a lower BIC and a higher C-index\n")
cat("  than line 3's three. Part 2 is the check of the change points themselves.\n")

## ---- 2. biomarker effects before and after the change points ------------------------------------------------------------------
hdr10("2. Biomarker effects before and after each change point (Cox, age adjusted by a smooth curve)")
cat(sprintf("  Surv(time2, Dementia) ~ ns(age, 4) + 9 covariates%s + biomarkers, each biomarker with its own effect in every segment between its\n",
            if (length(imps)) paste0(" + ", length(imps), " filled flags") else ""))
cat("  change points; compared with the same model in which every effect stays the same at every age.\n")
base_txt <- paste("Surv(time2, Dementia) ~ ns(age, df = 4) +", cov_txt)
f0 <- stats::as.formula(paste(base_txt, "+", paste(biomarkers, collapse = " + ")))
fit0 <- suppressWarnings(coxph(f0, data = dat)); c0 <- cindex(f0)
eff_out <- list(); hr_out <- list()
for (s in names(sets)[-1]) {
  terms <- character(0); seg_lab <- list()
  for (v in biomarkers) {
    cps <- sort(unique(sets[[s]][[v]]))
    if (!length(cps)) { terms <- c(terms, v); next }
    seg <- findInterval(dat$age, cps); edges <- c(45, cps, 80)
    for (j in 0:length(cps)) { nmj <- sprintf("%s__s%d", v, j); dat[[nmj]] <- dat[[v]] * (seg == j); terms <- c(terms, nmj)
      seg_lab[[nmj]] <- sprintf("%s, ages %g-%g", v, edges[j + 1], edges[j + 2] - 1) }
  }
  if (!length(seg_lab)) {           # this line found no change point: its model is the 'same at every age' model
    eff_out[[s]] <- data.frame(line = paste(s, "- none"), extra_parameters = 0L, lr_chisq = 0, p = NA_real_, bic_change = 0, c_index_test = c0)
    next }
  f1 <- stats::as.formula(paste(base_txt, "+", paste(terms, collapse = " + ")))
  fit1 <- suppressWarnings(coxph(f1, data = dat))
  lr <- 2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit0))); ddf <- length(coef(fit1)) - length(coef(fit0))
  eff_out[[s]] <- data.frame(line = s, extra_parameters = ddf, lr_chisq = lr, p = stats::pchisq(lr, ddf, lower.tail = FALSE),
                             bic_change = BIC(fit1) - BIC(fit0), c_index_test = cindex(f1))
  ci <- summary(fit1)$conf.int; nm1 <- names(seg_lab)
  hr_out[[s]] <- data.frame(line = s, segment = unlist(seg_lab), hr = ci[nm1, "exp(coef)"], lo = ci[nm1, "lower .95"], hi = ci[nm1, "upper .95"],
                            p = summary(fit1)$coefficients[nm1, "Pr(>|z|)"], row.names = NULL)
}
eff_out <- do.call(rbind, eff_out)
cat(sprintf("\n  %-40s %10s %16s %10s %12s %22s\n", "", "extra", "likelihood-ratio", "", "BIC change", "C-index (held-out 40%)"))
cat(sprintf("  %-40s %10s %16s %10s %12s %22s\n", "", "parameters", "chi-square", "p", "(< 0 better)", ""))
cat(sprintf("  %-40s %10s %16s %10s %12s %22.4f\n", "effects the same at every age", "-", "-", "-", "-", c0))
for (i in seq_len(nrow(eff_out))) cat(sprintf("  %-40s %10d %16.1f %10s %12.1f %22.4f\n", substr(eff_out$line[i], 1, 40), eff_out$extra_parameters[i], eff_out$lr_chisq[i],
                                             if (is.na(eff_out$p[i])) "-" else format.pval(eff_out$p[i], digits = 2, eps = 1e-4), eff_out$bic_change[i], eff_out$c_index_test[i]))
cat("  note: the change points were found in these same data, so these p-values are optimistic; read them as a description of how\n")
cat("  large the changes are in the survival model, not as an independent test.\n")
show_hr <- if (.c10$SHOW_ALL_HR) names(hr_out) else grep("^line 3", names(hr_out), value = TRUE)
if (!length(names(hr_out))) cat("\n  (no line has a change point, so there are no segment effects to show)\n")
for (s in show_hr) {
  h <- hr_out[[s]]; cat(sprintf("\n  %s: hazard ratio per SD in each segment (95%% interval)\n", s))
  for (i in seq_len(nrow(h))) cat(sprintf("    %-26s HR %6.3f  (%.3f-%.3f)  p %s\n", h$segment[i], h$hr[i], h$lo[i], h$hi[i], format.pval(h$p[i], digits = 2, eps = 1e-4)))
}
if (!.c10$SHOW_ALL_HR) cat("\n  (SHOW_ALL_HR <- TRUE prints the same table for every line)\n")

## ---- save --------------------------------------------------------------------------------------------------------------------
.c10$out_dir <- path.expand("~/dementia_v9"); dir.create(.c10$out_dir, showWarnings = FALSE)
.c10$stamp <- format(Sys.time(), "%Y%m%d-%H%M")
utils::write.csv(reg_out, file.path(.c10$out_dir, sprintf("cox_regimes_%s.csv", .c10$stamp)), row.names = FALSE)
utils::write.csv(eff_out, file.path(.c10$out_dir, sprintf("cox_effects_%s.csv", .c10$stamp)), row.names = FALSE)
utils::write.csv(do.call(rbind, hr_out), file.path(.c10$out_dir, sprintf("cox_hr_%s.csv", .c10$stamp)), row.names = FALSE)
COX_RESULT <- list(regimes = reg_out, effects = eff_out, hr = hr_out, sets = sets)
cat(sprintf("\n  saved: ~/dementia_v9/cox_{regimes,effects,hr}_%s.csv  (model results only)\n", .c10$stamp))
cat(sprintf("\nruntime %.1f min\n", as.numeric(difftime(Sys.time(), .c10$t0, units = "mins"))))
