## =============================================================================
## dementia_cp_11_extra_labs.R
## Dementia change-point paper -- adds more biomarkers to the existing analysis cohort.
##
## The cohort is NOT rebuilt. The script reads the cohort files written by dementia_cp_02_v9_cohort.R
## (~/dementia_v9/dementia_v9_{relax,analysis}_<CDR>.csv) and adds columns to them. Everything already in the
## files stays exactly as it is (same people, index dates, dementia cases, labs, covariates; the old part of
## every line is kept byte for byte). BigQuery recomputes each person's index date only to place the windows;
## if any index date differs from the file, the script stops before writing anything.
##
## New columns (value = latest plausible value in the window around the index date; several values at that
## same time are averaged):
##   dbp   diastolic blood pressure (mmHg)      window index -1 y .. +1 y (as SBP)        used if 30-150
##   hr    heart rate / pulse (beats per min)   window index -1 y .. +1 y                 used if 30-200
##   a1c   hemoglobin A1c (%)                   window index -3 y .. +1 y (as the lipids)  used if 3-20
##   glu   glucose, serum/plasma/blood (mg/dL; no point-of-care meters)                   used if 40-600
##   crea  creatinine, serum/plasma/blood (mg/dL)                                         used if 0.2-20
##   hgb   hemoglobin, blood (g/dL)                                                       used if 4-22; 40-220 = g/L, divided by 10
##   alb   albumin, serum/plasma (g/dL)                                                   used if 1-6.5; 10-65 = g/L, divided by 10
##   egfr  eGFR from crea, CKD-EPI 2021 (no race term), age at index, sex at birth (mL/min/1.73 m2)
## In this CDR most hemoglobin and albumin values of the main LOINC codes are in g/L (median 124 and 40; checked with the
## code list of section 1), so values in the g/L range are divided by 10; the two ranges do not overlap. Any other value
## outside the range (creatinine in umol/L, A1c in mmol/mol, ...) is not used; records in mmol/L are not used (as in 02).
## Second batch (2026-11-21; same rules, chosen before any of their values were seen): bun (urea nitrogen, mg/dL), na (sodium,
## mmol/L), ca (total calcium, mg/dL; mmol/L values x 4.008), alt, ast, alp (U/L), bili (total bilirubin, mg/dL; umol/L records
## dropped), wbc, plt (10^3/uL; values per uL / 1000), rdw (%), mcv (fL), tsh (mIU/L); all with window index -3 y .. +1 y.
## Codes: standard LOINC / SNOMED measurement concepts chosen by name (rules in NEW below). Section 1 lists every
## code used, with the number of EHR participants and the median value, and the largest similar codes NOT used.
## Missing value: imp_<lab> = 1, filled with the median within race x sex at birth (as 02 does for its labs).
## egfr: computed from the (measured or filled) creatinine; imp_egfr = imp_crea, and also 1 when sex at birth is
## not female or male (the equation needs it; then eGFR is filled like a missing value).
##
## Participant-level files stay in this workspace (app disk + workspace bucket). The console shows aggregates
## only; participant counts of 1-20 print as <=20.
## =============================================================================

if (!exists("CDR")) CDR <- ""                 # "" = the latest All of Us CDR attached (C2025Q4R6 here); or a full path
if (!exists("INVENTORY")) INVENTORY <- TRUE   # section 1: list the codes in use (one extra BigQuery scan); TRUE = every lab,
                                              # FALSE = none, or the labs to list, e.g. c("wbc", "mcv")
INV_TOP <- 5                                   # inventory: the 5 largest used codes per lab (+ the 3 largest left out)
INV_SET <- INVENTORY; rm(INVENTORY)           # applies to this run only (a value left over from an earlier run is not reused)
INV_ON <- isTRUE(INV_SET) || is.character(INV_SET)

t0 <- Sys.time()
options(width = 160, scipen = 999, warn = 1)
for (p in c("jsonlite", "bigrquery")) if (!requireNamespace(p, quietly = TRUE)) utils::install.packages(p, quiet = TRUE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
hdr <- function(s) cat("\n==================== ", s, " ====================\n", sep = "")
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
sh <- function(cmd, args = character(), err = TRUE) {
  if (!nzchar(Sys.which(cmd))) return(structure(character(), status = 127L))
  tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = err)),
           error = function(e) structure(paste("ERROR:", conditionMessage(e)), status = 1L))
}
ok <- function(x) { s <- attr(x, "status"); is.null(s) || identical(as.integer(s), 0L) }

## ---- 0. workspace, data, billing (as in 02) ---------------------------------------------
hdr("0. Setup")
ctx_file <- path.expand("~/.workbench/context.json")
ctx <- if (file.exists(ctx_file)) tryCatch(jsonlite::fromJSON(ctx_file, simplifyVector = FALSE), error = function(e) NULL)
BILLING <- Sys.getenv("GOOGLE_PROJECT")
if (!nzchar(BILLING)) BILLING <- as.character(ctx$workspace$googleProjectId %||% "")
if (!nzchar(BILLING)) BILLING <- Sys.getenv("GOOGLE_CLOUD_PROJECT")
if (!nzchar(BILLING)) { g <- sh("gcloud", c("config", "get-value", "project"), err = FALSE); g <- trimws(tail(g[nzchar(trimws(g))], 1)); if (ok(g) && length(g)) BILLING <- g }
if (!nzchar(BILLING)) stop("No billing project found (~/.workbench/context.json missing). Run this inside a Workbench 2.0 app.")
res <- list()
if (nzchar(Sys.which("wb"))) {
  res <- tryCatch(jsonlite::fromJSON(paste(sh("wb", c("resource", "list", "--format=json"), err = FALSE), collapse = "\n"),
                                     simplifyVector = FALSE), error = function(e) list())
  if (is.list(res) && !is.null(res$resources)) res <- res$resources
}
is_bq <- function(r) grepl("BQ_DATASET|BIG_?QUERY_DATASET", r$resourceType %||% "") && !is.null(r$projectId) && !is.null(r$datasetId)
BQ_PATHS <- vapply(Filter(is_bq, res), function(r) paste0(r$projectId, ".", r$datasetId), "")
resolve_res <- function(id) {
  for (flag in c("--id=", "--name=")) {
    p <- sh("wb", c("resource", "resolve", paste0(flag, id)), err = FALSE)
    p <- trimws(p[nzchar(trimws(p))])
    if (ok(p) && length(p)) return(tail(p, 1))
  }
  NA_character_
}
own <- function(r) { st <- toupper(as.character(r$stewardshipType %||% "")); !nzchar(st) || st == "CONTROLLED" }
BUCKETS <- unlist(lapply(Filter(function(r) grepl("GCS_BUCKET", r$resourceType %||% "") && own(r), res), function(r) {
  b <- if (!is.null(r$bucketName)) r$bucketName else resolve_res(as.character(r$id %||% r$name %||% ""))
  if (is.na(b) || !nzchar(b)) NULL else sub("/+$", "", if (startsWith(b, "gs://")) b else paste0("gs://", b))
}))
if (nzchar(Sys.getenv("WORKSPACE_BUCKET"))) BUCKETS <- unique(c(BUCKETS, sub("/+$", "", Sys.getenv("WORKSPACE_BUCKET"))))
BUCKETS <- BUCKETS[!grepl("^gs://(vwb|fc)-aou-datasets", BUCKETS)]
if (!nzchar(CDR)) {
  hit <- BQ_PATHS[grepl("\\.C[0-9]{4}Q[0-9]R[0-9]+$", BQ_PATHS)]
  hit <- hit[order(sub("^.*\\.", "", hit), decreasing = TRUE)]
  if (!length(hit)) stop("No All of Us CDR dataset (named like C2025Q4R6) is attached to this workspace.", call. = FALSE)
  CDR <- hit[1]
}
CDR_ID <- sub("^.*\\.", "", CDR)
cat("  script:             dementia_cp_11_extra_labs.R, version 1.6 (2026-11-21; MCV codes fixed, blood-count SNOMED codes added, shorter output)\n")
cat("  billing project:   ", BILLING, "\n")
cat("  CDR used:          ", CDR, "\n")
cat("  buckets:           ", if (length(BUCKETS)) paste(BUCKETS, collapse = ", ") else "<none found>", "\n")
tab_exists <- function(path, table) {
  pp <- strsplit(path, ".", fixed = TRUE)[[1]]
  tryCatch(isTRUE(bigrquery::bq_table_exists(bigrquery::bq_table(pp[1], pp[2], table))), error = function(e) FALSE)
}
HAS_CRIT <- tab_exists(CDR, "cb_criteria")
if (tab_exists(CDR, "cb_search_person")) {
  EHR_SQL <- sprintf("SELECT person_id FROM `%s.cb_search_person` WHERE has_ehr_data = 1", CDR)
} else {
  PR <- strsplit(CDR, ".", fixed = TRUE)[[1]]
  idx <- BQ_PATHS[startsWith(sub("^.*\\.", "", BQ_PATHS), paste0(CDR_ID, "_index"))]
  if (!length(idx)) {
    ds <- tryCatch(bigrquery::bq_project_datasets(PR[1]), error = function(e) list())
    nm <- vapply(ds, function(d) d$dataset %||% "", "")
    nm <- sort(nm[startsWith(nm, paste0(CDR_ID, "_index"))], decreasing = TRUE)
    if (length(nm)) idx <- paste0(PR[1], ".", nm[1])
  }
  EHR_SQL <- if (length(idx) && tab_exists(idx[1], "T_ENT_person")) sprintf("SELECT id AS person_id FROM `%s.T_ENT_person` WHERE has_ehr_data = TRUE", idx[1]) else
    sprintf("SELECT person_id FROM `%s.person`", CDR)
}

## ---- the cohort files of 02 ------------------------------------------------------------------
OUT_DIR <- file.path(path.expand("~"), "dementia_v9"); dir.create(OUT_DIR, showWarnings = FALSE)
BAK_DIR <- file.path(OUT_DIR, "before_extra_labs")
F <- c(relax = file.path(OUT_DIR, sprintf("dementia_v9_relax_%s.csv", CDR_ID)),
       analysis = file.path(OUT_DIR, sprintf("dementia_v9_analysis_%s.csv", CDR_ID)))
for (k in names(F)) if (!file.exists(F[[k]]) && length(BUCKETS)) {     # not on this app's disk: take the copy 02 put in the bucket
  r <- sh("gsutil", c("-q", "cp", shQuote(paste0(BUCKETS[1], "/dementia_v9/", CDR_ID, "/", basename(F[[k]]))), shQuote(OUT_DIR)))
  if (ok(r) && file.exists(F[[k]])) cat("  copied from the bucket:", basename(F[[k]]), "\n")
}
if (!all(file.exists(F))) stop("Cohort files of ", CDR_ID, " not found in ", OUT_DIR, " or the bucket: run dementia_cp_02_v9_cohort.R first.", call. = FALSE)

## ---- 1. the new labs: code rules ------------------------------------------------------------
## inc / exc: patterns on the lower-case concept name (standard LOINC or SNOMED measurement concepts).
## broad: the wider name pattern of the inventory, which also shows the largest codes the rules leave out.
NEW <- list(
  dbp  = list(label = "diastolic blood pressure", unit = "mmHg", win = 1, lo = 30, hi = 150,
              inc = "diastolic blood pressure|^non-invasive diastolic arterial pressure$",
              exc = paste0("exercise|ambulatory|24 hour|24h|mean|average|maximum|minimum|fetal|fetus|pulmonary|venous|percentile|z-score|target|goal|home|self|",
                           "central|ankle|thigh|leg|(^|[^-n])invasive|intra-?arterial|arterial line|night|daytime|sleep|awake|change|difference|ratio|/|index|",
                           "orthostatic|tilt|alarm|limit|systolic|\\bpre\\b|\\bpost\\b|peak|lowest|highest|maternal|anesthesia"),
              broad = "diastolic", broad_exc = "^$"),
  hr   = list(label = "heart rate / pulse", unit = "beats/min", win = 1, lo = 30, hi = 200,
              inc = "^(resting )?heart rate|^(peripheral |radial |apical )?pulse rate|^pulse$",
              exc = paste0("fetal|fetus|variability|exercise|maximum|minimum|mean|average|24 hour|target|peak|recovery|predicted|reserve|alarm|limit|",
                           "deceleration|baseline|change|ekg|ecg|electrocardiogram|holter|monitor|device|pacemaker|pacing|sleep|night|percentile|",
                           "intrinsic|ventricular|atrial|stress|training|maternal|\\bpost\\b|\\bpre\\b|/|ratio|rhythm|method|anesthesia|",
                           "(^|[^-n])invasive|intra[- ]?arterial|arterial line"),
              broad = "heart rate|pulse rate|^pulse$", broad_exc = "^$"),
  a1c  = list(label = "hemoglobin A1c", unit = "%", win = 3, lo = 3, hi = 20,
              inc = "^hemoglobin a1c/hemoglobin\\.total in blood|^hemoglobin a1c measurement|^glycated hemoglobin|^glycosylated hemoglobin|^haemoglobin a1c|^hba1c",
              exc = paste0("ifcc|mmol/mol|\\[moles|\\[mass/volume|estimated average glucose|presence|interpretation|goal|target|panel|percentile|fetal|",
                           "a1a|a1b|reference|cutoff|interval|average|future|dcct"),
              broad = "a1c|glycated h|glycosylated h|hba1c", broad_exc = "^$"),
  glu  = list(label = "glucose (serum, plasma, blood)", unit = "mg/dL", win = 3, lo = 40, hi = 600,
              inc = paste0("^glucose \\[mass/volume\\] in (serum or plasma|serum|plasma|blood|venous blood|arterial blood)|",
                           "^fasting glucose \\[mass/volume\\] in (serum or plasma|serum|plasma|blood)"),
              exc = paste0("glucometer|capillary|test strip|\\bpost\\b|minute|hour|challenge|tolerance|dose|meal|baseline|specimen|meter|monitor|",
                           "interstitial|presence|percentile|dialysate|\\bpre\\b"),
              broad = "glucose", broad_exc = "^$"),
  crea = list(label = "creatinine (serum, plasma, blood)", unit = "mg/dL", win = 3, lo = 0.2, hi = 20,
              inc = paste0("^creatinine \\[mass/volume\\] in (serum or plasma|serum|plasma|blood|venous blood|arterial blood)|",
                           "^creatinine measurement, serum|^serum creatinine( level| measurement)?$"),
              exc = "dialysis|\\bpost\\b|\\bpre\\b|baseline|clearance|urine|ratio|fluid|dialysate|hour|percentile|--",
              broad = "creatinine", broad_exc = "^$"),
  hgb  = list(label = "hemoglobin (blood)", unit = "g/dL", win = 3, lo = 4, hi = 22, alt = c(lo = 40, hi = 220, f = 0.1),
              inc = paste0("^hemoglobin \\[mass/volume\\] in (blood|arterial blood|venous blood|capillary blood|mixed venous blood)|",
                           "^(hemoglobin|haemoglobin) (estimation|determination|measurement)$"),
              exc = "a1c|fetal|plasma|serum|urine|free|glycated|presence|percentile|--|blood product",
              broad = "hemoglobin|haemoglobin", broad_exc = "a1c|corpuscular|glycated|glycosylated"),
  alb  = list(label = "albumin (serum, plasma)", unit = "g/dL", win = 3, lo = 1, hi = 6.5, alt = c(lo = 10, hi = 65, f = 0.1),
              inc = "^albumin \\[mass/volume\\] in (serum or plasma|serum|plasma|blood)|^serum albumin( measurement| level)?$",
              exc = paste0("urine|csf|fluid|ratio|micro|globulin|prealbumin|fraction|percent|presence|hour|dialysate|--|electrophoresis|cerebral|",
                           "ascites|pleural|peritoneal|synovial"),
              broad = "albumin", broad_exc = "^$"),
  ## ---- second batch (chosen 2026-11-21, before any of these values were looked at): routine panels (basic / comprehensive
  ## metabolic panel, blood count) with a plausible link to dementia, plus TSH; near-duplicates of existing ones left out
  ## (hematocrit ~ hemoglobin, total protein ~ albumin, total cholesterol ~ LDL + HDL + TG/5). No '/' in these exclusion lists:
  ## the names contain '[mass/volume]' etc.
  bun  = list(label = "urea nitrogen (serum, plasma, blood)", unit = "mg/dL", win = 3, lo = 2, hi = 150, batch = 2,
              inc = "^urea nitrogen \\[mass/volume\\] in (serum or plasma|serum|plasma|blood|venous blood|arterial blood)",
              exc = "urine|dialysis|dialysate|\\bpost\\b|\\bpre\\b|ratio|clearance|hour|fluid|--|percentile",
              broad = "urea", broad_exc = "^$"),
  na   = list(label = "sodium (serum, plasma, blood)", unit = "mmol/L", win = 3, lo = 110, hi = 170, keep_mmol = TRUE, batch = 2,
              inc = "^sodium \\[moles/volume\\] in (serum or plasma|serum|plasma|blood|venous blood|arterial blood)",
              exc = "urine|stool|sweat|fluid|dialysate|corrected|\\bpost\\b|\\bpre\\b|ratio|hour|--|percentile",
              broad = "sodium", broad_exc = "^$"),
  ca   = list(label = "calcium, total (serum, plasma)", unit = "mg/dL", win = 3, lo = 6, hi = 14, keep_mmol = TRUE, batch = 2,
              alt = c(lo = 1.5, hi = 3.5, f = 4.008), alt_note = "mmol/L, times 4.008",
              inc = "^calcium \\[(mass|moles)/volume\\] in (serum or plasma|serum|plasma|blood)",
              exc = "ionized|ionised|free|corrected|albumin|urine|fluid|\\bpost\\b|\\bpre\\b|ratio|hour|--|percentile",
              broad = "calcium", broad_exc = "channel|blocker|score|oxalate|phosphate crystal"),
  alt  = list(label = "ALT (serum, plasma)", unit = "U/L", win = 3, lo = 2, hi = 2000, batch = 2,
              inc = "^alanine aminotransferase \\[enzymatic activity/volume\\] in (serum or plasma|serum|plasma|blood)",
              exc = "urine|fluid|\\bpost\\b|\\bpre\\b|ratio|--|percentile",
              broad = "alanine aminotransferase|sgpt", broad_exc = "^$"),
  ast  = list(label = "AST (serum, plasma)", unit = "U/L", win = 3, lo = 2, hi = 2000, batch = 2,
              inc = "^aspartate aminotransferase \\[enzymatic activity/volume\\] in (serum or plasma|serum|plasma|blood)",
              exc = "urine|fluid|mitochondrial|\\bpost\\b|\\bpre\\b|ratio|--|percentile",
              broad = "aspartate aminotransferase|sgot", broad_exc = "^$"),
  alp  = list(label = "alkaline phosphatase (serum, plasma)", unit = "U/L", win = 3, lo = 10, hi = 2000, batch = 2,
              inc = "^alkaline phosphatase \\[enzymatic activity/volume\\] in (serum or plasma|serum|plasma|blood)",
              exc = "bone|liver|intestinal|placental|isoenzyme|heat|leukocyte|fraction|urine|fluid|\\bpost\\b|\\bpre\\b|ratio|--|percentile",
              broad = "alkaline phosphatase", broad_exc = "^$"),
  bili = list(label = "bilirubin, total (serum, plasma)", unit = "mg/dL", win = 3, lo = 0.05, hi = 30, drop_unit = 8749L, batch = 2,
              inc = "^bilirubin\\.total \\[mass/volume\\] in (serum or plasma|serum|plasma|blood)",
              exc = "urine|fluid|neonatal|cord|transcutaneous|presence|\\bpost\\b|\\bpre\\b|ratio|--|percentile",
              broad = "bilirubin", broad_exc = "^$"),
  wbc  = list(label = "leukocytes (blood)", unit = "10^3/uL", win = 3, lo = 0.5, hi = 100, batch = 2,
              alt = c(lo = 500, hi = 100000, f = 0.001), alt_note = "per uL, divided by 1000",
              inc = "^leukocytes \\[#/volume\\] (corrected for nucleated erythrocytes )?in (blood|venous blood|capillary blood)|^white blood cell count$",
              exc = "urine|fluid|csf|stool|synovial|pleural|peritoneal|ascites|product|estimate|fetus|cord|\\bpost\\b|\\bpre\\b|ratio|--|percentile",
              broad = "leukocytes \\[#|white blood cell|\\bwbc\\b", broad_exc = "^$"),
  plt  = list(label = "platelets (blood)", unit = "10^3/uL", win = 3, lo = 5, hi = 1500, batch = 2,
              alt = c(lo = 5000, hi = 1500000, f = 0.001), alt_note = "per uL, divided by 1000",
              inc = "^platelets \\[#/volume\\] in (blood|venous blood|capillary blood)|^platelet count$|^platelet count, blood, (manual|automated)$",
              exc = paste0("urine|fluid|plasma|citrate|giant|large|immature|reticulated|corrected|product|estimate|fetus|cord|",
                           "\\bpost\\b|\\bpre\\b|ratio|--|percentile"),
              broad = "platelets \\[#|platelet count", broad_exc = "^$"),
  rdw  = list(label = "red cell distribution width", unit = "%", win = 3, lo = 8, hi = 40, batch = 2,
              inc = "^erythrocyte distribution width \\[ratio\\]|^red cell distribution width determination$",
              exc = "entitic|fetus|cord|\\bpost\\b|\\bpre\\b|--|percentile",
              broad = "distribution width|\\brdw\\b", broad_exc = "platelet"),
  mcv  = list(label = "mean corpuscular volume", unit = "fL", win = 3, lo = 50, hi = 150, batch = 2,
              inc = "^mcv \\[entitic volume\\]|^erythrocyte mean corpuscular volume( \\[entitic volume\\]| determination)|^mcv - mean corpuscular volume$",
              exc = "reticulocyte|fetus|cord|\\bpost\\b|\\bpre\\b|--|percentile",
              broad = "mean corpuscular volume|\\bmcv\\b", broad_exc = "hemoglobin"),
  tsh  = list(label = "TSH (serum, plasma)", unit = "mIU/L", win = 3, lo = 0.005, hi = 100, batch = 2,
              inc = "^thyrotropin \\[units/volume\\] in (serum or plasma|serum|plasma|blood)|^(serum|plasma) tsh measurement$",
              exc = "receptor|antibod|immunoglobulin|stimulat|\\btrh\\b|\\bpost\\b|\\bpre\\b|ratio|minute|hour|dose|challenge|--|percentile",
              broad = "thyrotropin|\\btsh\\b", broad_exc = "^$")
)
NEW_LABS <- names(NEW)
ALL_NEW_COLS <- c(NEW_LABS, "egfr", paste0("imp_", c(NEW_LABS, "egfr")))   # every column this script can add
sq <- function(x) gsub("'", "\\\\'", x)
NEW_CS_SQL <- paste(vapply(NEW_LABS, function(k) sprintf(paste0(
  "SELECT '%1$s' AS lab, concept_id FROM `%2$s.concept` WHERE domain_id = 'Measurement' AND standard_concept = 'S' ",
  "AND vocabulary_id IN ('LOINC', 'SNOMED') AND REGEXP_CONTAINS(LOWER(concept_name), r'%3$s') ",
  "AND NOT REGEXP_CONTAINS(LOWER(concept_name), r'%4$s')"), k, CDR, sq(NEW[[k]]$inc), sq(NEW[[k]]$exc)), ""), collapse = "\nUNION ALL\n")
fl <- function(x) if (is.null(x) || is.na(x)) "CAST(NULL AS FLOAT64)" else sprintf("CAST(%s AS FLOAT64)", format(x))
## keep_mmol: records in mmol/L are kept (sodium; calcium); every other lab drops them, as 02 does. drop_unit: records in this unit
## are dropped (bilirubin in umol/L: those values overlap the mg/dL range)
RNG_SQL <- paste(vapply(NEW_LABS, function(k) sprintf(paste0("SELECT '%s' AS lab, %s AS lo, %s AS hi, %d AS win, %s AS lo2, %s AS hi2, %s AS f2, ",
                                                            "%s AS keep_mmol, %s AS drop_unit"), k,
                                                     fl(NEW[[k]]$lo), fl(NEW[[k]]$hi), as.integer(NEW[[k]]$win), fl(NEW[[k]]$alt["lo"]),
                                                     fl(NEW[[k]]$alt["hi"]), fl(NEW[[k]]$alt["f"]), if (isTRUE(NEW[[k]]$keep_mmol)) "TRUE" else "FALSE",
                                                     if (is.null(NEW[[k]]$drop_unit)) "CAST(NULL AS INT64)" else sprintf("CAST(%d AS INT64)", as.integer(NEW[[k]]$drop_unit))), ""),
                 collapse = " UNION ALL ")
rng_txt <- function(k) sprintf("values %g-%g used%s%s%s", NEW[[k]]$lo, NEW[[k]]$hi,
                               if (is.null(NEW[[k]]$alt)) "" else sprintf(", %g-%g taken as %s", NEW[[k]]$alt["lo"], NEW[[k]]$alt["hi"],
                                                                         NEW[[k]]$alt_note %||% "g/L, divided by 10"),
                               if (isTRUE(NEW[[k]]$keep_mmol)) "; mmol/L records kept" else "",
                               if (is.null(NEW[[k]]$drop_unit)) "" else sprintf("; unit %d (umol/L) dropped", as.integer(NEW[[k]]$drop_unit)))

## SBP and BMI codes of 02, used here only to recompute the index date (copied from 02; the check below compares every
## recomputed date with the cohort file)
IX_CS <- list(sbp = list(ids = c(21490853, 21492239, 3004249, 3028737), add = c(4152194, 3018586, 4232915, 3035856, 3009395, 4354252)),
              bmi = list(ids = c(3038553, 36304833), add = c(4245997, 44783982)))
ix_cs_sql <- function(k) {
  ids <- paste(IX_CS[[k]]$ids, collapse = ", "); exact <- paste(c(IX_CS[[k]]$ids, IX_CS[[k]]$add), collapse = ", ")
  if (HAS_CRIT) return(sprintf(paste0(
    "SELECT DISTINCT c.concept_id FROM `%1$s.cb_criteria` c JOIN (SELECT CAST(cr.id AS STRING) AS id ",
    "FROM `%1$s.cb_criteria` cr WHERE concept_id IN (%2$s) AND full_text LIKE '%%_rank1]%%') a ",
    "ON (c.path LIKE CONCAT('%%.', a.id, '.%%') OR c.path LIKE CONCAT('%%.', a.id) ",
    "OR c.path LIKE CONCAT(a.id, '.%%') OR c.path = a.id) WHERE is_standard = 1 AND is_selectable = 1 ",
    "UNION DISTINCT SELECT id FROM UNNEST([%3$s]) AS id"), CDR, ids, exact))
  sprintf(paste0("SELECT descendant_concept_id AS concept_id FROM `%1$s.concept_ancestor` WHERE ancestor_concept_id IN (%2$s) ",
                 "UNION DISTINCT SELECT id FROM UNNEST([%3$s]) AS id"), CDR, ids, exact)
}
IX_TAGGED <- paste(vapply(names(IX_CS), function(k) sprintf("SELECT '%s' AS v, concept_id FROM (%s)", k, ix_cs_sql(k)), ""), collapse = "\nUNION ALL\n")

dl <- function(tb, page_size) {
  args <- list(tb, bigint = "numeric", quiet = TRUE, page_size = page_size)
  if ("api" %in% names(formals(bigrquery::bq_table_download))) args$api <- "json"
  do.call(bigrquery::bq_table_download, args)
}
q <- function(sql, what, people = TRUE) {
  cat(sprintf("  %-40s ", what)); t1 <- Sys.time()
  out <- tryCatch({
    tb <- bigrquery::bq_project_query(BILLING, sql, quiet = TRUE)
    for (ps in c(5000, 1000)) {
      r <- tryCatch(dl(tb, ps), error = function(e) e)
      if (!inherits(r, "error") || !grepl("rows were requested", conditionMessage(r))) break
    }
    if (inherits(r, "error")) stop(conditionMessage(r), call. = FALSE)
    r
  }, error = function(e) { cat("FAILED\n"); stop(what, ": ", conditionMessage(e), call. = FALSE) })
  cat(sprintf("%s rows  (%.1f min)\n", if (people) msk(nrow(out)) else format(nrow(out), big.mark = ","),
              as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  as.data.frame(out)
}

hdr("1. Codes of the new labs")
cat("  rule: standard LOINC / SNOMED measurement concepts whose lower-case name matches 'uses' and not 'leaves out'\n")
INV_LABS <- if (is.character(INV_SET)) intersect(INV_SET, NEW_LABS) else if (INV_ON) NEW_LABS else character(0)
for (k in INV_LABS) cat(sprintf("  %-5s %-34s uses: %s\n        %-34s leaves out: %s\n", k, NEW[[k]]$label, NEW[[k]]$inc, "", NEW[[k]]$exc))
if (length(INV_LABS) < length(NEW_LABS)) cat("  (rules of the other labs: see the NEW list in the script)\n")
if (INV_ON) {
  SQL_INV <- sprintf("
WITH ehr AS (%1$s),
rng AS (%2$s),
sel AS (SELECT DISTINCT lab, concept_id FROM (%3$s)),
cand AS (%4$s),
m AS (
  SELECT cand.lab, m.measurement_concept_id AS cid, m.person_id, m.value_as_number AS val, m.unit_concept_id AS u
  FROM `%5$s.measurement` m JOIN cand ON m.measurement_concept_id = cand.concept_id JOIN ehr ON ehr.person_id = m.person_id
  WHERE m.measurement_datetime >= TIMESTAMP('2007-01-01') AND m.measurement_datetime < TIMESTAMP('2020-01-01')),
agg AS (
  SELECT m.lab, m.cid, COUNT(DISTINCT m.person_id) AS persons, APPROX_QUANTILES(m.val, 10)[OFFSET(1)] AS p10,
    APPROX_QUANTILES(m.val, 10)[OFFSET(5)] AS p50, APPROX_QUANTILES(m.val, 10)[OFFSET(9)] AS p90, AVG(IF(m.val IS NULL, 1.0, 0.0)) AS no_num,
    AVG(IF(m.val IS NULL, NULL, IF((m.val >= rng.lo AND m.val <= rng.hi) OR (m.val >= rng.lo2 AND m.val <= rng.hi2), 1.0, 0.0))) AS in_range
  FROM m JOIN rng ON m.lab = rng.lab GROUP BY m.lab, m.cid),
un AS (SELECT lab, cid, u, ROW_NUMBER() OVER (PARTITION BY lab, cid ORDER BY COUNT(*) DESC, u) AS rk FROM m GROUP BY lab, cid, u),
ag2 AS (SELECT agg.*, sel.concept_id IS NOT NULL AS used FROM agg LEFT JOIN sel ON sel.lab = agg.lab AND sel.concept_id = agg.cid),
rk AS (SELECT ag2.*, ROW_NUMBER() OVER (PARTITION BY lab, used ORDER BY persons DESC, cid) AS r_grp FROM ag2)
SELECT rk.lab, rk.cid AS concept_id, c.concept_name, c.vocabulary_id, c.standard_concept, rk.persons, rk.p10, rk.p50, rk.p90, rk.no_num,
  rk.in_range, rk.used,
  uc.concept_name AS unit
FROM rk JOIN `%5$s.concept` c ON c.concept_id = rk.cid
LEFT JOIN un ON un.lab = rk.lab AND un.cid = rk.cid AND un.rk = 1
LEFT JOIN `%5$s.concept` uc ON uc.concept_id = un.u
WHERE (rk.used AND rk.r_grp <= 25) OR (NOT rk.used AND rk.r_grp <= 6)",
    EHR_SQL, RNG_SQL, NEW_CS_SQL,
    paste(c(vapply(NEW_LABS, function(k) sprintf(paste0(
      "SELECT '%1$s' AS lab, concept_id FROM `%2$s.concept` WHERE domain_id = 'Measurement' ",
      "AND REGEXP_CONTAINS(LOWER(concept_name), r'%3$s') AND NOT REGEXP_CONTAINS(LOWER(concept_name), r'%4$s')"),
      k, CDR, sq(NEW[[k]]$broad), sq(NEW[[k]]$broad_exc)), ""), "SELECT lab, concept_id FROM sel"), collapse = "\nUNION DISTINCT\n"),
    CDR)
  INV <- q(SQL_INV, "codes in use (EHR participants, 2007-2019)", people = FALSE)
  NSEL <- q(sprintf("SELECT lab, COUNT(DISTINCT concept_id) AS n FROM (%s) GROUP BY lab", NEW_CS_SQL), "number of codes per new lab", people = FALSE)
  for (k in INV_LABS) {
    x <- INV[INV$lab == k, , drop = FALSE]; x <- x[order(!as.logical(x$used), -x$persons), , drop = FALSE]
    n_sel <- NSEL$n[NSEL$lab == k] %||% 0; xu <- x[as.logical(x$used), , drop = FALSE]; xn <- x[!as.logical(x$used), , drop = FALSE]
    cat(sprintf("\n  %s: %s (%s; %s; window index -%d y .. +1 y) -- %d code(s) match the rule, %s with records\n",
                k, NEW[[k]]$label, NEW[[k]]$unit, rng_txt(k), NEW[[k]]$win, n_sel, if (nrow(xu) >= 25) ">= 25" else nrow(xu)))
    if (!nrow(x)) { cat("      (no records)\n"); next }
    more <- max(0, nrow(xu) - INV_TOP); x <- rbind(utils::head(xu, INV_TOP), utils::head(xn, 3))
    g4 <- function(v) ifelse(is.na(v), "NA", formatC(v, format = "g", digits = 4))
    pc1 <- function(v) ifelse(is.na(v), "-", sprintf("%.1f%%", 100 * v))
    cat(sprintf("      %-5s %-9s %-54s %-6s %9s %8s %8s %8s %9s %9s  %s\n", "used", "concept", "name", "vocab", "persons", "p10", "median", "p90",
                "in range", "no value", "most common unit"))
    for (i in seq_len(nrow(x))) cat(sprintf("      %-5s %-9s %-54s %-6s %9s %8s %8s %8s %9s %9s  %s\n", if (isTRUE(as.logical(x$used[i]))) "yes" else "-",
                                            x$concept_id[i], substr(x$concept_name[i], 1, 54), substr(x$vocabulary_id[i], 1, 6), msk(x$persons[i]),
                                            g4(x$p10[i]), g4(x$p50[i]), g4(x$p90[i]), pc1(x$in_range[i]), pc1(x$no_num[i]), substr(x$unit[i] %||% "NA", 1, 22)))
    if (more > 0) cat(sprintf("      ... and %s more used code(s), each with fewer people than those above\n", if (nrow(xu) >= 25) paste0(more, "+") else more))
  }
  cat("\n  'used = -': the largest codes with a similar name that the rule leaves out. 'in range': share of numeric values inside the\n")
  cat("  range used (low = another unit). 'no value': share of records without a number. p10 / median / p90 of the numeric values.\n")
}

## ---- 2. extraction: index date (recomputed) and the new labs ---------------------------------
hdr("2. Extraction (BigQuery)")
lab_cols <- paste(vapply(NEW_LABS, function(k) sprintf("MAX(IF(lv.lab = '%1$s', lv.val, NULL)) AS %1$s", k), ""), collapse = ",\n  ")
SQL_NEW <- sprintf("
WITH ehr AS (%1$s),
ix_cs AS (%2$s),
ixm AS (
  SELECT lc.v AS lab, m.person_id, m.measurement_datetime AS dt
  FROM `%3$s.measurement` m JOIN ix_cs lc ON m.measurement_concept_id = lc.concept_id JOIN ehr ON ehr.person_id = m.person_id
  WHERE m.measurement_datetime >= TIMESTAMP('2006-12-01') AND COALESCE(m.unit_concept_id, 0) != 8753),
ix0 AS (
  SELECT person_id, DATE(MIN(IF(lab = 'sbp' AND DATE(dt) > DATE '2009-12-31', dt, NULL))) AS s,
                    DATE(MIN(IF(lab = 'bmi' AND DATE(dt) > DATE '2009-12-31', dt, NULL))) AS b
  FROM ixm GROUP BY person_id),
ix AS (
  SELECT person_id,
    CASE
      WHEN s IS NOT NULL AND b IS NOT NULL AND EXTRACT(YEAR FROM s) = EXTRACT(YEAR FROM b) THEN GREATEST(s, b)
      WHEN s IS NULL OR b IS NULL THEN COALESCE(s, b)
      WHEN ABS(DATE_DIFF(s, b, DAY)) > 365 THEN LEAST(s, b)
      ELSE GREATEST(s, b)
    END AS index_date
  FROM ix0),
ixk AS (SELECT * FROM ix WHERE index_date >= DATE '2010-01-01' AND index_date <= DATE '2018-12-31'),   -- each CTE is used once or twice: BigQuery recomputes it per use
rng AS (%4$s),
new_cs AS (SELECT DISTINCT lab, concept_id FROM (%5$s)),
w AS (
  SELECT lc.lab, m.person_id, m.measurement_datetime AS dt,
    IF(m.value_as_number >= rng.lo AND m.value_as_number <= rng.hi, m.value_as_number, m.value_as_number * rng.f2) AS val
  FROM `%3$s.measurement` m JOIN new_cs lc ON m.measurement_concept_id = lc.concept_id
  JOIN rng ON rng.lab = lc.lab
  JOIN ixk ON ixk.person_id = m.person_id
  WHERE m.measurement_datetime >= TIMESTAMP('2006-12-01') AND m.measurement_datetime < TIMESTAMP('2020-01-01')
    AND (rng.keep_mmol OR COALESCE(m.unit_concept_id, 0) != 8753) AND (rng.drop_unit IS NULL OR COALESCE(m.unit_concept_id, 0) != rng.drop_unit)
    AND m.value_as_number IS NOT NULL
    AND ((m.value_as_number >= rng.lo AND m.value_as_number <= rng.hi) OR (m.value_as_number >= rng.lo2 AND m.value_as_number <= rng.hi2))
    AND m.measurement_datetime >= TIMESTAMP(DATE_SUB(ixk.index_date, INTERVAL rng.win YEAR))
    AND m.measurement_datetime <= TIMESTAMP(DATE_ADD(ixk.index_date, INTERVAL 1 YEAR))),
w2 AS (SELECT lab, person_id, dt, val, MAX(dt) OVER (PARTITION BY lab, person_id) AS last_dt FROM w),
lv AS (SELECT lab, person_id, AVG(val) AS val FROM w2 WHERE dt = last_dt GROUP BY lab, person_id)
SELECT ixk.person_id, ixk.index_date,
  %6$s
FROM ixk LEFT JOIN lv ON lv.person_id = ixk.person_id
GROUP BY ixk.person_id, ixk.index_date", EHR_SQL, IX_TAGGED, CDR, RNG_SQL, NEW_CS_SQL, lab_cols)
NEWV <- q(SQL_NEW, sprintf("measurement (index date + %d new labs)", length(NEW_LABS)))
NEWV$index_date <- as.character(as.Date(NEWV$index_date))

## ---- 3. add the columns to the cohort files ---------------------------------------------------
hdr("3. Cohort files")
## the base of each file: the file as 02 wrote it. On a re-run (the file already has the new columns) the copy kept
## in before_extra_labs/ is the base, so the columns are replaced rather than added twice.
dir.create(BAK_DIR, showWarnings = FALSE)
base_lines <- list(); bak <- setNames(file.path(BAK_DIR, basename(F)), names(F))
hdr_of <- function(l) gsub('^"|"$', "", strsplit(l[1], ",", fixed = TRUE)[[1]])
for (k in names(F)) {
  cur <- readLines(F[[k]]); h <- hdr_of(cur)
  if (!any(ALL_NEW_COLS %in% h)) { base_lines[[k]] <- cur; file.copy(F[[k]], bak[[k]], overwrite = TRUE)
  } else {
    if (!file.exists(bak[[k]])) stop(basename(F[[k]]), " already has the new columns but its original is not in ", BAK_DIR,
                                     ": run dementia_cp_02_v9_cohort.R again first.", call. = FALSE)
    b <- readLines(bak[[k]])
    if (length(b) != length(cur) || !identical(hdr_of(b), setdiff(h, ALL_NEW_COLS)))
      stop(basename(F[[k]]), " does not match its copy in ", BAK_DIR, ": run dementia_cp_02_v9_cohort.R again first.", call. = FALSE)
    base_lines[[k]] <- b
  }
}
rd <- function(l) utils::read.csv(text = l, stringsAsFactors = FALSE)
REL <- rd(base_lines$relax); AN <- rd(base_lines$analysis)
cat(sprintf("  relax file:    %s people, %s dementia cases\n", msk(nrow(REL)), msk(sum(REL$Dementia))))
cat(sprintf("  analysis file: %s people, %s dementia cases\n", msk(nrow(AN)), msk(sum(AN$Dementia))))

## index dates: every person of the cohort must get the same index date as in the file
m <- match(REL$person_id, NEWV$person_id)
n_none <- sum(is.na(m)); n_diff <- sum(!is.na(m) & NEWV$index_date[m] != as.character(as.Date(REL$index_date)))
if (n_none + n_diff > 0) {
  pm <- part_mask(c(nrow(REL) - n_none - n_diff, n_none, n_diff))   # no count of 1-20, also not by subtraction
  cat(sprintf("  index date recomputed: NOT the same as the file -- same %s, no date %s, another date %s\n", pm[1], pm[2], pm[3]))
  stop("The recomputed index dates do not match the cohort file (different CDR or changed 02 code lists?). Nothing was written.", call. = FALSE)
}
cat(sprintf("  index date recomputed: the same as the file for all %s people\n", msk(nrow(REL))))
if (!all(AN$person_id %in% REL$person_id)) stop("The analysis file has people who are not in the relax file. Nothing was written.", call. = FALSE)

## values, missing flags, filling (median within race x sex at birth; overall median if a whole group has no value)
X <- NEWV[m, NEW_LABS, drop = FALSE]
none <- NEW_LABS[vapply(NEW_LABS, function(k) !any(is.finite(X[[k]])), logical(1))]
if (length(none)) cat("  NOTE: no value at all in this cohort for:", paste(none, collapse = ", "), "-> not added\n")
NEW_LABS <- setdiff(NEW_LABS, none)
for (k in NEW_LABS) X[[paste0("imp_", k)]] <- as.integer(is.na(X[[k]]))
grp <- paste(REL$race, REL$sex_at_birth, sep = " | ")
fill <- function(x, g) {
  med <- tapply(x, g, stats::median, na.rm = TRUE); v <- med[g]; v[is.na(v)] <- stats::median(x, na.rm = TRUE)
  ifelse(is.na(x), v, x)
}
for (k in NEW_LABS) X[[k]] <- fill(X[[k]], grp)
## eGFR, CKD-EPI 2021 (Inker et al., NEJM 2021), creatinine in mg/dL, age in years
NEW_COLS <- c(NEW_LABS, if ("crea" %in% NEW_LABS) "egfr", paste0("imp_", c(NEW_LABS, if ("crea" %in% NEW_LABS) "egfr")))
if ("crea" %in% NEW_LABS) {
fem <- REL$sex_at_birth %in% "Female"; mal <- REL$sex_at_birth %in% "Male"
kap <- ifelse(fem, 0.7, 0.9); alp <- ifelse(fem, -0.241, -0.302); sr <- X$crea / kap
eg <- 142 * pmin(sr, 1)^alp * pmax(sr, 1)^(-1.200) * 0.9938^REL$age * ifelse(fem, 1.012, 1)
eg[!(fem | mal)] <- NA
X$imp_egfr <- as.integer(X$imp_crea == 1 | !(fem | mal))
## sex at birth known: eGFR from the measured or filled creatinine. Not known: the median of the measured eGFRs
## (a race x sex-at-birth median does not exist for this group, since no one in it has an eGFR)
X$egfr <- ifelse(fem | mal, eg, stats::median(eg[X$imp_egfr == 0], na.rm = TRUE))
}
X <- X[, NEW_COLS, drop = FALSE]
if (anyNA(X)) stop("Some new values are still missing after filling. Nothing was written.", call. = FALSE)

## write: the old part of every line is kept as it is; the new columns are added at the end
fmt <- function(v) { s <- trimws(formatC(v, digits = 10, format = "fg")); s[is.na(v)] <- "NA"; s }
add_cols <- function(lines, rows) {
  body <- vapply(NEW_COLS, function(cn) fmt(X[[cn]][rows]), character(length(rows)))
  if (is.null(dim(body))) body <- matrix(body, nrow = length(rows))
  c(paste0(lines[1], ",", paste0('"', NEW_COLS, '"', collapse = ",")), paste0(lines[-1], ",", apply(body, 1, paste, collapse = ",")))
}
OUT <- list(relax = add_cols(base_lines$relax, seq_len(nrow(REL))), analysis = add_cols(base_lines$analysis, match(AN$person_id, REL$person_id)))
for (k in names(F)) {
  tmp <- paste0(F[[k]], ".tmp"); writeLines(OUT[[k]], tmp)
  chk <- utils::read.csv(tmp, stringsAsFactors = FALSE); orig <- rd(base_lines[[k]])
  same <- identical(names(chk)[seq_len(ncol(orig))], names(orig)) && nrow(chk) == nrow(orig) &&
    all(vapply(names(orig), function(cn) isTRUE(all.equal(chk[[cn]], orig[[cn]], check.attributes = FALSE)), logical(1)))
  if (!same) { unlink(tmp); stop("Check failed: the old columns of ", basename(F[[k]]), " would change. Nothing was written.", call. = FALSE) }
  file.rename(tmp, F[[k]])
}
cat(sprintf("  old columns: unchanged (checked after writing); new columns: %d labs + egfr, each with its imp_ flag\n", length(setdiff(NEW_COLS[!grepl("^imp_", NEW_COLS)], "egfr"))))
f_def <- file.path(OUT_DIR, "dementia_v9_extra_labs_definitions.txt")
writeLines(c(paste0("Extra labs added ", format(Sys.time(), "%Y-%m-%d %H:%M UTC", tz = "UTC"), " from ", CDR, " by dementia_cp_11_extra_labs.R"),
             vapply(NEW_LABS, function(k) sprintf("%s: %s, %s; latest value in index -%d y .. +1 y; %s", k, NEW[[k]]$label,
                                                  NEW[[k]]$unit, NEW[[k]]$win, rng_txt(k)), ""),
             "egfr: CKD-EPI 2021 creatinine equation (no race term), age at index, sex at birth",
             "imp_<lab> = 1: missing in the window (or, for egfr, sex at birth not female/male), filled with the median within race x sex at birth",
             "Code rules: see the header and the NEW list of dementia_cp_11_extra_labs.R"), f_def)
if (length(BUCKETS)) {
  dest <- paste0(BUCKETS[1], "/dementia_v9/", CDR_ID, "/")
  for (f in c(F, f_def)) { r <- sh("gsutil", c("-q", "cp", shQuote(f), shQuote(dest)))
    cat(sprintf("  %s -> %s%s\n", basename(f), dest, if (ok(r)) "" else "   COPY FAILED, file is on this app's disk only")) }
  for (f in bak) { r <- sh("gsutil", c("-q", "cp", shQuote(f), shQuote(paste0(dest, "before_extra_labs/"))))
    if (!ok(r)) cat("  copy of the original", basename(f), "to the bucket FAILED (it is in", BAK_DIR, ")\n") }
} else cat("  no bucket found; files are only on this app's disk:", OUT_DIR, "\n")

## ---- 4. report: the analysis set --------------------------------------------------------------
hdr("4. New labs in the analysis set")
A <- utils::read.csv(F[["analysis"]], stringsAsFactors = FALSE); N <- nrow(A); cs <- A$Dementia == 1
cat(sprintf("  N %s, dementia cases %s. 'measured' = a value in the window (not filled).\n\n", msk(N), msk(sum(cs))))
cat(sprintf("    %-5s %-34s %18s %16s %20s %22s\n", "", "", "measured: all", "cases", "median (IQR) measured", "1st - 99th percentile"))
show <- c("sbp", "bmi", "ldlc", "hdl", "trig", NEW_LABS, "egfr")
LBL <- c(sbp = "systolic blood pressure", bmi = "body mass index", ldlc = "LDL cholesterol", hdl = "HDL cholesterol", trig = "triglycerides",
         vapply(NEW, function(z) z$label, ""), egfr = "eGFR (CKD-EPI 2021)")
cov_tab <- NULL
for (l in show) {
  fl <- paste0("imp_", l); if (!fl %in% names(A)) next
  me <- A[[fl]] == 0; k1 <- sum(me); k2 <- sum(me & cs)
  s1 <- part_mask(c(k1, N - k1))[1]; s2 <- part_mask(c(k2, sum(cs) - k2))[1]
  pc <- function(s, a, b) if (grepl("<=|hidden", s)) s else sprintf("%5.1f%%", 100 * a / b)
  x <- A[[l]][me]; qq <- if (length(x) > 20) stats::quantile(x, c(0.01, 0.25, 0.5, 0.75, 0.99)) else rep(NA, 5)
  cat(sprintf("    %-5s %-34s %18s %16s %20s %22s\n", l, substr(LBL[[l]], 1, 34), pc(s1, k1, N), pc(s2, k2, sum(cs)),
              if (anyNA(qq)) "-" else sprintf("%.3g (%.3g-%.3g)", qq[3], qq[2], qq[4]), if (anyNA(qq)) "-" else sprintf("%.3g - %.3g", qq[1], qq[5])))
  cov_tab <- rbind(cov_tab, data.frame(lab = l, measured = k1 / N, measured_cases = if (sum(cs)) k2 / sum(cs) else NA))
}
cat("\n  Pairs with |r| >= 0.3 between measured values (people with both measured; Pearson):\n")
cm <- show[paste0("imp_", show) %in% names(A)]; pr <- NULL
for (i in seq_along(cm)) for (j in seq_len(i - 1L)) { a <- cm[j]; b <- cm[i]; u <- A[[paste0("imp_", a)]] == 0 & A[[paste0("imp_", b)]] == 0
  r <- if (sum(u) > 100) stats::cor(A[[a]][u], A[[b]][u]) else NA_real_
  if (isTRUE(abs(r) >= 0.3)) pr <- rbind(pr, data.frame(a = a, b = b, r = r)) }
if (is.null(pr)) cat("    none\n") else { pr <- pr[order(-abs(pr$r)), ]; txt <- sprintf("%s-%s %+.2f", pr$a, pr$b, pr$r)
  cat(paste0("    ", vapply(split(txt, ceiling(seq_along(txt) / 8)), paste, "", collapse = " | "), "\n"), sep = "") }
cat("\n  Next: dementia_cp_12_panel.R (which of these biomarkers enter the change-point analysis is decided there, by the rules\n")
cat("  written in its header, before any change point is looked at).\n")
cat(sprintf("\n  CDR %s | runtime %.1f min\n", CDR, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
