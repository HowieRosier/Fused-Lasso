## =============================================================================
## dementia_cp_02_v9_cohort.R
## Dementia change-point paper -- analysis cohort on the current All of Us data
## (v9, C2025Q4R6) with corrected variable definitions.
##
## Same design as the paper (Ziye's Dementia_Cohort.R):
##   EHR participants with an outpatient visit (9202); no dementia diagnosis before 2010;
##   index date = first SBP / BMI after 2009 (same combination rule); index year 2010-2018;
##   lab windows SBP index -1 y..+1 y, other labs -3 y..+1 y; missing labs -> median within
##   race x sex at birth; age 30-79, SBP 90-200, TC 130-320, HDL 20-100, BMI 18.5-40 at index;
##   an outpatient visit >= index + 7 days; time2 = years to the first dementia diagnosis,
##   else to the last outpatient visit; covariates = any record on or before the index date.
##   Analysis set = time2 > 0, age 45-79, LDL < 400, TG < 700 (the paper's Table 1 rule).
## Corrected:
##   1. dementia and the condition covariates include all descendant codes (concept_ancestor);
##   2. lipid labs use mass-concentration codes only (no percentile, presence, ratio or mmol/L codes),
##      values recorded in mmol/L are dropped, and codes that carry the same labs in v9 are added as
##      they are, without their sub-codes (SNOMED-coded SBP / BMI / HDL / TG / TC and similar; see the CS list);
##   3. the latest NON-missing value in the window is used (the old code took the latest
##      record even when its value was empty);
##   4. index dates on Feb 29 keep their lab windows (the old code lost them);
##   5. smoking = a current or former smoking status recorded on or before index (the old code
##      counted any smoking-status record, "never smoker" included);
##   6. antihypertensive drugs = ATC C02, C03, C07, C08, C09 (the old list held 13 diuretics);
##   7. social isolation is also searched in the observation table (its code is an observation).
##   Kept as they were (near zero in the paper): physical inactivity, alcohol, low education.
##
## Participant-level files go only to this workspace's bucket; the console shows aggregates
## only, and participant counts of 1-20 never print, directly or by subtraction.
## =============================================================================

## Settings. Edit here, or set them in the console before source().
if (!exists("CDR"))  CDR  <- ""    # "" = the latest All of Us CDR attached (C2025Q4R6 here); or a full path
if (!exists("SAVE")) SAVE <- TRUE  # write the cohort files to the workspace bucket
if (!exists("IMPUTE")) IMPUTE <- "median"  # missing labs: "median" (Ziye's code; the paper's cohort matches it) or "mean" (sensitivity)

t0 <- Sys.time()
options(width = 160, scipen = 999, warn = 1)
for (p in c("jsonlite", "bigrquery", "dplyr", "lubridate")) {
  if (!requireNamespace(p, quietly = TRUE)) utils::install.packages(p, quiet = TRUE)
}
suppressPackageStartupMessages({ library(dplyr); library(lubridate) })

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
hdr <- function(s) cat("\n==================== ", s, " ====================\n", sep = "")
## AoU dissemination rule: no participant count of 1-20 may leave the Workbench, including counts that could be
## recovered by subtracting printed numbers. msk() hides 1-20; part_mask() also hides the next-smallest cell of a
## table whose total is printed; chain_mask() hides a step that differs from the last printed step by 1-20.
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
chain_mask <- function(n) {
  out <- character(length(n)); last <- NA
  for (k in seq_along(n)) {
    if (SMALL(n[k])) out[k] <- "<=20"
    else if (!is.na(last) && SMALL(abs(last - n[k]))) out[k] <- "~above (1-20)"
    else { out[k] <- msk(n[k]); last <- n[k] }
  }
  out
}
exact_shown <- function(x) !grepl("<=|1-20|hidden", x)
sh <- function(cmd, args = character(), err = TRUE) {
  if (!nzchar(Sys.which(cmd))) return(structure(character(), status = 127L))
  tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = err)),
           error = function(e) structure(paste("ERROR:", conditionMessage(e)), status = 1L))
}
ok <- function(x) { s <- attr(x, "status"); is.null(s) || identical(as.integer(s), 0L) }

## ---- 0. workspace, data, billing ---------------------------------------------------
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

cat("  billing project:   ", BILLING, "\n")
cat("  attached datasets: ", if (length(BQ_PATHS)) paste(BQ_PATHS, collapse = ", ") else "<none found>", "\n")
cat("  buckets:           ", if (length(BUCKETS)) paste(BUCKETS, collapse = ", ") else "<none found>", "\n")
if (!nzchar(CDR)) {
  hit <- BQ_PATHS[grepl("\\.C[0-9]{4}Q[0-9]R[0-9]+$", BQ_PATHS)]
  hit <- hit[order(sub("^.*\\.", "", hit), decreasing = TRUE)]
  if (!length(hit)) stop("No All of Us CDR dataset (named like C2025Q4R6) is attached to this workspace.", call. = FALSE)
  CDR <- hit[1]
}
CDR_ID <- sub("^.*\\.", "", CDR)
PR <- strsplit(CDR, ".", fixed = TRUE)[[1]]
cat("  CDR used:          ", CDR, "\n")

tab_exists <- function(path, table) {
  pp <- strsplit(path, ".", fixed = TRUE)[[1]]
  tryCatch(isTRUE(bigrquery::bq_table_exists(bigrquery::bq_table(pp[1], pp[2], table))), error = function(e) FALSE)
}
HAS_CRIT <- tab_exists(CDR, "cb_criteria")
if (tab_exists(CDR, "cb_search_person")) {
  EHR_SQL <- sprintf("SELECT person_id FROM `%s.cb_search_person` WHERE has_ehr_data = 1", CDR)
  EHR_SRC <- "cb_search_person.has_ehr_data = 1 (as in the paper)"
} else {
  idx <- BQ_PATHS[startsWith(sub("^.*\\.", "", BQ_PATHS), paste0(CDR_ID, "_index"))]
  if (!length(idx)) {
    ds <- tryCatch(bigrquery::bq_project_datasets(PR[1]), error = function(e) list())
    nm <- vapply(ds, function(d) d$dataset %||% "", "")
    nm <- sort(nm[startsWith(nm, paste0(CDR_ID, "_index"))], decreasing = TRUE)
    if (length(nm)) idx <- paste0(PR[1], ".", nm[1])
  }
  if (length(idx) && tab_exists(idx[1], "T_ENT_person")) {
    EHR_SQL <- sprintf("SELECT id AS person_id FROM `%s.T_ENT_person` WHERE has_ehr_data = TRUE", idx[1])
    EHR_SRC <- paste0(idx[1], ".T_ENT_person.has_ehr_data")
  } else {
    EHR_SQL <- sprintf("SELECT person_id FROM `%s.person`", CDR)
    EHR_SRC <- "NONE: no has-EHR flag found, all persons kept"
  }
}
cat("  concept sets:       diagnoses with all descendants; labs via ", if (HAS_CRIT) "cb_criteria paths (as the old Dataset Builder)" else "concept_ancestor", " + added v9 codes as they are",
    "; antihypertensives via ATC\n", sep = "")
cat("  has-EHR filter:    ", EHR_SRC, "\n")
cat("  missing labs:       race x sex-at-birth", IMPUTE, "\n")

## ---- 1. concept sets ---------------------------------------------------------------
CS <- list(
  dementia   = list(d = "condition", ids = c(374888, 378419, 4182210)),
  diabetes   = list(d = "condition", ids = c(192279, 195771, 201254, 201820, 201826, 37016349, 37016354, 37017432, 37311673,
                                             376065, 376112, 376683, 376979, 378743, 4008576, 4024659, 4029423, 4034964, 4044391,
                                             4058243, 4082346, 4159742, 4174977, 4175440, 4193704, 4221495, 4226121, 42536605,
                                             4311629, 4311708, 435216, 43530656, 43530689, 43530690, 43531578, 442793, 443412,
                                             443727, 443730, 443731, 443732, 45757277, 45757363, 45757435)),
  hearing    = list(d = "condition", ids = c(374366, 377889, 4110815)),
  isolation  = list(d = "condition", ids = c(4309238)),
  tbi        = list(d = "condition", ids = c(4001336, 4132546)),
  depression = list(d = "condition", ids = c(4077577, 4098302, 4149320, 4149321, 4151170, 4152280, 4176002, 4282096, 4282316,
                                             42872722, 4307111, 4327337, 4336957, 433991, 435220, 43531624, 440383)),
  ## labs: ids = the paper's mass-concentration codes, expanded through cb_criteria as the old Dataset Builder did;
  ## add = codes that carry the same measurement in v9 (SNOMED-coded or position-specific SBP, SNOMED BMI / HDL /
  ## TG / TC, TG in blood or fasting), checked with dementia_cp_04_concept_check.R (values and units match) and used
  ## as they are: their sub-codes are not pulled in (SNOMED "Cholesterol measurement" has HDL, LDL and ratio codes
  ## under it). Not used: All of Us enrollment physical measurements (PPI), ambulatory or exercise BP, ratios,
  ## percentiles, molar codes.
  sbp        = list(d = "measurement", ids = c(21490853, 21492239, 3004249, 3028737),
                    add = c(4152194, 3018586, 4232915, 3035856, 3009395, 4354252)),
  bmi        = list(d = "measurement", ids = c(3038553, 36304833), add = c(4245997, 44783982)),
  tc         = list(d = "measurement", ids = c(3027114, 3049873), add = c(4299360)),
  hdl        = list(d = "measurement", ids = c(3007070, 3013473, 3053286), add = c(4042059)),
  ldl        = list(d = "measurement", ids = c(3009966, 3028288, 3028437, 3035899, 3053341, 4012479), add = c(3035009)),
  trig       = list(d = "measurement", ids = c(3022192, 4032789), add = c(4017787, 3022038, 3027997, 3007943)),
  smoker     = list(d = "observation", ids = c(43054909)),
  inactivity = list(d = "observation", ids = c(4195380, 44808782)),
  alcohol    = list(d = "observation", ids = c(4052460)),
  education  = list(d = "observation", ids = c(42528763))
)
ATC <- c("C02", "C03", "C07", "C08", "C09")
SMOKE_YES <- "current|former|ex-smoker|every day|some day|heavy|light|smoker|smokes|smoked|stopped smoking|quit smoking|tobacco user|^yes$"
SMOKE_NO  <- "never|non-smoker|nonsmoker|non smoker|unknown if ever|not a smoker|passive|second hand|secondhand|doesn.t|does not|don.t|not use"
LOW_EDU <- c("10th grade", "11th grade", "12th grade, no diploma", "2nd grade", "3rd grade", "5th grade", "6th grade",
             "7th grade", "8th grade", "9th grade", "GED or equivalent", "High school graduate",
             "Kindergarten", "No formal schooling", "Nursery school to 8th grade", "No matching concept",
             "Less than high school")

cs_sql <- function(k) {
  ids <- paste(CS[[k]]$ids, collapse = ", "); dom <- CS[[k]]$d
  exact <- paste(c(CS[[k]]$ids, CS[[k]]$add), collapse = ", ")  # seed and added codes themselves (added codes are not expanded)
  if (dom == "observation") return(sprintf("SELECT id AS concept_id FROM UNNEST([%s]) AS id", exact))
  if (dom == "measurement" && HAS_CRIT) return(sprintf(paste0(
    "SELECT DISTINCT c.concept_id FROM `%1$s.cb_criteria` c JOIN (SELECT CAST(cr.id AS STRING) AS id ",
    "FROM `%1$s.cb_criteria` cr WHERE concept_id IN (%2$s) AND full_text LIKE '%%_rank1]%%') a ",
    "ON (c.path LIKE CONCAT('%%.', a.id, '.%%') OR c.path LIKE CONCAT('%%.', a.id) ",
    "OR c.path LIKE CONCAT(a.id, '.%%') OR c.path = a.id) WHERE is_standard = 1 AND is_selectable = 1 ",
    "UNION DISTINCT SELECT id FROM UNNEST([%3$s]) AS id"), CDR, ids, exact))
  sprintf(paste0("SELECT descendant_concept_id AS concept_id FROM `%1$s.concept_ancestor` WHERE ancestor_concept_id IN (%2$s) ",
                 "UNION DISTINCT SELECT id FROM UNNEST([%3$s]) AS id"), CDR, ids, exact)
}
ANTI_SQL <- sprintf(paste0("SELECT DISTINCT ca.descendant_concept_id AS concept_id FROM `%1$s.concept_ancestor` ca ",
                           "JOIN `%1$s.concept` c ON ca.ancestor_concept_id = c.concept_id ",
                           "WHERE c.vocabulary_id = 'ATC' AND c.concept_code IN (%2$s)"), CDR, paste0("'", ATC, "'", collapse = ", "))
tagged <- function(keys) paste(vapply(keys, function(k) sprintf("SELECT '%s' AS v, concept_id FROM (%s)", k, cs_sql(k)), ""),
                               collapse = "\nUNION ALL\n")

## bq_table_download() fetches results in pages; a page over the API's size limit comes back short and
## bigrquery then stops ("... rows were requested, but only ... received"). Small fixed pages avoid that.
dl <- function(tb, page_size) {
  args <- list(tb, bigint = "numeric", quiet = TRUE, page_size = page_size)
  if ("api" %in% names(formals(bigrquery::bq_table_download))) args$api <- "json"
  do.call(bigrquery::bq_table_download, args)
}
q <- function(sql, what, people = TRUE) {
  cat(sprintf("  %-34s", what))
  t1 <- Sys.time()
  out <- tryCatch({
    tb <- bigrquery::bq_project_query(BILLING, sql, quiet = TRUE)
    for (ps in c(5000, 1000)) {
      res <- tryCatch(dl(tb, ps), error = function(e) e)
      if (!inherits(res, "error") || !grepl("rows were requested", conditionMessage(res))) break
    }
    if (inherits(res, "error")) stop(conditionMessage(res), call. = FALSE)
    res
  }, error = function(e) { cat("FAILED\n"); stop(what, ": ", conditionMessage(e), call. = FALSE) })
  cat(sprintf("%s rows  (%.1f min)\n", if (people) msk(nrow(out)) else format(nrow(out), big.mark = ","),
              as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  as.data.frame(out)
}

hdr("1. Concept sets")
SEEDS <- unique(unlist(lapply(CS, function(x) c(x$ids, x$add))))
NAMES <- q(sprintf("SELECT concept_id, concept_name, domain_id, vocabulary_id, standard_concept FROM `%s.concept` WHERE concept_id IN (%s)",
                   CDR, paste(SEEDS, collapse = ", ")), "seed concept names", people = FALSE)
SIZES <- q(paste(c(vapply(names(CS), function(k) sprintf("SELECT '%s' AS v, COUNT(DISTINCT concept_id) AS n FROM (%s)", k, cs_sql(k)), ""),
                   sprintf("SELECT 'anti' AS v, COUNT(DISTINCT concept_id) AS n FROM (%s)", ANTI_SQL)),
                 collapse = "\nUNION ALL\n"), "expanded concept-set sizes", people = FALSE)
ATCN <- q(sprintf("SELECT concept_code, concept_name FROM `%s.concept` WHERE vocabulary_id = 'ATC' AND concept_code IN (%s)",
                  CDR, paste0("'", ATC, "'", collapse = ", ")), "ATC classes", people = FALSE)
for (k in names(CS)) {
  all_ids <- c(CS[[k]]$ids, CS[[k]]$add); n_add <- length(CS[[k]]$add)
  nm <- NAMES[match(all_ids, NAMES$concept_id), ]
  miss <- sum(is.na(nm$concept_id))
  cat(sprintf("  %-11s %2d seed(s)%s -> %s concept(s)%s\n", k, length(CS[[k]]$ids),
              if (n_add) sprintf(" + %d added", n_add) else "",
              SIZES$n[SIZES$v == k] %||% "?", if (miss) sprintf("  [%d id(s) not in this CDR]", miss) else ""))
  show <- if (k %in% c("diabetes", "depression")) head(nm, 3) else nm
  for (i in seq_len(nrow(show))) if (!is.na(show$concept_id[i]))
    cat(sprintf("    %s %-9s %-55s %-11s %s\n", if (all_ids[i] %in% CS[[k]]$add) "+" else " ", show$concept_id[i],
                substr(show$concept_name[i], 1, 55), show$domain_id[i], show$vocabulary_id[i]))
  if (k %in% c("diabetes", "depression")) cat(sprintf("      ... %d more\n", nrow(nm) - 3))
}
cat(sprintf("  %-11s ATC %s -> %s drug concept(s)\n", "anti", paste(ATC, collapse = ", "), SIZES$n[SIZES$v == "anti"] %||% "?"))
for (i in seq_len(nrow(ATCN))) cat(sprintf("      %-9s %s\n", ATCN$concept_code[i], ATCN$concept_name[i]))

SMOKE_VALS <- q(sprintf("
WITH ehr AS (%1$s)
SELECT IFNULL(c.concept_name, '<no value>') AS value_name, COUNT(DISTINCT o.person_id) AS n
FROM `%2$s.observation` o LEFT JOIN `%2$s.concept` c ON o.value_as_concept_id = c.concept_id
WHERE o.observation_concept_id IN (%3$s) AND o.person_id IN (SELECT person_id FROM ehr)
GROUP BY 1 ORDER BY 2 DESC LIMIT 30", EHR_SQL, CDR, paste(CS$smoker$ids, collapse = ", ")), "smoking-status answers", people = FALSE)
cat("  smoking-status answers (persons with the answer at any time; counted as smoker?)\n")
is_smoker <- function(x) { x <- tolower(x); grepl(SMOKE_YES, x) & !grepl(SMOKE_NO, x) }
for (i in seq_len(nrow(SMOKE_VALS)))
  cat(sprintf("      %-55s %10s   %s\n", substr(SMOKE_VALS$value_name[i], 1, 55), msk(SMOKE_VALS$n[i]),
              if (is_smoker(SMOKE_VALS$value_name[i])) "smoker" else "-"))

## ---- 2. extraction: one query per OMOP table -----------------------------------------
hdr("2. Extraction (BigQuery)")
LAB <- c("sbp", "bmi", "tc", "hdl", "ldl", "trig")
lab_agg <- paste(vapply(LAB, function(l) sprintf(
  "ARRAY_AGG(IF(lab = '%1$s' AND in_win AND val IS NOT NULL, STRUCT(dt, val), NULL) IGNORE NULLS ORDER BY dt DESC LIMIT 1) AS a_%1$s", l), ""),
  collapse = ",\n  ")
lab_out <- paste(vapply(LAB, function(l) sprintf("a_%1$s[SAFE_OFFSET(0)].val AS %1$s", l), ""), collapse = ",\n  ")
## 8753 = mmol/L: none of these labs is analysed in molar units here
SQL_MEAS <- sprintf("
WITH ehr AS (%1$s),
lab_cs AS (SELECT v AS lab, concept_id FROM (%2$s)),
meas AS (
  SELECT lc.lab, m.person_id, m.measurement_datetime AS dt, m.value_as_number AS val
  FROM `%3$s.measurement` m JOIN lab_cs lc ON m.measurement_concept_id = lc.concept_id
  WHERE m.measurement_datetime >= TIMESTAMP('2006-12-01') AND COALESCE(m.unit_concept_id, 0) != 8753
    AND m.person_id IN (SELECT person_id FROM ehr)),
e AS (
  SELECT *,
    MIN(IF(lab = 'sbp' AND DATE(dt) > DATE '2009-12-31', dt, NULL)) OVER (PARTITION BY person_id) AS e_sbp,
    MIN(IF(lab = 'bmi' AND DATE(dt) > DATE '2009-12-31', dt, NULL)) OVER (PARTITION BY person_id) AS e_bmi
  FROM meas),
ix AS (
  SELECT lab, person_id, dt, val,
    CASE
      WHEN s IS NOT NULL AND b IS NOT NULL AND EXTRACT(YEAR FROM s) = EXTRACT(YEAR FROM b) THEN GREATEST(s, b)
      WHEN s IS NULL OR b IS NULL THEN COALESCE(s, b)
      WHEN ABS(DATE_DIFF(s, b, DAY)) > 365 THEN LEAST(s, b)
      ELSE GREATEST(s, b)
    END AS index_date
  FROM (SELECT *, DATE(e_sbp) AS s, DATE(e_bmi) AS b FROM e)),
w AS (
  SELECT *,
    dt >= TIMESTAMP(DATE_SUB(index_date, INTERVAL IF(lab = 'sbp', 1, 3) YEAR))
      AND dt <= TIMESTAMP(DATE_ADD(index_date, INTERVAL 1 YEAR)) AS in_win
  FROM ix WHERE index_date IS NOT NULL),
agg AS (
  SELECT person_id, ANY_VALUE(index_date) AS index_date,
  %4$s
  FROM w GROUP BY person_id)
SELECT person_id, index_date,
  %5$s
FROM agg", EHR_SQL, tagged(LAB), CDR, lab_agg, lab_out)

COND <- c("dementia", "diabetes", "hearing", "isolation", "tbi", "depression")
SQL_COND <- sprintf("
WITH ehr AS (%1$s), cs AS (%2$s)
SELECT co.person_id,
  MIN(IF(cs.v = 'dementia', co.condition_start_datetime, NULL)) AS dem_first,
  MAX(IF(cs.v = 'dementia', co.condition_start_datetime, NULL)) AS dem_last,
  LOGICAL_OR(cs.v = 'dementia' AND co.condition_start_datetime < TIMESTAMP('2010-01-01')) AS dem_before_2010,
  %3$s
FROM `%4$s.condition_occurrence` co JOIN cs ON co.condition_concept_id = cs.concept_id
WHERE co.person_id IN (SELECT person_id FROM ehr)
GROUP BY co.person_id", EHR_SQL, tagged(COND),
  paste(vapply(COND[-1], function(k) sprintf("MIN(IF(cs.v = '%1$s', co.condition_start_datetime, NULL)) AS %1$s_first", k), ""),
        collapse = ",\n  "), CDR)

SQL_DRUG <- sprintf("
WITH ehr AS (%1$s), cs AS (%2$s)
SELECT de.person_id, MIN(de.drug_exposure_start_datetime) AS anti_first
FROM `%3$s.drug_exposure` de
WHERE de.drug_concept_id IN (SELECT concept_id FROM cs) AND de.person_id IN (SELECT person_id FROM ehr)
GROUP BY de.person_id", EHR_SQL, ANTI_SQL, CDR)

SQL_OBS <- sprintf("
WITH ehr AS (%1$s), iso AS (%2$s)
SELECT o.person_id,
  MIN(IF(o.observation_concept_id IN (%3$s)
         AND REGEXP_CONTAINS(LOWER(IFNULL(c.concept_name, '')), r'%4$s')
         AND NOT REGEXP_CONTAINS(LOWER(IFNULL(c.concept_name, '')), r'%5$s'), o.observation_datetime, NULL)) AS smoker_first,
  MIN(IF(o.observation_concept_id IN (%6$s), o.observation_datetime, NULL)) AS inactivity_first,
  MIN(IF(iso.concept_id IS NOT NULL, o.observation_datetime, NULL)) AS isolation_obs_first
FROM `%7$s.observation` o
LEFT JOIN `%7$s.concept` c ON o.value_as_concept_id = c.concept_id
LEFT JOIN iso ON o.observation_concept_id = iso.concept_id
WHERE (o.observation_concept_id IN (%3$s, %6$s) OR iso.concept_id IS NOT NULL) AND o.person_id IN (SELECT person_id FROM ehr)
GROUP BY o.person_id", EHR_SQL, cs_sql("isolation"), paste(CS$smoker$ids, collapse = ", "), SMOKE_YES, SMOKE_NO,
  paste(CS$inactivity$ids, collapse = ", "), CDR)

SQL_OBS_RAW <- sprintf("
WITH ehr AS (%1$s)
SELECT o.person_id, o.observation_concept_id, o.observation_datetime, c.concept_name AS value_as_concept_name
FROM `%2$s.observation` o LEFT JOIN `%2$s.concept` c ON o.value_as_concept_id = c.concept_id
WHERE o.observation_concept_id IN (%3$s, %4$s) AND o.observation_datetime >= TIMESTAMP('2006-12-01')
  AND o.person_id IN (SELECT person_id FROM ehr)", EHR_SQL, CDR, paste(CS$alcohol$ids, collapse = ", "),
  paste(CS$education$ids, collapse = ", "))

SQL_VISIT <- sprintf("
WITH ehr AS (%1$s)
SELECT v.person_id, MAX(v.visit_start_datetime) AS visit_last
FROM `%2$s.visit_occurrence` v
WHERE v.visit_concept_id = 9202 AND v.person_id IN (SELECT person_id FROM ehr)
GROUP BY v.person_id", EHR_SQL, CDR)

SQL_PERSON <- sprintf("
WITH ehr AS (%1$s)
SELECT p.person_id, p.birth_datetime AS date_of_birth,
  g.concept_name AS gender, r.concept_name AS race, e.concept_name AS ethnicity, s.concept_name AS sex_at_birth
FROM `%2$s.person` p
LEFT JOIN `%2$s.concept` g ON p.gender_concept_id = g.concept_id
LEFT JOIN `%2$s.concept` r ON p.race_concept_id = r.concept_id
LEFT JOIN `%2$s.concept` e ON p.ethnicity_concept_id = e.concept_id
LEFT JOIN `%2$s.concept` s ON p.sex_at_birth_concept_id = s.concept_id
WHERE p.person_id IN (SELECT person_id FROM ehr)", EHR_SQL, CDR)

MEAS <- q(SQL_MEAS, "measurement (index date + 6 labs)")
CONDX <- q(SQL_COND, "condition_occurrence")
DRUG <- q(SQL_DRUG, "drug_exposure (antihypertensives)")
OBS <- q(SQL_OBS, "observation (smoking, inactivity, isolation)")
OBSR <- q(SQL_OBS_RAW, "observation (alcohol, education)")
VIS <- q(SQL_VISIT, "visit_occurrence (9202)")
PER <- q(SQL_PERSON, "person")

## ---- 3. cohort -----------------------------------------------------------------------
LABS_R <- c("sbp", "bmi", "tc", "hdl", "ldlc", "trig")
build <- function() {
  att <- list()
  note <- function(step, n) att[[length(att) + 1]] <<- data.frame(step = step, n = n)
  d <- data.frame(person_id = unique(VIS$person_id)) %>% left_join(CONDX, by = "person_id")
  note("outpatient visit (9202), has EHR", nrow(d))
  d <- d %>% filter(!coalesce(dem_before_2010, FALSE))
  note("- dementia dx before 2010-01-01", nrow(d))
  d <- d %>% mutate(Dementia = as.numeric(!is.na(dem_first))) %>% inner_join(MEAS, by = "person_id")
  note("index date (first SBP/BMI after 2009)", nrow(d))
  d <- d %>% filter(is.na(dem_last) | as.Date(dem_last) >= index_date)
  note("no dementia dx, or a dx on/after index", nrow(d))
  d <- d %>% rename(ldlc = ldl) %>% left_join(PER, by = "person_id") %>%
    filter(year(index_date) >= 2010, year(index_date) <= 2018)
  note("index year 2010-2018", nrow(d))
  for (l in LABS_R) d[[paste0("imp_", l)]] <- as.integer(is.na(d[[l]]))
  # missing labs <- median within race x sex at birth (Ziye's code; the paper's Table 1 matches it), or the mean.
  # The mean only uses values inside the ranges applied below (SBP 90-200, TC 130-320, HDL 20-100, BMI 18.5-40,
  # LDL < 400, TG < 700): otherwise a few implausible EHR values pull a group's mean outside those ranges and the
  # range filter then drops every imputed member of that group (seen on v9: ~49,000 people, mostly men).
  OK_RANGE <- list(sbp = c(90, 200), tc = c(130, 320), hdl = c(20, 100), bmi = c(18.5, 40), ldlc = c(1e-9, 400), trig = c(1e-9, 700))
  fill <- function(x, lab) {
    if (!identical(IMPUTE, "mean")) return(median(x, na.rm = TRUE))
    r <- OK_RANGE[[lab]]; mean(x[!is.na(x) & x >= r[1] & x <= r[2]])
  }
  d <- d %>% group_by(race, sex_at_birth) %>%
    mutate(across(all_of(LABS_R), ~ ifelse(is.na(.x), fill(.x, cur_column()), .x), .names = "{.col}")) %>% ungroup()
  d <- d %>% mutate(age = as.numeric(difftime(as.POSIXct(index_date), as_datetime(date_of_birth), units = "days")) / 365.25)
  d <- d %>% filter(age >= 30 & age <= 79, sbp >= 90 & sbp <= 200, tc >= 130 & tc <= 320,
                    hdl >= 20 & hdl <= 100, bmi >= 18.5 & bmi < 40)
  note("age 30-79 + SBP/TC/HDL/BMI ranges", nrow(d))
  d <- d %>% left_join(VIS, by = "person_id") %>%
    filter(!is.na(visit_last) & visit_last >= as.POSIXct(index_date) + days(7))
  note("outpatient visit >= index + 7 days", nrow(d))
  d <- d %>% mutate(time2 = ifelse(!is.na(dem_first), as.numeric(as.Date(dem_first) - index_date) / 365,
                                   as.numeric(as.Date(visit_last) - index_date) / 365)) %>%
    filter(time2 >= 0)
  note("time2 >= 0  (relax file, age 30-79)", nrow(d))
  # covariates: any record on or before the index date (midnight)
  d <- d %>% left_join(DRUG, by = "person_id") %>% left_join(OBS, by = "person_id")
  ix <- as.POSIXct(d$index_date)
  first_le <- function(x) as.integer(!is.na(x) & x <= ix)
  d$smoker       <- first_le(d$smoker_first)
  d$anti         <- first_le(d$anti_first)
  d$diabete      <- first_le(d$diabetes_first)
  d$hearing_loss <- first_le(d$hearing_first)
  d$inactivity   <- first_le(d$inactivity_first)
  d$isolation    <- as.integer(first_le(d$isolation_first) | first_le(d$isolation_obs_first))
  d$brain        <- first_le(d$tbi_first)
  d$depression   <- first_le(d$depression_first)
  win <- OBSR %>% inner_join(d %>% transmute(person_id, idx = as.POSIXct(index_date)), by = "person_id") %>%
    filter(observation_datetime <= idx, observation_datetime >= idx - years(3))
  alc <- win %>% filter(observation_concept_id %in% CS$alcohol$ids) %>% distinct(person_id)
  edu <- win %>% filter(observation_concept_id %in% CS$education$ids) %>% arrange(person_id, desc(observation_datetime)) %>%
    group_by(person_id) %>% slice_head(n = 1) %>% ungroup() %>%
    transmute(person_id, low_education = as.integer(value_as_concept_name %in% LOW_EDU))
  d <- d %>% mutate(alcohol = as.integer(person_id %in% alc$person_id)) %>%
    left_join(edu, by = "person_id") %>% mutate(low_education = coalesce(low_education, 0L))
  # recoding, as in the paper
  d$ethnicity[d$ethnicity %in% c("PMI: Skip", "What Race Ethnicity: Race Ethnicity None Of These",
                                 "PMI: Prefer Not To Answer", "No matching concept")] <- "No information"
  d$race[d$race %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander", "I prefer not to answer",
                       "None Indicated", "None of these", "PMI: Skip", "Middle Eastern or North African",
                       "More than one population")] <- "Others"
  d$gender[d$gender %in% c("Gender Identity: Additional Options", "Gender Identity: Non Binary", "I prefer not to answer",
                           "Gender Identity: Transgender", "Not man only, not woman only, prefer not to answer, or skipped",
                           "PMI: Skip", "No matching concept")] <- "Others"
  v30 <- d %>% filter(is.finite(time2), time2 > 0, is.finite(age), age >= 30, age < 80)
  v45 <- v30 %>% filter(age >= 45)
  an  <- v45 %>% filter(ldlc < 400, trig < 700)
  note("V30: time2 > 0, age 30-79", nrow(v30))
  note("V45: time2 > 0, age 45-79", nrow(v45))
  note("analysis set: V45 + LDL<400, TG<700", nrow(an))
  list(att = bind_rows(att), cohort = d, an = an,
       ev = c(relax = sum(d$Dementia), v30 = sum(v30$Dementia), an = sum(an$Dementia)))
}

## ---- 4. report (the paper's v8 numbers are shown for orientation only) --------------
PAPER <- list(prev = c(smoker = 7233, anti = 15778, diabete = 12841, hearing_loss = 7410, inactivity = 21, isolation = 0,
                       alcohol = 0, brain = 841, depression = 15997, low_education = 0),
              mean = c(sbp = 126.22, bmi = 28.32, ldlc = 105.35, hdl = 52.29, trig = 116.12, tc = 185.52, age = 59.05),
              sd = c(sbp = 15.03, bmi = 4.18, ldlc = 23.79, hdl = 11.79, trig = 53.69, tc = 25.33, age = 8.11),
              gender = c(Female = 59318, Male = 44184, Others = 1423),
              race = c(Asian = 1791, `Black or African American` = 15199, Others = 18021, White = 69914),
              ethnicity = c(`Hispanic or Latino` = 11385, `No information` = 3324, `Not Hispanic or Latino` = 90216),
              regimes = c(`[45,47)` = 6086, `[47,55)` = 30771, `[55,68)` = 51817, `[68,74)` = 12051, `[74,79]` = 4200),
              regime_events = c(`[45,47)` = 123, `[47,55)` = 674, `[55,68)` = 1552, `[68,74)` = 649, `[74,79]` = 344))
REF_N <- c("time2 >= 0  (relax file, age 30-79)" = 141693, "V30: time2 > 0, age 30-79" = 141617,
           "V45: time2 > 0, age 45-79" = 104925, "analysis set: V45 + LDL<400, TG<700" = 104925)
REF_EV <- c(3990, 3915, 3342)
EV_LAB <- c("relax file: dementia events", "V30: dementia events", "analysis set: dementia events")
EV_DEN <- c("time2 >= 0  (relax file, age 30-79)", "V30: time2 > 0, age 30-79", "analysis set: V45 + LDL<400, TG<700")
line <- function(label, shown, paper = NA, exact = NA, digits = 0, share = NULL) {
  ps <- if (is.na(paper)) "" else if (digits == 0) msk(paper) else formatC(paper, format = "f", digits = digits)
  flag <- ""
  if (!is.na(paper) && !is.na(exact) && exact_shown(shown)) {
    flag <- if (!is.null(share)) sprintf("%5.1f%% vs %5.1f%%", 100 * exact / share[1], 100 * paper / share[2])
            else if (isTRUE(round(exact, digits) == round(paper, digits))) "same"
            else if (paper != 0) sprintf("%+.1f%%", 100 * (exact - paper) / paper) else ""
  }
  cat(sprintf("    %-44s %15s %13s   %s\n", label, shown, ps, flag))
}
counts_header <- function(share = FALSE) cat(sprintf("    %-44s %15s %13s   %s\n", "", sprintf("%s (this run)", CDR_ID), "paper (v8)",
                                                     if (share) "share: this run vs paper" else ""))

R <- build()
hdr("3. Attrition and events")
a <- R$att; sn <- setNames(chain_mask(a$n), a$step); se <- chain_mask(R$ev)
counts_header()
for (i in seq_len(nrow(a))) line(a$step[i], sn[[i]], if (a$step[i] %in% names(REF_N)) REF_N[[a$step[i]]] else NA, a$n[i])
for (i in 1:3) line(paste(EV_LAB[i], "(rate)"), se[i], REF_EV[i], R$ev[i],
                    share = if (exact_shown(sn[[EV_DEN[i]]])) c(a$n[a$step == EV_DEN[i]], REF_N[[EV_DEN[i]]]) else NULL)

v <- R$an; N <- nrow(v); P <- 104925
sh_ok <- exact_shown(sn[["analysis set: V45 + LDL<400, TG<700"]])
SH <- if (sh_ok) c(N, P) else NULL
hdr("4. Analysis set (Table 1)")
counts_header(share = TRUE)
line("median follow-up, years (time2)", formatC(median(v$time2), format = "f", digits = 1), 9.6, median(v$time2), 1)
for (b in names(PAPER$mean)) {
  line(sprintf("%s mean", b), formatC(mean(v[[b]], na.rm = TRUE), format = "f", digits = 2), PAPER$mean[[b]], mean(v[[b]], na.rm = TRUE), 2)
  line(sprintf("%s SD", b), formatC(sd(v[[b]], na.rm = TRUE), format = "f", digits = 2), PAPER$sd[[b]], sd(v[[b]], na.rm = TRUE), 2)
}
for (k in names(PAPER$prev)) {
  x1 <- sum(v[[k]] == 1, na.rm = TRUE)
  line(sprintf("%s = 1 (n)", k), part_mask(c(x1, N - x1))[1], PAPER$prev[[k]], x1, share = SH)
}
for (g in c("gender", "race", "ethnicity")) {
  val <- ifelse(is.na(v[[g]]), "<missing>", as.character(v[[g]]))
  tb <- table(val); pm <- setNames(part_mask(as.numeric(tb)), names(tb))
  for (lv in names(PAPER[[g]])) line(sprintf("%s: %s", g, lv), if (lv %in% names(tb)) pm[[lv]] else "0",
                                    PAPER[[g]][[lv]], if (lv %in% names(tb)) as.numeric(tb[[lv]]) else 0, share = SH)
  for (lv in setdiff(names(tb), names(PAPER[[g]]))) line(sprintf("%s: %s  (not in paper)", g, substr(lv, 1, 22)), pm[[lv]])
}
reg <- cut(v$age, c(45, 47, 55, 68, 74, Inf), right = FALSE, labels = names(PAPER$regimes))
tn <- as.numeric(table(reg)); te <- as.numeric(tapply(v$Dementia, reg, sum)); te[is.na(te)] <- 0
pn <- part_mask(tn); pe <- part_mask(te)
for (i in seq_along(tn)) {
  line(sprintf("age %s: N", names(PAPER$regimes)[i]), pn[i], PAPER$regimes[[i]], tn[i], share = SH)
  line(sprintf("age %s: events (rate)", names(PAPER$regimes)[i]), pe[i], PAPER$regime_events[[i]], te[i],
       share = if (exact_shown(pn[i])) c(tn[i], PAPER$regimes[[i]]) else NULL)
}
ty <- table(factor(year(v$index_date), levels = 2010:2018))
cat(sprintf("    %-44s %s\n", "index year 2010 ... 2018 (N)", paste(part_mask(as.numeric(ty)), collapse = " | ")))
cat("\n    labs in the analysis set: imputed share, and 1st / 50th / 99th percentile after imputation\n")
for (l in LABS_R) {
  k <- sum(v[[paste0("imp_", l)]] == 1); s <- part_mask(c(k, N - k))[1]
  qq <- stats::quantile(v[[l]], c(0.01, 0.5, 0.99), na.rm = TRUE)
  cat(sprintf("      %-5s imputed %8s %s   p1 %7.1f   p50 %7.1f   p99 %7.1f\n", l, s,
              if (exact_shown(s)) sprintf("(%4.1f%%)", 100 * k / N) else "       ", qq[1], qq[2], qq[3]))
}

## ---- 5. save ------------------------------------------------------------------------
hdr("5. Save")
if (SAVE) {
  out_dir <- file.path(path.expand("~"), "dementia_v9"); dir.create(out_dir, showWarnings = FALSE)
  keep <- c("person_id", "Dementia", "index_date", "ldlc", "sbp", "tc", "hdl", "trig", "bmi", "age", "time2", "gender", "race",
            "ethnicity", "sex_at_birth", "date_of_birth", "smoker", "anti", "diabete", "hearing_loss", "inactivity", "isolation",
            "alcohol", "brain", "depression", "low_education", paste0("imp_", LABS_R), "dem_first", "visit_last")
  tag <- if (identical(IMPUTE, "mean")) "dementia_v9_meanimp" else "dementia_v9"
  f_relax <- file.path(out_dir, sprintf("%s_relax_%s.csv", tag, CDR_ID))
  f_an    <- file.path(out_dir, sprintf("%s_analysis_%s.csv", tag, CDR_ID))
  f_def   <- file.path(out_dir, sprintf("%s_definitions.txt", tag))
  utils::write.csv(R$cohort[, intersect(keep, names(R$cohort))], f_relax, row.names = FALSE)
  utils::write.csv(R$an[, intersect(keep, names(R$an))], f_an, row.names = FALSE)
  writeLines(c(
    paste0("Dementia change-point analysis cohort, built ", format(Sys.time(), "%Y-%m-%d %H:%M UTC", tz = "UTC"), " from ", CDR),
    paste0("relax file: ", basename(f_relax), " (time2 >= 0, ages 30-79)"),
    paste0("analysis set: ", basename(f_an), " (time2 > 0, ages 45-79, LDL < 400, TG < 700)"),
    sprintf("imp_<lab> = 1: lab missing in its window, filled with the race x sex-at-birth %s", IMPUTE),
    "Definitions: see the header of dementia_cp_02_v9_cohort.R"), f_def)
  if (length(BUCKETS)) {
    dest <- paste0(BUCKETS[1], "/dementia_v9/", CDR_ID, "/")
    for (f in c(f_relax, f_an, f_def)) {
      r <- sh("gsutil", c("-q", "cp", shQuote(f), shQuote(dest)))
      cat(sprintf("  %s -> %s%s\n", basename(f), dest, if (ok(r)) "" else "   COPY FAILED, file is on this app's disk only"))
    }
  } else {
    cat("  no bucket found; files are only on this app's disk:", out_dir, "\n")
  }
}
cat(sprintf("\n  CDR %s | runtime %.1f min\n", CDR, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
