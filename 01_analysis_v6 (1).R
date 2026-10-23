# =============================================================================
# JGSS-2026-346: MAIN ANALYSIS
# =============================================================================
#
# Reads the two compact files written by 00_wrangle_hrs.R and
# 00_wrangle_share.R, produces every main-text result, and saves the workspace
# 02_sensitivity_v6.R runs on.
#
# WHAT IS ESTIMATED
#
#  ADL and IADL onset are discrete-time survival models on person-wave data,
#  estimated as pooled logistic regressions with wave dummies for the baseline
#  hazard and standard errors clustered on the person. Both limitations are
#  treated as absorbing.
#
#  Verbal memory is a linear mixed model of the 0 to 20 word-recall score with
#  a random intercept, a random slope for time, and baseline age interacted
#  with time (Part 7b). Part 7b also fits the specification without that
#  interaction, because the comparison between the two is a supplementary
#  table.
#
#  A threshold definition of memory decline, a drop of one baseline standard
#  deviation, is estimated alongside ADL and IADL in Part 5 with the same
#  discrete-time machinery and is reported as a sensitivity analysis.
#
#  Two decompositions, for two different questions. Part 8 is a
#  Karlson-Holm-Breen decomposition with education as the exposure. Part 9 is
#  a g-computation of the United States minus Europe gap. Part 9 is not a KHB:
#  country fixed effects are impossible when country is the exposure, and
#  pooling the two surveys breaks the wave dummies, which are disjoint by
#  survey.
#
#  Horizons are fixed and common to both surveys. HRS observes 2, 4, 6, 8 and
#  10 years and SHARE 2, 4, 7 and 9, so reading each survey at its own end of
#  follow-up would compare the United States at 10 years against Europe at 9,
#  and memory at 8 years against 7. For memory that asymmetry alone reverses
#  the ranking.
#
#  The healthy profile for the Model 3 predictions is the three
#  doctor-diagnosed conditions. A wider profile matched 8.0% of the US sample
#  against 16.1% of the Nordic one, which put the cross-national comparison on
#  thin support in the region the claim is about. Self-rated health, mobility
#  limitation and ever-smoker are left out of the profile: the first two are
#  reported very differently across countries and the third cannot be set to
#  zero for a 65-year-old. All of them stay in the model.
#
#  Model 1, age and wave only, carries the total effect of education on the
#  already-restricted sample and is what the figures report. Model 3 at a
#  healthy profile answers the separate cross-national question.
#
# PARTS
#   0 to 3   settings, load, risk sets, person-wave construction
#   4        descriptives and incidence rates
#   5        discrete-time models, all three outcomes
#   6        whether the gradient varies with age
#   7        predicted cumulative incidence
#   7b       linear mixed models of word recall
#   8        KHB decomposition, education as the exposure
#   9        g-computation of the United States minus Europe gap
#   10       horizon sensitivity and the education crossover
#   11       saves main_v6.RData
#   12       Figure 1, Figure S1, Figure S2 and the two forest plots
#
# OUTPUTS
#   V4_01_table1.csv               Table 1, weighted
#   V4_02_tableS_selection.csv     sample construction
#   V4_03_incidence.csv            Table 2, weighted and unweighted
#   V4_04_models.csv               M1 / M1c / M2 / M3, all outcomes
#   V5_05b_all_coefficients.csv    every coefficient, M1 / M2 / M3
#   V4_05_models_floor.csv         the same with the cognitive floor
#   V4_06_age_interaction.csv      education x age
#   V5_07_cuminc_series.csv        cumulative incidence, wave by wave
#   V5_07_cuminc_horizon.csv       cumulative incidence at fixed horizons
#   S5_E2_memory_lmm_age.csv       mixed-model coefficients, both specifications
#   S5_E2_memory_predicted_age.csv predicted recall trajectories
#   S5_E2_memory_gaps.csv          low minus high recall gaps
#   V5_08_khb_education.csv        KHB, education as exposure
#   V5_09_gap_decomposition.csv    US-Europe gap, g-computation
#   V5_10_horizon_sensitivity.csv  does the ranking depend on the horizon
#   V5_10_crossover.csv            US high-educated vs European low-educated
#   main_v6.RData                  everything, for the sensitivity and figures
#
# Output file names carry the version in which each output first appeared and
# are left alone so that earlier result files still match. They do not track
# the current version number.
#
# RUNTIME: Parts 0 to 7 take roughly half an hour. Part 8 is the bootstrap and
# dominates; set N_BOOT low once to time a single cell before committing to a
# full run.
# =============================================================================
# =============================================================================
# PART 0: SETTINGS
# =============================================================================

library(tidyverse)
library(haven)
library(sandwich)
library(lmtest)
library(MASS)
library(lme4)
select <- dplyr::select

# ---- paths -----------------------------------------------------------------
# Relative to the working directory, which should be the repository root. The
# two inputs are what the wrangling scripts write. Copy config_example.R to
# config.R to point them somewhere else; config.R is untracked and is sourced
# here if it exists. Every table and figure goes to OUT_DIR.

DERIVED_DIR <- "data/derived"
OUT_DIR     <- "output"

HRS_FILE   <- file.path(DERIVED_DIR, "hrs_compact.rds")
SHARE_FILE <- file.path(DERIVED_DIR, "share_compact.rds")

if (file.exists("config.R")) source("config.R")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SD_MULT <- 1.0     # memory event and floor, in SD of the pooled baseline
REF_AGE <- 65
N_SIM   <- 1000
USE_COUNTRY_FE <- TRUE

HRS_WAVES   <- data.frame(wave = c(12, 13, 14, 15, 16), time = c(2, 4, 6, 8, 10))
SHARE_WAVES <- data.frame(wave = c(6, 7, 8, 9),         time = c(2, 4, 7, 9))

# memory stops in 2020 in both surveys: harmonized word recall is released for
# RAND HRS waves 3-15 and SHARE waves 4-9, and the submitted paper used the
# same window
MEM_HRS_WAVES   <- HRS_WAVES[HRS_WAVES$wave <= 15, ]
MEM_SHARE_WAVES <- SHARE_WAVES[SHARE_WAVES$wave <= 8, ]

MED_RAW <- c("poor_health", "partnered", "has_mobility",
             "has_hypertension", "has_diabetes", "obese", "ever_smoker")
MED_M2  <- c("poor_health", "partnered_obs", "has_mobility")
MED_M3  <- c(MED_M2, "has_hypertension", "has_diabetes", "obese", "ever_smoker")

# The healthy profile: the four objective markers only. poor_health,
# partnered and ever_smoker stay IN the model as mediators but are held at
# observed means rather than used to define "healthy". Self-rated health is
# reported very differently across countries, ever-smoker is a lifetime
# exposure that cannot be set to zero for a 65-year-old, and partnership is a
# social resource rather than a risk factor.
# Mobility limitation came out too: US 33.4% against 14.5% in the Nordics is
# the same instrumentation worry as the IADL index, and building the profile
# on a difficulty scale bakes that into the cross-national claim. What is left
# are three "has a doctor ever told you" items, which travel better.
HEALTHY_PROFILE <- c(has_hypertension = 0, has_diabetes = 0, obese = 0)

# The horizon each outcome is reported at: the longest both surveys reach.
# Memory caps at 6 because it runs to 7 years in SHARE and 8 in HRS.
HORIZON     <- list(ADL = 8, IADL = 8, Memory = 6)

HORIZON_DEFAULT <- list(ADL = 8, IADL = 8, Memory = 6)
# HORIZON is saved into main_*.RData, so a load() from an older session can
# replace the list with a bare value and every HORIZON[["Memory"]] below then
# fails. Repair it here rather than failing deep inside a model or figure
# block, and say so, since a silent repair would hide a stale workspace.
if (!is.list(HORIZON) ||
    !all(c("ADL", "IADL", "Memory") %in% names(HORIZON))) {
  warning("HORIZON was not the expected list; restoring the default. ",
          "Check for a load() of an older RData in this session.",
          call. = FALSE)
  HORIZON <- HORIZON_DEFAULT
}

# How far the memory trajectories are PREDICTED, as opposed to reported. HRS
# observes recall at 0, 2, 4, 6 and 8 years and SHARE at 0, 2, 4 and 7, so a
# value at 8 is inside the range for the HRS and one year beyond SHARE's last
# memory interview. Six years is what the main text reports, because it is the
# only horizon inside both surveys' observed range; eight is produced for the
# supplement and must be labelled as extrapolated for the European regions.
MEM_MAX <- 8

HORIZON_ALL <- list(ADL = c(2, 4, 6, 8), IADL = c(2, 4, 6, 8),
                    Memory = c(2, 4, 6))

N_BOOT        <- 500                  # cluster bootstrap draws for the KHB
KHB_OUTCOMES  <- c("ADL", "IADL")     # memory's total effect is ~0, so there
# is nothing to decompose; add "Memory"
# to produce the supplementary numbers

REGIONS_ALL <- c("US", "Nordic", "Western", "Mediterranean", "Eastern")

set.seed(20261009)


# =============================================================================
# PART 1: LOAD
# =============================================================================

cat("=== Loading ===\n")

strip_labels <- function(d) {
  d[] <- lapply(d, function(x) {
    if (inherits(x, c("haven_labelled", "labelled"))) x <- haven::zap_labels(x)
    x
  })
  d
}
make_id <- function(x) if (is.numeric(x)) sprintf("%.0f", x) else as.character(x)

stopifnot(file.exists(HRS_FILE), file.exists(SHARE_FILE))
hrs   <- strip_labels(readRDS(HRS_FILE))
share <- strip_labels(readRDS(SHARE_FILE))

# The files must be the surveys they are named after, because reading one
# survey's data as the other's is silent and ruinous.
stopifnot(identical(unique(hrs$study), "HRS"),
          identical(unique(share$study), "SHARE"))
cat("HRS   from", HRS_FILE, "\n  written", format(file.info(HRS_FILE)$mtime), "\n")
cat("SHARE from", SHARE_FILE, "\n  written", format(file.info(SHARE_FILE)$mtime), "\n")

prep <- function(dat, id_col, wt_col, survey_name, adl_cands, iadl_cands) {
  out <- dat
  out$id     <- make_id(out[[id_col]])
  out$weight <- as.numeric(out[[wt_col]])
  out$survey <- survey_name
  if (!"country" %in% names(out)) out$country <- survey_name
  out$country <- as.character(out$country)
  
  raw_ed <- out$education_3cat
  out$education_3cat <- factor(as.character(raw_ed),
                               levels = c("High", "Medium", "Low"))
  if (all(is.na(out$education_3cat)) && !all(is.na(raw_ed)))
    stop("education_3cat in ", survey_name, " is not High/Medium/Low")
  
  a <- adl_cands[adl_cands  %in% names(out)][1]
  i <- iadl_cands[iadl_cands %in% names(out)][1]
  if (is.na(a) || is.na(i)) stop("baseline ADL or IADL column not found in ", survey_name)
  out$adl_baseline  <- out[[a]]
  out$iadl_baseline <- out[[i]]
  out
}

# The first candidate is what the wrangle scripts write, and it is the one
# that should always be used: adl_baseline is the verified five-item index and
# iadl_baseline is the five-item index built by hand from its components. The
# second candidate in each pair is the survey's own pre-computed index, kept
# only so this still runs against an older compact file. Do NOT rely on it:
# r5iadla covers three items against r11iadl5a's five, so analysing those two
# together compares instruments rather than health. build_person_wave carries
# a scale guard that catches it.
hrs_b   <- prep(hrs,   "hhidpn",  "wtresp_w11", "HRS",
                c("adl_baseline", "r11adl5a"), c("iadl_baseline", "r11iadl5a"))
share_b <- prep(share, "mergeid", "wtresp_w5",  "SHARE",
                c("adl_baseline", "r5adla"),   c("iadl_baseline", "r5iadla"))

cat("HRS", nrow(hrs_b), "rows; SHARE", nrow(share_b), "rows\n")


# =============================================================================
# PART 2: RISK SETS
# =============================================================================

cat("\n=== Risk sets ===\n")

PERSON_COLS <- c("id", "survey", "country", "weight", "education_3cat",
                 "age_baseline", "female", "region", "tr20_baseline",
                 "adl_baseline", "iadl_baseline", "death_year", MED_RAW)

pick <- function(d) {
  have    <- PERSON_COLS[PERSON_COLS %in% names(d)]
  missing <- setdiff(PERSON_COLS, have)
  if (length(missing))
    cat("  NOTE, missing in", unique(d$survey)[1], ":",
        paste(missing, collapse = ", "), "\n")
  d[, have, drop = FALSE]
}

person_all <- bind_rows(pick(hrs_b), pick(share_b))

# free of ADL and IADL at baseline, with education
phys_base <- person_all %>%
  filter(!is.na(adl_baseline), !is.na(iadl_baseline),
         adl_baseline == 0, iadl_baseline == 0,
         !is.na(education_3cat))

# partnership: a missing indicator, which should now be near zero after the
# SHARE recode fix. build_rhs drops it where it does not vary.
phys_base$partnered_miss <- as.integer(is.na(phys_base$partnered))
phys_base$partnered_obs  <- ifelse(is.na(phys_base$partnered), 0L,
                                   as.integer(phys_base$partnered))

cat("Partnership missing, by region (%):\n")
phys_base %>% group_by(region) %>%
  summarise(pct = round(100 * mean(partnered_miss), 2), .groups = "drop") %>%
  print()

# the memory floor, computed BEFORE it is applied so it is not circular
BL_SD   <- sd(phys_base$tr20_baseline, na.rm = TRUE)
BL_MEAN <- mean(phys_base$tr20_baseline, na.rm = TRUE)
MEM_CUT <- SD_MULT * BL_SD

cat(sprintf("\nPooled baseline recall: mean %.2f, SD %.2f\n", BL_MEAN, BL_SD))
cat(sprintf("Memory floor and decline threshold: %.2f points (%.1f SD)\n",
            MEM_CUT, SD_MULT))

mem_base <- phys_base %>%
  filter(!is.na(tr20_baseline), tr20_baseline >= MEM_CUT)

cat("\nfree of ADL and IADL, with education :", nrow(phys_base), "\n")
cat("also above the memory floor          :", nrow(mem_base), "\n")

cat("\nWho the memory floor removes (%), by region and education:\n")
phys_base %>%
  mutate(removed = !(id %in% mem_base$id)) %>%
  group_by(region, education_3cat) %>%
  summarise(pct = round(100 * mean(removed), 1), .groups = "drop") %>%
  pivot_wider(names_from = education_3cat, values_from = pct) %>%
  print()


# =============================================================================
# PART 3: PERSON-WAVE DATA
# =============================================================================

cat("\n=== Person-wave data ===\n")

get_tr20_col <- function(dat, wave) {
  cand <- c(paste0("tr20_w", wave), paste0("tr20_w", wave, "_p"),
            paste0("tr20_w", wave, "_w"))
  hit <- cand[cand %in% names(dat)]
  if (!length(hit)) NA_character_ else hit[1]
}

get_wave_value <- function(wide, wave, kind) {
  if (kind == "mem") {
    cl <- get_tr20_col(wide, wave)
    return(if (is.na(cl)) NULL else wide[[cl]])
  }
  cand <- c(paste0(kind, "_w", wave), paste0("r", wave, kind, "5a"),
            paste0("r", wave, kind, "a"))
  hit <- cand[cand %in% names(wide)]
  if (!length(hit)) NULL else wide[[hit[1]]]
}

build_person_wave <- function(wide, wave_tab, kind) {
  blocks <- list()
  for (k in seq_len(nrow(wave_tab))) {
    val <- get_wave_value(wide, wave_tab$wave[k], kind)
    if (is.null(val)) { cat("    no", kind, "for wave", wave_tab$wave[k], "\n"); next }
    
    # SCALE GUARD. ADL and IADL are both five-item counts analysed as >= 1,
    # so every column reaching the models must top out at 5. This exists
    # because the Gateway SHARE IADL index carries three items (money, meds,
    # phone) against RAND's five, which made every between-region IADL
    # comparison a comparison of instruments rather than of health. The
    # wrangle scripts now build the index by hand from the five named
    # components; this catches it if they ever drift apart again. See
    # check_harmonization.R for the full version.
    if (kind %in% c("adl", "iadl")) {
      mx <- suppressWarnings(max(val, na.rm = TRUE))
      if (is.finite(mx)) {
        if (mx > 5)
          stop(sprintf("%s wave %d has maximum %s: not a five-item count",
                       kind, wave_tab$wave[k], mx))
        if (mx < 5)
          warning(sprintf(paste("%s wave %d tops out at %s, not 5. Check the",
                                "item set before using any between-region",
                                "comparison of this outcome."),
                          kind, wave_tab$wave[k], mx), call. = FALSE)
      }
    }
    
    blocks[[length(blocks) + 1]] <- data.frame(
      id = wide$id, wave = wave_tab$wave[k], time = wave_tab$time[k],
      value = val, stringsAsFactors = FALSE)
  }
  bind_rows(blocks)
}

make_event <- function(long, kind, baseline_lookup) {
  if (kind == "mem") {
    bl <- baseline_lookup[long$id]
    long$event <- ifelse(is.na(long$value) | is.na(bl), NA_integer_,
                         ifelse(long$value <= bl - MEM_CUT, 1L, 0L))
  } else {
    long$event <- ifelse(is.na(long$value), NA_integer_,
                         ifelse(long$value >= 1, 1L, 0L))
  }
  long
}

# keep each person up to and including their first event; a wave with a
# missing outcome is skipped but does not end follow-up
trim_to_first_event <- function(long) {
  long %>%
    filter(!is.na(event)) %>%
    arrange(id, time) %>%
    group_by(id) %>%
    mutate(prior = cumsum(dplyr::lag(event, default = 0))) %>%
    filter(prior == 0) %>% ungroup() %>% select(-prior)
}

build_long <- function(kind, base) {
  wt_h <- if (kind == "mem") MEM_HRS_WAVES   else HRS_WAVES
  wt_s <- if (kind == "mem") MEM_SHARE_WAVES else SHARE_WAVES
  bl   <- setNames(base$tr20_baseline, base$id)
  
  out <- bind_rows(build_person_wave(hrs_b,   wt_h, kind),
                   build_person_wave(share_b, wt_s, kind)) %>%
    make_event(kind, bl) %>%
    trim_to_first_event() %>%
    inner_join(base, by = "id") %>%
    mutate(age = age_baseline + time)
  
  cat("  ", kind, ": ", nrow(out), " person-waves, ",
      length(unique(out$id)), " people, ", sum(out$event), " events\n", sep = "")
  out
}

LONG <- list(ADL    = build_long("adl",  phys_base),
             IADL   = build_long("iadl", phys_base),
             Memory = build_long("mem",  mem_base))

# the same outcomes on the memory-restricted sample, for the Table S check
LONG_FLOOR <- list(ADL  = build_long("adl",  mem_base),
                   IADL = build_long("iadl", mem_base))

person_years <- function(long) {
  long %>% arrange(id, time) %>% group_by(id) %>%
    mutate(gap = time - dplyr::lag(time, default = 0)) %>% ungroup()
}

cat("\nPerson-years and events, for the Table 2 footnote:\n")
for (nm in names(LONG)) {
  py <- person_years(LONG[[nm]]) %>% summarise(s = sum(gap)) %>% pull(s)
  cat("  ", nm, ": ", format(py, big.mark = ","), " person-years, ",
      format(sum(LONG[[nm]]$event), big.mark = ","), " events, ",
      format(length(unique(LONG[[nm]]$id)), big.mark = ","), " people\n", sep = "")
}


# =============================================================================
# PART 4: DESCRIPTIVES
# =============================================================================

cat("\n=== Table 1 ===\n")

analytic <- phys_base %>% filter(id %in% LONG[["ADL"]]$id)

wmean <- function(x, w) weighted.mean(x, w, na.rm = TRUE)
wpct  <- function(x, w) 100 * weighted.mean(x, w, na.rm = TRUE)
wsd   <- function(x, w) {
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  sqrt(sum(w * (x - weighted.mean(x, w))^2) / sum(w))
}

# Kish's effective sample size. Unequal weights inflate the variance of a
# weighted estimate, so the usual sqrt(p(1-p)/n) understates it. n_eff is the
# size of an equally weighted sample carrying the same information, and every
# interval below uses it in place of n.
n_eff <- function(w) sum(w)^2 / sum(w^2)

# A weighted percentage with a 95% interval on the effective sample size.
wpct_ci <- function(x, w) {
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  p  <- weighted.mean(x, w); ne <- n_eff(w)
  se <- sqrt(max(p * (1 - p), 0) / ne)
  sprintf("%.1f (%.1f\u2013%.1f)", 100 * p,
          100 * max(0, p - 1.96 * se), 100 * min(1, p + 1.96 * se))
}

# A weighted mean with its SD and a 95% interval, again on n_eff.
wmean_ci <- function(x, w) {
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  m  <- weighted.mean(x, w); s <- wsd(x, w); ne <- n_eff(w)
  sprintf("%.1f (%.1f) [%.1f\u2013%.1f]", m, s,
          m - 1.96 * s / sqrt(ne), m + 1.96 * s / sqrt(ne))
}

t1_col <- function(d) {
  dw <- d %>% filter(!is.na(weight), weight > 0); w <- dw$weight
  c("N (unweighted)"            = sprintf("%d", nrow(d)),
    "N in weighted estimates"   = sprintf("%d", nrow(dw)),
    "Effective N (Kish)"        = sprintf("%.0f", n_eff(w)),
    "Age, mean (SD) [95% CI]"   = wmean_ci(dw$age_baseline, w),
    "Female, % (95% CI)"              = wpct_ci(dw$female, w),
    "Education high, % (95% CI)"      = wpct_ci(dw$education_3cat == "High", w),
    "Education medium, % (95% CI)"    = wpct_ci(dw$education_3cat == "Medium", w),
    "Education low, % (95% CI)"       = wpct_ci(dw$education_3cat == "Low", w),
    "Poor or fair health, % (95% CI)" = wpct_ci(dw$poor_health, w),
    "Partnered, % (95% CI)"           = wpct_ci(dw$partnered_obs, w),
    "Mobility limitation, % (95% CI)" = wpct_ci(dw$has_mobility, w),
    "Hypertension, % (95% CI)"        = wpct_ci(dw$has_hypertension, w),
    "Diabetes, % (95% CI)"            = wpct_ci(dw$has_diabetes, w),
    "Obese, % (95% CI)"               = wpct_ci(dw$obese, w),
    "Ever smoker, % (95% CI)"         = wpct_ci(dw$ever_smoker, w),
    "Word recall, mean (SD) [95% CI]" = wmean_ci(dw$tr20_baseline, w))
}

tab1 <- sapply(REGIONS_ALL, function(r) t1_col(analytic %>% filter(region == r)))
tab1 <- cbind(tab1, Total = t1_col(analytic))
tab1 <- data.frame(Characteristic = rownames(tab1), tab1,
                   row.names = NULL, check.names = FALSE)
print(tab1, right = FALSE)
write.csv(tab1, file.path(OUT_DIR, "V4_01_table1.csv"), row.names = FALSE)

# does the submitted sentence about the US still hold?
cat("\nIs the US lowest on low education and highest on the rest?\n")
chk <- analytic %>% filter(!is.na(weight), weight > 0) %>%
  group_by(region) %>%
  summarise(low_edu  = wpct(education_3cat == "Low", weight),
            hyperten = wpct(has_hypertension, weight),
            diabetes = wpct(has_diabetes, weight),
            obese    = wpct(obese, weight),
            mobility = wpct(has_mobility, weight), .groups = "drop")
for (v in c("low_edu", "hyperten", "diabetes", "obese", "mobility")) {
  us <- chk[[v]][chk$region == "US"]; oth <- chk[[v]][chk$region != "US"]
  ok <- if (v == "low_edu") us < min(oth) else us > max(oth)
  cat(sprintf("  %-9s US %5.1f   others %5.1f to %5.1f   %s\n",
              v, us, min(oth), max(oth), ifelse(ok, "holds", "DOES NOT HOLD")))
}

# sample construction, Table S
sel_col <- function(d) c(
  "N"                      = sprintf("%d", nrow(d)),
  "Age, mean (SD)"         = sprintf("%.1f (%.1f)", mean(d$age_baseline, na.rm = TRUE),
                                     sd(d$age_baseline, na.rm = TRUE)),
  "Female, %"              = sprintf("%.1f", 100 * mean(d$female, na.rm = TRUE)),
  "Education low, %"       = sprintf("%.1f", 100 * mean(d$education_3cat == "Low", na.rm = TRUE)),
  "Poor or fair health, %" = sprintf("%.1f", 100 * mean(d$poor_health, na.rm = TRUE)),
  "Hypertension, %"        = sprintf("%.1f", 100 * mean(d$has_hypertension, na.rm = TRUE)),
  "Word recall, mean (SD)" = sprintf("%.1f (%.1f)", mean(d$tr20_baseline, na.rm = TRUE),
                                     sd(d$tr20_baseline, na.rm = TRUE)))

tabS <- cbind(
  `1. Age-eligible, free of ADL and IADL limitation at baseline` = sel_col(person_all),
  `2. ... with education observed`                               = sel_col(phys_base),
  `3. ... with at least one follow-up observation (functional analytic sample)` = sel_col(analytic),
  `4. Of (2), baseline word recall above the memory floor`       = sel_col(mem_base),
  `5. ... with at least one follow-up observation (memory analytic sample)`     =
    sel_col(mem_base %>% filter(id %in% LONG[["Memory"]]$id)),
  `Excluded at step 3: no follow-up observed`                    =
    sel_col(phys_base %>% filter(!(id %in% analytic$id))))
tabS <- data.frame(Characteristic = rownames(tabS), tabS,
                   row.names = NULL, check.names = FALSE)
cat("\n=== Sample construction ===\n")
print(tabS, right = FALSE)
write.csv(tabS, file.path(OUT_DIR, "V4_02_tableS_selection.csv"), row.names = FALSE)


# ---- Table 2, incidence ----------------------------------------------------
cat("\n=== Table 2, incidence ===\n")

pois_ci <- function(ev, py) {
  if (ev == 0) return(c(0, 0, 1000 * qgamma(.975, 1) / py))
  c(1000 * ev / py, 1000 * qgamma(.025, ev) / py, 1000 * qgamma(.975, ev + 1) / py)
}

incidence_table <- function(long, label, weighted) {
  d <- person_years(long)
  if (weighted) d <- d %>% filter(!is.na(weight), weight > 0)
  d$.w <- if (weighted) d$weight else 1
  
  d %>%
    group_by(region, female, education_3cat) %>%
    summarise(n_people = n_distinct(id),
              events   = sum(event * .w),
              pyears   = sum(gap * .w),
              events_n = sum(event),
              pyears_n = sum(gap), .groups = "drop") %>%
    rowwise() %>%
    mutate(rate1000 = 1000 * events / pyears,
           ci_low   = if (weighted) NA_real_ else pois_ci(events_n, pyears_n)[2],
           ci_high  = if (weighted) NA_real_ else pois_ci(events_n, pyears_n)[3]) %>%
    ungroup() %>%
    mutate(outcome = label, sex = ifelse(female == 1, "Female", "Male"),
           weighted = weighted) %>%
    select(outcome, region, sex, education_3cat, n_people, events, pyears,
           rate1000, ci_low, ci_high, weighted)
}

incidence <- bind_rows(
  lapply(names(LONG), function(nm) incidence_table(LONG[[nm]], nm, FALSE)),
  lapply(names(LONG), function(nm) incidence_table(LONG[[nm]], nm, TRUE)))

cat("\nUnweighted rates with exact Poisson intervals (Table 2):\n")
incidence %>% filter(!weighted) %>%
  mutate(val = sprintf("%.1f (%.1f-%.1f)", rate1000, ci_low, ci_high)) %>%
  select(outcome, region, sex, education_3cat, val) %>%
  pivot_wider(names_from = education_3cat, values_from = val) %>%
  arrange(outcome, sex, match(region, REGIONS_ALL)) %>%
  print(n = 40, width = 200)

write.csv(incidence, file.path(OUT_DIR, "V4_03_incidence.csv"), row.names = FALSE)


# =============================================================================
# PART 5: MODELS
# =============================================================================

cat("\n=== Models ===\n")

build_rhs <- function(model_type, region_name, dat = NULL) {
  terms <- c("education_3cat", "age", "factor(wave)")
  if (USE_COUNTRY_FE && region_name != "US") terms <- c(terms, "factor(country)")
  if (is.null(dat) || length(unique(na.omit(dat$partnered_miss))) > 1)
    terms <- c(terms, "partnered_miss")
  if (model_type %in% c("M2", "M3")) terms <- c(terms, MED_M2)
  if (model_type == "M3")            terms <- c(terms, setdiff(MED_M3, MED_M2))
  paste(terms, collapse = " + ")
}

# M1c is M1 refitted on the people with complete M3 data, so that the M1-to-M3
# comparison is adjustment rather than a change of sample
fit_dt <- function(long, region_name, sex_val, outcome_label, model_type,
                   weighted = FALSE) {
  
  sex_lab <- ifelse(sex_val == 1, "Female", "Male")
  d <- long %>% filter(region == region_name, female == sex_val)
  if (!nrow(d)) return(NULL)
  
  spec <- if (model_type == "M1c") "M1" else model_type
  cc   <- if (model_type == "M1c") "M3" else spec
  
  # RESTRICT FIRST, BUILD THE FITTING FORMULA SECOND.
  #
  # build_rhs() includes partnered_miss only when it varies, and if it is
  # asked that question on the unfiltered cell it can answer yes and then
  # lose its only non-zero rows to the complete-case restriction that the
  # mediators impose. The term is then constant, glm returns NA for its
  # coefficient, and the anyNA() guard below discards the whole cell.
  #
  # That is what removed Mediterranean men from the verbal memory column of
  # Table 3: Models 1 and 2 estimated and Models 1c and 3 did not, because
  # only those two take complete cases on the four Model 3 mediators. In that
  # cell partnered_miss has a coefficient of 0.88 with a standard error of
  # 0.94, which is a handful of person-waves, and they all sat in the rows
  # the mediators dropped.
  f_cc <- as.formula(paste("event ~", build_rhs(cc, region_name, d)))
  d <- d[complete.cases(model.frame(f_cc, data = d, na.action = na.pass)), ]
  if (weighted) d <- d %>% filter(!is.na(weight), weight > 0)
  if (!nrow(d) || sum(d$event) < 10) return(NULL)
  f <- as.formula(paste("event ~", build_rhs(spec, region_name, d)))
  
  # Survey weights go in normalised to mean 1. In a binomial GLM `weights` is
  # a PRIOR weight, read as a number of trials, and the raw HRS weights run to
  # tens of thousands: summed over a cell that is a prior weight near 1e9, the
  # fitted probabilities are driven to 0 and 1, and the coefficients diverge
  # to +/- several hundred, which exp() turns into 0 or Inf. Scaling all the
  # weights by a constant multiplies the log-likelihood by that constant and
  # leaves the maximum exactly where it is, so the estimates are unchanged.
  # The standard errors come from vcovCL, so the scale does not affect them.
  m <- if (weighted) {
    d$.w <- d$weight / mean(d$weight)
    suppressWarnings(glm(f, data = d, family = binomial, weights = .w))
  } else glm(f, data = d, family = binomial)
  
  # A diverged or rank-deficient fit must not be reported as an estimate. The
  # warning names the offending term, so that a cell dropped here is visible
  # rather than silently absent from the output.
  bad <- names(coef(m))[is.na(coef(m))]
  if (!m$converged || length(bad) ||
      max(abs(coef(m)), na.rm = TRUE) > 30) {
    warning("dropped a non-finite fit: ", outcome_label, " ", region_name, " ",
            sex_lab, " ", model_type, if (weighted) " weighted" else "",
            " -- ", if (length(bad)) paste("aliased:", paste(bad, collapse = ", "))
            else if (!m$converged) "did not converge"
            else "coefficients diverged", call. = FALSE)
    return(NULL)
  }
  
  V <- vcovCL(m, cluster = d$id)
  est <- coeftest(m, vcov = V); ci <- coefci(m, vcov = V)
  keep <- grep("education", rownames(est))
  if (!length(keep)) return(NULL)
  
  data.frame(outcome = outcome_label, region = region_name, sex = sex_lab,
             model = model_type, weighted = weighted,
             term = rownames(est)[keep],
             OR = exp(est[keep, "Estimate"]),
             CI_low = exp(ci[keep, 1]), CI_high = exp(ci[keep, 2]),
             p = est[keep, "Pr(>|z|)"],
             n_rows = nrow(d), n_people = length(unique(d$id)),
             n_events = sum(d$event), row.names = NULL)
}

run_models <- function(L, weighted = FALSE) {
  out <- data.frame()
  for (nm in names(L)) {
    for (reg in REGIONS_ALL) for (sx in c(0, 1))
      for (mt in c("M1", "M1c", "M2", "M3")) {
        r <- fit_dt(L[[nm]], reg, sx, nm, mt, weighted)
        if (!is.null(r)) out <- bind_rows(out, r)
      }
    cat("  ", nm, "done\n")
  }
  as_tibble(out) %>%
    mutate(stars = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*",
                             TRUE ~ ""),
           OR_CI = sprintf("%.2f (%.2f-%.2f)%s", OR, CI_low, CI_high, stars))
}

models <- run_models(LONG)

cat("\n--- Low vs High, fully adjusted (M3) ---\n")
models %>% filter(term == "education_3catLow", model == "M3") %>%
  select(outcome, region, sex, OR_CI, n_people, n_events) %>%
  arrange(outcome, sex, match(region, REGIONS_ALL)) %>% print(n = 40)

cat("\n--- All specifications ---\n")
models %>% filter(term == "education_3catLow") %>%
  mutate(model = factor(model, levels = c("M1", "M1c", "M2", "M3"))) %>%
  select(outcome, region, sex, model, OR_CI) %>%
  pivot_wider(names_from = model, values_from = OR_CI) %>%
  arrange(outcome, sex, match(region, REGIONS_ALL)) %>%
  print(n = 40, width = 200)

write.csv(models, file.path(OUT_DIR, "V4_04_models.csv"), row.names = FALSE)

# weighted models, for the supplement
cat("\n--- Weighted with baseline survey weights (Table S) ---\n")
models_w <- run_models(LONG, weighted = TRUE)
models_w %>% filter(term == "education_3catLow", model == "M3") %>%
  select(outcome, region, sex, OR_CI) %>%
  arrange(outcome, sex, match(region, REGIONS_ALL)) %>% print(n = 40)

# the cognitive floor applied to the functional outcomes, for the supplement
cat("\n--- ADL and IADL with the cognitive floor applied (Table S) ---\n")
models_floor <- run_models(LONG_FLOOR)
models %>% filter(term == "education_3catLow", model == "M3",
                  outcome %in% c("ADL", "IADL")) %>%
  select(outcome, region, sex, main = OR_CI) %>%
  left_join(models_floor %>%
              filter(term == "education_3catLow", model == "M3") %>%
              select(outcome, region, sex, with_floor = OR_CI),
            by = c("outcome", "region", "sex")) %>%
  arrange(outcome, sex, match(region, REGIONS_ALL)) %>% print(n = 30)

write.csv(bind_rows(models_w %>% mutate(set = "weighted"),
                    models_floor %>% mutate(set = "cognitive floor")),
          file.path(OUT_DIR, "V4_05_models_floor.csv"), row.names = FALSE)


# ---- every coefficient, not only education, for the supplement -------------
# The supplementary tables in the first submission reported the full model, so
# the same is produced here: one row per term per cell, for M1 and M3.

cat("\n=== Full coefficient tables ===\n")

all_terms <- function(long, region_name, sex_val, outcome_label, model_type) {
  d <- long %>% filter(region == region_name, female == sex_val)
  if (!nrow(d)) return(NULL)
  # restrict first, then build the formula; see the note in fit_dt
  f <- as.formula(paste("event ~", build_rhs(model_type, region_name, d)))
  d <- d[complete.cases(model.frame(f, data = d, na.action = na.pass)), ]
  if (!nrow(d) || sum(d$event) < 10) return(NULL)
  f <- as.formula(paste("event ~", build_rhs(model_type, region_name, d)))
  m <- glm(f, data = d, family = binomial)
  if (!m$converged || anyNA(coef(m))) return(NULL)
  V <- vcovCL(m, cluster = d$id)
  est <- coeftest(m, vcov = V); cf <- coefci(m, vcov = V)
  data.frame(outcome = outcome_label, region = region_name,
             sex = ifelse(sex_val == 1, "Female", "Male"),
             model = model_type, term = rownames(est),
             beta = est[, 1], se = est[, 2],
             OR = exp(est[, 1]), CI_low = exp(cf[, 1]), CI_high = exp(cf[, 2]),
             p = est[, 4], n_rows = nrow(d), n_people = length(unique(d$id)),
             n_events = sum(d$event), row.names = NULL)
}

coefs <- data.frame()
for (nm in names(LONG)) {
  for (reg in REGIONS_ALL) for (sx in c(0, 1)) for (mt in c("M1", "M2", "M3")) {
    r <- all_terms(LONG[[nm]], reg, sx, nm, mt)
    if (!is.null(r)) coefs <- bind_rows(coefs, r)
  }
  cat("  ", nm, "done\n")
}
coefs <- as_tibble(coefs) %>%
  mutate(stars = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*",
                           TRUE ~ ""),
         OR_CI = sprintf("%.2f (%.2f-%.2f)%s", OR, CI_low, CI_high, stars))

cat("  ", nrow(coefs), "coefficient rows across",
    nrow(distinct(coefs, outcome, region, sex, model)), "fitted models\n")
write.csv(coefs, file.path(OUT_DIR, "V5_05b_all_coefficients.csv"),
          row.names = FALSE)


# =============================================================================
# PART 6: DOES THE GRADIENT VARY WITH AGE?
# =============================================================================
# One coefficient per contrast per cell, rather than predictions far from the
# sample mean under a linear age term.

cat("\n=== Education x age interaction ===\n")

age_int <- data.frame()
for (nm in names(LONG)) {
  for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
    d <- LONG[[nm]] %>% filter(region == reg, female == sx)
    if (!nrow(d)) next
    # age_c replaces age in the RHS so the main effect and the interaction use
    # the same column; otherwise age and I(age - 65) are collinear, R expands
    # the interaction over all three education levels and aliases one out.
    d$age_c <- d$age - REF_AGE
    # restrict first, then rebuild the formula; see the note in fit_dt
    rhs <- sub("\\bage\\b", "age_c", build_rhs("M3", reg, d))
    f   <- as.formula(paste("event ~", rhs, "+ education_3cat:age_c"))
    d <- d[complete.cases(model.frame(f, data = d, na.action = na.pass)), ]
    if (sum(d$event) < 20) next
    rhs <- sub("\\bage\\b", "age_c", build_rhs("M3", reg, d))
    f   <- as.formula(paste("event ~", rhs, "+ education_3cat:age_c"))
    m <- try(glm(f, data = d, family = binomial), silent = TRUE)
    if (inherits(m, "try-error")) next
    
    # An aliased coefficient means some term is perfectly collinear in this
    # cell, which after taking complete cases is usually a country dummy or a
    # mediator with no variation left. This loop needs only the
    # education-by-age terms, so report what aliased, drop the cell when those
    # are the aliased ones, and otherwise carry on. Anything printed here
    # belongs in the supplement note for the age-interaction table.
    bad <- names(coef(m))[is.na(coef(m))]
    if (length(bad)) {
      cat(sprintf("    aliased in %-7s %-14s %-6s : %s\n", nm, reg,
                  ifelse(sx == 1, "Female", "Male"),
                  paste(bad, collapse = ", ")))
      if (any(grepl("education.*age_c", bad))) next
    }
    
    V   <- try(vcovCL(m, cluster = d$id), silent = TRUE)
    if (inherits(V, "try-error")) V <- vcov(m)
    est <- try(coeftest(m, vcov = V), silent = TRUE)
    if (inherits(est, "try-error")) {
      cat(sprintf("    coeftest failed in %s %s %s, cell skipped\n",
                  nm, reg, ifelse(sx == 1, "Female", "Male")))
      next
    }
    k <- grep("education.*age_c", rownames(est))
    if (!length(k)) next
    age_int <- bind_rows(age_int, data.frame(
      outcome = nm, region = reg, sex = ifelse(sx == 1, "Female", "Male"),
      term = rownames(est)[k], estimate = est[k, "Estimate"],
      se = est[k, "Std. Error"], p = est[k, "Pr(>|z|)"], row.names = NULL))
  }
  cat("  ", nm, "done\n")
}

age_int <- as_tibble(age_int)
cat("\nLow vs High education, change in log-odds per year of age:\n")
age_int %>% filter(grepl("Low", term)) %>%
  mutate(val = sprintf("%+.4f (p = %.3f)", estimate, p)) %>%
  select(outcome, region, sex, val) %>%
  pivot_wider(names_from = sex, values_from = val) %>%
  arrange(outcome, match(region, REGIONS_ALL)) %>% print(n = 20)

cat("\nSignificant at p < .05:", sum(age_int$p[grepl("Low", age_int$term)] < .05),
    "of", sum(grepl("Low", age_int$term)), "\n")

write.csv(age_int, file.path(OUT_DIR, "V4_06_age_interaction.csv"), row.names = FALSE)



# =============================================================================
# PART 7: CUMULATIVE INCIDENCE
# =============================================================================
# Three things this part settles.
#
#  M1 as well as M3. M1 is age and wave only and carries the TOTAL effect of
#  education on the already-restricted sample, so M1 is what the figures
#  report and M3 at a healthy profile answers the separate cross-national
#  question.
#
#  A narrowed healthy profile. The seven-variable profile matched only 8.0% of
#  the US sample against 16.1% in the Nordics, so the comparison sat on thin
#  support in exactly the region the claim is about. HEALTHY_PROFILE is now
#  the three doctor-diagnosed conditions; every other mediator stays IN the
#  model and is held at its observed mean rather than at a healthy value.
#
#  Fixed common horizons. HRS observes 2,4,6,8,10 years and SHARE 2,4,7,9, so
#  "end of follow-up" compared the US at 10 years with Europe at 9, and memory
#  at 8 against 7. For memory that asymmetry alone reversed the ranking. Each
#  interval hazard is converted to a per-year rate and cumulated to a horizon
#  both surveys reach. This assumes a constant hazard within each observed
#  interval, which a discrete-time model with unequal spacing already assumes.
#
#  Country fixed effects stay in the regional models, so the odds ratios are
#  within-country contrasts. Each country is predicted separately and the
#  curves averaged using summed baseline survey weights, which are calibrated
#  to national populations.

cat("\n=== Cumulative incidence ===\n")

# cum is one draw's cumulative incidence at `times`; returns its value at t,
# interpolating inside the interval that contains t at a constant hazard.
cum_at <- function(cum, times, t) {
  S <- 1 - cum
  prev_t <- 0; prev_S <- 1; out <- 1
  for (i in seq_along(times)) {
    hj <- 1 - S[i] / prev_S
    if (t >= times[i]) {
      out <- out * (1 - hj)
      prev_t <- times[i]; prev_S <- S[i]
      if (abs(t - times[i]) < 1e-9) return(1 - out)
    } else {
      out <- out * (1 - hj) ^ ((t - prev_t) / (times[i] - prev_t))
      return(1 - out)
    }
  }
  NA_real_                       # t lies beyond observed follow-up here
}

# the country-averaged cumulative incidence, one row per simulation draw
cum_draws <- function(m, draws, d, cw, waves, ed, model_type, profile,
                      mediator_means = NULL) {
  cum_mix <- matrix(0, nrow = nrow(draws), ncol = length(waves))
  for (k in seq_len(nrow(cw))) {
    nd <- d[rep(1, length(waves)), ]
    nd$wave <- waves
    nd$age  <- REF_AGE
    nd$education_3cat <- factor(ed, levels = levels(d$education_3cat))
    nd$country <- cw$country[k]
    
    # a control, not a mediator, so it is set whatever the model
    if ("partnered_miss" %in% names(nd)) nd$partnered_miss <- 0
    
    if (model_type != "M1") {
      # every mediator starts at its observed mean, in this cell unless an
      # external set of means is supplied (used by the Part 9 counterfactual)
      mm <- if (is.null(mediator_means))
        sapply(MED_M3, function(v) mean(d[[v]], na.rm = TRUE)) else mediator_means
      for (v in MED_M3) if (v %in% names(nd)) nd[[v]] <- mm[[v]]
      # the healthy profile then overrides the four objective markers
      if (profile == "healthy")
        for (v in names(HEALTHY_PROFILE))
          if (v %in% names(nd)) nd[[v]] <- HEALTHY_PROFILE[[v]]
    }
    
    # A rank-deficient cell leaves some coefficients NA, usually a country
    # dummy that is collinear with the wave dummies because that country is
    # observed at only some waves. coef(m) keeps those entries, vcovCL drops
    # them, so the design has to be cut to the estimable set to match `draws`.
    keep <- !is.na(coef(m))
    X <- model.matrix(delete.response(terms(m)), data = nd, xlev = m$xlevels)
    X <- X[, keep, drop = FALSE]
    if (!identical(colnames(X), names(coef(m))[keep]))
      stop("design mismatch in cum_draws")
    h   <- plogis(draws %*% t(X))
    cum_mix <- cum_mix + cw$w[k] * (1 - t(apply(1 - h, 1, cumprod)))
  }
  cum_mix
}

# the fitted model, its draws and the country weights for one cell
cell_fit <- function(long, region_name, sex_val, model_type) {
  d <- long %>% filter(region == region_name, female == sex_val)
  if (!nrow(d)) return(NULL)
  # restrict first, then build the formula; see the note in fit_dt. Here the
  # old order was survivable rather than fatal, because cell_fit has no
  # anyNA() guard and cum_draws drops NA coefficients with its keep vector,
  # so an aliased constant term was silently ignored. That is why the
  # cumulative incidence curves have a Mediterranean memory cell for men
  # while Table 3 does not. Fixing the order removes the aliased term
  # instead of carrying it, and changes no estimate: a column that is
  # constant at zero contributes nothing either way.
  f <- as.formula(paste("event ~", build_rhs(model_type, region_name, d)))
  d <- d[complete.cases(model.frame(f, data = d, na.action = na.pass)), ]
  if (!nrow(d) || sum(d$event) < 10) return(NULL)
  f <- as.formula(paste("event ~", build_rhs(model_type, region_name, d)))
  
  m  <- glm(f, data = d, family = binomial)
  cw <- d %>% distinct(id, country, weight) %>%
    filter(!is.na(weight), weight > 0) %>%
    group_by(country) %>% summarise(w = sum(weight), .groups = "drop")
  if (!nrow(cw)) cw <- d %>% distinct(id, country) %>% count(country, name = "w")
  cw$w <- cw$w / sum(cw$w)
  
  waves <- sort(unique(d$wave))
  
  # Same rank-deficiency handling as cum_draws. mvrnorm fails with
  # "incompatible arguments" if the mean vector keeps the NA coefficients
  # while the covariance matrix has already dropped them.
  keep <- !is.na(coef(m))
  if (!all(keep))
    cat(sprintf("    rank-deficient (%s %s %s), dropped: %s\n",
                region_name, ifelse(sex_val == 1, "Female", "Male"),
                model_type, paste(names(coef(m))[!keep], collapse = ", ")))
  
  list(m = m, d = d, cw = cw, waves = waves,
       times = d$time[match(waves, d$wave)],
       draws = MASS::mvrnorm(N_SIM, coef(m)[keep],
                             vcovCL(m, cluster = d$id)))
}

cum_incidence <- function(long, region_name, sex_val, outcome_label,
                          model_type = "M3", profile = "none") {
  
  fit <- cell_fit(long, region_name, sex_val, model_type)
  if (is.null(fit)) return(NULL)
  sex_lab <- ifelse(sex_val == 1, "Female", "Male")
  
  series <- data.frame(); points <- data.frame()
  for (ed in levels(fit$d$education_3cat)) {
    
    cum_mix <- cum_draws(fit$m, fit$draws, fit$d, fit$cw, fit$waves, ed,
                         model_type, profile)
    
    base <- data.frame(outcome = outcome_label, region = region_name,
                       sex = sex_lab, education = ed, model = model_type,
                       profile = profile)
    
    series <- bind_rows(series, cbind(base, data.frame(
      wave = fit$waves, time = fit$times,
      cum_inc = apply(cum_mix, 2, mean),
      ci_low  = apply(cum_mix, 2, quantile, .025),
      ci_high = apply(cum_mix, 2, quantile, .975),
      n_countries = nrow(fit$cw))))
    
    # every horizon that this cell can reach, so Part 10 gets its sensitivity
    for (tt in HORIZON_ALL[[outcome_label]]) {
      if (max(fit$times) < tt) next
      vals <- apply(cum_mix, 1, cum_at, times = fit$times, t = tt)
      points <- bind_rows(points, cbind(base, data.frame(
        horizon = tt, main = isTRUE(tt == HORIZON[[outcome_label]]),
        cum_inc = mean(vals),
        ci_low  = unname(quantile(vals, .025)),
        ci_high = unname(quantile(vals, .975)),
        n_countries = nrow(fit$cw))))
    }
  }
  list(series = series, points = points)
}

RUNS <- list(c(model = "M1", profile = "none"),
             c(model = "M3", profile = "healthy"),
             c(model = "M3", profile = "average"))

cuminc <- data.frame(); cumpoint <- data.frame()
for (rn in RUNS) {
  for (nm in names(LONG)) for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
    r <- cum_incidence(LONG[[nm]], reg, sx, nm, rn[["model"]], rn[["profile"]])
    if (is.null(r)) next
    cuminc   <- bind_rows(cuminc, r$series)
    cumpoint <- bind_rows(cumpoint, r$points)
  }
  cat("  ", rn[["model"]], rn[["profile"]], "done\n")
}
cuminc <- as_tibble(cuminc); cumpoint <- as_tibble(cumpoint)

cat("\nCells:", nrow(distinct(cuminc, outcome, region, sex, model, profile)),
    "of", 30 * length(RUNS), "\n")

for (mp in list(c("M1", "none"), c("M3", "healthy"))) {
  cat("\n=== ", mp[1], " / ", mp[2], ", at the main horizon ===\n\n", sep = "")
  cumpoint %>% filter(model == mp[1], profile == mp[2], main) %>%
    mutate(val = sprintf("%.1f (%.1f-%.1f)", 100 * cum_inc, 100 * ci_low,
                         100 * ci_high)) %>%
    select(outcome, horizon, sex, region, education, val) %>%
    pivot_wider(names_from = education, values_from = val) %>%
    arrange(outcome, sex, match(region, REGIONS_ALL)) %>%
    print(n = 40, width = 200)
}

# support for the profile, by region: a profile almost nobody in the US
# matches would make the cross-national prediction an extrapolation
cat("\n=== Share of the baseline sample matching the healthy profile ===\n")
print(round(100 * sapply(split(phys_base, phys_base$region), function(d) {
  ok <- rep(TRUE, nrow(d))
  for (v in names(HEALTHY_PROFILE)) ok <- ok & (d[[v]] == HEALTHY_PROFILE[[v]])
  mean(ok, na.rm = TRUE)
}), 1))

write.csv(cuminc,   file.path(OUT_DIR, "V5_07_cuminc_series.csv"),  row.names = FALSE)
write.csv(cumpoint, file.path(OUT_DIR, "V5_07_cuminc_horizon.csv"), row.names = FALSE)


# =============================================================================
# PART 7b: MEMORY, LINEAR MIXED MODELS
# =============================================================================
# Moved here from Part E of the sensitivity script, because the mixed model is
# now the main memory result. Unchanged except for this header.
#
# It fits two specifications on the same cells: "as published", with baseline
# age as a level term only, and "age x time", which lets the rate of decline
# vary with baseline age. The reported results are "age x time". The other is
# kept so the supplement can report the comparison without refitting.
#
# Writes S5_E2_memory_lmm_age.csv       every fixed-effect coefficient
#        S5_E2_memory_predicted_age.csv predicted trajectories by education
#        S5_E2_memory_gaps.csv          the education gap at 0, 6 and 8 years


HORIZON_DEFAULT <- list(ADL = 8, IADL = 8, Memory = 6)
# HORIZON is saved into main_*.RData, so a load() from an older session can
# replace the list with a bare value and every HORIZON[["Memory"]] below then
# fails. Repair it here rather than failing deep inside a model or figure
# block, and say so, since a silent repair would hide a stale workspace.
if (!is.list(HORIZON) ||
    !all(c("ADL", "IADL", "Memory") %in% names(HORIZON))) {
  warning("HORIZON was not the expected list; restoring the default. ",
          "Check for a load() of an older RData in this session.",
          call. = FALSE)
  HORIZON <- HORIZON_DEFAULT
}

# How far the memory trajectories are PREDICTED, as opposed to reported. HRS
# observes recall at 0, 2, 4, 6 and 8 years and SHARE at 0, 2, 4 and 7, so a
# value at 8 is inside the range for the HRS and one year beyond SHARE's last
# memory interview. Six years is what the main text reports, because it is the
# only horizon inside both surveys' observed range; eight is produced for the
# supplement and must be labelled as extrapolated for the European regions.
MEM_MAX <- 8

stopifnot(all(sapply(c("hrs_b", "share_b", "mem_base", "MEM_HRS_WAVES",
                       "MEM_SHARE_WAVES", "REF_AGE", "REGIONS_ALL",
                       "MED_M3", "HORIZON"), exists)),
          is.function(build_person_wave))
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(lme4)})

OUT_DIR <- if (exists("OUT_DIR")) OUT_DIR else "output"

LMM_PROFILE <- c(age_c = 0, poor_health = 0, partnered_obs = 1,
                 has_mobility = 0, has_hypertension = 0, has_diabetes = 0,
                 obese = 0, ever_smoker = 0)
LMM_PROFILE_LAB <- "partnered, no poor health, no mobility limitation, no hypertension, no diabetes, not obese, never smoked, age 65"

BASE <- list(
  M1 = "education_3cat * time + age_c",
  M2 = "education_3cat * time + age_c + poor_health + partnered_obs + has_mobility",
  M3 = "education_3cat * time + age_c + poor_health + partnered_obs + has_mobility + has_hypertension + has_diabetes + obese + ever_smoker")

# the only change: let the slope vary with baseline age
SPECS <- list(
  "as published"   = BASE,
  "age x time"     = lapply(BASE, function(r) sub("age_c", "age_c + age_c:time",
                                                  r, fixed = TRUE)))

mem_cont <- bind_rows(build_person_wave(hrs_b,   MEM_HRS_WAVES,   "mem"),
                      build_person_wave(share_b, MEM_SHARE_WAVES, "mem")) %>%
  filter(!is.na(value)) %>%
  inner_join(mem_base, by = "id") %>%
  mutate(age_c = age_baseline - REF_AGE)

cat(nrow(mem_cont), "observations on", length(unique(mem_cont$id)), "people\n")

co <- data.frame(); pr_out <- data.frame(); gaps <- data.frame()

for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
  sex_lab <- ifelse(sx == 1, "Female", "Male")
  d0 <- mem_cont %>% filter(region == reg, female == sx)
  if (nrow(d0) < 200) next
  
  vars3 <- c("value", "education_3cat", "time", "age_c", MED_M3)
  if (USE_COUNTRY_FE && reg != "US") vars3 <- c(vars3, "country")
  d <- d0[complete.cases(d0[, vars3, drop = FALSE]), ]
  if (nrow(d) < 200) next
  cat(sprintf("  %-14s %-6s %d obs on %d people\n", reg, sex_lab,
              nrow(d), length(unique(d$id))))
  
  for (sp in names(SPECS)) for (mt in names(SPECS[[sp]])) {
    rhs <- SPECS[[sp]][[mt]]
    if (USE_COUNTRY_FE && reg != "US") rhs <- paste(rhs, "+ factor(country)")
    f <- as.formula(paste("value ~", rhs, "+ (1 + time | id)"))
    m <- try(lmer(f, data = d, REML = FALSE,
                  control = lmerControl(optimizer = "bobyqa")), silent = TRUE)
    if (inherits(m, "try-error")) {
      m <- try(lmer(update(f, . ~ . - (1 + time | id) + (1 | id)),
                    data = d, REML = FALSE), silent = TRUE)
      if (!inherits(m, "try-error"))
        message("    ", reg, " ", sex_lab, " ", mt, " [", sp,
                "]: random slope failed, intercept only")
    }
    if (inherits(m, "try-error")) {
      message("    dropped ", reg, " ", sex_lab, " ", mt, " [", sp, "]"); next
    }
    
    # keep EVERY coefficient, not only the education terms, so the full
    # model can be tabulated without another run
    cf <- summary(m)$coefficients
    co <- bind_rows(co, data.frame(
      spec = sp, model = mt, region = reg, sex = sex_lab,
      term = rownames(cf),
      part = ifelse(grepl(":time|time:", rownames(cf)), "slope", "level"),
      estimate = cf[, 1], se = cf[, 2],
      ci_low = cf[, 1] - 1.96 * cf[, 2], ci_high = cf[, 1] + 1.96 * cf[, 2],
      p = 2 * pnorm(-abs(cf[, 3])),
      n_obs = nrow(d), n_people = length(unique(d$id)), row.names = NULL))
    
    nd <- expand.grid(education_3cat = levels(d$education_3cat),
                      time = seq(0, MEM_MAX, by = 0.5),
                      stringsAsFactors = FALSE)
    nd$education_3cat <- factor(nd$education_3cat,
                                levels = levels(d$education_3cat))
    for (v in names(LMM_PROFILE)) nd[[v]] <- LMM_PROFILE[[v]]
    if (USE_COUNTRY_FE && reg != "US")
      nd$country <- d$country[which.max(table(d$country))]
    
    pp <- try({
      tt   <- delete.response(terms(m, fixed.only = TRUE))
      xlev <- stats::.getXlevels(tt, m@frame)
      X    <- model.matrix(tt, data = nd, xlev = xlev)
      fe   <- fixef(m)
      stopifnot(all(names(fe) %in% colnames(X)))
      X    <- X[, names(fe), drop = FALSE]
      V    <- as.matrix(vcov(m))
      fit  <- as.vector(X %*% fe)
      list(fit = fit, se = sqrt(pmax(0, rowSums((X %*% V) * X))))
    }, silent = TRUE)
    if (!inherits(pp, "try-error"))
      pr_out <- bind_rows(pr_out, data.frame(
        spec = sp, model = mt, region = reg, sex = sex_lab,
        education = as.character(nd$education_3cat), time = nd$time,
        recall = pp$fit, se = pp$se,
        ci_low = pp$fit - 1.96 * pp$se, ci_high = pp$fit + 1.96 * pp$se,
        profile = LMM_PROFILE_LAB, row.names = NULL))
    
    # ---- education contrasts at baseline and at the horizon --------------
    # The gap at the end is level + T * slope, a contrast of two correlated
    # parameters, so its standard error needs the covariance between them.
    # Build the two design rows, difference them, and apply the same delta
    # method used for the predictions: for a contrast c, var = c' V c.
    ct <- try({
      tt   <- delete.response(terms(m, fixed.only = TRUE))
      xlev <- stats::.getXlevels(tt, m@frame)
      fe   <- fixef(m); V <- as.matrix(vcov(m))
      base <- nd[1, , drop = FALSE]
      out  <- data.frame()
      for (ed in setdiff(levels(d$education_3cat), "High")) {
        for (TT in c(0, HORIZON[["Memory"]], MEM_MAX)) {
          r1 <- r0 <- base
          r1$education_3cat <- factor(ed,      levels = levels(d$education_3cat))
          r0$education_3cat <- factor("High",  levels = levels(d$education_3cat))
          r1$time <- r0$time <- TT
          X1 <- model.matrix(tt, data = r1, xlev = xlev)[, names(fe), drop = FALSE]
          X0 <- model.matrix(tt, data = r0, xlev = xlev)[, names(fe), drop = FALSE]
          cc <- as.vector(X1 - X0)
          est <- sum(cc * fe)
          se  <- sqrt(max(0, as.numeric(t(cc) %*% V %*% cc)))
          out <- rbind(out, data.frame(
            spec = sp, model = mt, region = reg, sex = sex_lab,
            education = ed, time = TT, estimate = est, se = se,
            ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
            p = 2 * pnorm(-abs(est / se)), row.names = NULL))
        }
      }
      out
    }, silent = TRUE)
    if (!inherits(ct, "try-error")) gaps <- bind_rows(gaps, ct)
  }
}

write.csv(co,     file.path(OUT_DIR, "S5_E2_memory_lmm_age.csv"),       row.names = FALSE)
write.csv(pr_out, file.path(OUT_DIR, "S5_E2_memory_predicted_age.csv"), row.names = FALSE)
write.csv(gaps,   file.path(OUT_DIR, "S5_E2_memory_gaps.csv"),          row.names = FALSE)

# ---- the comparison this was run for --------------------------------------
cmp <- co %>%
  filter(grepl("^education_3cat", term)) %>%
  select(spec, model, region, sex, term, part, estimate, p) %>%
  pivot_wider(names_from = spec, values_from = c(estimate, p))

for (pt in c("level", "slope")) {
  cat("\n===== ", pt, ", low and medium vs high =====\n", sep = "")
  cmp %>% filter(part == pt, model %in% c("M1", "M3")) %>%
    mutate(edu   = ifelse(grepl("Low", term), "Low", "Medium"),
           shift = `estimate_age x time` - `estimate_as published`,
           flip  = (`p_as published` < .05) != (`p_age x time` < .05)) %>%
    arrange(model, sex, match(region, REGIONS_ALL), edu) %>%
    select(model, region, sex, edu, `estimate_as published`,
           `estimate_age x time`, shift, flip) %>%
    as_tibble() %>% print(n = 100, width = 200)
}

cat("\nCells where the verdict changes:",
    sum((cmp$`p_as published` < .05) != (cmp$`p_age x time` < .05)),
    "of", nrow(cmp), "\n")
cat("Median absolute shift in the education-by-time terms:",
    sprintf("%.4f words/year\n",
            median(abs(cmp$`estimate_age x time`[cmp$part == "slope"] -
                         cmp$`estimate_as published`[cmp$part == "slope"]))))
cat("Median absolute shift in the level terms:",
    sprintf("%.4f words\n",
            median(abs(cmp$`estimate_age x time`[cmp$part == "level"] -
                         cmp$`estimate_as published`[cmp$part == "level"]))))

# ---- the education gap at baseline and at the horizon ----------------------
if (nrow(gaps)) {
  cat("\n===== education gap in words, baseline and horizon =====\n")
  gaps %>%
    filter(spec == "age x time", model %in% c("M1", "M3")) %>%
    mutate(val = sprintf("%+.2f (%+.2f to %+.2f)%s", estimate, ci_low, ci_high,
                         ifelse(p < .05, "*", ""))) %>%
    select(model, region, sex, education, time, val) %>%
    tidyr::pivot_wider(names_from = c(education, time), values_from = val) %>%
    arrange(model, sex, match(region, REGIONS_ALL)) %>%
    as_tibble() %>% print(n = 40, width = 250)
  cat("\nA gap at time 0 must reproduce the education main effect in the\n")
  cat("coefficient table; if it does not, the contrast is built wrong.\n")
}

### to estimate how much is the mean value of scores

mem_cont %>%
  group_by(region, female, education_3cat, wave) %>%
  summarise(n = n(), mean_tr20 = mean(value), sd = sd(value), .groups = "drop") %>%
  mutate(sex = ifelse(female == 1, "Female", "Male")) %>%
  write.csv("T_00_memory_means_by_sex.csv", row.names = FALSE)
mem_base %>% group_by(region, female, education_3cat) %>%
  summarise(n = n(), mean_tr20 = mean(tr20_baseline),
            sd = sd(tr20_baseline), .groups = "drop") %>%
  write.csv("T_00_memory_baseline_by_sex.csv", row.names = FALSE)
# =============================================================================
# PART 8: KHB, EDUCATION AS THE EXPOSURE
# =============================================================================

# The method. Residualise each mediator on the exposure and the controls, then
# fit the outcome model with the RESIDUALS in place of the mediators. That
# model's exposure coefficient is the total effect, on the same scale as the
# full model's, so the difference between them is a clean indirect effect.
#
# Here the matrices are built once and the draws
# call glm.fit() and .lm.fit() on row indice

cat("\n=== KHB, education as the exposure ===\n")

khb_cell <- function(long, region_name, sex_val, outcome_label, med_set,
                     spec_label) {
  
  d <- long %>% filter(region == region_name, female == sex_val)
  if (!nrow(d)) return(NULL)
  
  # the full model's complete cases, so reduced and full share a sample
  f_full <- as.formula(paste("event ~", build_rhs("M3", region_name, d)))
  d <- d[complete.cases(model.frame(f_full, data = d, na.action = na.pass)), ]
  if (!nrow(d) || sum(d$event) < 20) return(NULL)
  
  # X: exposure plus controls. M: the mediators under test.
  f_base <- as.formula(paste("~", build_rhs("M1", region_name, d)))
  X <- model.matrix(f_base, data = d)
  M <- as.matrix(d[, med_set, drop = FALSE])
  y <- d$event
  k <- grep("education_3catLow", colnames(X))
  if (!length(k)) return(NULL)
  
  # one KHB pass on a given set of rows
  pass <- function(idx) {
    Xi <- X[idx, , drop = FALSE]; Mi <- M[idx, , drop = FALSE]; yi <- y[idx]
    R  <- .lm.fit(Xi, Mi)$residuals               # mediators net of X
    if (is.null(dim(R))) R <- matrix(R, ncol = 1)
    colnames(R) <- paste0("r_", colnames(Mi))
    
    fit_tot <- glm.fit(cbind(Xi, R),  yi, family = binomial())
    fit_dir <- glm.fit(cbind(Xi, Mi), yi, family = binomial())
    if (!fit_tot$converged || !fit_dir$converged ||
        anyNA(coef(fit_tot)) || anyNA(coef(fit_dir))) return(c(NA, NA))
    c(total = unname(coef(fit_tot)[k]), direct = unname(coef(fit_dir)[k]))
  }
  
  pt <- pass(seq_len(nrow(d)))
  if (anyNA(pt)) return(NULL)
  
  # cluster bootstrap: resample people, take all of their rows
  ids  <- unique(d$id)
  rows <- split(seq_len(nrow(d)), d$id)
  B    <- matrix(NA_real_, nrow = N_BOOT, ncol = 2)
  for (b in seq_len(N_BOOT)) {
    idx <- unlist(rows[sample(ids, length(ids), replace = TRUE)],
                  use.names = FALSE)
    B[b, ] <- tryCatch(pass(idx), error = function(e) c(NA, NA))
  }
  ok <- stats::complete.cases(B)
  Bi <- B[ok, 1] - B[ok, 2]                       # bootstrap indirect effects
  
  q <- function(v) unname(quantile(v, c(.025, .975), na.rm = TRUE))
  # Bootstrap p for the indirect effect, twice the smaller tail, in the
  # (1 + r) / (B + 1) form. The naive proportion returns exactly 0 whenever no
  # resample crosses zero, which would be reported as p < .001 when 500 draws
  # cannot resolve anything below 2/501.
  nb    <- length(Bi)
  p_ind <- min(1, 2 * (1 + min(sum(Bi <= 0), sum(Bi >= 0))) / (nb + 1))
  data.frame(
    outcome = outcome_label, region = region_name,
    sex = ifelse(sex_val == 1, "Female", "Male"), spec = spec_label,
    term = c("total", "direct", "indirect"),
    estimate = c(pt[["total"]], pt[["direct"]], pt[["total"]] - pt[["direct"]]),
    ci_low  = c(q(B[ok, 1])[1], q(B[ok, 2])[1], q(Bi)[1]),
    ci_high = c(q(B[ok, 1])[2], q(B[ok, 2])[2], q(Bi)[2]),
    pct_mediated = c(NA, NA,
                     100 * (pt[["total"]] - pt[["direct"]]) / pt[["total"]]),
    p_boot = c(NA, NA, p_ind),
    n_boot_ok = sum(ok), n_events = sum(y), n_people = length(ids),
    row.names = NULL)
}

khb <- data.frame()
for (nm in KHB_OUTCOMES) {
  for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
    for (sp in list(list(MED_M2, "M1 to M2"), list(MED_M3, "M1 to M3"))) {
      r <- khb_cell(LONG[[nm]], reg, sx, nm, sp[[1]], sp[[2]])
      if (!is.null(r)) khb <- bind_rows(khb, r)
    }
  }
  cat("  ", nm, "done\n")
}
khb <- as_tibble(khb)

# Benjamini-Hochberg across the ten region-sex cells within each outcome and
# spec, so the indirect effects are not read off ten uncorrected intervals
fdr <- khb %>% filter(term == "indirect") %>%
  group_by(outcome, spec) %>%
  mutate(p_fdr = p.adjust(p_boot, method = "BH")) %>%
  ungroup() %>%
  select(outcome, region, sex, spec, p_fdr)

khb <- khb %>%
  left_join(fdr, by = c("outcome", "region", "sex", "spec")) %>%
  mutate(p_fdr        = ifelse(term == "indirect", p_fdr, NA_real_),
         sig_indirect = ifelse(term == "indirect", p_fdr < .05, NA))

# A share mediated is indirect / total, so it is uninterpretable wherever the
# total effect is indistinguishable from zero: the denominator is noise and the
# percentage runs off to several hundred. Memory is excluded for this reason;
# the same test has to be applied cell by cell to the functional outcomes,
# because IADL among Nordic women has a null total effect too.
khb <- khb %>%
  group_by(outcome, region, sex, spec) %>%
  mutate(total_sig = {
    i <- which(term == "total")
    !(ci_low[i] <= 0 & ci_high[i] >= 0)
  },
  # fully mediated: the direct effect's interval covers zero
  fully_mediated = {
    i <- which(term == "direct")
    ci_low[i] <= 0 & ci_high[i] >= 0
  }) %>%
  ungroup() %>%
  mutate(pct_mediated = ifelse(term == "indirect" & total_sig, pct_mediated,
                               NA_real_))

cat("\nCells where the total effect is null, so no share is reported:\n")
khb %>% filter(term == "indirect", !total_sig) %>%
  select(outcome, region, sex, spec, estimate, ci_low, ci_high) %>%
  print(n = 20)

cat("\nIndirect effect and share mediated, M1 to M3:\n")
khb %>% filter(spec == "M1 to M3", term == "indirect") %>%
  mutate(val = sprintf("%.3f (%.3f to %.3f)  %s%s",
                       estimate, ci_low, ci_high,
                       ifelse(is.na(pct_mediated), "total null",
                              sprintf("%.0f%%", pct_mediated)),
                       ifelse(fully_mediated %in% TRUE, " (full)", ""))) %>%
  select(outcome, region, sex, val) %>%
  pivot_wider(names_from = sex, values_from = val) %>%
  arrange(outcome, match(region, REGIONS_ALL)) %>% print(n = 30, width = 200)

write.csv(khb, file.path(OUT_DIR, "V5_08_khb_education.csv"), row.names = FALSE)


# =============================================================================
# PART 9: WHAT EXPLAINS THE US-EUROPE GAP (Not in the manuscript but is for informative purposes)
# =============================================================================
# This is the cross-national counterpart of Part 8, and it is NOT a KHB.
#
# A pooled KHB with region as the exposure cannot be fitted here. Country fixed
# effects are impossible when country is the exposure, and pooling HRS with
# SHARE breaks factor(wave): HRS waves are 12 to 16 and SHARE waves 6 to 9, so
# the wave dummies are disjoint by survey and separate the US from Europe
# perfectly.
#
# Instead this is g-computation. Take the American model and ask what American
# incidence would be if Americans carried the European region's mediator
# distribution. The share of the gap that closes is the share attributable to
# the risk factor burden. It is on the probability scale, so no rescaling
# problem arises, and the intervals come from the same simulation draws.
#
#   A  US, at US mediator means
#   B  US model, at region r's mediator means   (the counterfactual)
#   C  region r, at its own means
#   gap = A - C,  explained = A - B,  residual = B - C
#
# The counterfactual uses American coefficients. Doing it the other way round,
# giving Europeans the American risk profile, is reported alongside, because
# which population supplies the coefficients is a choice and the two can differ.

cat("\n=== US-Europe gap, decomposed by g-computation ===\n")

med_means <- function(d) sapply(MED_M3, function(v) mean(d[[v]], na.rm = TRUE))

gap_decomp <- function(long, outcome_label, sex_val) {
  tt  <- HORIZON[[outcome_label]]
  us  <- cell_fit(long, "US", sex_val, "M3")
  if (is.null(us)) return(NULL)
  sex_lab <- ifelse(sex_val == 1, "Female", "Male")
  out <- data.frame()
  
  for (reg in setdiff(REGIONS_ALL, "US")) {
    eu <- cell_fit(long, reg, sex_val, "M3")
    if (is.null(eu)) next
    m_us <- med_means(us$d); m_eu <- med_means(eu$d)
    
    for (ed in levels(us$d$education_3cat)) {
      at <- function(fit, means) {
        cm <- cum_draws(fit$m, fit$draws, fit$d, fit$cw, fit$waves, ed,
                        "M3", "average", mediator_means = means)
        apply(cm, 1, cum_at, times = fit$times, t = tt)
      }
      if (max(us$times) < tt || max(eu$times) < tt) next
      
      A  <- at(us, m_us)      # US as it is
      B  <- at(us, m_eu)      # Americans with European risk factors
      Cc <- at(eu, m_eu)      # the European region as it is
      D  <- at(eu, m_us)      # Europeans with American risk factors
      
      gap <- A - Cc
      out <- bind_rows(out, data.frame(
        outcome = outcome_label, sex = sex_lab, region = reg, education = ed,
        horizon = tt,
        us_inc = mean(A), eu_inc = mean(Cc), gap = mean(gap),
        gap_low = unname(quantile(gap, .025)),
        gap_high = unname(quantile(gap, .975)),
        explained_us_coef  = mean(A - B),
        pct_expl_us_coef   = 100 * mean(A - B) / mean(gap),
        explained_eu_coef  = mean(D - Cc),
        pct_expl_eu_coef   = 100 * mean(D - Cc) / mean(gap),
        row.names = NULL))
    }
  }
  out
}

gapdec <- data.frame()
for (nm in names(LONG)) for (sx in c(0, 1)) {
  gapdec <- bind_rows(gapdec, gap_decomp(LONG[[nm]], nm, sx))
  cat("  ", nm, ifelse(sx == 1, "women", "men"), "done\n")
}
gapdec <- as_tibble(gapdec)

cat("\nUS minus region, percentage points, and the share the risk factor\n")
cat("distribution accounts for (American coefficients):\n\n")
gapdec %>% filter(education == "Low") %>%
  mutate(val = sprintf("%+.1f (%+.1f to %+.1f)  %.0f%% expl.",
                       100 * gap, 100 * gap_low, 100 * gap_high,
                       pct_expl_us_coef)) %>%
  select(outcome, sex, region, val) %>%
  pivot_wider(names_from = sex, values_from = val) %>%
  arrange(outcome, match(region, REGIONS_ALL)) %>% print(n = 30, width = 200)

write.csv(gapdec, file.path(OUT_DIR, "V5_09_gap_decomposition.csv"),
          row.names = FALSE)


# =============================================================================
# PART 10: DOES THE CROSS-NATIONAL FINDING DEPEND ON THE HORIZON?
# =============================================================================
# HRS and SHARE stop at different times, so the question is whether the
# ranking is an artefact of where the clock stops. It is not, for
# ADL and IADL, where the gap grows steadily with the horizon. It is also not
# the reason memory shows no gap: the memory gap sits near zero at every
# horizon rather than being squeezed out by a short one.

cat("\n=== Horizon sensitivity ===\n\n")

horizon_tab <- cumpoint %>%
  filter(model == "M3", profile == "healthy") %>%
  select(outcome, region, sex, education, horizon, cum_inc) %>%
  pivot_wider(names_from = region, values_from = cum_inc) %>%
  filter(if_all(all_of(REGIONS_ALL), ~ !is.na(.x))) %>%
  rowwise() %>%
  mutate(eu_mean = mean(c_across(all_of(setdiff(REGIONS_ALL, "US")))),
         gap_pp  = 100 * (US - eu_mean),
         us_rank = rank(-c_across(all_of(REGIONS_ALL)))[1]) %>%
  ungroup()

horizon_tab %>%
  group_by(outcome, horizon) %>%
  summarise(mean_gap_pp = round(mean(gap_pp), 1),
            us_worst_in = sprintf("%d of %d", sum(us_rank == 1), n()),
            .groups = "drop") %>%
  arrange(outcome, horizon) %>% print(n = 30)

write.csv(horizon_tab, file.path(OUT_DIR, "V5_10_horizon_sensitivity.csv"),
          row.names = FALSE)

cat("\n=== US high-educated vs European low-educated ===\n\n")

crossover <- cumpoint %>%
  filter(main, education %in% c("High", "Low")) %>%
  select(outcome, horizon, model, profile, region, sex, education, cum_inc) %>%
  { inner_join(
    filter(., region == "US", education == "High") %>%
      select(outcome, horizon, model, profile, sex, us_high = cum_inc),
    filter(., region != "US", education == "Low") %>%
      select(outcome, horizon, model, profile, sex, region, eu_low = cum_inc),
    by = c("outcome", "horizon", "model", "profile", "sex")) } %>%
  mutate(us_worse = us_high > eu_low)

crossover %>%
  group_by(outcome, model, profile, horizon) %>%
  summarise(holds = sprintf("%d of %d", sum(us_worse), n()),
            us_high = sprintf("%.1f", 100 * mean(us_high)),
            eu_low_range = sprintf("%.1f to %.1f", 100 * min(eu_low),
                                   100 * max(eu_low)), .groups = "drop") %>%
  arrange(outcome, model, profile) %>% print(n = 30, width = 200)

write.csv(crossover, file.path(OUT_DIR, "V5_10_crossover.csv"), row.names = FALSE)


# =============================================================================
# PART 10b: SHARED FIGURE CONSTANTS
# =============================================================================
# These live above the save so that 02_sensitivity_v6.R can draw its own
# figures in the same house style when it is run in a fresh session, where it
# has only main_v6.RData to work from. Part 12 below and Part G of the
# sensitivity script both use them.

FIG_DIR <- OUT_DIR
EDU_COL <- c(High = "#2166AC", Medium = "#8C8C88", Low = "#B2182B")
EDU_ORD <- c("High", "Medium", "Low")          # legend and drawing order
OUT_LAB <- c(ADL = "ADL", IADL = "IADL", Memory = "Verbal memory decline")
SEX_ORD <- c("Male", "Female")                 # men first, so men on top

# where the outcome is a rotated ROW strip it has only one panel's height to
# fit in, so the long label is broken over two lines
CURVE_LAB <- c(ADL = "ADL", IADL = "IADL",
               Memory = "Verbal memory\ndecline")


theme_ms <- function(base = 15) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#e3e3e0", linewidth = .3),
          strip.placement    = "outside",
          strip.text.x       = element_text(face = "bold", size = 17),
          strip.text.y.left  = element_text(face = "bold", size = 17,
                                            angle = 90),
          axis.text          = element_text(size = 15, colour = "#1f2328"),
          axis.title.x       = element_text(size = 16, margin = margin(t = 10)),
          axis.title.y       = element_blank(),
          legend.position    = "bottom",
          legend.title        = element_blank(),
          legend.text        = element_text(size = 15),
          panel.spacing.x    = unit(1.1, "lines"),
          panel.spacing.y    = unit(1.4, "lines"))
}


# =============================================================================
# PART 11: SAVE
# =============================================================================

save(mem_cont, MEM_MAX,
     LONG, LONG_FLOOR, phys_base, mem_base, analytic, person_all,
     models, models_w, models_floor, coefs, age_int, incidence,
     cuminc, cumpoint, khb, gapdec, horizon_tab, crossover,
     hrs_b, share_b, HRS_WAVES, SHARE_WAVES, MEM_HRS_WAVES, MEM_SHARE_WAVES,
     MED_M2, MED_M3, MED_RAW, HEALTHY_PROFILE, HORIZON, HORIZON_ALL, REGIONS_ALL,
     BL_SD, BL_MEAN, MEM_CUT, SD_MULT, REF_AGE, N_SIM, N_BOOT, USE_COUNTRY_FE,
     build_rhs, fit_dt, run_models, cum_incidence, cum_draws, cum_at, cell_fit,
     khb_cell, gap_decomp, person_years, incidence_table, build_long,
     build_person_wave, make_event, trim_to_first_event, get_tr20_col,
     get_wave_value, pois_ci, med_means,
     FIG_DIR, EDU_COL, EDU_ORD, OUT_LAB, SEX_ORD, CURVE_LAB, theme_ms,
     file = file.path(OUT_DIR, "main_v6.RData"))



# =============================================================================
# =============================================================================
# PART 12: FIGURES
# =============================================================================
# Figure 1   three-panel forest plot at six years: predicted ADL onset,
#            predicted IADL onset and predicted word recall, Model 3.
# Figure S1  predicted word recall from baseline to six years, Models 1 and 3.
# Figure S3  the old Figure 1, cumulative incidence at the common horizon,
#            Model 1, kept for the supplement.
# Figure S4  the old Figure 2, the same for Model 3.
#
# Figure S2  cumulative incidence of ADL and IADL onset and predicted word
#            recall, all Model 3, all to six years.

# =============================================================================
# Figure 1. Three-panel forest plot at the six-year horizon.
#   Panel 1: predicted cumulative incidence of ADL onset
#   Panel 2: predicted cumulative incidence of IADL onset
#   Panel 3: predicted word recall from the linear mixed model
# Rows are region x sex, shared across panels. Within each row one point per
# education level, with 95% intervals. Everything is Model 3 at six years.
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)   # install.packages("patchwork") if needed

# ---- paths ------------------------------------------------------------------
# point these at your own output files
F_CUMINC <- "V5_07_cuminc_horizon.csv"
F_MEMORY <- "S5_E2_memory_predicted_age.csv"
OUT      <- "Figure1_forest3.png"

# ---- settings ---------------------------------------------------------------
MODEL   <- "M3"
PROF    <- "average"      # cumulative-incidence profile used for M3
FIG_HORIZON <- 6

REGIONS <- c("US", "Nordic", "Western", "Mediterranean", "Eastern")
EDU     <- c("High", "Medium", "Low")
COL     <- c(High = "#2166AC", Medium = "#8C8C88", Low = "#B2182B")
INK     <- "#1f2328"; INK2 <- "#57606a"; GRID <- "#e3e3e0"

# ---- data -------------------------------------------------------------------
cum <- read.csv(F_CUMINC) %>%
  filter(model == MODEL, profile == PROF, horizon == FIG_HORIZON,
         outcome %in% c("ADL", "IADL")) %>%
  mutate(est = 100 * cum_inc, lo = 100 * ci_low, hi = 100 * ci_high)

mem <- read.csv(F_MEMORY) %>%
  filter(spec == "age x time", model == MODEL, time == FIG_HORIZON) %>%
  mutate(outcome = "Memory", est = recall, lo = ci_low, hi = ci_high)

dat <- bind_rows(
  cum %>% select(outcome, region, sex, education, est, lo, hi),
  mem %>% select(outcome, region, sex, education, est, lo, hi))

# one row per region x sex, US men at the top
strata <- expand.grid(sex = c("Male", "Female"), region = REGIONS,
                      stringsAsFactors = FALSE)[, c("region", "sex")]
strata$label <- paste0(strata$region, ", ",
                       ifelse(strata$sex == "Male", "men", "women"))
strata$y <- rev(seq_len(nrow(strata)))

dat <- dat %>%
  inner_join(strata, by = c("region", "sex")) %>%
  mutate(education = factor(education, levels = EDU),
         # spread the three education points inside the row
         yy = y + unname(c(High = 0.26, Medium = 0,
                           Low = -0.26)[as.character(education)]))

# grey band behind every other region
bands <- strata %>%
  group_by(region) %>%
  summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop") %>%
  mutate(ord = match(region, REGIONS)) %>%
  filter(ord %% 2 == 1)

# ADL and IADL share one x range so the level difference between them shows
cum_max <- 1.03 * max(cum$hi)


# ---- one panel --------------------------------------------------------------
BASE <- 14          # single knob: raise this and everything follows

panel <- function(oc, xlab, xlim = NULL, show_y = FALSE) {
  d <- dat %>% filter(outcome == oc)
  p <- ggplot(d) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = "#f6f6f4") +
    geom_linerange(aes(y = yy, xmin = lo, xmax = hi, colour = education),
                   linewidth = 0.8) +
    geom_point(aes(y = yy, x = est, colour = education), size = 2.6) +
    scale_colour_manual(values = COL, breaks = EDU,
                        labels = paste(EDU, "education"), name = NULL) +
    scale_y_continuous(breaks = strata$y, labels = strata$label,
                       limits = c(0.1, nrow(strata) + 0.9), expand = c(0, 0)) +
    labs(x = xlab, y = NULL) +
    theme_minimal(base_size = BASE) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
          axis.ticks         = element_blank(),
          axis.title.x       = element_text(colour = INK, size = BASE + 2,
                                            margin = margin(t = 10)),
          axis.text.x        = element_text(colour = INK2, size = BASE),
          legend.position    = "bottom")
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  if (!show_y) p <- p + theme(axis.text.y = element_blank())
  else p <- p + theme(axis.text.y = element_text(colour = INK, hjust = 1,
                                                 size = BASE + 1))
  p
}

p1 <- panel("ADL",    "Predicted ADL (%)",   c(0, cum_max), show_y = TRUE)
p2 <- panel("IADL",   "Predicted IADL (%)",  c(0, cum_max))
p3 <- panel("Memory", "Predicted word recall")

fig <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text     = element_text(size = BASE + 2, colour = INK),
        legend.key.width = unit(1.4, "lines"))

ggsave(OUT, fig, width = 16, height = 9.5, dpi = 300, bg = "white")
cat("written", OUT, "\n")


# =============================================================================
# Figure S1. Predicted word recall from baseline to six years, Models 1 and 3
# =============================================================================
# Rows are sex, columns are region, colour is education. Model 1 is drawn
# solid and Model 3 dashed, so the reader can see that adjustment moves the
# level of the trajectories and not their separation.

local({
  d <- read.csv(file.path(OUT_DIR, "S5_E2_memory_predicted_age.csv"))
  d <- subset(d, spec == "age x time" & model %in% c("M1", "M3"))
  d$education <- factor(d$education, levels = c("High", "Medium", "Low"))
  d$sex       <- factor(d$sex,       levels = c("Male", "Female"),
                        labels = c("Men", "Women"))
  d$region    <- factor(d$region,    levels = REGIONS_ALL)
  d$model     <- factor(d$model,     levels = c("M1", "M3"))
  
  p <- ggplot(d, aes(time, recall, colour = education,
                     linetype = model, group = interaction(education, model))) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = education),
                alpha = .12, colour = NA) +
    geom_line(linewidth = .8) +
    facet_grid(sex ~ region) +
    scale_colour_manual(values = c(High = "#2166AC", Medium = "#8C8C88",
                                   Low = "#B2182B"),
                        labels = paste(c("High", "Medium", "Low"),
                                       "education"), name = NULL) +
    scale_fill_manual(values = c(High = "#2166AC", Medium = "#8C8C88",
                                 Low = "#B2182B"), guide = "none") +
    scale_linetype_manual(values = c(M1 = "solid", M3 = "22"),
                          labels = c("Model 1", "Model 3"), name = NULL) +
    scale_x_continuous(breaks = c(0, 2, 4, 6)) +
    labs(x = "Years since baseline", y = "Predicted word recall") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          strip.text       = element_text(face = "bold"))
  
  ggsave(file.path(OUT_DIR, "FigureS1_memory_trajectories.png"), p,
         width = 15, height = 8, dpi = 300, bg = "white")
  cat("written FigureS1_memory_trajectories.png\n")
})


# =============================================================================
# Figure S2. Model 3 to six years: cumulative incidence of ADL and IADL onset,
#            and predicted word recall
# =============================================================================
# Rows are sex by outcome, columns are region, colour is education. The two
# functional rows are cumulative incidence in per cent; the memory row is
# predicted recall in words, so the rows carry free y scales.
#
# SHARE observes at 2, 4, 7 and 9 years, so a plain truncation would end the
# European curves at four. They are interpolated to exactly six instead.

local({
  TMAX <- HORIZON[["Memory"]]
  
  ser <- read.csv(file.path(OUT_DIR, "V5_07_cuminc_series.csv")) %>%
    filter(model == "M3", profile == "average",
           outcome %in% c("ADL", "IADL")) %>%
    mutate(across(c(cum_inc, ci_low, ci_high), ~ 100 * .x))
  
  # anchor at the origin, then cut or interpolate to exactly TMAX
  to_six <- function(d) {
    t <- c(0, d$time); y <- c(0, d$cum_inc)
    lo <- c(0, d$ci_low); hi <- c(0, d$ci_high)
    k <- t <= TMAX
    out <- data.frame(time = t[k], value = y[k], lo = lo[k], hi = hi[k])
    if (max(out$time) < TMAX)
      out <- rbind(out, data.frame(
        time  = TMAX,
        value = approx(t, y,  xout = TMAX)$y,
        lo    = approx(t, lo, xout = TMAX)$y,
        hi    = approx(t, hi, xout = TMAX)$y))
    out
  }
  
  cum <- ser %>%
    arrange(time) %>%
    group_by(outcome, region, sex, education) %>%
    group_modify(~ to_six(.x)) %>%
    ungroup()
  
  mem <- read.csv(file.path(OUT_DIR, "S5_E2_memory_predicted_age.csv")) %>%
    filter(spec == "age x time", model == "M3", time <= TMAX) %>%
    transmute(outcome = "Memory", region, sex, education, time,
              value = recall, lo = ci_low, hi = ci_high)
  
  d <- bind_rows(cum, mem) %>%
    mutate(
      education = factor(education, levels = c("High", "Medium", "Low")),
      region    = factor(region, levels = REGIONS_ALL),
      row = factor(
        paste(ifelse(sex == "Male", "Men", "Women"), outcome),
        levels = c("Men ADL", "Men IADL", "Men Memory",
                   "Women ADL", "Women IADL", "Women Memory")))
  
  labs_row <- c(
    "Men ADL"        = "Men\nADL onset (%)",
    "Men IADL"       = "Men\nIADL onset (%)",
    "Men Memory"     = "Men\nRecall (words)",
    "Women ADL"      = "Women\nADL onset (%)",
    "Women IADL"     = "Women\nIADL onset (%)",
    "Women Memory"   = "Women\nRecall (words)")
  
  p <- ggplot(d, aes(time, value, colour = education, fill = education)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .14, colour = NA) +
    geom_line(linewidth = .8) +
    facet_grid(row ~ region, scales = "free_y", switch = "y",
               labeller = labeller(row = labs_row)) +
    scale_colour_manual(values = c(High = "#2166AC", Medium = "#8C8C88",
                                   Low = "#B2182B"),
                        labels = paste(c("High", "Medium", "Low"),
                                       "education"), name = NULL) +
    scale_fill_manual(values = c(High = "#2166AC", Medium = "#8C8C88",
                                 Low = "#B2182B"), guide = "none") +
    scale_x_continuous(breaks = c(0, 2, 4, 6)) +
    labs(x = "Years since baseline", y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor   = element_blank(),
          strip.placement    = "outside",
          strip.text.y.left  = element_text(angle = 90, face = "bold"),
          strip.text.x       = element_text(face = "bold", size = 14),
          legend.position    = "bottom")
  
  ggsave(file.path(OUT_DIR, "FigureS2_M3_six_years.png"), p,
         width = 14, height = 12, dpi = 300, bg = "white")
  cat("written FigureS2_M3_six_years.png\n")
})

# =============================================================================
# Figures S3 and S4, the cumulative-incidence forest plots
# =============================================================================
#
# Both come from cumpoint, the point estimates at each horizon, which is
# already in memory.


cat("\n=== Figures ===\n")

# ---- Figures 1 and 2: forest plots -----------------------------------------

FOREST_SPEC <- list(M1 = list(profile = "none",    file = "F1_forest_M1.png"),
                    M3 = list(profile = "average", file = "F2_forest_M3.png"))

facet_os <- function(row_var, col_var) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    facet_grid(rows = vars({{ row_var }}), cols = vars({{ col_var }}),
               switch = "y", axes = "all_x")
  else
    facet_grid(rows = vars({{ row_var }}), cols = vars({{ col_var }}),
               switch = "y")
}

forest_fig <- function(model_type, prof) {
  d <- cumpoint %>%
    filter(model == model_type, profile == prof, main) %>%
    mutate(across(c(cum_inc, ci_low, ci_high), ~ 100 * .x),
           region    = factor(region, levels = rev(REGIONS_ALL)),
           sex       = factor(sex, levels = SEX_ORD),
           education = factor(education, levels = rev(EDU_ORD)),
           panel     = factor(OUT_LAB[outcome], levels = OUT_LAB),
           lab = sprintf("%.1f (%.1f-%.1f)", cum_inc, ci_low, ci_high))
  
  # labels sit in one aligned column just right of the widest interval in the
  # whole figure, with room beyond that so they are not clipped
  hi    <- max(d$ci_high)
  lab_x <- hi * 1.10
  x_max <- hi * 2.25
  
  ggplot(d, aes(y = region, colour = education)) +
    # faint rules between regions, so the three education estimates read as
    # one block; a discrete axis puts the five regions at 1 to 5
    geom_hline(yintercept = seq(1.5, length(REGIONS_ALL) - 0.5),
               colour = "#e3e3e0", linewidth = .3) +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high),
                   position = position_dodge(width = .7), linewidth = 1.1) +
    geom_point(aes(x = cum_inc), position = position_dodge(width = .7),
               size = 2.6) +
    geom_text(aes(x = lab_x, label = lab),
              position = position_dodge(width = .7),
              hjust = 0, size = 4.1, show.legend = FALSE) +
    geom_blank(aes(x = x_max)) +
    facet_os(panel, sex) +
    scale_colour_manual(values = EDU_COL, breaks = EDU_ORD,
                        labels = paste(EDU_ORD, "education")) +
    scale_x_continuous(breaks = seq(0, floor(hi / 10) * 10, by = 10),
                       expand = expansion(mult = c(.02, 0))) +
    labs(x = "Predicted cumulative incidence (%)") +
    theme_ms()
}

for (mt in names(FOREST_SPEC)) {
  sp <- FOREST_SPEC[[mt]]
  ggsave(file.path(FIG_DIR, sp$file), forest_fig(mt, sp$profile),
         width = 11.5, height = 14, dpi = 300)
  cat("  ", sp$file, " (", mt, ", profile ", sp$profile, ")\n", sep = "")
}
