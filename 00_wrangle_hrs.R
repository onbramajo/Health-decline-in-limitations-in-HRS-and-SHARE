# =============================================================================
# HRS COMPACT DATASET - WAVE 11 BASELINE (2012) THROUGH WAVE 16 (2022)
# =============================================================================
#
# Reads the RAND HRS Longitudinal File and writes the compact person-level
# file the analysis scripts use: one row per respondent interviewed at wave
# 11, aged 50 to 100, and free of ADL and IADL limitations at that wave.
#
# WHAT IS KEPT, AND WHY
#
#  1. Cognition mode is explicit. From wave 14 RAND splits word recall into
#     r14tr20p (face-to-face and telephone) and r14tr20w (web). Both are kept,
#     along with a mode indicator, so the choice is visible and testable and
#     web respondents are not silently dropped.
#
#  2. Immediate and delayed recall at every wave, not only at baseline and
#     endpoint, so the separate-components check can run on person-wave data.
#
#  3. ADL and IADL component items, so a six-item sensitivity analysis adding
#     toileting is possible. r11adl5a was verified against the data to contain
#     bathing, dressing, eating, bed transfer and walking across a room.
#     Toileting is not in it.
#
#  4. Education as completed years as well as the three categories, so a
#     years-based sensitivity can be run against the ISCED crosswalk.
#
#  Endpoint-only variables are not kept. The analysis is discrete-time on
#  person-wave data and does not use them.
#
# NOTE ON MOBILITY
#   r11mobila is a 0-5 index; the SHARE equivalent r5mobilb is 0-4. The item
#   sets differ: HRS has walking several blocks and walking one block, SHARE
#   has walking 100 metres. It is used only as a binary covariate and is never
#   compared across surveys.
#
# NOTE ON select()
#   select <- dplyr::select below. MASS, which the analysis scripts load,
#   masks select().
#
# =============================================================================

rm(list = ls())

library(tidyverse)
library(haven)

select <- dplyr::select     # guard against MASS and friends

# ---- paths -----------------------------------------------------------------
# Relative to the working directory, which should be the repository root.
# Put the raw RAND HRS file in data/raw/, or copy config_example.R to
# config.R and point RAW_HRS somewhere else. config.R is untracked and is
# sourced here if it exists.

RAW_DIR     <- "data/raw"
DERIVED_DIR <- "data/derived"

RAW_HRS  <- file.path(RAW_DIR, "randhrs1992_2022v1.dta")
OUT_FILE <- file.path(DERIVED_DIR, "hrs_compact.rds")

if (file.exists("config.R")) source("config.R")
dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(RAW_HRS))


# =============================================================================
# LOAD
# =============================================================================

cat("Loading RAND HRS...\n")
HRS <- read_dta(RAW_HRS)
cat("  N =", nrow(HRS), " Variables =", ncol(HRS), "\n\n")


# =============================================================================
# BUILD THE FIVE-ITEM IADL INDEX
# =============================================================================
# Done here, on the raw file, because the baseline exclusion below has to use
# it. If this ran after the filter, the sample would still be selected on
# RAND's index and people already limited on shopping or meals would enter the
# risk set and then count as incident cases at the first follow-up, which
# would replace one artefact with another.

IADL_ITEMS <- c("moneya", "medsa", "shopa", "mealsa", "phonea")

for (w in 11:16) {
  it <- paste0("r", w, IADL_ITEMS)
  stopifnot(all(it %in% names(HRS)))        # stop rather than build a short index
  HRS[[paste0("iadl5_w", w)]] <-
    HRS[[it[1]]] + HRS[[it[2]]] + HRS[[it[3]]] + HRS[[it[4]]] + HRS[[it[5]]]
}

cat("=== FIVE-ITEM IADL INDEX, RAW FILE ===\n")
for (w in 11:16) {
  v <- HRS[[paste0("iadl5_w", w)]]
  cat(sprintf("  iadl5_w%-3d max = %s  mean = %.3f  any = %.1f%%  NA = %.1f%%\n",
              w, max(v, na.rm = TRUE), mean(v, na.rm = TRUE),
              100 * mean(v >= 1, na.rm = TRUE), 100 * mean(is.na(v))))
}

# Does the hand-built index agree with RAND's where both are observed? It
# should, on every case. If this ever prints FALSE the item set has changed.
tgt <- as.numeric(HRS$r11iadl5a)
s5  <- as.numeric(HRS$iadl5_w11)
ok  <- !is.na(s5) & !is.na(tgt)
cat(sprintf("\n  agrees with r11iadl5a on %d jointly observed cases: %s\n",
            sum(ok), all(s5[ok] == tgt[ok])))
cat(sprintf("  RAND fills the index for %d cases where an item is unobserved;\n",
            sum(is.na(s5) & !is.na(tgt))))
cat("  those are treated as missing here, which is SHARE's rule.\n\n")


# =============================================================================
# MODE CHECK FOR COGNITION
# =============================================================================
# Run this before anything else. If a large share of waves 14-15 is web-only,
# using the p variables alone restricts the sample, and that needs stating.

cat("=== COGNITION MODE, WAVES 14 AND 15 ===\n")
for (w in c(14, 15)) {
  p <- HRS[[paste0("r", w, "tr20p")]]
  wv <- HRS[[paste0("r", w, "tr20w")]]
  cat(sprintf("  wave %d:  ftf/tel = %d,  web = %d,  both = %d\n",
              w, sum(!is.na(p)), sum(!is.na(wv)), sum(!is.na(p) & !is.na(wv))))
}
cat("\nIf 'both' is 0, each person is in one mode only, and using the\n")
cat("ftf/tel variable alone drops the web respondents from those waves.\n\n")


# =============================================================================
# BASELINE SAMPLE
# =============================================================================
# Both IADL definitions are printed so the cost of the switch is visible.

cat("=== BASELINE (wave 11, 2012) ===\n")
cat("  interviewed:          ", sum(HRS$r11iwstat == 1, na.rm = TRUE), "\n")
cat("  aged 50-100:          ", HRS %>% filter(r11iwstat == 1, r11agey_e >= 50,
                                               r11agey_e <= 100) %>% nrow(), "\n")
cat("  free of ADL and IADL, RAND 5-item index: ",
    HRS %>% filter(r11iwstat == 1, r11agey_e >= 50, r11agey_e <= 100,
                   r11adl5a == 0, r11iadl5a == 0) %>% nrow(), "\n")
cat("  free of ADL and IADL, hand-built index:  ",
    HRS %>% filter(r11iwstat == 1, r11agey_e >= 50, r11agey_e <= 100,
                   r11adl5a == 0, iadl5_w11 == 0) %>% nrow(), "\n")
cat("  (the second is the sample used; it is smaller because a missing\n")
cat("   component item now excludes rather than counts as no difficulty)\n\n")


# =============================================================================
# BUILD
# =============================================================================

hrs_compact <- HRS %>%
  filter(r11iwstat == 1) %>%
  filter(r11agey_e >= 50 & r11agey_e <= 100) %>%
  filter(r11adl5a == 0 & iadl5_w11 == 0) %>%
  
  select(
    # --- id ---------------------------------------------------------------
    hhidpn,
    
    # --- weights ----------------------------------------------------------
    wtresp_w11 = r11wtresp, wthh_w11 = r11wthh,
    
    # --- time-invariant demographics --------------------------------------
    birth_year = rabyear, birth_month = rabmonth, gender = ragender,
    race = raracem, hispanic = rahispan,
    education_years = raedyrs, education_degree = raedegrm,
    
    # --- baseline demographics --------------------------------------------
    age_baseline = r11agey_e, mstat_baseline = r11mstat,
    
    # --- mortality --------------------------------------------------------
    death_year = radyear, death_month = radmonth,
    
    # --- interview status, every wave --------------------------------------
    iwstat_w11 = r11iwstat, iwstat_w12 = r12iwstat, iwstat_w13 = r13iwstat,
    iwstat_w14 = r14iwstat, iwstat_w15 = r15iwstat, iwstat_w16 = r16iwstat,
    
    # --- baseline health ---------------------------------------------------
    srh_baseline = r11shlt,
    hibp_baseline = r11hibpe, diab_baseline = r11diabe, cancer_baseline = r11cancre,
    lung_baseline = r11lunge, heart_baseline = r11hearte, stroke_baseline = r11stroke,
    arthritis_baseline = r11arthre, psych_baseline = r11psyche,
    cesd_baseline = r11cesd,
    smoking_ever = r11smokev, smoking_now = r11smoken,
    bmi_baseline = r11bmi,
    
    # --- ADL: RAND's summary index, verified to match SHARE's --------------
    # adl5a = bathe, dress, eat, BED, walk across room (no toilet)
    adl_baseline = r11adl5a,
    adl_w12 = r12adl5a, adl_w13 = r13adl5a, adl_w14 = r14adl5a,
    adl_w15 = r15adl5a, adl_w16 = r16adl5a,
    
    # --- IADL: the hand-built five-item index (see the header note) --------
    iadl_baseline = iadl5_w11,
    iadl_w12 = iadl5_w12, iadl_w13 = iadl5_w13, iadl_w14 = iadl5_w14,
    iadl_w15 = iadl5_w15, iadl_w16 = iadl5_w16,
    
    # RAND's own index, kept only to document the discrepancy. NOT analysed.
    iadl_rand_baseline = r11iadl5a,
    
    # --- ADL component items at baseline, for the 6-item sensitivity -------
    bath_baseline = r11batha, dress_baseline = r11dressa, eat_baseline = r11eata,
    bed_baseline = r11beda, toilet_baseline = r11toilta, walkr_baseline = r11walkra,
    
    # --- IADL component items at baseline, for item-level checks -----------
    money_baseline = r11moneya, meds_baseline = r11medsa,
    shop_baseline = r11shopa, meals_baseline = r11mealsa,
    phone_baseline = r11phonea, map_baseline = r11mapa,
    
    # --- mobility (covariate only; see header note) -------------------------
    mobility_baseline = r11mobila,
    
    # --- cognition: baseline and every follow-up wave with data ------------
    imrc_baseline = r11imrc, dlrc_baseline = r11dlrc, tr20_baseline = r11tr20,
    ser7_baseline = r11ser7, bwc20_baseline = r11bwc20, cogtot_baseline = r11cogtot,
    
    imrc_w12 = r12imrc, dlrc_w12 = r12dlrc, tr20_w12 = r12tr20,
    imrc_w13 = r13imrc, dlrc_w13 = r13dlrc, tr20_w13 = r13tr20,
    
    # waves 14-15: p = face-to-face/telephone, w = web. Both kept.
    imrc_w14_p = r14imrcp, dlrc_w14_p = r14dlrcp, tr20_w14_p = r14tr20p,
    imrc_w14_w = r14imrcw, dlrc_w14_w = r14dlrcw, tr20_w14_w = r14tr20w,
    imrc_w15_p = r15imrcp, dlrc_w15_p = r15dlrcp, tr20_w15_p = r15tr20p,
    imrc_w15_w = r15imrcw, dlrc_w15_w = r15dlrcw, tr20_w15_w = r15tr20w,
    
    ser7_w15 = r15ser7p,
    
    # --- self-rated health at follow-up ------------------------------------
    srh_w12 = r12shlt, srh_w13 = r13shlt, srh_w14 = r14shlt,
    srh_w15 = r15shlt, srh_w16 = r16shlt,
    
    # --- wealth -------------------------------------------------------------
    wealth_total = h11atoth, wealth_net = h11atotb, income_hh = h11itot
  ) %>%
  
  mutate(
    # --- metadata -----------------------------------------------------------
    study = "HRS", region = "US", country_name = "United States",
    baseline_wave = 11, baseline_year = 2012,
    
    # --- demographics -------------------------------------------------------
    female = ifelse(gender == 2, 1, 0),
    
    race_eth = case_when(
      hispanic == 1 ~ "Hispanic",
      race == 1 & hispanic == 0 ~ "White",
      race == 2 & hispanic == 0 ~ "Black",
      TRUE ~ "Other"),
    
    # Education: years mapped onto ISCED bands.
    # <12 yrs  -> ISCED 0-2 (less than upper secondary)
    # 12-15    -> ISCED 3-4 (upper secondary / post-secondary non-tertiary)
    # 16+      -> ISCED 5-6 (tertiary)
    # SHARE uses native ISCED. Describe this crosswalk in the Methods.
    education_3cat = case_when(
      education_years < 12                             ~ "Low",
      education_years >= 12 & education_years < 16      ~ "Medium",
      education_years >= 16                             ~ "High",
      TRUE                                              ~ NA_character_),
    
    partnered = case_when(
      mstat_baseline %in% c(1, 2, 3)          ~ 1,
      mstat_baseline %in% c(4, 5, 6, 7, 8)    ~ 0,
      TRUE                                    ~ NA_real_),
    
    # --- health covariates ---------------------------------------------------
    poor_health      = ifelse(srh_baseline %in% c(4, 5), 1, 0),
    has_hypertension = ifelse(hibp_baseline == 1, 1, 0),
    has_diabetes     = ifelse(diab_baseline == 1, 1, 0),
    has_mobility     = ifelse(mobility_baseline >= 1, 1, 0),
    obese            = ifelse(bmi_baseline >= 30, 1, 0),
    ever_smoker      = ifelse(smoking_ever == 1, 1, 0),
    current_smoker   = ifelse(smoking_now == 1, 1, 0),
    depressed        = ifelse(cesd_baseline >= 4, 1, 0),
    
    chronic_count = rowSums(cbind(
      hibp_baseline == 1, diab_baseline == 1, cancer_baseline == 1,
      lung_baseline == 1, heart_baseline == 1, stroke_baseline == 1,
      arthritis_baseline == 1), na.rm = TRUE),
    
    # --- six-item ADL at baseline (adds toileting) ---------------------------
    # For the sensitivity check only. Everyone has adl_baseline == 0 by
    # construction, so this will be 0 or 1 depending on toileting alone.
    adl6_baseline = rowSums(cbind(bath_baseline, dress_baseline, eat_baseline,
                                  bed_baseline, toilet_baseline, walkr_baseline)),
    
    # --- cognition: choose the mode -------------------------------------------
    # PRIMARY: face-to-face/telephone only, consistent with waves 11-13.
    tr20_w14 = tr20_w14_p,
    tr20_w15 = tr20_w15_p,
    imrc_w14 = imrc_w14_p, dlrc_w14 = dlrc_w14_p,
    imrc_w15 = imrc_w15_p, dlrc_w15 = dlrc_w15_p,
    
    # ALTERNATIVE: both modes pooled, with an indicator, for a sensitivity
    tr20_w14_any = coalesce(tr20_w14_p, tr20_w14_w),
    tr20_w15_any = coalesce(tr20_w15_p, tr20_w15_w),
    web_w14 = ifelse(!is.na(tr20_w14_w) & is.na(tr20_w14_p), 1, 0),
    web_w15 = ifelse(!is.na(tr20_w15_w) & is.na(tr20_w15_p), 1, 0),
    
    # --- mortality and observation -------------------------------------------
    died_before_w16 = ifelse(!is.na(death_year) & death_year <= 2021, 1, 0),
    observed_w15 = iwstat_w15 == 1,
    observed_w16 = iwstat_w16 == 1,
    
    # --- wealth ---------------------------------------------------------------
    wealth_quartile = ntile(wealth_total, 4)
  ) %>%
  
  filter(!is.na(education_3cat))

cat("Built: N =", nrow(hrs_compact), " Variables =", ncol(hrs_compact), "\n\n")


# =============================================================================
# CHECKS
# =============================================================================

cat("=== SAMPLE ===\n")
cat(sprintf("  N: %d\n", nrow(hrs_compact)))
cat(sprintf("  Mean age: %.1f (SD %.1f)\n",
            mean(hrs_compact$age_baseline), sd(hrs_compact$age_baseline)))
cat(sprintf("  Female: %.1f%%\n", mean(hrs_compact$female) * 100))

cat("\n=== EDUCATION ===\n")
print(hrs_compact %>% count(education_3cat) %>% mutate(pct = round(100 * n / sum(n), 1)))

cat("\n=== SIX-ITEM ADL AT BASELINE ===\n")
cat("(everyone is 0 on the five-item index; this shows how many have\n")
cat(" toileting difficulty only)\n")
print(table(hrs_compact$adl6_baseline, useNA = "ifany"))

cat("\n=== COGNITION COVERAGE BY WAVE ===\n")
for (v in c("tr20_baseline", "tr20_w12", "tr20_w13", "tr20_w14", "tr20_w15")) {
  cat(sprintf("  %-14s n = %5d  mean = %.2f\n", v,
              sum(!is.na(hrs_compact[[v]])),
              mean(hrs_compact[[v]], na.rm = TRUE)))
}
cat("\n  with both modes pooled:\n")
for (v in c("tr20_w14_any", "tr20_w15_any")) {
  cat(sprintf("  %-14s n = %5d  mean = %.2f\n", v,
              sum(!is.na(hrs_compact[[v]])),
              mean(hrs_compact[[v]], na.rm = TRUE)))
}
cat("\n  web-only respondents: w14 =", sum(hrs_compact$web_w14),
    " w15 =", sum(hrs_compact$web_w15), "\n")

# Covers ADL and IADL, and prints the maximum of each index. The maximum is
# the number that reveals an item-set mismatch between the two surveys.
cat("\n=== ADL/IADL COVERAGE BY WAVE ===\n")
for (v in c("adl_baseline", "adl_w12", "adl_w13", "adl_w14", "adl_w15", "adl_w16",
            "iadl_baseline", "iadl_w12", "iadl_w13", "iadl_w14", "iadl_w15",
            "iadl_w16")) {
  x <- hrs_compact[[v]]
  cat(sprintf("  %-14s n = %5d  max = %s  any limitation = %.1f%%\n", v,
              sum(!is.na(x)), max(x, na.rm = TRUE),
              100 * mean(x >= 1, na.rm = TRUE)))
}
cat("\n  Every max must be 5, here and in the SHARE script. Anything else\n")
cat("  means the item sets have drifted apart again.\n")

cat("\n=== IADL COMPONENT ITEMS AT BASELINE ===\n")
cat("(all zero by construction; shown so the item set is on the record)\n")
for (v in c("money_baseline", "meds_baseline", "shop_baseline",
            "meals_baseline", "phone_baseline", "map_baseline")) {
  x <- hrs_compact[[v]]
  cat(sprintf("  %-16s n = %5d  mean = %.4f\n", v, sum(!is.na(x)),
              mean(x, na.rm = TRUE)))
}

cat("\n=== FOLLOW-UP ===\n")
cat(sprintf("  Observed w15: %.1f%%\n", 100 * mean(hrs_compact$observed_w15)))
cat(sprintf("  Observed w16: %.1f%%\n", 100 * mean(hrs_compact$observed_w16)))
cat(sprintf("  Died by 2021: %.1f%%\n", 100 * mean(hrs_compact$died_before_w16)))

print(hrs_compact %>%
        group_by(web_w15) %>%
        summarise(n = n(),
                  age = mean(age_baseline),
                  pct_high_edu = mean(education_3cat == "High") * 100,
                  tr20_base = mean(tr20_baseline, na.rm = TRUE)))


# =============================================================================
# SAVE
# =============================================================================

# Save guard. rm(list = ls()) at the head only helps if the head is actually
# executed, so these two lines check that HRS data is going to the HRS file
# when the script is run piecemeal rather than from the top.
stopifnot(identical(unique(hrs_compact$study), "HRS"))
stopifnot(grepl("^hrs_compact", basename(OUT_FILE)))

saveRDS(hrs_compact, OUT_FILE)
cat("\nSaved:", OUT_FILE, "\n")