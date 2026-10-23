# =============================================================================
# SHARE COMPACT DATASET - WAVE 5 BASELINE (2013) THROUGH WAVE 9 (2022)
# =============================================================================
#
# Reads the Gateway to Global Aging harmonised SHARE file and writes the
# compact person-level file the analysis scripts use: one row per respondent
# interviewed at wave 5, aged 50 to 100, and free of ADL and IADL limitations
# at that wave.
#
# WHAT IS KEPT, AND WHY
#
#  1. Immediate and delayed recall at every wave, not only at baseline and
#     endpoint, so the separate-components check can run on person-wave data.
#
#  2. ADL component items, for the six-item sensitivity analysis. r5adla was
#     verified against the data to contain bathing, dressing, eating, bed
#     transfer and walking across a room, matching HRS r11adl5a exactly.
#     Toileting is in neither index.
#
#  3. The five-item IADL index, built here by hand rather than taken from the
#     harmonised file. See the note below; this is not optional.
#
#  Endpoint-only variables are not kept. The analysis is discrete-time on
#  person-wave data and does not use them.
#
# NOTE ON IADL  --  READ THIS BEFORE CHANGING ANYTHING HERE
#   The pre-computed summary indices in the two surveys do not cover the same
#   items, and the mismatch is large enough to invalidate any cross-national
#   comparison of IADL levels:
#
#     Gateway SHARE  r{w}iadla    3 items: money, meds, phone
#     RAND HRS       r{w}iadl5a   5 items: money, meds, shopping, meals, phone
#
#   Confirmed against the data. r5iadla is reproduced exactly, including its
#   missingness pattern, by r5moneya + r5medsa + r5phonea, and by no other
#   triple. The two items it omits, shopping for groceries and preparing
#   meals, are the two most commonly endorsed of the five in SHARE wave 5
#   (7.1% and 4.9%, against 4.8%, 2.6% and 2.5% for the three it keeps). Using
#   the two summary indices as they come would measure American difficulty on
#   five items against European difficulty on the three rarest of those five,
#   and European IADL incidence would come out at roughly half European ADL
#   incidence, which is backwards: IADL tasks are more demanding and should
#   onset earlier.
#
#   All five component items are present in SHARE, so the index is built here
#   from the five named items with one explicit rule: plain addition, and a
#   wave is missing if any item is unobserved. That is what r5iadla already
#   does, so SHARE's behaviour is unchanged and the HRS script matches it.
#   RAND's index, by contrast, is filled when some items are unobserved, which
#   silently counts an unobserved item as no difficulty.
#
#   r5iadla is retained as iadl_gateway_baseline so the discrepancy can be
#   documented rather than just asserted.
#
# NOTE ON MOBILITY
#   r5mobilb is a 0-4 index; the HRS equivalent r11mobila is 0-5. The item
#   sets differ: SHARE has walking 100 metres, HRS has walking several blocks
#   and walking one block. It is used only as a binary covariate and is never
#   compared across surveys.
#
# NOTE ON WAVE SPACING
#   Waves 5, 6 and 7 are two years apart, 7 to 8 is three years (2017 to 2020,
#   the COVID telephone wave), and 8 to 9 is two. The analysis handles this
#   with wave dummies and a complementary log-log sensitivity check.
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
# Put the raw SHARE file in data/raw/, or copy config_example.R to config.R
# and point IN_FILE somewhere else. config.R is untracked and is sourced here
# if it exists. The working directory is never changed from inside the script:
# changing it to a missing directory fails quietly, and the output then lands
# somewhere other than the path the console reports.

RAW_DIR     <- "data/raw"
DERIVED_DIR <- "data/derived"

IN_FILE  <- file.path(RAW_DIR, "GH_SHARE_g.dta")
OUT_FILE <- file.path(DERIVED_DIR, "share_compact.rds")

if (file.exists("config.R")) source("config.R")
dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(IN_FILE))


# =============================================================================
# LOAD
# =============================================================================

cat("Loading Gateway Harmonized SHARE...\n")
GH_SHARE_g <- read_dta(IN_FILE)
cat("  N =", nrow(GH_SHARE_g), " Variables =", ncol(GH_SHARE_g), "\n\n")

cat("=== WAVE 5 COUNTRIES ===\n")
print(GH_SHARE_g %>% filter(r5iwstat == 1) %>% count(country) %>% arrange(country), n = 30)


# =============================================================================
# BUILD THE FIVE-ITEM IADL INDEX
# =============================================================================
# Done here, on the raw file, because the baseline exclusion below has to use
# it. If this ran after the filter, the sample would still be selected on the
# three-item index and people already limited on shopping or meals would enter
# the risk set and then count as incident cases at the first follow-up, which
# would replace one artefact with another pointing the other way.

IADL_ITEMS <- c("moneya", "medsa", "shopa", "mealsa", "phonea")

for (w in 5:9) {
  it <- paste0("r", w, IADL_ITEMS)
  stopifnot(all(it %in% names(GH_SHARE_g)))   # stop rather than build a short index
  GH_SHARE_g[[paste0("iadl5_w", w)]] <-
    GH_SHARE_g[[it[1]]] + GH_SHARE_g[[it[2]]] + GH_SHARE_g[[it[3]]] +
    GH_SHARE_g[[it[4]]] + GH_SHARE_g[[it[5]]]
}

cat("\n=== FIVE-ITEM IADL INDEX, RAW FILE ===\n")
for (w in 5:9) {
  v <- GH_SHARE_g[[paste0("iadl5_w", w)]]
  cat(sprintf("  iadl5_w%-3d max = %s  mean = %.3f  any = %.1f%%  NA = %.1f%%\n",
              w, max(v, na.rm = TRUE), mean(v, na.rm = TRUE),
              100 * mean(v >= 1, na.rm = TRUE), 100 * mean(is.na(v))))
}

# The three-item index must still be reproducible from its own components.
# If this prints FALSE, Gateway has changed the index and everything above
# needs rechecking.
tri <- as.numeric(GH_SHARE_g$r5moneya + GH_SHARE_g$r5medsa + GH_SHARE_g$r5phonea)
cat(sprintf("\n  r5iadla is money + meds + phone: %s\n",
            isTRUE(all.equal(tri, as.numeric(GH_SHARE_g$r5iadla)))))
cat(sprintf("  three-item any-difficulty %.1f%%  vs five-item %.1f%% at wave 5\n",
            100 * mean(GH_SHARE_g$r5iadla >= 1, na.rm = TRUE),
            100 * mean(GH_SHARE_g$iadl5_w5 >= 1, na.rm = TRUE)))
cat("  The gap between those two numbers is the whole problem.\n\n")


# =============================================================================
# BUILD
# =============================================================================
# Country codes: 11 Austria, 12 Germany, 13 Sweden, 14 Netherlands,
# 15 Spain, 16 Italy, 17 France, 18 Denmark, 20 Switzerland, 23 Belgium,
# 25 Israel (excluded), 28 Czechia, 31 Luxembourg, 34 Slovenia, 35 Estonia

cat("=== BASELINE (wave 5, 2013) ===\n")
cat("  free of ADL and IADL, Gateway 3-item index: ",
    GH_SHARE_g %>% filter(r5iwstat == 1, r5agey >= 50, r5agey <= 100,
                          country != 25, r5adla == 0, r5iadla == 0) %>% nrow(), "\n")
cat("  free of ADL and IADL, hand-built 5-item:    ",
    GH_SHARE_g %>% filter(r5iwstat == 1, r5agey >= 50, r5agey <= 100,
                          country != 25, r5adla == 0, iadl5_w5 == 0) %>% nrow(), "\n")
cat("  (the second is the sample used; it is smaller because people already\n")
cat("   limited on shopping or meals are now correctly excluded)\n\n")

share_compact <- GH_SHARE_g %>%
  filter(r5iwstat == 1) %>%
  filter(r5agey >= 50 & r5agey <= 100) %>%
  filter(r5adla == 0 & iadl5_w5 == 0) %>%
  filter(country != 25) %>%
  
  select(
    # --- ids ----------------------------------------------------------------
    mergeid, country,
    
    # --- weights -------------------------------------------------------------
    wtresp_w5 = r5wtresp, wthh_w5 = r5wthh,
    
    # --- time-invariant demographics -----------------------------------------
    birth_year = rabyear, gender = ragender,
    education_isced = raedisced, education_years = raedyrs,
    
    # --- baseline demographics ------------------------------------------------
    age_baseline = r5agey, mstat_baseline = r5mstat, partnered_baseline = r5mpart,
    
    # --- mortality -------------------------------------------------------------
    death_year = radyear, death_month = radmonth,
    
    # --- interview status, every wave -------------------------------------------
    iwstat_w5 = r5iwstat, iwstat_w6 = r6iwstat, iwstat_w7 = r7iwstat,
    iwstat_w8 = r8iwstat, iwstat_w9 = r9iwstat,
    
    # --- baseline health ---------------------------------------------------------
    srh_baseline = r5shlt,
    hibp_baseline = r5hibpe, diab_baseline = r5diabe, cancer_baseline = r5cancre,
    lung_baseline = r5lunge, heart_baseline = r5hearte, stroke_baseline = r5stroke,
    arthritis_baseline = r5arthre, alzdem_baseline = r5alzdeme,
    depression_baseline = r5eurod,
    smoking_ever = r5smokev, smoking_now = r5smoken,
    bmi_baseline = r5bmi, obese_baseline = r5obese,
    
    # --- ADL: Gateway's summary index, verified to match HRS adl5a ---------
    adl_baseline = r5adla,
    adl_w6 = r6adla, adl_w7 = r7adla, adl_w8 = r8adla, adl_w9 = r9adla,
    
    # --- IADL: the hand-built five-item index (see the header note) --------
    iadl_baseline = iadl5_w5,
    iadl_w6 = iadl5_w6, iadl_w7 = iadl5_w7, iadl_w8 = iadl5_w8,
    iadl_w9 = iadl5_w9,
    
    # Gateway's own three-item index, kept only to document the discrepancy.
    # NOT analysed.
    iadl_gateway_baseline = r5iadla,
    
    # --- ADL component items at baseline, for the 6-item sensitivity -------------
    bath_baseline = r5batha, dress_baseline = r5dressa, eat_baseline = r5eata,
    bed_baseline = r5beda, toilet_baseline = r5toilta, walkr_baseline = r5walkra,
    
    # --- IADL component items at baseline, for item-level checks -----------
    money_baseline = r5moneya, meds_baseline = r5medsa,
    shop_baseline = r5shopa, meals_baseline = r5mealsa,
    phone_baseline = r5phonea, map_baseline = r5mapa,
    
    # --- mobility (covariate only; see header note) --------------------------------
    mobility_baseline = r5mobilb,
    
    # --- cognition, baseline and every follow-up wave ------------------------------
    imrc_baseline = r5imrc, dlrc_baseline = r5dlrc, tr20_baseline = r5tr20,
    ser7_baseline = r5ser7, orient_baseline = r5orient, verbf_baseline = r5verbf,
    numeracy_baseline = r5numer_s,
    
    imrc_w6 = r6imrc, dlrc_w6 = r6dlrc, tr20_w6 = r6tr20,
    imrc_w7 = r7imrc, dlrc_w7 = r7dlrc, tr20_w7 = r7tr20,
    imrc_w8 = r8imrc, dlrc_w8 = r8dlrc, tr20_w8 = r8tr20,
    imrc_w9 = r9imrc, dlrc_w9 = r9dlrc, tr20_w9 = r9tr20,
    
    ser7_w8 = r8ser7, ser7_w9 = r9ser7,
    
    # --- self-rated health at follow-up ---------------------------------------------
    srh_w6 = r6shlt, srh_w7 = r7shlt, srh_w8 = r8shlt, srh_w9 = r9shlt,
    
    # --- wealth -----------------------------------------------------------------------
    wealth_total = h5atoth, wealth_net = h5atotb, income_hh = h5ittot
  ) %>%
  
  mutate(
    # --- metadata ---------------------------------------------------------------------
    study = "SHARE", baseline_wave = 5, baseline_year = 2013,
    
    region = case_when(
      country %in% c(13, 18)                  ~ "Nordic",
      country %in% c(15, 16, 17)              ~ "Mediterranean",
      country %in% c(28, 34, 35)              ~ "Eastern",
      country %in% c(11, 12, 14, 20, 23, 31)  ~ "Western",
      TRUE                                    ~ "Other"),
    
    country_name = case_when(
      country == 11 ~ "Austria",        country == 12 ~ "Germany",
      country == 13 ~ "Sweden",         country == 14 ~ "Netherlands",
      country == 15 ~ "Spain",          country == 16 ~ "Italy",
      country == 17 ~ "France",         country == 18 ~ "Denmark",
      country == 20 ~ "Switzerland",    country == 23 ~ "Belgium",
      country == 28 ~ "Czech Republic", country == 31 ~ "Luxembourg",
      country == 34 ~ "Slovenia",       country == 35 ~ "Estonia",
      TRUE          ~ "Unknown"),
    
    # --- demographics -------------------------------------------------------------------
    female = ifelse(gender == 2, 1, 0),
    
    # Native ISCED. HRS maps years onto the same bands; see the HRS script.
    education_3cat = case_when(
      education_isced %in% c(0, 1, 2) ~ "Low",
      education_isced %in% c(3, 4)    ~ "Medium",
      education_isced %in% c(5, 6)    ~ "High",
      TRUE                            ~ NA_character_),
    
    # Years-based version, for a sensitivity against the ISCED crosswalk
    education_3cat_yrs = case_when(
      education_years < 12                          ~ "Low",
      education_years >= 12 & education_years < 16  ~ "Medium",
      education_years >= 16                         ~ "High",
      TRUE                                          ~ NA_character_),
    
    partnered = case_when(
      mstat_baseline %in% c(1, 2, 3) | partnered_baseline == 1 ~ 1,
      mstat_baseline %in% c(4, 5, 6, 7, 8)                     ~ 0,
      TRUE                                                     ~ NA_real_),
    
    # --- health covariates ----------------------------------------------------------------
    poor_health      = ifelse(srh_baseline %in% c(4, 5), 1, 0),
    has_hypertension = ifelse(hibp_baseline == 1, 1, 0),
    has_diabetes     = ifelse(diab_baseline == 1, 1, 0),
    has_mobility     = ifelse(mobility_baseline >= 1, 1, 0),
    obese            = ifelse(bmi_baseline >= 30, 1, 0),
    ever_smoker      = ifelse(smoking_ever == 1, 1, 0),
    current_smoker   = ifelse(smoking_now == 1, 1, 0),
    depressed        = ifelse(depression_baseline >= 4, 1, 0),
    
    chronic_count = rowSums(cbind(
      hibp_baseline == 1, diab_baseline == 1, cancer_baseline == 1,
      lung_baseline == 1, heart_baseline == 1, stroke_baseline == 1,
      arthritis_baseline == 1), na.rm = TRUE),
    
    # --- six-item ADL at baseline (adds toileting) ------------------------------------------
    adl6_baseline = rowSums(cbind(bath_baseline, dress_baseline, eat_baseline,
                                  bed_baseline, toilet_baseline, walkr_baseline)),
    
    # --- mortality and observation ------------------------------------------------------------
    died_before_w9 = ifelse(!is.na(death_year) & death_year <= 2021, 1, 0),
    observed_w8 = iwstat_w8 == 1,
    observed_w9 = iwstat_w9 == 1,
    
    # --- wealth ---------------------------------------------------------------------------------
    wealth_quartile = ntile(wealth_total, 4)
  ) %>%
  
  filter(region != "Other") %>%
  filter(!is.na(education_3cat))

cat("\nBuilt: N =", nrow(share_compact), " Variables =", ncol(share_compact), "\n\n")


# =============================================================================
# CHECKS
# =============================================================================

cat("=== SAMPLE BY REGION ===\n")
print(share_compact %>% count(region) %>% mutate(pct = round(100 * n / sum(n), 1)))

cat("\n=== COUNTRIES ===\n")
print(share_compact %>% group_by(region) %>%
        summarise(countries = paste(sort(unique(country_name)), collapse = ", "),
                  n = n(), .groups = "drop"))

cat("\n=== DEMOGRAPHICS ===\n")
cat(sprintf("  Mean age: %.1f (SD %.1f)\n",
            mean(share_compact$age_baseline), sd(share_compact$age_baseline)))
cat(sprintf("  Female: %.1f%%\n", mean(share_compact$female) * 100))

cat("\n=== EDUCATION: ISCED vs YEARS ===\n")
cat("(if these diverge a lot, the crosswalk sensitivity is worth running)\n\n")
print(table(share_compact$education_3cat, share_compact$education_3cat_yrs,
            useNA = "ifany", dnn = c("ISCED", "Years")))

cat("\n=== SIX-ITEM ADL AT BASELINE ===\n")
cat("(everyone is 0 on the five-item index; this shows toileting-only cases)\n")
print(table(share_compact$adl6_baseline, useNA = "ifany"))

cat("\n=== COGNITION COVERAGE BY WAVE ===\n")
for (v in c("tr20_baseline", "tr20_w6", "tr20_w7", "tr20_w8", "tr20_w9")) {
  cat(sprintf("  %-14s n = %5d  mean = %.2f\n", v,
              sum(!is.na(share_compact[[v]])),
              mean(share_compact[[v]], na.rm = TRUE)))
}

# Covers ADL and IADL, and prints the maximum of each index. The maximum is
# the number that reveals an item-set mismatch between the two surveys.
cat("\n=== ADL/IADL COVERAGE BY WAVE ===\n")
for (v in c("adl_baseline", "adl_w6", "adl_w7", "adl_w8", "adl_w9",
            "iadl_baseline", "iadl_w6", "iadl_w7", "iadl_w8", "iadl_w9")) {
  x <- share_compact[[v]]
  cat(sprintf("  %-14s n = %5d  max = %s  any limitation = %.1f%%\n", v,
              sum(!is.na(x)), max(x, na.rm = TRUE),
              100 * mean(x >= 1, na.rm = TRUE)))
}
cat("\n  Every max must be 5, here and in the HRS script. Anything else\n")
cat("  means the item sets have drifted apart again.\n")

cat("\n=== IADL COMPONENT ITEMS AT BASELINE ===\n")
cat("(all zero by construction; shown so the item set is on the record)\n")
for (v in c("money_baseline", "meds_baseline", "shop_baseline",
            "meals_baseline", "phone_baseline", "map_baseline")) {
  x <- share_compact[[v]]
  cat(sprintf("  %-16s n = %5d  mean = %.4f\n", v, sum(!is.na(x)),
              mean(x, na.rm = TRUE)))
}

cat("\n=== FOLLOW-UP ===\n")
cat(sprintf("  Observed w8: %.1f%%\n", 100 * mean(share_compact$observed_w8)))
cat(sprintf("  Observed w9: %.1f%%\n", 100 * mean(share_compact$observed_w9)))
cat(sprintf("  Died by 2021: %.1f%%\n", 100 * mean(share_compact$died_before_w9)))

cat("\n=== BASELINE MEMORY BY REGION AND EDUCATION ===\n")
print(share_compact %>%
        group_by(region, education_3cat) %>%
        summarise(mean_tr20 = round(mean(tr20_baseline, na.rm = TRUE), 2),
                  .groups = "drop") %>%
        pivot_wider(names_from = education_3cat, values_from = mean_tr20))

cat("\n=== PARTNERSHIP SANITY CHECK ===\n")
print(table(share_compact$mstat_baseline, share_compact$partnered,
            useNA = "ifany"))
cat(sprintf("  missing: %.2f%% (expect about 0.10)\n",
            100 * mean(is.na(share_compact$partnered))))
cat(sprintf("  partnered: %.1f%% (expect about 77)\n",
            100 * mean(share_compact$partnered, na.rm = TRUE)))


# =============================================================================
# SAVE
# =============================================================================

# Save guard. rm(list = ls()) at the head only helps if the head is actually
# executed, so this line checks that SHARE data is going to the SHARE file
# when the script is run piecemeal rather than from the top.
stopifnot(identical(unique(share_compact$study), "SHARE"))

saveRDS(share_compact, OUT_FILE)
cat("\nSaved:", OUT_FILE, "\n")
cat("Written at:", format(file.info(OUT_FILE)$mtime), "\n")