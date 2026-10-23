# =============================================================================
# JGSS-2026-346 REVISION: SENSITIVITY AND ROBUSTNESS, for the supplement
# =============================================================================
#
# Run after 01_analysis_v6.R, in the same working directory. Reads
# main_v6.RData from OUT_DIR.
#
#
#   A  ascertainment frequency
#   B  complementary log-log link              functional outcomes
#   C  attrition weights                       functional outcomes
#   D  multiple comparisons                    all three outcomes
#   E  checks on the memory mixed model
#        E1 specification, with and without baseline age interacted with time
#        E2 random-effects structure
#        E3 each survey's first three follow-up interviews
#        E4 baseline survey weights
#   F  confirmed events                        not run, kept for reference
#
#   B and C apply to the functional outcomes only. The memory mixed model is
#   fitted by maximum likelihood on every observed score, so it is valid under
#   missing at random without censoring weights, and a complementary log-log
#   check has no counterpart in a model with an identity link.
#
#   The threshold definition of memory decline is estimated alongside ADL and
#   IADL in Part 5 of 01_analysis_v6.R, because it uses the same discrete-time
#   machinery, so it has no part of its own here.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(lme4)

OUT_DIR <- if (exists("OUT_DIR")) OUT_DIR else "output"
if (!exists("LONG")) load(file.path(OUT_DIR, "main_v6.RData"))
stopifnot(exists("LONG"), exists("build_rhs"), exists("hrs_b"))

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# Every object this script writes, declared empty up front, so the save() at
# the foot cannot fail on a missing name if a part bails out.
asc_rates <- data.frame(); asc_grad <- data.frame(); asc_cum <- data.frame()
asc_curve <- data.frame(); cll <- data.frame(); ipcw_fit <- data.frame()
multi <- data.frame()
mem_spec <- data.frame(); mem_ranef <- data.frame()
mem_three <- data.frame(); mem_weight <- data.frame()

# ---- figure constants ------------------------------------------------------

if (!exists("FIG_DIR"))  FIG_DIR  <- OUT_DIR
if (!exists("EDU_COL"))
  EDU_COL <- c(High = "#2166AC", Medium = "#8C8C88", Low = "#B2182B")
if (!exists("EDU_ORD"))  EDU_ORD  <- c("High", "Medium", "Low")
if (!exists("SEX_ORD"))  SEX_ORD  <- c("Male", "Female")   # men on top
if (!exists("OUT_LAB"))
  OUT_LAB <- c(ADL = "ADL", IADL = "IADL", Memory = "Verbal memory decline")
if (!exists("CURVE_LAB"))
  CURVE_LAB <- c(ADL = "ADL", IADL = "IADL",
                 Memory = "Verbal memory\ndecline")
if (!exists("theme_ms")) theme_ms <- function(base = 15) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "#e3e3e0",
                                                 linewidth = .3),
      strip.placement    = "outside",
      strip.text.x       = ggplot2::element_text(face = "bold", size = 17),
      strip.text.y.left  = ggplot2::element_text(face = "bold", size = 17,
                                                 angle = 90),
      axis.text          = ggplot2::element_text(size = 15,
                                                 colour = "#1f2328"),
      axis.title.x       = ggplot2::element_text(size = 16,
                                                 margin = ggplot2::margin(t = 10)),
      axis.title.y       = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.title       = ggplot2::element_blank(),
      legend.text        = ggplot2::element_text(size = 15),
      panel.spacing.x    = grid::unit(1.1, "lines"),
      panel.spacing.y    = grid::unit(1.4, "lines"))
}


# =============================================================================
# SHARED HELPERS
# =============================================================================

# Fit one cell, with the formula built AFTER the sample is fixed.
#
# The order matters. build_rhs() includes partnered_miss only when it varies,
# and asked that question on the unfiltered cell it can answer yes and then
# lose its only non-zero rows to the complete-case restriction the mediators
# impose. The term is then constant, glm returns NA for its coefficient, and
# the anyNA() guard below discards the whole cell. So the sample is restricted
# with a formula naming every candidate variable, and the fitting formula is
# built against the restricted sample, where the variation guard sees the data
# the model will actually use.
#
# cc_type is the model whose complete cases define the sample, which differs
# from model_type only for an M1-on-M3-cases refit. Returns NULL on any
# failure and reports the reason, so that a dropped cell is visible.
cell_glm <- function(d, reg, model_type, cc_type = model_type,
                     family = binomial(), wcol = NULL, label = "") {
  if (!nrow(d)) return(NULL)
  f_cc <- as.formula(paste("event ~", build_rhs(cc_type, reg, d)))
  d <- d[complete.cases(model.frame(f_cc, data = d, na.action = na.pass)), ]
  if (!is.null(wcol)) d <- d[!is.na(d[[wcol]]) & d[[wcol]] > 0, ]
  if (!nrow(d) || sum(d$event) < 10) return(NULL)
  
  f <- as.formula(paste("event ~", build_rhs(model_type, reg, d)))
  m <- if (is.null(wcol)) {
    try(glm(f, data = d, family = family), silent = TRUE)
  } else {
    d$.w <- d[[wcol]] / mean(d[[wcol]])
    try(suppressWarnings(glm(f, data = d, family = family, weights = .w)),
        silent = TRUE)
  }
  if (inherits(m, "try-error")) {
    message("  dropped ", label, " ", reg, ": glm failed"); return(NULL)
  }
  bad <- names(coef(m))[is.na(coef(m))]
  if (!m$converged || length(bad) || max(abs(coef(m)), na.rm = TRUE) > 30) {
    message("  dropped ", label, " ", reg, ": ",
            if (length(bad)) paste("aliased:", paste(bad, collapse = ", "))
            else if (!m$converged) "did not converge" else "coefficients diverged")
    return(NULL)
  }
  list(m = m, d = d)
}

# the education rows of one fitted model, in the usual shape
edu_row <- function(m, d, label, outcome, region, sex_lab, extra = NULL) {
  V   <- try(vcovCL(m, cluster = d$id), silent = TRUE)
  if (inherits(V, "try-error")) V <- vcov(m)
  est <- coeftest(m, vcov = V); ci <- coefci(m, vcov = V)
  k   <- grep("education_3cat", rownames(est))
  if (!length(k)) return(NULL)
  out <- data.frame(
    spec = label, outcome = outcome, region = region, sex = sex_lab,
    term = rownames(est)[k],
    OR = exp(est[k, 1]), CI_low = exp(ci[k, 1]), CI_high = exp(ci[k, 2]),
    p = est[k, 4], n_rows = nrow(d), n_events = sum(d$event),
    row.names = NULL)
  if (!is.null(extra)) out <- cbind(out, extra)
  out
}

# loop the ten region-sex cells of one outcome
fit_all_cells <- function(long, label, outcome, model_type = "M3",
                          family = binomial(), wcol = NULL) {
  out <- data.frame()
  for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
    sex_lab <- ifelse(sx == 1, "Female", "Male")
    fit <- cell_glm(long %>% filter(region == reg, female == sx),
                    reg, model_type, family = family, wcol = wcol,
                    label = paste(outcome, sex_lab))
    if (is.null(fit)) next
    r <- edu_row(fit$m, fit$d, label, outcome, reg, sex_lab)
    if (!is.null(r)) out <- bind_rows(out, r)
  }
  out
}

# the main-analysis estimates, to compare everything against
ref_tab <- function(model_type = "M3") {
  models %>%
    filter(model == model_type, !weighted) %>%
    transmute(outcome, region, sex, term, OR_main = OR, p_main = p,
              CI_main = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high))
}

compare_to_main <- function(x, label, model_type = "M3") {
  if (!nrow(x)) return(invisible(NULL))
  cmp <- x %>% filter(term == "education_3catLow") %>%
    left_join(ref_tab(model_type) %>% filter(term == "education_3catLow"),
              by = c("outcome", "region", "sex", "term")) %>%
    mutate(this = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high),
           ratio = OR / OR_main)
  cat("\n---", label, ": low vs high, against the main", model_type, "---\n")
  cmp %>% select(outcome, region, sex, main = CI_main, this) %>%
    arrange(outcome, sex, match(region, REGIONS_ALL)) %>%
    as_tibble() %>% print(n = 40)
  cat(sprintf("median ratio to main: %.3f   |  significance flips in %d of %d cells\n",
              median(cmp$ratio, na.rm = TRUE),
              sum((cmp$p < .05) != (cmp$p_main < .05), na.rm = TRUE),
              nrow(cmp)))
  invisible(cmp)
}



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
# memory interview. 

MEM_MAX <- 8

# The horizon the memory gap is read at, used by every check in Part E.
TEND <- HORIZON[["Memory"]]

# The reported memory specification, shared by E2, E3 and E4 so the three
# checks cannot drift apart from one another or from Part 7b.
LMM_RHS <- list(
  M1 = "education_3cat * time + age_c + age_c:time",
  M3 = "education_3cat * time + age_c + age_c:time + poor_health + partnered_obs + has_mobility + has_hypertension + has_diabetes + obese + ever_smoker")


mem_gap_at <- function(m, T, lev) {
  b  <- if (inherits(m, "merMod")) lme4::fixef(m) else coef(m)
  b  <- b[!is.na(b)]
  V  <- as.matrix(vcov(m))[names(b), names(b), drop = FALSE]
  a  <- paste0("education_3cat", lev)
  bt <- paste0("education_3cat", lev, ":time")
  if (!a %in% names(b)) return(NULL)
  cc <- setNames(rep(0, length(b)), names(b))
  cc[a] <- 1
  if (bt %in% names(b)) cc[bt] <- T
  est <- sum(cc * b)
  se  <- sqrt(max(0, as.numeric(t(cc) %*% V %*% cc)))
  data.frame(education = lev, time = T, estimate = est, se = se,
             ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
             p = 2 * pnorm(-abs(est / se)), row.names = NULL)
}

# One memory cell's complete cases, shared by E2, E3 and E4.
mem_cell <- function(reg, sx, extra = character(0)) {
  d0 <- mem_cont %>% filter(region == reg, female == sx)
  v  <- c("value", "education_3cat", "time", "age_c", MED_M3, extra)
  if (USE_COUNTRY_FE && reg != "US") v <- c(v, "country")
  v  <- intersect(v, names(d0))
  d0[complete.cases(d0[, v, drop = FALSE]), ]
}

mem_rhs <- function(mt, reg) {
  r <- LMM_RHS[[mt]]
  if (USE_COUNTRY_FE && reg != "US") r <- paste(r, "+ factor(country)")
  r
}

# =============================================================================
# PART A: ASCERTAINMENT FREQUENCY
# =============================================================================
# The two surveys do not look for events equally often. The HRS interviews at
# 2, 4, 6, 8 and 10 years after baseline and SHARE at 2, 4, 7 and 9. Because a
# reported limitation can be absent at the next interview even when the
# underlying limitation is treated as absorbing, the number of occasions on
# which a person is asked affects how many first reports are recorded. Wave
# dummies do not fix this: they let the hazard differ between waves but not the
# number of opportunities to be caught.
#
# This cannot confound the education gradient, since everyone within a country
# shares a schedule, so the gradient panel below is a check rather than a
# result. What it bears on is the cross-national comparison of levels.
#
# THE SPECIFICATION
#
#   survey   schedule   interviews   window
#   ----------------------------------------
#   HRS      2, 4, 6    3            6 years
#   SHARE    2, 4, 7    3            7 years
#
# Each survey is restricted to its first three follow-up interviews. The two
# panels are then observed the same number of times, and Europe has the longer
# observation window, so neither the frequency of observation nor the length of
# follow-up can favour the United States. The same restriction applies to all
# three outcomes, including verbal memory, whose series are shorter.
#
# Results are reported at six years, the last American observation, so the
# cumulative incidence comparison requires no extrapolation of the American
# curve. Note that within those six years the HRS contributes three looks and
# SHARE two, since the European third interview falls at seven years and is
# interpolated back; the quantity matched on both the number of looks and the
# exposure is therefore the incidence rate per person-year, which is the column
# to read for the ascertainment question. The cumulative incidence column is
# reported alongside it and is, if anything, conservative: it still allows the
# United States the extra look inside the horizon.

cat("\n=============== PART A: ascertainment frequency ===============\n")

try({
  HRS_FULL     <- HRS_WAVES
  SHARE_FULL   <- SHARE_WAVES
  MEM_HRS_FULL <- MEM_HRS_WAVES
  
  # each spec says which waves to KEEP, and the horizon to report at
  ASC_SPECS <- list(
    "matched looks" = list(
      hrs = c(12, 13, 14), share = c(6, 7, 8), mem_hrs = c(12, 13, 14),
      horizon = 6)      # the last American observation under this restriction
  )
  
  rate_tab <- function(L, label) {
    bind_rows(lapply(names(L), function(nm)
      person_years(L[[nm]]) %>%
        group_by(region, female, education_3cat) %>%
        summarise(ev = sum(event), py = sum(gap), .groups = "drop") %>%
        mutate(outcome = nm, rate = 1000 * ev / py, spec = label)))
  }
  
  LONG_V <- list()
  for (v in names(ASC_SPECS)) {
    sp <- ASC_SPECS[[v]]
    HRS_WAVES     <- HRS_FULL[HRS_FULL$wave %in% sp$hrs, ]
    SHARE_WAVES   <- SHARE_FULL[SHARE_FULL$wave %in% sp$share, ]
    MEM_HRS_WAVES <- MEM_HRS_FULL[MEM_HRS_FULL$wave %in% sp$mem_hrs, ]
    cat(sprintf("\n%s\n  HRS   %-10s window %gy, %d looks\n  SHARE %-10s window %gy, %d looks\n",
                v, paste(HRS_WAVES$time, collapse = ","), max(HRS_WAVES$time),
                nrow(HRS_WAVES), paste(SHARE_WAVES$time, collapse = ","),
                max(SHARE_WAVES$time), nrow(SHARE_WAVES)))
    LONG_V[[v]] <- list(ADL    = build_long("adl",  phys_base),
                        IADL   = build_long("iadl", phys_base),
                        Memory = build_long("mem",  mem_base))
  }
  HRS_WAVES <- HRS_FULL; SHARE_WAVES <- SHARE_FULL
  MEM_HRS_WAVES <- MEM_HRS_FULL          # restore before anything downstream
  
  # ---- A1. crude rates, the quantity the restriction is built for ------
  asc_rates <- bind_rows(c(list(rate_tab(LONG, "all waves")),
                           lapply(names(ASC_SPECS),
                                  function(v) rate_tab(LONG_V[[v]], v))))
  
  cat("\nPooled incidence per 1,000 person-years and the US / Europe ratio:\n")
  asc_rates %>%
    mutate(side = ifelse(region == "US", "US", "Europe")) %>%
    group_by(outcome, spec, side) %>%
    summarise(rate = 1000 * sum(ev) / sum(py), .groups = "drop") %>%
    pivot_wider(names_from = side, values_from = rate) %>%
    mutate(ratio = round(US / Europe, 2),
           across(c(US, Europe), ~ round(.x, 1))) %>%
    arrange(outcome, match(spec, c("all waves", names(ASC_SPECS)))) %>%
    as_tibble() %>% print(n = 30)
  
  # ---- A2. the education gradient, which should not move ---------------
  cat("\n--- Education gradient, M1, low vs high ---\n")
  ag <- data.frame()
  for (v in names(ASC_SPECS)) {
    for (nm in names(LONG_V[[v]])) {
      r <- fit_all_cells(LONG_V[[v]][[nm]], v, nm, model_type = "M1")
      if (nrow(r)) ag <- bind_rows(ag, r %>% mutate(spec = v))
    }
  }
  if (nrow(ag)) {
    m1ref <- models %>%
      filter(model == "M1", !weighted, term == "education_3catLow") %>%
      transmute(outcome, region, sex, OR_main = OR, p_main = p,
                main = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high))
    asc_grad <- ag %>%
      filter(term == "education_3catLow") %>%
      left_join(m1ref, by = c("outcome", "region", "sex")) %>%
      mutate(matched = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high))
    
    asc_grad %>% select(outcome, region, sex, spec, main, matched) %>%
      pivot_wider(names_from = spec, values_from = matched) %>%
      arrange(outcome, sex, match(region, REGIONS_ALL)) %>%
      as_tibble() %>% print(n = 40)
    for (v in names(ASC_SPECS)) {
      q <- asc_grad %>% filter(spec == v)
      if (!nrow(q)) next
      cat(sprintf("%-18s median ratio to main %.3f  |  significance flips %d of %d\n",
                  v, median(q$OR / q$OR_main, na.rm = TRUE),
                  sum((q$p < .05) != (q$p_main < .05), na.rm = TRUE), nrow(q)))
    }
  }
  
  # ---- A3. the cross-national comparison -------------------------------
  cat("\n--- Predicted cumulative incidence at six years ---\n")
  # Model 3 only, at the observed mediator means. Five regions and three
  # education groups already fill each panel, so a second model would put six
  # ribbons in it; add an entry to this list to draw both.
  ASC_MODELS <- list(c(model = "M3", profile = "average"))
  
  ac <- data.frame(); acv <- data.frame()
  for (v in names(ASC_SPECS)) {
    hz <- ASC_SPECS[[v]]$horizon
    for (sp in ASC_MODELS) {
      mt <- unname(sp["model"]); pf <- unname(sp["profile"])
      pts <- data.frame()
      for (nm in names(LONG_V[[v]])) for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
        r <- try(cum_incidence(LONG_V[[v]][[nm]], reg, sx, nm, mt, pf),
                 silent = TRUE)
        if (inherits(r, "try-error") || is.null(r)) next
        pts <- bind_rows(pts, r$points)
      }
      if (!nrow(pts)) { message("  no cells for ", v, " ", mt); next }
      pts <- as_tibble(pts)
      cat(sprintf("  %s, %s / %s: %d cells\n", v, mt, pf,
                  nrow(distinct(pts, outcome, region, sex, education))))
      
      acv <- bind_rows(acv, pts %>%
                         filter(horizon <= hz, !is.na(cum_inc)) %>%
                         select(outcome, region, sex, education, model, profile, horizon,
                                cum_inc, ci_low, ci_high) %>%
                         mutate(spec = v))
      
      # the main-analysis comparison must come from the SAME model and
      # profile, or the restriction would be confounded with the adjustment
      base <- cumpoint %>%
        filter(model == mt, profile == pf, horizon == hz) %>%
        select(outcome, region, sex, education, cum_main = cum_inc)
      ac <- bind_rows(ac, pts %>%
                        filter(horizon == hz) %>%
                        select(outcome, region, sex, education, model, profile, horizon,
                               cum_inc) %>%
                        left_join(base, by = c("outcome", "region", "sex", "education")) %>%
                        mutate(spec = v))
    }
  }
  asc_cum   <- ac
  asc_curve <- acv
  
  if (nrow(asc_cum)) {
    cat("\nMean US minus European mean, percentage points:\n")
    asc_cum %>%
      pivot_longer(c(cum_inc, cum_main), names_to = "which",
                   values_to = "val") %>%
      mutate(spec = ifelse(which == "cum_main", "all waves", spec)) %>%
      distinct(outcome, sex, education, region, spec, model, horizon,
               .keep_all = TRUE) %>%
      select(-which) %>%
      pivot_wider(names_from = region, values_from = val) %>%
      rowwise() %>%
      mutate(eu = mean(c_across(all_of(setdiff(REGIONS_ALL, "US"))),
                       na.rm = TRUE),
             gap_pp = 100 * (US - eu)) %>% ungroup() %>%
      group_by(outcome, model, horizon, spec) %>%
      summarise(gap = round(mean(gap_pp, na.rm = TRUE), 2), .groups = "drop") %>%
      pivot_wider(names_from = spec, values_from = gap) %>%
      as_tibble() %>% print(n = 20)
  }
  
  write.csv(asc_rates, file.path(OUT_DIR, "S5_A_ascertainment_rates.csv"),
            row.names = FALSE)
  write.csv(asc_grad,  file.path(OUT_DIR, "S5_A_ascertainment_gradient.csv"),
            row.names = FALSE)
  write.csv(asc_cum,   file.path(OUT_DIR, "S5_A_ascertainment_cuminc.csv"),
            row.names = FALSE)
  write.csv(asc_curve, file.path(OUT_DIR, "S5_A_ascertainment_curve.csv"),
            row.names = FALSE)
  
  # ---- A4. Figure S1 ---------------------------------------------------

  if (nrow(asc_curve)) {
    asc_fig <- asc_curve %>%
      filter(spec == "matched looks") %>%
      mutate(across(c(cum_inc, ci_low, ci_high), ~ 100 * .x),
             time      = horizon,
             region    = factor(region, levels = REGIONS_ALL),
             sex       = factor(sex, levels = SEX_ORD),
             education = factor(education, levels = rev(EDU_ORD)),
             rowlab    = factor(CURVE_LAB[outcome], levels = CURVE_LAB),
             modlab    = factor(paste("Model", sub("^M", "", model)),
                                levels = paste("Model",
                                               sub("^M", "",
                                                   vapply(ASC_MODELS, `[[`, "", "model")))))
    
    # cumulative incidence is 0 at baseline, so anchor each curve at the origin
    asc_fig <- asc_fig %>%
      distinct(outcome, region, sex, education, rowlab, modlab) %>%
      mutate(time = 0, cum_inc = 0, ci_low = 0, ci_high = 0) %>%
      bind_rows(asc_fig) %>%
      arrange(outcome, region, sex, education, modlab, time)
    
    # linetype only earns its place when there is more than one model to
    # tell apart; with one it would silently dash the only curve
    two_models <- nlevels(droplevels(asc_fig$modlab)) > 1
    g <- if (two_models) {
      ggplot(asc_fig, aes(time, cum_inc, colour = education, fill = education,
                          linetype = modlab,
                          group = interaction(education, modlab))) +
        scale_linetype_manual(values = c("Model 1" = "solid",
                                         "Model 3" = "22"))
    } else {
      ggplot(asc_fig, aes(time, cum_inc, colour = education, fill = education))
    }
    g <- g +
      geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = .14,
                  colour = NA, show.legend = FALSE) +
      geom_line(linewidth = .85) +
      facet_grid(sex + rowlab ~ region, switch = "y") +
      scale_colour_manual(values = EDU_COL, breaks = EDU_ORD,
                          labels = paste(EDU_ORD, "education")) +
      scale_fill_manual(values = EDU_COL, guide = "none") +
      scale_x_continuous(breaks = seq(0, max(asc_fig$time), by = 2)) +
      labs(x = "Years since baseline", y = "Cumulative incidence (%)") +
      theme_ms() +
      theme(panel.grid.major.y = element_line(colour = "#e3e3e0",
                                              linewidth = .3),
            strip.text.x       = element_text(face = "bold", size = 15),
            strip.text.y.left  = element_text(face = "bold", size = 12.5,
                                              angle = 90),
            axis.title.y       = element_text(size = 15, angle = 90,
                                              margin = margin(r = 10)),
            panel.spacing.y    = unit(.6, "lines"))
    
    ggsave(file.path(OUT_DIR, "FigureS5_cuminc_matched.png"), g,
           width = 14, height = 10, dpi = 300)
    cat(sprintf("\n  FigureS5_cuminc_matched.png written (matched looks, 2/4/6 years, %s)\n",
                paste(vapply(ASC_MODELS, `[[`, "", "model"), collapse = " and ")))
  }
})


# =============================================================================
# PART B: COMPLEMENTARY LOG-LOG LINK
# =============================================================================
# The two surveys space their interviews differently, so the intervals over
# which a hazard is measured are not equal. Under a logit link the coefficients
# depend on interval length; under a complementary log-log link they estimate a
# proportional-hazards model in continuous time and are invariant to it. This
# is the model-based counterpart to Part A, which fixes the same problem by
# design rather than by link function. Coefficients are hazard ratios, not odds
# ratios, so they are expected to be slightly smaller than the main estimates;
# what matters is whether the pattern and the significance hold.

cat("\n=============== PART B: complementary log-log ===============\n")

try({
  for (nm in names(LONG)) {
    r <- fit_all_cells(LONG[[nm]], "complementary log-log", nm,
                       model_type = "M3",
                       family = binomial(link = "cloglog"))
    if (nrow(r)) cll <- bind_rows(cll, r)
    cat("  ", nm, "done\n")
  }
  
  if (nrow(cll)) {
    compare_to_main(cll, "Complementary log-log")
    cat("Values are hazard ratios here and odds ratios in the main models, so\n")
    cat("a ratio slightly below 1 is expected; read the significance pattern.\n")
    write.csv(cll, file.path(OUT_DIR, "S5_B_cloglog.csv"), row.names = FALSE)
  }
})


# =============================================================================
# PART C: ATTRITION
# =============================================================================
# People who leave the study are not a random subset, and the
# main models assume the data are missing at random given the covariates. This
# builds inverse probability of censoring weights from baseline
# characteristics and refits.
#
# At each wave, among those still at risk after the previous wave, the
# probability of being observed is modelled on baseline characteristics. The
# weight is the running product of the inverse of those probabilities,
# stabilised by the marginal probability and truncated at the 1st and 99th
# percentiles so that no single person dominates.

cat("\n=============== PART C: attrition weights ===============\n")

try({
  for (nm in names(LONG)) {
    L <- LONG[[nm]]
    
    # Every wave each person was still AT RISK for, which is every wave of
    # their survey up to and including their first event. Capping at the last
    # wave they were actually seen at would make everyone look fully observed
    # and leave nothing for the weights to correct.
    stop_t <- L %>% group_by(id) %>%
      summarise(stop_t = if (any(event == 1)) min(time[event == 1]) else Inf,
                .groups = "drop")
    
    wave_tab <- L %>% distinct(survey, wave, time) %>% rename(sv = survey)
    
    grid <- L %>% distinct(id, survey) %>%
      left_join(stop_t, by = "id") %>%
      left_join(wave_tab, by = c("survey" = "sv"),
                relationship = "many-to-many") %>%
      filter(time <= stop_t) %>%
      select(id, wave, time)
    
    obs <- grid %>%
      left_join(L %>% transmute(id, wave, seen = 1L), by = c("id", "wave")) %>%
      mutate(seen = ifelse(is.na(seen), 0L, seen)) %>%
      left_join(phys_base %>% select(id, region, female, education_3cat,
                                     age_baseline, all_of(MED_M3),
                                     partnered_miss),
                by = "id")
    
    obs <- obs[complete.cases(obs %>%
                                select(seen, region, female, education_3cat, age_baseline,
                                       all_of(MED_M3))), ]
    if (!nrow(obs)) next
    
    mc <- try(glm(seen ~ education_3cat + age_baseline + factor(wave) + region +
                    poor_health + partnered_obs + has_mobility +
                    has_hypertension + has_diabetes + obese + ever_smoker,
                  data = obs, family = binomial), silent = TRUE)
    mn <- try(glm(seen ~ factor(wave), data = obs, family = binomial),
              silent = TRUE)
    if (inherits(mc, "try-error") || inherits(mn, "try-error")) next
    
    obs$p_cond <- fitted(mc); obs$p_marg <- fitted(mn)
    w <- obs %>% arrange(id, time) %>% group_by(id) %>%
      mutate(ipcw = cumprod(p_marg / p_cond)) %>% ungroup() %>%
      select(id, wave, ipcw)
    
    lim <- quantile(w$ipcw, c(.01, .99), na.rm = TRUE)
    w$ipcw <- pmin(pmax(w$ipcw, lim[1]), lim[2])
    cat(sprintf("  %-7s IPCW range %.2f to %.2f, mean %.2f\n", nm,
                min(w$ipcw), max(w$ipcw), mean(w$ipcw)))
    
    Lw <- L %>% left_join(w, by = c("id", "wave"))
    r <- fit_all_cells(Lw, "IPCW", nm, model_type = "M3", wcol = "ipcw")
    if (nrow(r)) ipcw_fit <- bind_rows(ipcw_fit, r)
    cat("  ", nm, "done\n")
  }
  
  if (nrow(ipcw_fit)) {
    compare_to_main(ipcw_fit, "IPCW for attrition")
    write.csv(ipcw_fit, file.path(OUT_DIR, "S5_C_ipcw.csv"), row.names = FALSE)
  }
})


# =============================================================================
# PART D: MULTIPLE COMPARISONS
# =============================================================================
# The education contrast is estimated in thirty region-by-sex-by-outcome cells,
# so some significant findings would be expected by chance alone. This applies
# the Benjamini-Hochberg false discovery rate correction WITHIN each outcome,
# across that outcome's ten region-by-sex strata.
#
# Benjamini-Hochberg depends on the whole ranked list, so a family spanning
# the three outcomes makes the adjusted values for ADL and IADL depend on
# which model is used for memory. The three outcomes are measured on different
# scales and address distinct hypotheses, so they are corrected separately.
# The across-thirty column is still computed; the Results quote the
# within-outcome one.

cat("\n=============== PART D: multiple comparisons ===============\n")

try({
  phys <- models %>%
    filter(!weighted, term == "education_3catLow",
           model %in% c("M1", "M3"), outcome %in% c("ADL", "IADL")) %>%
    select(model, outcome, region, sex, estimate = OR, p)
  
  mem <- read.csv(file.path(OUT_DIR, "S5_E2_memory_gaps.csv")) %>%
    filter(spec == "age x time", education == "Low",
           time == HORIZON[["Memory"]], model %in% c("M1", "M3")) %>%
    mutate(outcome = "Memory") %>%
    select(model, outcome, region, sex, estimate, p)
  
  multi <- bind_rows(phys, mem) %>%
    group_by(model) %>%
    mutate(p_bh_all = p.adjust(p, "BH")) %>%
    group_by(model, outcome) %>%
    mutate(p_bh_outcome = p.adjust(p, "BH")) %>%
    ungroup()
  
  cat("\nCells significant at 0.05, low vs high education:\n")
  multi %>% group_by(model, outcome) %>%
    summarise(cells = n(),
              raw = sum(p < .05),
              bh_within_outcome = sum(p_bh_outcome < .05),
              bh_all_30 = sum(p_bh_all < .05),
              .groups = "drop") %>%
    as_tibble() %>% print(n = 20)
  
  cat("\nHow close the within-outcome correction comes to changing anything:\n")
  multi %>% group_by(model, outcome) %>%
    summarise(worst_survivor = suppressWarnings(
      max(p_bh_outcome[p_bh_outcome < .05], na.rm = TRUE)),
      best_null = suppressWarnings(
        min(p_bh_outcome[p_bh_outcome >= .05], na.rm = TRUE)),
      flips = sum((p < .05) != (p_bh_outcome < .05)),
      .groups = "drop") %>%
    as_tibble() %>% print(n = 20)
  
  cat("\nQuote p_bh_outcome in the Results. The memory rows are the mixed\n")
  cat("model's gap at six years, not the threshold outcome.\n")
  
  write.csv(multi, file.path(OUT_DIR, "S5_D_multiplicity.csv"),
            row.names = FALSE)
})


# =============================================================================
# PART E: CHECKS ON THE MEMORY MIXED MODEL
# =============================================================================
# E1  specification, with and without baseline age interacted with time
# E2  random-effects structure, with and without the random slope
# E3  each survey restricted to its first three follow-up interviews
# E4  baseline survey weights
#
# All four read the reported specification from LMM_RHS above, so none of them
# can drift from Part 7b of the analysis script or from each other.


# ---- E1: specification ------------------------------------------------------
# No fitting required: Part 7b writes both specifications, so the comparison is
# a read and a join.

cat("\n=============== PART E1: specification ===============\n")

try({
  mem_spec <- read.csv(file.path(OUT_DIR, "S5_E2_memory_gaps.csv")) %>%
    filter(education == "Low", time == TEND) %>%
    select(spec, model, region, sex, estimate, se) %>%
    tidyr::pivot_wider(names_from = spec, values_from = c(estimate, se)) %>%
    mutate(difference = `estimate_age x time` - `estimate_as published`)
  print(as.data.frame(mem_spec %>% filter(model == "M3")),
        row.names = FALSE, digits = 3)
  cat("\nlargest absolute difference:",
      round(max(abs(mem_spec$difference), na.rm = TRUE), 3), "words\n")
  write.csv(mem_spec, file.path(OUT_DIR, "S5_E_memory_spec.csv"),
            row.names = FALSE)
})


# ---- E2: random-effects structure -------------------------------------------
# Refits with a random intercept only. This also settles the singular fits: a
# slope variance estimated at or near zero is what produces them, and they
# disappear without the slope.

cat("\n=============== PART E2: random effects ===============\n")

try({
  RAN <- c("intercept and slope" = "(1 + time | id)",
           "intercept only"      = "(1 | id)")
  co <- data.frame()
  for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
    sex_lab <- ifelse(sx == 1, "Female", "Male")
    d <- mem_cell(reg, sx)
    if (nrow(d) < 200) next
    
    for (rn in names(RAN)) for (mt in names(LMM_RHS)) {
      m <- try(lmer(as.formula(paste("value ~", mem_rhs(mt, reg), "+",
                                     RAN[[rn]])),
                    data = d, REML = FALSE,
                    control = lmerControl(optimizer = "bobyqa")),
               silent = TRUE)
      if (inherits(m, "try-error")) {
        message("    dropped ", reg, " ", sex_lab, " ", mt, " [", rn, "]")
        next
      }
      fe <- fixef(m); se <- sqrt(diag(as.matrix(vcov(m))))
      co <- bind_rows(co, data.frame(
        ranef = rn, model = mt, region = reg, sex = sex_lab,
        term = names(fe), estimate = as.vector(fe), se = as.vector(se),
        p = 2 * pnorm(-abs(fe / se)), singular = isSingular(m),
        n_obs = nrow(d), n_people = length(unique(d$id)), row.names = NULL))
    }
  }
  mem_ranef <- co
  write.csv(mem_ranef, file.path(OUT_DIR, "S5_E_memory_ranef.csv"),
            row.names = FALSE)
  
  if (nrow(mem_ranef)) {
    g <- mem_ranef %>%
      filter(model == "M3",
             term %in% c("education_3catLow", "education_3catLow:time")) %>%
      mutate(part = ifelse(grepl(":time", term), "slope", "level")) %>%
      select(ranef, region, sex, part, estimate) %>%
      tidyr::pivot_wider(names_from = part, values_from = estimate) %>%
      mutate(gap = level + TEND * slope) %>%
      select(region, sex, ranef, gap) %>%
      tidyr::pivot_wider(names_from = ranef, values_from = gap) %>%
      mutate(difference = `intercept only` - `intercept and slope`)
    print(as.data.frame(g), row.names = FALSE, digits = 3)
    cat("\nlargest absolute difference in the gap at", TEND, "years:",
        round(max(abs(g$difference), na.rm = TRUE), 3), "words\n")
    ints <- mem_ranef %>% filter(term == "(Intercept)")
    cat("singular fits, random slope:",
        sum(ints$singular[ints$ranef == "intercept and slope"]),
        "| intercept only:",
        sum(ints$singular[ints$ranef == "intercept only"]), "\n")
  }
})


# ---- E3: each survey's first three follow-up interviews ---------------------
# The counterpart for the memory model of the restricted-follow-up check on the
# functional outcomes. HRS observes recall at 0, 2, 4, 6 and 8 years and SHARE
# at 0, 2, 4 and 7, so SHARE already contributes three follow-ups and the
# restriction removes a wave from the HRS alone; the European estimates are
# unchanged by construction.

cat("\n=============== PART E3: first three interviews ===============\n")

try({
  three <- data.frame()
  for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
    sex_lab <- ifelse(sx == 1, "Female", "Male")
    d <- mem_cell(reg, sx)
    if (nrow(d) < 200) next
    keep_to <- if (reg == "US") 6 else 7     # the third follow-up interview
    d <- d %>% filter(time <= keep_to)
    cat(sprintf("  %-14s %-6s %d obs on %d people\n", reg, sex_lab,
                nrow(d), length(unique(d$id))))
    
    for (mt in names(LMM_RHS)) {
      m <- try(lmer(as.formula(paste("value ~", mem_rhs(mt, reg),
                                     "+ (1 + time | id)")),
                    data = d, REML = FALSE,
                    control = lmerControl(optimizer = "bobyqa")), silent = TRUE)
      if (inherits(m, "try-error")) {
        message("    dropped ", reg, " ", sex_lab, " ", mt); next
      }
      for (lev in setdiff(levels(droplevels(d$education_3cat)), "High")) {
        g <- mem_gap_at(m, TEND, lev)
        if (!is.null(g))
          three <- bind_rows(three, cbind(
            data.frame(analysis = "first three interviews", model = mt,
                       region = reg, sex = sex_lab, n_obs = nrow(d),
                       n_people = length(unique(d$id))), g))
      }
    }
  }
  mem_three <- three
  write.csv(mem_three, file.path(OUT_DIR, "S5_E_memory_three.csv"),
            row.names = FALSE)
})


# ---- E4: baseline survey weights --------------------------------------------
# A weighted linear model with standard errors clustered on the person is the
# same estimator as a weighted GEE with an independence working correlation,
# and it solves in closed form. Weights are normalised to mean one within each
# stratum, which leaves the point estimates unchanged.
#
# Each model is fitted unweighted and weighted, so the comparison isolates the
# weights rather than the change of estimator: comparing a weighted linear
# model against the unweighted mixed model would change both at once.

cat("\n=============== PART E4: survey weights ===============\n")

W_VAR <- "weight"

try({
  stopifnot(W_VAR %in% names(mem_cont))
  out <- data.frame()
  for (reg in REGIONS_ALL) for (sx in c(0, 1)) {
    sex_lab <- ifelse(sx == 1, "Female", "Male")
    d <- mem_cell(reg, sx, extra = W_VAR)
    d <- d[!is.na(d[[W_VAR]]) & d[[W_VAR]] > 0, ]
    if (nrow(d) < 200) next
    d$.w <- d[[W_VAR]] / mean(d[[W_VAR]])
    cat(sprintf("  %-14s %-6s %d obs on %d people\n", reg, sex_lab,
                nrow(d), length(unique(d$id))))
    
    for (mt in names(LMM_RHS)) {
      f <- as.formula(paste("value ~", mem_rhs(mt, reg)))
      for (nm in c("unweighted", "weighted")) {
        m <- if (nm == "weighted") lm(f, data = d, weights = d$.w)
        else                  lm(f, data = d)
        b <- coef(m); b <- b[!is.na(b)]
        V <- vcovCL(m, cluster = d$id)[names(b), names(b), drop = FALSE]
        for (lev in setdiff(levels(droplevels(d$education_3cat)), "High")) {
          a  <- paste0("education_3cat", lev)
          bt <- paste0(a, ":time")
          if (!a %in% names(b)) next
          cc <- setNames(rep(0, length(b)), names(b))
          cc[a] <- 1
          if (bt %in% names(b)) cc[bt] <- TEND
          est <- sum(cc * b)
          se  <- sqrt(max(0, as.numeric(t(cc) %*% V %*% cc)))
          out <- bind_rows(out, data.frame(
            analysis = paste("linear model,", nm), model = mt,
            region = reg, sex = sex_lab, n_obs = nrow(d),
            n_people = length(unique(d$id)), education = lev, time = TEND,
            estimate = est, se = se, ci_low = est - 1.96 * se,
            ci_high = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)),
            row.names = NULL))
        }
      }
    }
  }
  mem_weight <- out
  write.csv(mem_weight, file.path(OUT_DIR, "S5_E_memory_weighted.csv"),
            row.names = FALSE)
  
  if (nrow(mem_weight)) {
    cmp <- mem_weight %>%
      filter(education == "Low", model == "M3") %>%
      select(region, sex, analysis, estimate) %>%
      tidyr::pivot_wider(names_from = analysis, values_from = estimate) %>%
      mutate(difference = `linear model, weighted` -
               `linear model, unweighted`)
    print(as.data.frame(cmp), row.names = FALSE, digits = 3)
    cat("\nlargest absolute difference:",
        round(max(abs(cmp$difference), na.rm = TRUE), 3), "words\n")
  }
})



# =============================================================================
# SAVE


save(asc_rates, asc_grad, asc_cum, asc_curve,
     cll, ipcw_fit, multi,
     mem_spec, mem_ranef, mem_three, mem_weight,
     LMM_RHS, TEND, MEM_MAX,
     file = file.path(OUT_DIR, "sens_v6.RData"))

cat("\n=== Sensitivity analyses done ===\n")
for (o in list(c("A ascertainment",      nrow(asc_cum)),
               c("B complementary l-l",  nrow(cll)),
               c("C IPCW",               nrow(ipcw_fit)),
               c("D multiplicity",       nrow(multi)),
               c("E1 memory spec",       nrow(mem_spec)),
               c("E2 memory ranef",      nrow(mem_ranef)),
               c("E3 three interviews",  nrow(mem_three)),
               c("E4 survey weights",    nrow(mem_weight))))
  cat(sprintf("  %-22s %s\n", o[1],
              ifelse(as.numeric(o[2]) > 0, paste(o[2], "rows"),
                     "EMPTY, see above")))