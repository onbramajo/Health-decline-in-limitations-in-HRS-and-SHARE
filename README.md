# Health-decline-in-limitations-in-HRS-and-SHARE
Code to replicate study of Callaway Bramajo and Strozza on Health Decline in Functional Limitations in US and Europe 

Educational inequalities in functional limitation and verbal memory in the United States and Europe

Analysis code for the manuscript Educational inequalities in the onset of functional limitations and verbal memory decline in the United States and Europe (JGSS-2026-346, under review at The Journals of Gerontology, Series B).

The code reproduces every table and figure in the manuscript and its supplementary material from the Health and Retirement Study (HRS) and the Survey of Health, Ageing and Retirement in Europe (SHARE).


Data availability

This repository contains no data. HRS and SHARE are available free of charge to registered users from their own repositories, and neither permits redistribution, so the two input files have to be built from the raw releases before anything here will run.

Survey	Release used	Where to get it
HRS	RAND HRS Longitudinal File, waves 1 to 16 (1992 to 2022)	https://hrsdata.isr.umich.edu
SHARE	SHARE release 9.0.0, all datasets	https://share-eric.eu

Analysis uses HRS waves 11 to 16 and SHARE waves 5 to 9, with waves 11 and 5 as the respective baselines. Verbal memory uses HRS waves 11 to 15 and SHARE waves 5 to 8, so that both surveys end in the same calendar year.

Repository contents
00_wrangle_hrs.R        raw RAND HRS   ->  data/derived/hrs_compact.rds
00_wrangle_share.R      raw SHARE      ->  data/derived/share_compact.rds
01_analysis_v6.R        the two .rds   ->  output/, main_v6.RData
02_sensitivity_v6.R     main_v6.RData  ->  output/, sens_v6.RData
config_example.R        template for an optional local path override

Run them in that order; each depends on the output of the one before it. All four use paths relative to the repository root, so set the working directory there and nothing else needs editing:

data/raw/        the raw survey releases, which you download yourself
data/derived/    the two compact files the wrangling scripts write
output/          every table, figure and workspace

data/ and output/ are in .gitignore, so no data and no results are ever committed. If your raw files sit somewhere else, copy config_example.R to config.R and set the paths there. Each script sources config.R after its own defaults, so whatever you put in it wins, and config.R is untracked.

1. 00_wrangle_hrs.R

Reads randhrs1992_2022v1.dta and writes hrs_compact.rds, one row per respondent. Keeps the baseline sample as those interviewed at wave 11, aged 50 to 100, with no ADL and no IADL limitation at that wave, and carries ADL, IADL, word recall and the covariates forward through wave 16.

ADL is RAND's five-item index r{w}adl5a: bathing, dressing, eating, bed transfer, walking across a room. Toileting is not in it.
IADL is built here by hand from the five component items r{w}moneya, medsa, shopa, mealsa, phonea, by plain addition, with a wave counted as missing if any item is unobserved. RAND's own r{w}iadl5a fills the index when items are unobserved, which treats an unobserved item as no difficulty. The RAND version is kept as iadl_rand_baseline for comparison.
Word recall is the 0 to 20 total r{w}tr20, immediate plus delayed. From wave 14 RAND splits it by mode into r{w}tr20p (face to face and telephone) and r{w}tr20w (web); both are kept, pooled, and flagged with a mode indicator, so web respondents are not dropped.
Education is coded from completed years: under 12 Low, 12 to 15 Medium, 16 or more High. Years are kept as well, for the sensitivity check against the ISCED crosswalk.
2. 00_wrangle_share.R

Reads the Gateway to Global Aging harmonised file GH_SHARE_g.dta and writes share_compact.rds, on the same pattern: interviewed at wave 5, aged 50 to 100, free of ADL and IADL limitation, followed through wave 9. Israel is excluded.

ADL is r{w}adla, verified against the data to contain the same five items as the HRS index.
IADL is built by hand from the same five component items as in the HRS script, and this is not optional. The harmonised r{w}iadla has only three items, money, medication and telephone, against RAND's five, and the two it omits, shopping and meal preparation, are the two most commonly endorsed in SHARE wave 5. Taking both summary indices as they come compares American difficulty on five items with European difficulty on the three rarest of those five. r5iadla is kept as iadl_gateway_baseline so the discrepancy is documented rather than asserted. The script header has the full note.
Word recall is r{w}tr20, as in the HRS.
Education is coded from native ISCED: 0 to 2 Low, 3 to 4 Medium, 5 to 6 High. A years-based version is kept alongside it.
Region is coded from the country code: Nordic (Sweden, Denmark), Western (Austria, Germany, Netherlands, Switzerland, Belgium, Luxembourg), Mediterranean (Spain, Italy, France), Eastern (Czech Republic, Slovenia, Estonia).
Wave spacing is uneven: waves 5 to 6 and 6 to 7 are two years, 7 to 8 is three (2017 to 2020, the COVID telephone wave), 8 to 9 is two. The analysis handles this with wave dummies.
3. 01_analysis_v6.R

Main analysis. Reads the two compact files and runs in twelve parts:

Part	What it does
0 to 3	settings, load, risk sets, person-wave construction
4	descriptives and incidence rates
5	discrete-time models for ADL, IADL and the threshold memory outcome
6	whether the gradient varies with age
7	predicted cumulative incidence
7b	linear mixed models of word recall, the main memory analysis
8	Karlson-Holm-Breen decomposition, education as the exposure
9	g-computation of the US minus Europe gap
10	horizon sensitivity and the education crossover
11	saves main_v6.RData
12	Figure 1, Figure S1, Figure S2 and the two forest plots

ADL and IADL onset are modelled as discrete-time survival on person-wave data, estimated as pooled logistic regressions with wave dummies and standard errors clustered on the person. Verbal memory is modelled with linear mixed models of the recall score, with a random intercept, a random slope for time and baseline age interacted with time. The threshold definition of memory decline is estimated alongside ADL and IADL and reported as a sensitivity analysis.

4. 02_sensitivity_v6.R

Robustness checks, in the order the Methods presents them:

Part	Check	Applies to
A	ascertainment frequency, each survey cut to its first three follow-up interviews	all outcomes
B	complementary log-log link	ADL, IADL
C	inverse probability of censoring weights	ADL, IADL
D	Benjamini-Hochberg correction, within each outcome	all outcomes
E1	memory specification, with and without baseline age interacted with time	memory
E2	memory random-effects structure	memory
E3	memory, first three follow-up interviews per survey	memory
E4	memory, baseline survey weights	memory
F	confirmed-events definition, not run, kept for reference	ADL, IADL

Parts B and C have no counterpart for the memory model: it is fitted by maximum likelihood on every observed score, so it is valid under missing at random without censoring weights, and a complementary log-log check means nothing in a model with an identity link.

Outputs

Both analysis scripts write to OUT_DIR, which defaults to output/.

Tables and model output (CSV). V4_01_table1, V4_02_tableS_selection, V4_03_incidence, V4_04_models, V4_05_models_floor, V4_06_age_interaction, V5_05b_all_coefficients, V5_07_cuminc_series, V5_07_cuminc_horizon, V5_08_khb_education, V5_09_gap_decomposition, V5_10_horizon_sensitivity, V5_10_crossover, S5_E2_memory_lmm_age, S5_E2_memory_predicted_age, S5_E2_memory_gaps, T_00_memory_means_by_sex, T_00_memory_baseline_by_sex, and from the sensitivity script S5_A_*, S5_B_cloglog, S5_C_ipcw, S5_D_multiplicity and S5_E_memory_*.

Figures (PNG). Figure1_forest3, FigureS1_memory_trajectories, FigureS2_M3_six_years, FigureS5_cuminc_matched, F1_forest_M1, F2_forest_M3.

Workspaces. main_v6.RData and sens_v6.RData.

File names carry the version in which each output first appeared and are left alone so that earlier result files still match; they do not track the current version number.

Reproducing the analysis

Put the two raw survey files in data/raw/, set the working directory to the repository root, and run:

r
install.packages(c("tidyverse", "haven", "sandwich", "lmtest", "MASS",
                   "lme4", "patchwork"))
source("00_wrangle_hrs.R")       # writes data/derived/hrs_compact.rds
source("00_wrangle_share.R")     # writes data/derived/share_compact.rds
source("01_analysis_v6.R")       # writes output/ and main_v6.RData
source("02_sensitivity_v6.R")    # reads main_v6.RData, writes output/


Part 8 draws 500 person-level bootstrap resamples per cell and Part 9 runs a simulation, so the main script takes a while and its bootstrap intervals move slightly between runs unless a seed is set.


Notes on the analysis:

Horizons. HORIZON sets eight years for ADL and IADL and six for verbal memory. The manuscript reports six years for all three, because that is the only horizon inside both surveys' observed range for every outcome. The horizons are points at which the fitted cumulative incidence is read, not a restriction of follow-up: the models use every observed wave.

Interpolation. SHARE interviews at two, four, seven and nine years, so neither six nor eight years falls on an interview and the European values at those horizons are interpolated between adjacent waves. Memory to eight years is extrapolated for the European regions, one year beyond the last memory interview, which is why six years is what the manuscript reports.

Unweighted estimation. The main models are unweighted because the target is an association within region and sex rather than a population rate. Weighted models are reported as a sensitivity analysis.

Why KHB. A logistic model fixes the residual variance, so adding mediators rescales every coefficient whether or not it confounds anything. The Karlson-Holm-Breen decomposition puts the reduced and full models on one scale so the difference between them is a clean indirect effect.

Contact

Octavio Bramajo. Aarhus University.
onbramajo@gmail.com
