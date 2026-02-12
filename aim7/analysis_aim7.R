# =============================================================================
# AIM 7 STATISTICAL ANALYSIS: Influence of Session Duration on sRPE
# =============================================================================
# REVISED VERSION - Simplified Analysis Without Residualization
#
# Research Questions:
#   1. Is sRPE contaminated by duration?
#      - If duration predicts sRPE after controlling for intensity, this
#        suggests duration influences perceived exertion independent of
#        actual physiological intensity.
#
#   2. Does sRPE predict physiological load beyond what duration explains?
#      - If adding sRPE to a duration-only model meaningfully improves fit,
#        sRPE contains intensity information not redundant with duration.
#
#   3. Does sRPE-TL predict physiological load better than duration alone?
#      - Practical question: Is collecting sRPE worth the effort compared
#        to simply recording session duration?
#
# Interpretation Logic:
#   Results from analyses 2 and 3 are interpreted together:
#   - If sRPE predicts beyond duration AND sRPE-TL outperforms duration alone,
#     then sRPE-TL's value derives at least partly from genuine intensity info.
#   - If sRPE adds little beyond duration BUT sRPE-TL still outperforms duration,
#     then sRPE-TL's apparent validity comes from duration appearing twice.
#
# Contexts:
#   - Practices: duration_mins as exposure
#   - Games: toi_total_mins (individual time-on-ice) as exposure
#
# =============================================================================

setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim7")

# =============================================================================
# 1. LOAD PACKAGES
# =============================================================================
library(lme4)
library(lmerTest)
library(performance)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("=============================================================================\n")
cat("AIM 7: Influence of Session Duration on sRPE\n")
cat("=============================================================================\n\n")

data <- read.csv("../aim1/aim1_analysis_data.csv")

cat("Data loaded:", nrow(data), "observations\n")
cat("  Practices:", sum(data$context == "practice"), "\n")
cat("  Games:", sum(data$context == "game"), "\n")

# =============================================================================
# 3. PREPARE DATA
# =============================================================================
cat("\n=== PREPARING DATA ===\n")

# --- PRACTICES ---
practices <- data %>%
  filter(context == "practice") %>%
  filter(!is.na(sRPE), !is.na(duration_mins), !is.na(TRIMP), !is.na(TRIMP_min), !is.na(HR80_mins)) %>%
  group_by(athlete_name) %>%
  mutate(
    n_obs = n(),
    # Within/between decomposition
    TRIMP_min_within = TRIMP_min - mean(TRIMP_min, na.rm = TRUE),
    TRIMP_min_between = mean(TRIMP_min, na.rm = TRUE),
    duration_within = duration_mins - mean(duration_mins, na.rm = TRUE),
    duration_between = mean(duration_mins, na.rm = TRUE),
    sRPE_within = sRPE - mean(sRPE, na.rm = TRUE),
    sRPE_between = mean(sRPE, na.rm = TRUE),
    TRIMP_within = TRIMP - mean(TRIMP, na.rm = TRUE),
    TRIMP_between = mean(TRIMP, na.rm = TRUE),
    HR80_within = HR80_mins - mean(HR80_mins, na.rm = TRUE),
    HR80_between = mean(HR80_mins, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(n_obs >= 5) %>%  # Minimum 5 observations per athlete
  mutate(
    # Z-score for standardized coefficients
    TRIMP_min_within_z = as.numeric(scale(TRIMP_min_within)),
    TRIMP_min_between_z = as.numeric(scale(TRIMP_min_between)),
    duration_within_z = as.numeric(scale(duration_within)),
    duration_between_z = as.numeric(scale(duration_between)),
    sRPE_within_z = as.numeric(scale(sRPE_within)),
    sRPE_between_z = as.numeric(scale(sRPE_between)),
    # sRPE-TL
    sRPE_TL = sRPE * duration_mins
  ) %>%
  group_by(athlete_name) %>%
  mutate(
    sRPE_TL_within = sRPE_TL - mean(sRPE_TL, na.rm = TRUE),
    sRPE_TL_between = mean(sRPE_TL, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    sRPE_TL_within_z = as.numeric(scale(sRPE_TL_within)),
    sRPE_TL_between_z = as.numeric(scale(sRPE_TL_between))
  )

# --- GAMES ---
games <- data %>%
  filter(context == "game", !is.na(toi_total_mins)) %>%
  filter(!is.na(sRPE), !is.na(TRIMP), !is.na(TRIMP_min), !is.na(HR80_mins)) %>%
  group_by(athlete_name) %>%
  mutate(
    n_obs = n(),
    TRIMP_min_within = TRIMP_min - mean(TRIMP_min, na.rm = TRUE),
    TRIMP_min_between = mean(TRIMP_min, na.rm = TRUE),
    toi_within = toi_total_mins - mean(toi_total_mins, na.rm = TRUE),
    toi_between = mean(toi_total_mins, na.rm = TRUE),
    sRPE_within = sRPE - mean(sRPE, na.rm = TRUE),
    sRPE_between = mean(sRPE, na.rm = TRUE),
    TRIMP_within = TRIMP - mean(TRIMP, na.rm = TRUE),
    TRIMP_between = mean(TRIMP, na.rm = TRUE),
    HR80_within = HR80_mins - mean(HR80_mins, na.rm = TRUE),
    HR80_between = mean(HR80_mins, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(n_obs >= 5) %>%  # Minimum 5 observations per athlete
  mutate(
    TRIMP_min_within_z = as.numeric(scale(TRIMP_min_within)),
    TRIMP_min_between_z = as.numeric(scale(TRIMP_min_between)),
    toi_within_z = as.numeric(scale(toi_within)),
    toi_between_z = as.numeric(scale(toi_between)),
    sRPE_within_z = as.numeric(scale(sRPE_within)),
    sRPE_between_z = as.numeric(scale(sRPE_between)),
    sRPE_TL = sRPE * toi_total_mins
  ) %>%
  group_by(athlete_name) %>%
  mutate(
    sRPE_TL_within = sRPE_TL - mean(sRPE_TL, na.rm = TRUE),
    sRPE_TL_between = mean(sRPE_TL, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    sRPE_TL_within_z = as.numeric(scale(sRPE_TL_within)),
    sRPE_TL_between_z = as.numeric(scale(sRPE_TL_between))
  )

cat("After filtering (>=5 obs per athlete):\n")
cat("  Practices:", nrow(practices), "obs from", n_distinct(practices$athlete_name), "athletes\n")
cat("  Games:", nrow(games), "obs from", n_distinct(games$athlete_name), "athletes\n")

# =============================================================================
# 4. HELPER FUNCTIONS
# =============================================================================
safe_lmer <- function(formula, data) {
  tryCatch(
    lmer(formula, data = data, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))),
    error = function(e) NULL
  )
}

get_r2 <- function(model) {
  if (is.null(model)) return(NA)
  r2 <- tryCatch(performance::r2(model), error = function(e) NULL)
  if (!is.null(r2)) r2$R2_marginal else NA
}

get_beta <- function(model, term) {
  if (is.null(model)) return(c(beta = NA, SE = NA, CI_lower = NA, CI_upper = NA, p = NA))
  fe <- summary(model)$coefficients
  if (term %in% rownames(fe)) {
    beta <- fe[term, "Estimate"]
    se <- fe[term, "Std. Error"]
    ci <- tryCatch(confint(model, parm = term, method = "Wald"), error = function(e) c(NA, NA))
    c(beta = beta, 
      SE = se, 
      CI_lower = ci[1],
      CI_upper = ci[2],
      p = fe[term, "Pr(>|t|)"])
  } else {
    c(beta = NA, SE = NA, CI_lower = NA, CI_upper = NA, p = NA)
  }
}

# =============================================================================
# 5. ANALYSIS 1: IS sRPE CONTAMINATED BY DURATION?
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# ANALYSIS 1: Is sRPE Contaminated by Duration?                           #\n")
cat("# (If sRPE purely reflects intensity, duration effect should be ~0)       #\n")
cat("###########################################################################\n")

# --- PRACTICES ---
cat("\n=== PRACTICES ===\n")

m_srpe_both_p <- safe_lmer(
  sRPE ~ TRIMP_min_within_z + TRIMP_min_between_z + 
    duration_within_z + duration_between_z + (1 | athlete_name),
  data = practices
)

beta_int_p <- get_beta(m_srpe_both_p, "TRIMP_min_within_z")
beta_dur_p <- get_beta(m_srpe_both_p, "duration_within_z")

cat("Model: sRPE ~ TRIMPmin + Duration + (1|Athlete)\n\n")
cat("Within-athlete effects (standardized):\n")
cat("  Intensity (TRIMPmin): β =", round(beta_int_p["beta"], 3),
    ", 95% CI [", round(beta_int_p["CI_lower"], 3), ",", round(beta_int_p["CI_upper"], 3), "]",
    ", p =", ifelse(beta_int_p["p"] < 0.001, "< .001", round(beta_int_p["p"], 3)), "\n")
cat("  Duration:             β =", round(beta_dur_p["beta"], 3),
    ", 95% CI [", round(beta_dur_p["CI_lower"], 3), ",", round(beta_dur_p["CI_upper"], 3), "]",
    ", p =", ifelse(beta_dur_p["p"] < 0.001, "< .001", round(beta_dur_p["p"], 3)), "\n")

# Relative magnitude
ratio_p <- abs(beta_dur_p["beta"]) / abs(beta_int_p["beta"])
cat("\nDuration/Intensity ratio:", round(ratio_p, 2), "\n")
if (ratio_p > 0.5) {
  cat("--> Duration effect is substantial relative to intensity.\n")
  cat("    sRPE appears contaminated by duration.\n")
} else {
  cat("--> Duration effect is small relative to intensity.\n")
  cat("    sRPE primarily reflects physiological intensity.\n")
}

# --- GAMES ---
cat("\n=== GAMES ===\n")

m_srpe_both_g <- safe_lmer(
  sRPE ~ TRIMP_min_within_z + TRIMP_min_between_z + 
    toi_within_z + toi_between_z + (1 | athlete_name),
  data = games
)

beta_int_g <- get_beta(m_srpe_both_g, "TRIMP_min_within_z")
beta_toi_g <- get_beta(m_srpe_both_g, "toi_within_z")

cat("Model: sRPE ~ TRIMPmin + TOI + (1|Athlete)\n\n")
cat("Within-athlete effects (standardized):\n")
cat("  Intensity (TRIMPmin): β =", round(beta_int_g["beta"], 3),
    ", 95% CI [", round(beta_int_g["CI_lower"], 3), ",", round(beta_int_g["CI_upper"], 3), "]",
    ", p =", ifelse(beta_int_g["p"] < 0.001, "< .001", round(beta_int_g["p"], 3)), "\n")
cat("  Time-on-Ice:          β =", round(beta_toi_g["beta"], 3),
    ", 95% CI [", round(beta_toi_g["CI_lower"], 3), ",", round(beta_toi_g["CI_upper"], 3), "]",
    ", p =", ifelse(beta_toi_g["p"] < 0.001, "< .001", round(beta_toi_g["p"], 3)), "\n")

ratio_g <- abs(beta_toi_g["beta"]) / abs(beta_int_g["beta"])
cat("\nTOI/Intensity ratio:", round(ratio_g, 2), "\n")

# Store Analysis 1 results
analysis1_results <- data.frame(
  context = c("Practice", "Practice", "Game", "Game"),
  predictor = c("Intensity (TRIMPmin)", "Duration", "Intensity (TRIMPmin)", "TOI"),
  beta = c(beta_int_p["beta"], beta_dur_p["beta"], beta_int_g["beta"], beta_toi_g["beta"]),
  SE = c(beta_int_p["SE"], beta_dur_p["SE"], beta_int_g["SE"], beta_toi_g["SE"]),
  CI_lower = c(beta_int_p["CI_lower"], beta_dur_p["CI_lower"], beta_int_g["CI_lower"], beta_toi_g["CI_lower"]),
  CI_upper = c(beta_int_p["CI_upper"], beta_dur_p["CI_upper"], beta_int_g["CI_upper"], beta_toi_g["CI_upper"]),
  p = c(beta_int_p["p"], beta_dur_p["p"], beta_int_g["p"], beta_toi_g["p"])
)

# =============================================================================
# 6. ANALYSIS 2: DOES sRPE PREDICT TRIMP BEYOND DURATION?
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# ANALYSIS 2: Does sRPE Predict Physiological Load Beyond Duration?       #\n")
cat("# (Comparing duration-only model vs. duration + sRPE model)               #\n")
cat("###########################################################################\n")

# --- PRACTICES ---
cat("\n=== PRACTICES ===\n")

# Duration-only model predicting TRIMP
m_trimp_dur_p <- safe_lmer(
  TRIMP ~ duration_within_z + duration_between_z + (1 | athlete_name),
  data = practices
)

# Duration + sRPE model predicting TRIMP
m_trimp_both_p <- safe_lmer(
  TRIMP ~ duration_within_z + duration_between_z + 
    sRPE_within_z + sRPE_between_z + (1 | athlete_name),
  data = practices
)

r2_dur_p <- get_r2(m_trimp_dur_p)
r2_both_p <- get_r2(m_trimp_both_p)
r2_increment_p <- r2_both_p - r2_dur_p

beta_srpe_trimp_p <- get_beta(m_trimp_both_p, "sRPE_within_z")

cat("Predicting TRIMP:\n")
cat("  Duration only:     R² =", round(r2_dur_p, 4), "\n")
cat("  Duration + sRPE:   R² =", round(r2_both_p, 4), "\n")
cat("  Increment (ΔR²):       ", round(r2_increment_p, 4), "\n")
cat("\nsRPE effect (controlling for duration):\n")
cat("  β =", round(beta_srpe_trimp_p["beta"], 3),
    ", p =", ifelse(beta_srpe_trimp_p["p"] < 0.001, "< .001", round(beta_srpe_trimp_p["p"], 3)), "\n")

# Same for HR80
m_hr80_dur_p <- safe_lmer(
  HR80_mins ~ duration_within_z + duration_between_z + (1 | athlete_name),
  data = practices
)

m_hr80_both_p <- safe_lmer(
  HR80_mins ~ duration_within_z + duration_between_z + 
    sRPE_within_z + sRPE_between_z + (1 | athlete_name),
  data = practices
)

r2_hr80_dur_p <- get_r2(m_hr80_dur_p)
r2_hr80_both_p <- get_r2(m_hr80_both_p)
r2_hr80_increment_p <- r2_hr80_both_p - r2_hr80_dur_p

beta_srpe_hr80_p <- get_beta(m_hr80_both_p, "sRPE_within_z")

cat("\nPredicting HR80 (minutes above 80% HRmax):\n")
cat("  Duration only:     R² =", round(r2_hr80_dur_p, 4), "\n")
cat("  Duration + sRPE:   R² =", round(r2_hr80_both_p, 4), "\n")
cat("  Increment (ΔR²):       ", round(r2_hr80_increment_p, 4), "\n")
cat("\nsRPE effect (controlling for duration):\n")
cat("  β =", round(beta_srpe_hr80_p["beta"], 3),
    ", p =", ifelse(beta_srpe_hr80_p["p"] < 0.001, "< .001", round(beta_srpe_hr80_p["p"], 3)), "\n")

# --- GAMES ---
cat("\n=== GAMES ===\n")

m_trimp_toi_g <- safe_lmer(
  TRIMP ~ toi_within_z + toi_between_z + (1 | athlete_name),
  data = games
)

m_trimp_both_g <- safe_lmer(
  TRIMP ~ toi_within_z + toi_between_z + 
    sRPE_within_z + sRPE_between_z + (1 | athlete_name),
  data = games
)

r2_toi_g <- get_r2(m_trimp_toi_g)
r2_both_g <- get_r2(m_trimp_both_g)
r2_increment_g <- r2_both_g - r2_toi_g

beta_srpe_trimp_g <- get_beta(m_trimp_both_g, "sRPE_within_z")

cat("Predicting TRIMP:\n")
cat("  TOI only:          R² =", round(r2_toi_g, 4), "\n")
cat("  TOI + sRPE:        R² =", round(r2_both_g, 4), "\n")
cat("  Increment (ΔR²):       ", round(r2_increment_g, 4), "\n")
cat("\nsRPE effect (controlling for TOI):\n")
cat("  β =", round(beta_srpe_trimp_g["beta"], 3),
    ", p =", ifelse(beta_srpe_trimp_g["p"] < 0.001, "< .001", round(beta_srpe_trimp_g["p"], 3)), "\n")

# HR80 for games
m_hr80_toi_g <- safe_lmer(
  HR80_mins ~ toi_within_z + toi_between_z + (1 | athlete_name),
  data = games
)

m_hr80_both_g <- safe_lmer(
  HR80_mins ~ toi_within_z + toi_between_z + 
    sRPE_within_z + sRPE_between_z + (1 | athlete_name),
  data = games
)

r2_hr80_toi_g <- get_r2(m_hr80_toi_g)
r2_hr80_both_g <- get_r2(m_hr80_both_g)
r2_hr80_increment_g <- r2_hr80_both_g - r2_hr80_toi_g

beta_srpe_hr80_g <- get_beta(m_hr80_both_g, "sRPE_within_z")

cat("\nPredicting HR80:\n")
cat("  TOI only:          R² =", round(r2_hr80_toi_g, 4), "\n")
cat("  TOI + sRPE:        R² =", round(r2_hr80_both_g, 4), "\n")
cat("  Increment (ΔR²):       ", round(r2_hr80_increment_g, 4), "\n")
cat("\nsRPE effect (controlling for TOI):\n")
cat("  β =", round(beta_srpe_hr80_g["beta"], 3),
    ", p =", ifelse(beta_srpe_hr80_g["p"] < 0.001, "< .001", round(beta_srpe_hr80_g["p"], 3)), "\n")

# Store Analysis 2 results
analysis2_results <- data.frame(
  context = c("Practice", "Practice", "Game", "Game"),
  outcome = c("TRIMP", "HR80", "TRIMP", "HR80"),
  r2_duration_only = c(r2_dur_p, r2_hr80_dur_p, r2_toi_g, r2_hr80_toi_g),
  r2_duration_plus_srpe = c(r2_both_p, r2_hr80_both_p, r2_both_g, r2_hr80_both_g),
  r2_increment = c(r2_increment_p, r2_hr80_increment_p, r2_increment_g, r2_hr80_increment_g),
  srpe_beta = c(beta_srpe_trimp_p["beta"], beta_srpe_hr80_p["beta"], 
                beta_srpe_trimp_g["beta"], beta_srpe_hr80_g["beta"]),
  srpe_p = c(beta_srpe_trimp_p["p"], beta_srpe_hr80_p["p"],
             beta_srpe_trimp_g["p"], beta_srpe_hr80_g["p"])
)

# =============================================================================
# 7. ANALYSIS 3: DOES sRPE-TL PREDICT BETTER THAN DURATION ALONE?
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# ANALYSIS 3: Does sRPE-TL Predict Better Than Duration Alone?            #\n")
cat("# (Practical question: Is collecting sRPE worth the effort?)              #\n")
cat("###########################################################################\n")

# --- PRACTICES ---
cat("\n=== PRACTICES ===\n")

# Duration-only model
m_p_trimp_dur <- safe_lmer(
  TRIMP ~ duration_within_z + duration_between_z + (1 | athlete_name),
  data = practices
)

# sRPE-TL model
m_p_trimp_tl <- safe_lmer(
  TRIMP ~ sRPE_TL_within_z + sRPE_TL_between_z + (1 | athlete_name),
  data = practices
)

r2_p_trimp_dur <- get_r2(m_p_trimp_dur)
r2_p_trimp_tl <- get_r2(m_p_trimp_tl)
tl_improvement_trimp_p <- r2_p_trimp_tl - r2_p_trimp_dur

cat("Predicting TRIMP:\n")
cat("  Duration alone:    R² =", round(r2_p_trimp_dur, 4), "\n")
cat("  sRPE-TL:           R² =", round(r2_p_trimp_tl, 4), "\n")
cat("  Improvement (ΔR²):     ", round(tl_improvement_trimp_p, 4), 
    " (", round(tl_improvement_trimp_p/r2_p_trimp_dur*100, 1), "% relative improvement)\n")

# HR80
m_p_hr80_dur <- safe_lmer(
  HR80_mins ~ duration_within_z + duration_between_z + (1 | athlete_name),
  data = practices
)

m_p_hr80_tl <- safe_lmer(
  HR80_mins ~ sRPE_TL_within_z + sRPE_TL_between_z + (1 | athlete_name),
  data = practices
)

r2_p_hr80_dur <- get_r2(m_p_hr80_dur)
r2_p_hr80_tl <- get_r2(m_p_hr80_tl)
tl_improvement_hr80_p <- r2_p_hr80_tl - r2_p_hr80_dur

cat("\nPredicting HR80:\n")
cat("  Duration alone:    R² =", round(r2_p_hr80_dur, 4), "\n")
cat("  sRPE-TL:           R² =", round(r2_p_hr80_tl, 4), "\n")
cat("  Improvement (ΔR²):     ", round(tl_improvement_hr80_p, 4),
    " (", round(tl_improvement_hr80_p/r2_p_hr80_dur*100, 1), "% relative improvement)\n")

# --- GAMES ---
cat("\n=== GAMES ===\n")

m_g_trimp_toi <- safe_lmer(
  TRIMP ~ toi_within_z + toi_between_z + (1 | athlete_name),
  data = games
)

m_g_trimp_tl <- safe_lmer(
  TRIMP ~ sRPE_TL_within_z + sRPE_TL_between_z + (1 | athlete_name),
  data = games
)

r2_g_trimp_toi <- get_r2(m_g_trimp_toi)
r2_g_trimp_tl <- get_r2(m_g_trimp_tl)
tl_improvement_trimp_g <- r2_g_trimp_tl - r2_g_trimp_toi

cat("Predicting TRIMP:\n")
cat("  TOI alone:         R² =", round(r2_g_trimp_toi, 4), "\n")
cat("  sRPE-TL:           R² =", round(r2_g_trimp_tl, 4), "\n")
cat("  Improvement (ΔR²):     ", round(tl_improvement_trimp_g, 4),
    " (", round(tl_improvement_trimp_g/r2_g_trimp_toi*100, 1), "% relative improvement)\n")

m_g_hr80_toi <- safe_lmer(
  HR80_mins ~ toi_within_z + toi_between_z + (1 | athlete_name),
  data = games
)

m_g_hr80_tl <- safe_lmer(
  HR80_mins ~ sRPE_TL_within_z + sRPE_TL_between_z + (1 | athlete_name),
  data = games
)

r2_g_hr80_toi <- get_r2(m_g_hr80_toi)
r2_g_hr80_tl <- get_r2(m_g_hr80_tl)
tl_improvement_hr80_g <- r2_g_hr80_tl - r2_g_hr80_toi

cat("\nPredicting HR80:\n")
cat("  TOI alone:         R² =", round(r2_g_hr80_toi, 4), "\n")
cat("  sRPE-TL:           R² =", round(r2_g_hr80_tl, 4), "\n")
cat("  Improvement (ΔR²):     ", round(tl_improvement_hr80_g, 4),
    " (", round(tl_improvement_hr80_g/r2_g_hr80_toi*100, 1), "% relative improvement)\n")

# Store Analysis 3 results
analysis3_results <- data.frame(
  context = c("Practice", "Practice", "Game", "Game"),
  outcome = c("TRIMP", "HR80", "TRIMP", "HR80"),
  r2_duration_only = c(r2_p_trimp_dur, r2_p_hr80_dur, r2_g_trimp_toi, r2_g_hr80_toi),
  r2_srpe_tl = c(r2_p_trimp_tl, r2_p_hr80_tl, r2_g_trimp_tl, r2_g_hr80_tl),
  improvement = c(tl_improvement_trimp_p, tl_improvement_hr80_p, 
                  tl_improvement_trimp_g, tl_improvement_hr80_g),
  relative_improvement_pct = c(
    tl_improvement_trimp_p/r2_p_trimp_dur*100,
    tl_improvement_hr80_p/r2_p_hr80_dur*100,
    tl_improvement_trimp_g/r2_g_trimp_toi*100,
    tl_improvement_hr80_g/r2_g_hr80_toi*100
  )
)

# =============================================================================
# 8. INTEGRATED INTERPRETATION
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# INTEGRATED INTERPRETATION                                               #\n")
cat("###########################################################################\n")

cat("\n=== PRACTICES ===\n")
cat("Analysis 2: Does sRPE add to duration?\n")
cat("  TRIMP: sRPE adds ΔR² =", round(r2_increment_p, 4), "\n")
cat("  HR80:  sRPE adds ΔR² =", round(r2_hr80_increment_p, 4), "\n")

cat("\nAnalysis 3: Does sRPE-TL outperform duration?\n")
cat("  TRIMP: sRPE-TL improvement =", round(tl_improvement_trimp_p, 4), "\n")
cat("  HR80:  sRPE-TL improvement =", round(tl_improvement_hr80_p, 4), "\n")

# Interpretation logic
srpe_adds_p <- (r2_increment_p > 0.01) || (r2_hr80_increment_p > 0.01)
tl_better_p <- (tl_improvement_trimp_p > 0.01) || (tl_improvement_hr80_p > 0.01)

cat("\nInterpretation:\n")
if (srpe_adds_p && tl_better_p) {
  cat("  sRPE adds meaningful information beyond duration AND sRPE-TL outperforms duration.\n")
  cat("  --> sRPE-TL's value derives at least partly from genuine intensity information.\n")
} else if (!srpe_adds_p && tl_better_p) {
  cat("  sRPE adds little beyond duration, BUT sRPE-TL still outperforms duration.\n")
  cat("  --> sRPE-TL's apparent validity may come from duration appearing twice\n")
  cat("      in the calculation (once in sRPE ratings, once in multiplication).\n")
} else if (srpe_adds_p && !tl_better_p) {
  cat("  sRPE adds information beyond duration, but sRPE-TL doesn't outperform duration.\n")
  cat("  --> Unexpected pattern; requires further investigation.\n")
} else {
  cat("  sRPE adds little beyond duration AND sRPE-TL doesn't outperform duration.\n")
  cat("  --> sRPE may not provide meaningful value over simple duration tracking.\n")
}

cat("\n=== GAMES ===\n")
cat("Analysis 2: Does sRPE add to TOI?\n")
cat("  TRIMP: sRPE adds ΔR² =", round(r2_increment_g, 4), "\n")
cat("  HR80:  sRPE adds ΔR² =", round(r2_hr80_increment_g, 4), "\n")

cat("\nAnalysis 3: Does sRPE-TL outperform TOI?\n")
cat("  TRIMP: sRPE-TL improvement =", round(tl_improvement_trimp_g, 4), "\n")
cat("  HR80:  sRPE-TL improvement =", round(tl_improvement_hr80_g, 4), "\n")

srpe_adds_g <- (r2_increment_g > 0.01) || (r2_hr80_increment_g > 0.01)
tl_better_g <- (tl_improvement_trimp_g > 0.01) || (tl_improvement_hr80_g > 0.01)

cat("\nInterpretation:\n")
if (srpe_adds_g && tl_better_g) {
  cat("  sRPE adds meaningful information beyond TOI AND sRPE-TL outperforms TOI.\n")
  cat("  --> sRPE-TL's value derives at least partly from genuine intensity information.\n")
} else if (!srpe_adds_g && tl_better_g) {
  cat("  sRPE adds little beyond TOI, BUT sRPE-TL still outperforms TOI.\n")
  cat("  --> sRPE-TL's apparent validity may come from duration appearing twice.\n")
} else if (srpe_adds_g && !tl_better_g) {
  cat("  sRPE adds information beyond TOI, but sRPE-TL doesn't outperform TOI.\n")
  cat("  --> Unexpected pattern; requires further investigation.\n")
} else {
  cat("  sRPE adds little beyond TOI AND sRPE-TL doesn't outperform TOI.\n")
  cat("  --> sRPE may not provide meaningful value over simple TOI tracking.\n")
}

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================
cat("\n=== SAVING RESULTS ===\n")

write.csv(analysis1_results, "results_analysis1_contamination.csv", row.names = FALSE)
write.csv(analysis2_results, "results_analysis2_srpe_beyond_duration.csv", row.names = FALSE)
write.csv(analysis3_results, "results_analysis3_value_added.csv", row.names = FALSE)

# Summary table
summary_aim7 <- data.frame(
  finding = c(
    "N_practices", "N_games",
    "N_athletes_practices", "N_athletes_games",
    "Practice: Intensity β on sRPE", 
    "Practice: Duration β on sRPE",
    "Practice: Duration/Intensity ratio",
    "Game: Intensity β on sRPE", 
    "Game: TOI β on sRPE", 
    "Game: TOI/Intensity ratio",
    "Practice: sRPE adds to duration (TRIMP ΔR²)",
    "Practice: sRPE adds to duration (HR80 ΔR²)",
    "Game: sRPE adds to TOI (TRIMP ΔR²)",
    "Game: sRPE adds to TOI (HR80 ΔR²)",
    "Practice: sRPE-TL improvement over Duration (TRIMP)",
    "Practice: sRPE-TL improvement over Duration (HR80)",
    "Game: sRPE-TL improvement over TOI (TRIMP)",
    "Game: sRPE-TL improvement over TOI (HR80)"
  ),
  value = c(
    nrow(practices), nrow(games),
    n_distinct(practices$athlete_name), n_distinct(games$athlete_name),
    beta_int_p["beta"], beta_dur_p["beta"], ratio_p,
    beta_int_g["beta"], beta_toi_g["beta"], ratio_g,
    r2_increment_p, r2_hr80_increment_p,
    r2_increment_g, r2_hr80_increment_g,
    tl_improvement_trimp_p, tl_improvement_hr80_p,
    tl_improvement_trimp_g, tl_improvement_hr80_g
  )
)

write.csv(summary_aim7, "summary_aim7.csv", row.names = FALSE)

# Save models
save(m_srpe_both_p, m_srpe_both_g,
     m_trimp_dur_p, m_trimp_both_p,
     m_hr80_dur_p, m_hr80_both_p,
     m_trimp_toi_g, m_trimp_both_g,
     m_hr80_toi_g, m_hr80_both_g,
     m_p_trimp_dur, m_p_trimp_tl,
     m_p_hr80_dur, m_p_hr80_tl,
     m_g_trimp_toi, m_g_trimp_tl,
     m_g_hr80_toi, m_g_hr80_tl,
     file = "aim7_models.RData")

cat("Files saved:\n")
cat("  - results_analysis1_contamination.csv\n")
cat("  - results_analysis2_srpe_beyond_duration.csv\n")
cat("  - results_analysis3_value_added.csv\n")
cat("  - summary_aim7.csv\n")
cat("  - aim7_models.RData\n")

# =============================================================================
# 10. FINAL SUMMARY
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# AIM 7 FINAL SUMMARY                                                     #\n")
cat("###########################################################################\n")

cat("\n=== ANALYSIS 1: Is sRPE contaminated by duration? ===\n")
cat("PRACTICES:\n")
cat("  Duration β =", round(beta_dur_p["beta"], 3), 
    ", Intensity β =", round(beta_int_p["beta"], 3),
    ", Ratio =", round(ratio_p, 2), "\n")
cat("GAMES:\n")
cat("  TOI β =", round(beta_toi_g["beta"], 3),
    ", Intensity β =", round(beta_int_g["beta"], 3),
    ", Ratio =", round(ratio_g, 2), "\n")

if (ratio_p > 0.5 || ratio_g > 0.5) {
  cat("--> sRPE is substantially influenced by duration.\n")
} else {
  cat("--> sRPE primarily reflects intensity with minimal duration contamination.\n")
}

cat("\n=== ANALYSIS 2: Does sRPE add information beyond duration? ===\n")
cat("PRACTICES:\n")
cat("  TRIMP: ΔR² =", round(r2_increment_p, 4), "\n")
cat("  HR80:  ΔR² =", round(r2_hr80_increment_p, 4), "\n")
cat("GAMES:\n")
cat("  TRIMP: ΔR² =", round(r2_increment_g, 4), "\n")
cat("  HR80:  ΔR² =", round(r2_hr80_increment_g, 4), "\n")

cat("\n=== ANALYSIS 3: Does sRPE-TL outperform duration alone? ===\n")
cat("PRACTICES:\n")
cat("  TRIMP: Improvement =", round(tl_improvement_trimp_p, 4),
    " (", round(tl_improvement_trimp_p/r2_p_trimp_dur*100, 1), "%)\n")
cat("  HR80:  Improvement =", round(tl_improvement_hr80_p, 4),
    " (", round(tl_improvement_hr80_p/r2_p_hr80_dur*100, 1), "%)\n")
cat("GAMES:\n")
cat("  TRIMP: Improvement =", round(tl_improvement_trimp_g, 4),
    " (", round(tl_improvement_trimp_g/r2_g_trimp_toi*100, 1), "%)\n")
cat("  HR80:  Improvement =", round(tl_improvement_hr80_g, 4),
    " (", round(tl_improvement_hr80_g/r2_g_hr80_toi*100, 1), "%)\n")

cat("\n=== BOTTOM LINE ===\n")
cat("Interpreting Analyses 2 and 3 together reveals whether sRPE-TL's validity\n")
cat("derives from genuine intensity information or from duration appearing twice.\n")

cat("\n=============================================================================\n")
cat("AIM 7 ANALYSIS COMPLETE\n")
cat("=============================================================================\n")