# =============================================================================
# AIM 6 STATISTICAL ANALYSIS: Wellbeing, Perceived Exertion, and Performance
# =============================================================================
# REVISED VERSION - Descriptive Associations (No Causal/Pathway Language)
#
# Research Questions:
#   1. Does morning wellbeing predict same-day sRPE?
#   2. Does this association persist after accounting for session duration/TOI?
#   3. Does morning wellbeing predict post-game performance?
#
# Rationale for Question 2:
#   If wellbeing predicts sRPE but this association disappears after accounting
#   for duration, then wellbeing is essentially providing information about 
#   session exposure that could be obtained directly from duration data.
#   If the association persists, wellbeing captures something about how athletes
#   experience sessions that is not reducible to how long they trained.
#
# =============================================================================

setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim6")

# =============================================================================
# 1. LOAD PACKAGES
# =============================================================================
library(lme4)
library(lmerTest)
library(ordinal)
library(performance)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("=============================================================================\n")
cat("AIM 6: Wellbeing, Perceived Exertion, and Performance\n")
cat("=============================================================================\n\n")

data <- read.csv("../merged_data_full.csv")
data$session_date <- as.Date(data$session_date)

# =============================================================================
# 3. PREPARE WELLBEING COMPOSITE
# =============================================================================
cat("=== PREPARING DATA ===\n")

# Wellbeing items:
# - confidence, excitement, sleep_quality: higher = better
# - fatigue, stress: higher = WORSE (reverse score)

# Step 1: Filter to complete cases
analysis_data <- data %>%
  filter(!is.na(sRPE),
         !is.na(confidence),
         !is.na(excitement),
         !is.na(fatigue),
         !is.na(sleep_quality),
         !is.na(stress),
         full_participation == 1) %>%
  mutate(
    # Reverse score fatigue and stress (so higher = better for all items)
    fatigue_r = 10 - fatigue,
    stress_r = 10 - stress,
    
    # Create unified exposure variable (duration for practices, TOI for games)
    exposure_mins = case_when(
      context == "practice" ~ duration_mins,
      context == "game" ~ toi_total_mins,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(exposure_mins))

# Step 2: Within-athlete z-scoring of wellbeing items
# This removes between-athlete scaling differences and ensures the composite
# reflects deviations from each athlete's typical state
analysis_data <- analysis_data %>%
  group_by(athlete_name) %>%
  mutate(
    # Z-score each item WITHIN athlete
    confidence_z = as.numeric(scale(confidence)),
    excitement_z = as.numeric(scale(excitement)),
    fatigue_r_z = as.numeric(scale(fatigue_r)),
    sleep_quality_z = as.numeric(scale(sleep_quality)),
    stress_r_z = as.numeric(scale(stress_r))
  ) %>%
  ungroup()

# Handle athletes with no variance (constant values) - replace NaN with 0
analysis_data <- analysis_data %>%
  mutate(across(c(confidence_z, excitement_z, fatigue_r_z, sleep_quality_z, stress_r_z),
                ~ifelse(is.nan(.), 0, .)))

# Step 3: Create composite (average of within-athlete z-scores)
analysis_data <- analysis_data %>%
  mutate(
    wellbeing_composite = (confidence_z + excitement_z + fatigue_r_z + 
                             sleep_quality_z + stress_r_z) / 5
  )

cat("Observations with complete data:", nrow(analysis_data), "\n")
cat("  Practices:", sum(analysis_data$context == "practice"), "\n")
cat("  Games:", sum(analysis_data$context == "game"), "\n")
cat("  Athletes:", n_distinct(analysis_data$athlete_name), "\n")

# =============================================================================
# 4. WITHIN/BETWEEN DECOMPOSITION
# =============================================================================

# Separate by context for proper decomposition
practices <- analysis_data %>%
  filter(context == "practice") %>%
  group_by(athlete_name) %>%
  mutate(
    n_obs = n(),
    wellbeing_mean = mean(wellbeing_composite, na.rm = TRUE),
    wellbeing_within = wellbeing_composite - wellbeing_mean,
    duration_mean = mean(exposure_mins, na.rm = TRUE),
    duration_within = exposure_mins - duration_mean
  ) %>%
  ungroup() %>%
  filter(n_obs >= 5) %>%  # Minimum 5 observations per athlete
  mutate(
    # Standardize for interpretability
    wellbeing_within_z = as.numeric(scale(wellbeing_within)),
    wellbeing_between_z = as.numeric(scale(wellbeing_mean)),
    duration_within_z = as.numeric(scale(duration_within)),
    duration_between_z = as.numeric(scale(duration_mean))
  )

games <- analysis_data %>%
  filter(context == "game") %>%
  group_by(athlete_name) %>%
  mutate(
    n_obs = n(),
    wellbeing_mean = mean(wellbeing_composite, na.rm = TRUE),
    wellbeing_within = wellbeing_composite - wellbeing_mean,
    toi_mean = mean(exposure_mins, na.rm = TRUE),
    toi_within = exposure_mins - toi_mean
  ) %>%
  ungroup() %>%
  filter(n_obs >= 5) %>%  # Minimum 5 observations per athlete
  mutate(
    wellbeing_within_z = as.numeric(scale(wellbeing_within)),
    wellbeing_between_z = as.numeric(scale(wellbeing_mean)),
    toi_within_z = as.numeric(scale(toi_within)),
    toi_between_z = as.numeric(scale(toi_mean))
  )

cat("\nAfter filtering (>=5 obs per athlete):\n")
cat("  Practices:", nrow(practices), "obs from", n_distinct(practices$athlete_name), "athletes\n")
cat("  Games:", nrow(games), "obs from", n_distinct(games$athlete_name), "athletes\n")

# =============================================================================
# 5. HELPER FUNCTIONS
# =============================================================================
safe_lmer <- function(formula, data) {
  tryCatch(
    lmer(formula, data = data, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))),
    error = function(e) NULL
  )
}

extract_effect <- function(model, term) {
  if (is.null(model)) return(c(beta = NA, SE = NA, CI_lower = NA, CI_upper = NA, p = NA))
  fe <- summary(model)$coefficients
  if (term %in% rownames(fe)) {
    beta <- fe[term, "Estimate"]
    se <- fe[term, "Std. Error"]
    ci <- confint(model, parm = term, method = "Wald")
    c(beta = beta, 
      SE = se, 
      CI_lower = ci[1],
      CI_upper = ci[2],
      p = fe[term, "Pr(>|t|)"])
  } else {
    c(beta = NA, SE = NA, CI_lower = NA, CI_upper = NA, p = NA)
  }
}

get_r2 <- function(model) {
  if (is.null(model)) return(NA)
  r2 <- tryCatch(performance::r2(model), error = function(e) NULL)
  if (!is.null(r2)) r2$R2_marginal else NA
}

# =============================================================================
# 6. DESCRIPTIVES
# =============================================================================
cat("\n=== DESCRIPTIVES ===\n")

desc_by_context <- analysis_data %>%
  group_by(context) %>%
  summarise(
    n = n(),
    n_athletes = n_distinct(athlete_name),
    wellbeing_mean = mean(wellbeing_composite),
    wellbeing_sd = sd(wellbeing_composite),
    sRPE_mean = mean(sRPE),
    sRPE_sd = sd(sRPE),
    exposure_mean = mean(exposure_mins),
    exposure_sd = sd(exposure_mins),
    .groups = "drop"
  )

cat("\nBy Context:\n")
print(desc_by_context %>% mutate(across(where(is.numeric), ~round(., 2))))

# =============================================================================
# 7. WELLBEING-sRPE ASSOCIATIONS - PRACTICES
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# PRACTICES: Wellbeing and sRPE Associations                              #\n")
cat("###########################################################################\n")

# --- A. Unadjusted association: Does wellbeing predict sRPE? ---
cat("\n=== A. Wellbeing → sRPE (Unadjusted) ===\n")

m_p_srpe <- safe_lmer(
  sRPE ~ wellbeing_within_z + wellbeing_between_z + (1 | athlete_name),
  data = practices
)

eff_p_srpe <- extract_effect(m_p_srpe, "wellbeing_within_z")
eff_p_srpe_bw <- extract_effect(m_p_srpe, "wellbeing_between_z")

cat("Within-athlete: β =", round(eff_p_srpe["beta"], 3), 
    ", 95% CI [", round(eff_p_srpe["CI_lower"], 3), ",", round(eff_p_srpe["CI_upper"], 3), "]",
    ", p =", round(eff_p_srpe["p"], 3), "\n")
cat("Between-athlete: β =", round(eff_p_srpe_bw["beta"], 3),
    ", p =", round(eff_p_srpe_bw["p"], 3), "\n")

# --- B. Adjusted association: Does wellbeing predict sRPE after accounting for duration? ---
cat("\n=== B. Wellbeing → sRPE (Adjusted for Duration) ===\n")

m_p_srpe_adj <- safe_lmer(
  sRPE ~ wellbeing_within_z + wellbeing_between_z + 
    duration_within_z + duration_between_z + (1 | athlete_name),
  data = practices
)

eff_p_srpe_adj <- extract_effect(m_p_srpe_adj, "wellbeing_within_z")
eff_p_dur <- extract_effect(m_p_srpe_adj, "duration_within_z")

cat("Within-athlete (wellbeing): β =", round(eff_p_srpe_adj["beta"], 3), 
    ", 95% CI [", round(eff_p_srpe_adj["CI_lower"], 3), ",", round(eff_p_srpe_adj["CI_upper"], 3), "]",
    ", p =", round(eff_p_srpe_adj["p"], 3), "\n")
cat("Within-athlete (duration): β =", round(eff_p_dur["beta"], 3), 
    ", p =", ifelse(eff_p_dur["p"] < 0.001, "< .001", round(eff_p_dur["p"], 3)), "\n")

# Calculate attenuation
attenuation_p <- (eff_p_srpe["beta"] - eff_p_srpe_adj["beta"]) / eff_p_srpe["beta"] * 100
cat("\nAttenuation after including duration:", round(attenuation_p, 1), "%\n")

# --- C. Supplementary: Does wellbeing predict duration? ---
cat("\n=== C. Wellbeing → Duration (Supplementary) ===\n")

m_p_dur <- safe_lmer(
  exposure_mins ~ wellbeing_within_z + wellbeing_between_z + (1 | athlete_name),
  data = practices
)

eff_p_wb_dur <- extract_effect(m_p_dur, "wellbeing_within_z")

cat("Within-athlete: β =", round(eff_p_wb_dur["beta"], 3), 
    ", 95% CI [", round(eff_p_wb_dur["CI_lower"], 3), ",", round(eff_p_wb_dur["CI_upper"], 3), "]",
    ", p =", round(eff_p_wb_dur["p"], 3), "\n")

# --- Practice Summary ---
cat("\n=== PRACTICES SUMMARY ===\n")
cat("Unadjusted wellbeing-sRPE association: β =", round(eff_p_srpe["beta"], 3), "\n")
cat("Adjusted wellbeing-sRPE association:   β =", round(eff_p_srpe_adj["beta"], 3), "\n")
cat("Attenuation:", round(attenuation_p, 1), "%\n")
if (abs(attenuation_p) > 50) {
  cat("--> Substantial attenuation: duration accounts for much of the shared variance.\n")
} else if (eff_p_srpe_adj["p"] < 0.05) {
  cat("--> Association persists: wellbeing captures aspects of perceived exertion\n")
  cat("    not reducible to session duration.\n")
}

# =============================================================================
# 8. WELLBEING-sRPE ASSOCIATIONS - GAMES
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# GAMES: Wellbeing and sRPE Associations                                  #\n")
cat("###########################################################################\n")

# --- A. Unadjusted association ---
cat("\n=== A. Wellbeing → sRPE (Unadjusted) ===\n")

m_g_srpe <- safe_lmer(
  sRPE ~ wellbeing_within_z + wellbeing_between_z + (1 | athlete_name),
  data = games
)

eff_g_srpe <- extract_effect(m_g_srpe, "wellbeing_within_z")
eff_g_srpe_bw <- extract_effect(m_g_srpe, "wellbeing_between_z")

cat("Within-athlete: β =", round(eff_g_srpe["beta"], 3), 
    ", 95% CI [", round(eff_g_srpe["CI_lower"], 3), ",", round(eff_g_srpe["CI_upper"], 3), "]",
    ", p =", round(eff_g_srpe["p"], 3), "\n")
cat("Between-athlete: β =", round(eff_g_srpe_bw["beta"], 3),
    ", p =", round(eff_g_srpe_bw["p"], 3), "\n")

# --- B. Adjusted association ---
cat("\n=== B. Wellbeing → sRPE (Adjusted for Time-on-Ice) ===\n")

m_g_srpe_adj <- safe_lmer(
  sRPE ~ wellbeing_within_z + wellbeing_between_z + 
    toi_within_z + toi_between_z + (1 | athlete_name),
  data = games
)

eff_g_srpe_adj <- extract_effect(m_g_srpe_adj, "wellbeing_within_z")
eff_g_toi <- extract_effect(m_g_srpe_adj, "toi_within_z")

cat("Within-athlete (wellbeing): β =", round(eff_g_srpe_adj["beta"], 3), 
    ", 95% CI [", round(eff_g_srpe_adj["CI_lower"], 3), ",", round(eff_g_srpe_adj["CI_upper"], 3), "]",
    ", p =", round(eff_g_srpe_adj["p"], 3), "\n")
cat("Within-athlete (TOI): β =", round(eff_g_toi["beta"], 3), 
    ", p =", ifelse(eff_g_toi["p"] < 0.001, "< .001", round(eff_g_toi["p"], 3)), "\n")

# Calculate attenuation
attenuation_g <- (eff_g_srpe["beta"] - eff_g_srpe_adj["beta"]) / eff_g_srpe["beta"] * 100
cat("\nAttenuation after including TOI:", round(attenuation_g, 1), "%\n")

# --- C. Supplementary: Does wellbeing predict TOI? ---
cat("\n=== C. Wellbeing → Time-on-Ice (Supplementary) ===\n")

m_g_toi <- safe_lmer(
  exposure_mins ~ wellbeing_within_z + wellbeing_between_z + (1 | athlete_name),
  data = games
)

eff_g_wb_toi <- extract_effect(m_g_toi, "wellbeing_within_z")

cat("Within-athlete: β =", round(eff_g_wb_toi["beta"], 3), 
    ", 95% CI [", round(eff_g_wb_toi["CI_lower"], 3), ",", round(eff_g_wb_toi["CI_upper"], 3), "]",
    ", p =", round(eff_g_wb_toi["p"], 3), "\n")

# --- Games Summary ---
cat("\n=== GAMES SUMMARY ===\n")
cat("Unadjusted wellbeing-sRPE association: β =", round(eff_g_srpe["beta"], 3), "\n")
cat("Adjusted wellbeing-sRPE association:   β =", round(eff_g_srpe_adj["beta"], 3), "\n")
cat("Attenuation:", round(attenuation_g, 1), "%\n")
cat("Note: In games, TOI is largely determined by coaching decisions and game situation,\n")
cat("      not athlete wellbeing. Persistence of the wellbeing-sRPE association after\n")
cat("      controlling for TOI provides stronger evidence that wellbeing relates to\n")
cat("      how athletes perceive exertion independent of exposure.\n")

# =============================================================================
# 9. PERFORMANCE RATINGS (Games Only)
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# PERFORMANCE RATINGS                                                     #\n")
cat("###########################################################################\n")

games_perf <- games %>%
  filter(!is.na(self_performance) | !is.na(coach_performance)) %>%
  mutate(
    self_perf_ord = factor(self_performance, ordered = TRUE),
    coach_perf_ord = factor(coach_performance, ordered = TRUE)
  )

cat("\nGames with performance data:", nrow(games_perf), "\n")

# Check category distributions
cat("\nSelf-Performance Distribution:\n")
print(table(games_perf$self_performance, useNA = "ifany"))

cat("\nCoach-Performance Distribution:\n")
print(table(games_perf$coach_performance, useNA = "ifany"))

# Self-rated performance
cat("\n=== Self-Rated Performance (CLMM) ===\n")
m_perf_self <- tryCatch(
  clmm(self_perf_ord ~ wellbeing_within_z + wellbeing_between_z + (1 | athlete_name),
       data = games_perf %>% filter(!is.na(self_performance))),
  error = function(e) {
    cat("Model failed:", e$message, "\n")
    NULL
  }
)

if (!is.null(m_perf_self)) {
  s <- summary(m_perf_self)
  coefs <- s$coefficients
  wb_row <- coefs["wellbeing_within_z", ]
  cat("Wellbeing → Self Performance:\n")
  cat("  β =", round(wb_row["Estimate"], 3),
      ", SE =", round(wb_row["Std. Error"], 3),
      ", OR =", round(exp(wb_row["Estimate"]), 2),
      ", p =", round(wb_row["Pr(>|z|)"], 3), "\n")
  cat("  Interpretation: OR =", round(exp(wb_row["Estimate"]), 2), 
      "means a 1 SD increase in within-athlete wellbeing is associated with\n")
  cat("  ", round(exp(wb_row["Estimate"]), 2), "times the odds of a higher performance rating.\n")
}

# Coach-rated performance
cat("\n=== Coach-Rated Performance (CLMM) ===\n")
m_perf_coach <- tryCatch(
  clmm(coach_perf_ord ~ wellbeing_within_z + wellbeing_between_z + (1 | athlete_name),
       data = games_perf %>% filter(!is.na(coach_performance))),
  error = function(e) {
    cat("Model failed:", e$message, "\n")
    NULL
  }
)

if (!is.null(m_perf_coach)) {
  s <- summary(m_perf_coach)
  coefs <- s$coefficients
  wb_row <- coefs["wellbeing_within_z", ]
  cat("Wellbeing → Coach Performance:\n")
  cat("  β =", round(wb_row["Estimate"], 3),
      ", SE =", round(wb_row["Std. Error"], 3),
      ", OR =", round(exp(wb_row["Estimate"]), 2),
      ", p =", round(wb_row["Pr(>|z|)"], 3), "\n")
}

# =============================================================================
# 10. ITEM-LEVEL EXPLORATORY ANALYSES
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# EXPLORATORY: Item-Level Effects                                         #\n")
cat("###########################################################################\n")

# Prepare item-level within-athlete deviations for practices
practices_items <- practices %>%
  group_by(athlete_name) %>%
  mutate(
    confidence_within = confidence - mean(confidence, na.rm = TRUE),
    excitement_within = excitement - mean(excitement, na.rm = TRUE),
    fatigue_within = fatigue - mean(fatigue, na.rm = TRUE),  # Note: NOT reversed for item-level
    sleep_within = sleep_quality - mean(sleep_quality, na.rm = TRUE),
    stress_within = stress - mean(stress, na.rm = TRUE)  # Note: NOT reversed for item-level
  ) %>%
  ungroup() %>%
  mutate(
    confidence_within_z = as.numeric(scale(confidence_within)),
    excitement_within_z = as.numeric(scale(excitement_within)),
    fatigue_within_z = as.numeric(scale(fatigue_within)),
    sleep_within_z = as.numeric(scale(sleep_within)),
    stress_within_z = as.numeric(scale(stress_within))
  )

items <- c("confidence", "excitement", "fatigue", "sleep", "stress")
item_results <- list()

for (item in items) {
  within_var <- paste0(item, "_within_z")
  
  # sRPE model (unadjusted)
  f1 <- as.formula(paste0("sRPE ~ ", within_var, " + (1 | athlete_name)"))
  m1 <- safe_lmer(f1, data = practices_items)
  
  # sRPE model (adjusted for duration)
  f2 <- as.formula(paste0("sRPE ~ ", within_var, " + duration_within_z + (1 | athlete_name)"))
  m2 <- safe_lmer(f2, data = practices_items)
  
  # Duration model
  f3 <- as.formula(paste0("exposure_mins ~ ", within_var, " + (1 | athlete_name)"))
  m3 <- safe_lmer(f3, data = practices_items)
  
  eff1 <- extract_effect(m1, within_var)
  eff2 <- extract_effect(m2, within_var)
  eff3 <- extract_effect(m3, within_var)
  
  item_results[[item]] <- data.frame(
    item = item,
    srpe_unadj_beta = eff1["beta"],
    srpe_unadj_p = eff1["p"],
    srpe_adj_beta = eff2["beta"],
    srpe_adj_p = eff2["p"],
    duration_beta = eff3["beta"],
    duration_p = eff3["p"]
  )
}

results_item_level <- bind_rows(item_results)
rownames(results_item_level) <- NULL

cat("\nItem-Level Effects (Practices):\n")
cat("Note: Fatigue and stress are in ORIGINAL direction (higher = worse)\n")
cat("Columns: sRPE (unadj), sRPE (adj for duration), Duration\n\n")
print(results_item_level %>% mutate(across(where(is.numeric), ~round(., 3))))

# =============================================================================
# 11. SAVE RESULTS
# =============================================================================
cat("\n=== SAVING RESULTS ===\n")

# Main results - wellbeing-sRPE associations
results_associations <- data.frame(
  context = c("Practice", "Practice", "Practice", "Game", "Game", "Game"),
  model = c("Wellbeing→sRPE (unadjusted)", 
            "Wellbeing→sRPE (adjusted for duration)",
            "Wellbeing→Duration",
            "Wellbeing→sRPE (unadjusted)",
            "Wellbeing→sRPE (adjusted for TOI)",
            "Wellbeing→TOI"),
  beta = c(eff_p_srpe["beta"], eff_p_srpe_adj["beta"], eff_p_wb_dur["beta"],
           eff_g_srpe["beta"], eff_g_srpe_adj["beta"], eff_g_wb_toi["beta"]),
  SE = c(eff_p_srpe["SE"], eff_p_srpe_adj["SE"], eff_p_wb_dur["SE"],
         eff_g_srpe["SE"], eff_g_srpe_adj["SE"], eff_g_wb_toi["SE"]),
  CI_lower = c(eff_p_srpe["CI_lower"], eff_p_srpe_adj["CI_lower"], eff_p_wb_dur["CI_lower"],
               eff_g_srpe["CI_lower"], eff_g_srpe_adj["CI_lower"], eff_g_wb_toi["CI_lower"]),
  CI_upper = c(eff_p_srpe["CI_upper"], eff_p_srpe_adj["CI_upper"], eff_p_wb_dur["CI_upper"],
               eff_g_srpe["CI_upper"], eff_g_srpe_adj["CI_upper"], eff_g_wb_toi["CI_upper"]),
  p = c(eff_p_srpe["p"], eff_p_srpe_adj["p"], eff_p_wb_dur["p"],
        eff_g_srpe["p"], eff_g_srpe_adj["p"], eff_g_wb_toi["p"])
)

# Performance results
results_performance <- data.frame(
  outcome = c("Self_Performance", "Coach_Performance"),
  wellbeing_beta = c(
    if (!is.null(m_perf_self)) summary(m_perf_self)$coefficients["wellbeing_within_z", "Estimate"] else NA,
    if (!is.null(m_perf_coach)) summary(m_perf_coach)$coefficients["wellbeing_within_z", "Estimate"] else NA
  ),
  wellbeing_SE = c(
    if (!is.null(m_perf_self)) summary(m_perf_self)$coefficients["wellbeing_within_z", "Std. Error"] else NA,
    if (!is.null(m_perf_coach)) summary(m_perf_coach)$coefficients["wellbeing_within_z", "Std. Error"] else NA
  ),
  wellbeing_OR = c(
    if (!is.null(m_perf_self)) exp(summary(m_perf_self)$coefficients["wellbeing_within_z", "Estimate"]) else NA,
    if (!is.null(m_perf_coach)) exp(summary(m_perf_coach)$coefficients["wellbeing_within_z", "Estimate"]) else NA
  ),
  wellbeing_p = c(
    if (!is.null(m_perf_self)) summary(m_perf_self)$coefficients["wellbeing_within_z", "Pr(>|z|)"] else NA,
    if (!is.null(m_perf_coach)) summary(m_perf_coach)$coefficients["wellbeing_within_z", "Pr(>|z|)"] else NA
  )
)

# Summary
summary_aim6 <- data.frame(
  finding = c(
    "N_practices", "N_games", "N_athletes_practices", "N_athletes_games",
    "Practice: Wellbeing→sRPE (unadj) β",
    "Practice: Wellbeing→sRPE (unadj) p",
    "Practice: Wellbeing→sRPE (adj) β",
    "Practice: Wellbeing→sRPE (adj) p",
    "Practice: Attenuation %",
    "Practice: Wellbeing→Duration p",
    "Game: Wellbeing→sRPE (unadj) β",
    "Game: Wellbeing→sRPE (unadj) p",
    "Game: Wellbeing→sRPE (adj) β",
    "Game: Wellbeing→sRPE (adj) p", 
    "Game: Attenuation %",
    "Game: Wellbeing→TOI p",
    "Self_Performance OR",
    "Self_Performance p",
    "Coach_Performance OR",
    "Coach_Performance p"
  ),
  value = c(
    nrow(practices), nrow(games), 
    n_distinct(practices$athlete_name), n_distinct(games$athlete_name),
    eff_p_srpe["beta"], eff_p_srpe["p"],
    eff_p_srpe_adj["beta"], eff_p_srpe_adj["p"],
    attenuation_p,
    eff_p_wb_dur["p"],
    eff_g_srpe["beta"], eff_g_srpe["p"],
    eff_g_srpe_adj["beta"], eff_g_srpe_adj["p"],
    attenuation_g,
    eff_g_wb_toi["p"],
    results_performance$wellbeing_OR[1],
    results_performance$wellbeing_p[1],
    results_performance$wellbeing_OR[2],
    results_performance$wellbeing_p[2]
  )
)

# Save all
write.csv(desc_by_context, "descriptives_by_context.csv", row.names = FALSE)
write.csv(results_associations, "results_associations.csv", row.names = FALSE)
write.csv(results_performance, "results_performance.csv", row.names = FALSE)
write.csv(results_item_level, "results_item_level.csv", row.names = FALSE)
write.csv(summary_aim6, "summary_aim6.csv", row.names = FALSE)

save(m_p_srpe, m_p_srpe_adj, m_p_dur,
     m_g_srpe, m_g_srpe_adj, m_g_toi,
     m_perf_self, m_perf_coach,
     file = "aim6_models.RData")

cat("\nFiles saved:\n")
cat("  - descriptives_by_context.csv\n")
cat("  - results_associations.csv\n")
cat("  - results_performance.csv\n")
cat("  - results_item_level.csv\n")
cat("  - summary_aim6.csv\n")
cat("  - aim6_models.RData\n")

# =============================================================================
# 12. FINAL SUMMARY
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# AIM 6 FINAL SUMMARY                                                     #\n")
cat("###########################################################################\n")

cat("\n=== PRACTICES ===\n")
cat("Wellbeing-sRPE association:\n")
cat("  Unadjusted: β =", round(eff_p_srpe["beta"], 3), 
    ", 95% CI [", round(eff_p_srpe["CI_lower"], 3), ",", round(eff_p_srpe["CI_upper"], 3), "]",
    ", p =", round(eff_p_srpe["p"], 3), "\n")
cat("  Adjusted:   β =", round(eff_p_srpe_adj["beta"], 3), 
    ", 95% CI [", round(eff_p_srpe_adj["CI_lower"], 3), ",", round(eff_p_srpe_adj["CI_upper"], 3), "]",
    ", p =", round(eff_p_srpe_adj["p"], 3), "\n")
cat("  Attenuation:", round(attenuation_p, 1), "%\n")
cat("Wellbeing-Duration association: β =", round(eff_p_wb_dur["beta"], 3), 
    ", p =", round(eff_p_wb_dur["p"], 3), "\n")

cat("\n=== GAMES ===\n")
cat("Wellbeing-sRPE association:\n")
cat("  Unadjusted: β =", round(eff_g_srpe["beta"], 3), 
    ", 95% CI [", round(eff_g_srpe["CI_lower"], 3), ",", round(eff_g_srpe["CI_upper"], 3), "]",
    ", p =", round(eff_g_srpe["p"], 3), "\n")
cat("  Adjusted:   β =", round(eff_g_srpe_adj["beta"], 3), 
    ", 95% CI [", round(eff_g_srpe_adj["CI_lower"], 3), ",", round(eff_g_srpe_adj["CI_upper"], 3), "]",
    ", p =", round(eff_g_srpe_adj["p"], 3), "\n")
cat("  Attenuation:", round(attenuation_g, 1), "%\n")
cat("Wellbeing-TOI association: β =", round(eff_g_wb_toi["beta"], 3), 
    ", p =", round(eff_g_wb_toi["p"], 3), "\n")

cat("\n=== PERFORMANCE ===\n")
if (!is.null(m_perf_self)) {
  cat("Self-rated: OR =", round(results_performance$wellbeing_OR[1], 2),
      ", p =", round(results_performance$wellbeing_p[1], 3), "\n")
}
if (!is.null(m_perf_coach)) {
  cat("Coach-rated: OR =", round(results_performance$wellbeing_OR[2], 2),
      ", p =", round(results_performance$wellbeing_p[2], 3), "\n")
}

cat("\n=== INTERPRETATION ===\n")
if (eff_p_srpe["p"] < 0.05 || eff_g_srpe["p"] < 0.05) {
  cat("Morning wellbeing predicts same-day sRPE.\n")
  if (abs(attenuation_p) > 50 && abs(attenuation_g) > 50) {
    cat("This association is substantially attenuated after accounting for duration/TOI,\n")
    cat("indicating that duration accounts for much of the shared variance between\n")
    cat("wellbeing and sRPE.\n")
  } else if (eff_p_srpe_adj["p"] < 0.05 || eff_g_srpe_adj["p"] < 0.05) {
    cat("This association persists after accounting for duration/TOI, indicating that\n")
    cat("wellbeing captures aspects of perceived exertion not reducible to session exposure.\n")
  }
} else {
  cat("Morning wellbeing does not significantly predict same-day sRPE in this sample.\n")
}

cat("\n=============================================================================\n")
cat("AIM 6 ANALYSIS COMPLETE\n")
cat("=============================================================================\n")