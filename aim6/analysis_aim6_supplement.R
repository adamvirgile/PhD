# =============================================================================
# AIM 6 SUPPLEMENTARY: Duration/TOI Control - Full Analysis
# =============================================================================
# Checking wellbeing → sRPE relationship controlling for:
#   - Practices: session duration
#   - Games: TOI (time-on-ice)
# =============================================================================

setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim6")

library(lme4)
library(lmerTest)
library(performance)
library(dplyr)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================
cat("=============================================================================\n")
cat("AIM 6: Duration/TOI Control - Full Analysis\n")
cat("=============================================================================\n\n")

data <- read.csv("../merged_data_full.csv")
data$session_date <- as.Date(data$session_date)

# Prepare wellbeing composite
analysis_data <- data %>%
  filter(!is.na(sRPE),
         !is.na(confidence),
         !is.na(excitement),
         !is.na(fatigue),
         !is.na(sleep_quality),
         !is.na(stress),
         full_participation == 1) %>%
  mutate(
    fatigue_r = 10 - fatigue,
    stress_r = 10 - stress,
    confidence_z = scale(confidence)[,1],
    excitement_z = scale(excitement)[,1],
    fatigue_r_z = scale(fatigue_r)[,1],
    sleep_quality_z = scale(sleep_quality)[,1],
    stress_r_z = scale(stress_r)[,1],
    wellbeing_composite = (confidence_z + excitement_z + fatigue_r_z + 
                             sleep_quality_z + stress_r_z) / 5
  )

# =============================================================================
# 2. PRACTICES - Duration Control
# =============================================================================
cat("=== PRACTICES ===\n")

practices <- analysis_data %>% 
  filter(context == "practice") %>%
  group_by(athlete_name) %>%
  mutate(
    wellbeing_mean = mean(wellbeing_composite, na.rm = TRUE),
    wellbeing_within = wellbeing_composite - wellbeing_mean,
    duration_mean = mean(duration_mins, na.rm = TRUE),
    duration_within = duration_mins - duration_mean
  ) %>%
  ungroup() %>%
  mutate(
    wellbeing_within_z = scale(wellbeing_within)[,1],
    wellbeing_between_z = scale(wellbeing_mean)[,1],
    duration_within_z = scale(duration_within)[,1],
    duration_between_z = scale(duration_mean)[,1]
  )

cat("Practice observations:", nrow(practices), "\n")

# Without duration
m_practice_wb <- lmer(
  sRPE ~ wellbeing_within_z + wellbeing_between_z + (1 | athlete_name),
  data = practices
)

# With duration
m_practice_both <- lmer(
  sRPE ~ wellbeing_within_z + wellbeing_between_z + 
    duration_within_z + duration_between_z + (1 | athlete_name),
  data = practices
)

fe_p1 <- summary(m_practice_wb)$coefficients
fe_p2 <- summary(m_practice_both)$coefficients

cat("\nPractices - Wellbeing effect:\n")
cat("  Without duration: β =", round(fe_p1["wellbeing_within_z", "Estimate"], 4),
    ", p =", round(fe_p1["wellbeing_within_z", "Pr(>|t|)"], 4), "\n")
cat("  With duration:    β =", round(fe_p2["wellbeing_within_z", "Estimate"], 4),
    ", p =", round(fe_p2["wellbeing_within_z", "Pr(>|t|)"], 4), "\n")
cat("  Duration effect:  β =", round(fe_p2["duration_within_z", "Estimate"], 4),
    ", p < .001\n")

# =============================================================================
# 3. GAMES - TOI Control
# =============================================================================
cat("\n=== GAMES ===\n")

games <- analysis_data %>% 
  filter(context == "game", !is.na(toi_total_mins)) %>%
  group_by(athlete_name) %>%
  mutate(
    wellbeing_mean = mean(wellbeing_composite, na.rm = TRUE),
    wellbeing_within = wellbeing_composite - wellbeing_mean,
    toi_mean = mean(toi_total_mins, na.rm = TRUE),
    toi_within = toi_total_mins - toi_mean
  ) %>%
  ungroup() %>%
  mutate(
    wellbeing_within_z = scale(wellbeing_within)[,1],
    wellbeing_between_z = scale(wellbeing_mean)[,1],
    toi_within_z = scale(toi_within)[,1],
    toi_between_z = scale(toi_mean)[,1]
  )

cat("Game observations with TOI:", nrow(games), "\n")

# Check wellbeing-TOI correlation
cor_wb_toi <- cor(games$wellbeing_within, games$toi_within, use = "complete.obs")
cat("Wellbeing-TOI within-athlete correlation:", round(cor_wb_toi, 3), "\n")

# Without TOI
m_games_wb <- lmer(
  sRPE ~ wellbeing_within_z + wellbeing_between_z + (1 | athlete_name),
  data = games
)

# With TOI
m_games_both <- lmer(
  sRPE ~ wellbeing_within_z + wellbeing_between_z + 
    toi_within_z + toi_between_z + (1 | athlete_name),
  data = games
)

fe_g1 <- summary(m_games_wb)$coefficients
fe_g2 <- summary(m_games_both)$coefficients

cat("\nGames - Wellbeing effect:\n")
cat("  Without TOI: β =", round(fe_g1["wellbeing_within_z", "Estimate"], 4),
    ", p =", round(fe_g1["wellbeing_within_z", "Pr(>|t|)"], 4), "\n")
cat("  With TOI:    β =", round(fe_g2["wellbeing_within_z", "Estimate"], 4),
    ", p =", round(fe_g2["wellbeing_within_z", "Pr(>|t|)"], 4), "\n")
cat("  TOI effect:  β =", round(fe_g2["toi_within_z", "Estimate"], 4),
    ", p =", round(fe_g2["toi_within_z", "Pr(>|t|)"], 4), "\n")

# =============================================================================
# 4. COMBINED MODEL - Context-Specific Duration Control
# =============================================================================
cat("\n=== COMBINED MODEL ===\n")

# Create unified duration variable (duration for practices, TOI for games)
combined <- analysis_data %>%
  mutate(
    exposure_mins = case_when(
      context == "practice" ~ duration_mins,
      context == "game" ~ toi_total_mins,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(exposure_mins)) %>%
  group_by(athlete_name) %>%
  mutate(
    wellbeing_mean = mean(wellbeing_composite, na.rm = TRUE),
    wellbeing_within = wellbeing_composite - wellbeing_mean,
    exposure_mean = mean(exposure_mins, na.rm = TRUE),
    exposure_within = exposure_mins - exposure_mean
  ) %>%
  ungroup() %>%
  mutate(
    wellbeing_within_z = scale(wellbeing_within)[,1],
    wellbeing_between_z = scale(wellbeing_mean)[,1],
    exposure_within_z = scale(exposure_within)[,1],
    exposure_between_z = scale(exposure_mean)[,1],
    context_f = factor(context, levels = c("practice", "game"))
  )

cat("Combined observations:", nrow(combined), "\n")

# Model with context and exposure
m_combined <- lmer(
  sRPE ~ wellbeing_within_z * context_f + wellbeing_between_z + 
    exposure_within_z * context_f + exposure_between_z +
    (1 | athlete_name),
  data = combined
)

cat("\n--- Combined Model Summary ---\n")
print(summary(m_combined))

# =============================================================================
# 5. SUMMARY TABLE
# =============================================================================
cat("\n=== SUMMARY TABLE ===\n")

results_full <- data.frame(
  context = c("Practice", "Practice", "Game", "Game"),
  model = c("Wellbeing_only", "With_Duration", "Wellbeing_only", "With_TOI"),
  wellbeing_beta = c(
    fe_p1["wellbeing_within_z", "Estimate"],
    fe_p2["wellbeing_within_z", "Estimate"],
    fe_g1["wellbeing_within_z", "Estimate"],
    fe_g2["wellbeing_within_z", "Estimate"]
  ),
  wellbeing_p = c(
    fe_p1["wellbeing_within_z", "Pr(>|t|)"],
    fe_p2["wellbeing_within_z", "Pr(>|t|)"],
    fe_g1["wellbeing_within_z", "Pr(>|t|)"],
    fe_g2["wellbeing_within_z", "Pr(>|t|)"]
  ),
  duration_toi_beta = c(
    NA,
    fe_p2["duration_within_z", "Estimate"],
    NA,
    fe_g2["toi_within_z", "Estimate"]
  ),
  duration_toi_p = c(
    NA,
    fe_p2["duration_within_z", "Pr(>|t|)"],
    NA,
    fe_g2["toi_within_z", "Pr(>|t|)"]
  )
)

cat("\nResults by Context:\n")
print(results_full %>% mutate(across(where(is.numeric), ~round(., 4))))

write.csv(results_full, "results_duration_toi_control.csv", row.names = FALSE)

# Save models
save(m_practice_wb, m_practice_both, m_games_wb, m_games_both, m_combined,
     file = "aim6_duration_toi_models.RData")

cat("\nFiles saved:\n")
cat("  - results_duration_toi_control.csv\n")
cat("  - aim6_duration_toi_models.RData\n")

# =============================================================================
# 6. FINAL SUMMARY
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# DURATION/TOI CONTROL - FINAL SUMMARY                                    #\n")
cat("###########################################################################\n")

cat("\n=== PRACTICES ===\n")
cat("Wellbeing → sRPE:\n")
cat("  Without duration control: β =", round(fe_p1["wellbeing_within_z", "Estimate"], 3),
    ", p =", round(fe_p1["wellbeing_within_z", "Pr(>|t|)"], 3), "\n")
cat("  With duration control:    β =", round(fe_p2["wellbeing_within_z", "Estimate"], 3),
    ", p =", round(fe_p2["wellbeing_within_z", "Pr(>|t|)"], 3), "\n")

cat("\n=== GAMES ===\n")
cat("Wellbeing-TOI correlation: r =", round(cor_wb_toi, 3), "\n")
cat("Wellbeing → sRPE:\n")
cat("  Without TOI control: β =", round(fe_g1["wellbeing_within_z", "Estimate"], 3),
    ", p =", round(fe_g1["wellbeing_within_z", "Pr(>|t|)"], 3), "\n")
cat("  With TOI control:    β =", round(fe_g2["wellbeing_within_z", "Estimate"], 3),
    ", p =", round(fe_g2["wellbeing_within_z", "Pr(>|t|)"], 3), "\n")

cat("\n=== CONCLUSION ===\n")
if (fe_g2["wellbeing_within_z", "Pr(>|t|)"] < 0.05) {
  cat("Wellbeing effect PERSISTS in games after controlling for TOI.\n")
} else {
  cat("Wellbeing effect is NOT SIGNIFICANT in games after controlling for TOI.\n")
}

cat("\n=============================================================================\n")