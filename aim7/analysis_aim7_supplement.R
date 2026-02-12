# =============================================================================
# AIM 7 SUPPLEMENTARY: Decomposing sRPE-TL Value
# =============================================================================
# Location: C:\Users\avirgile\Dropbox\Adam Only\School\PhD\UVM Coursework\Dissertation\data\aim7
#
# Key Question:
#   sRPE-TL predicts TRIMP better than duration alone.
#   But sRPE is contaminated by duration (β ≈ 0.68).
#   
#   So is sRPE-TL's improvement due to:
#   (A) sRPE capturing real intensity? (good)
#   (B) Duration being counted twice? (redundant)
#
# Method:
#   1. Residualize sRPE for duration → sRPE_intensity (pure intensity signal)
#   2. Create sRPE_intensity_TL = sRPE_intensity × Duration
#   3. Compare predictive validity:
#      - Duration alone → TRIMP
#      - sRPE-TL → TRIMP (includes duration contamination)
#      - sRPE_intensity_TL → TRIMP (intensity component only)
#
#   If sRPE_intensity_TL still improves prediction, the value comes from
#   real intensity information, not double-counting duration.
# =============================================================================

setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim7")

library(lme4)
library(lmerTest)
library(performance)
library(dplyr)

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("=============================================================================\n")
cat("AIM 7: Decomposing sRPE-TL — Intensity vs Duration Components\n")
cat("=============================================================================\n\n")

data <- read.csv("../aim1/aim1_analysis_data.csv")

# =============================================================================
# 2. PRACTICES
# =============================================================================
cat("###########################################################################\n")
cat("# PRACTICES                                                               #\n")
cat("###########################################################################\n")

practices <- data %>%
  filter(context == "practice") %>%
  filter(!is.na(sRPE), !is.na(duration_mins), !is.na(TRIMP), !is.na(HR80_mins))

cat("\nN =", nrow(practices), "observations\n")

# --- Step 1: Residualize sRPE for duration ---
cat("\n=== Step 1: Extract intensity component from sRPE ===\n")

# Model: sRPE ~ Duration (to get residuals)
m_srpe_dur_p <- lmer(
  sRPE ~ duration_mins + (1 | athlete_name),
  data = practices
)

# sRPE_intensity = residuals (what's left after removing duration effect)
practices$sRPE_intensity <- residuals(m_srpe_dur_p)

cat("sRPE decomposed into:\n")
cat("  - Duration component (removed)\n")
cat("  - Intensity component (residuals): SD =", round(sd(practices$sRPE_intensity), 3), "\n")

# Verify: sRPE_intensity should NOT correlate with duration
cor_check <- cor(practices$sRPE_intensity, practices$duration_mins)
cat("  - Correlation(sRPE_intensity, duration) =", round(cor_check, 4), "(should be ~0)\n")

# --- Step 2: Create load metrics ---
cat("\n=== Step 2: Create load metrics ===\n")

practices <- practices %>%
  mutate(
    # Original sRPE-TL (includes duration contamination)
    sRPE_TL = sRPE * duration_mins,
    
    # Intensity-only sRPE-TL (duration effect removed from sRPE)
    sRPE_intensity_TL = sRPE_intensity * duration_mins
  )

cat("Created:\n")
cat("  - sRPE_TL = sRPE × Duration (original, includes duration contamination)\n")
cat("  - sRPE_intensity_TL = sRPE_intensity × Duration (intensity component only)\n")

# --- Step 3: Within/between decomposition ---
practices <- practices %>%
  group_by(athlete_name) %>%
  mutate(
    duration_within = duration_mins - mean(duration_mins),
    duration_between = mean(duration_mins),
    sRPE_TL_within = sRPE_TL - mean(sRPE_TL),
    sRPE_TL_between = mean(sRPE_TL),
    sRPE_int_TL_within = sRPE_intensity_TL - mean(sRPE_intensity_TL),
    sRPE_int_TL_between = mean(sRPE_intensity_TL)
  ) %>%
  ungroup() %>%
  mutate(
    duration_within_z = scale(duration_within)[,1],
    duration_between_z = scale(duration_between)[,1],
    sRPE_TL_within_z = scale(sRPE_TL_within)[,1],
    sRPE_TL_between_z = scale(sRPE_TL_between)[,1],
    sRPE_int_TL_within_z = scale(sRPE_int_TL_within)[,1],
    sRPE_int_TL_between_z = scale(sRPE_int_TL_between)[,1]
  )

# --- Step 4: Compare predictive validity for TRIMP ---
cat("\n=== Step 3: Compare predictive validity for TRIMP ===\n")

# Model 1: Duration alone
m1_p <- lmer(TRIMP ~ duration_within_z + duration_between_z + (1 | athlete_name),
             data = practices)

# Model 2: sRPE-TL (original, includes duration contamination)
m2_p <- lmer(TRIMP ~ sRPE_TL_within_z + sRPE_TL_between_z + (1 | athlete_name),
             data = practices)

# Model 3: sRPE_intensity_TL (intensity component only)
m3_p <- lmer(TRIMP ~ sRPE_int_TL_within_z + sRPE_int_TL_between_z + (1 | athlete_name),
             data = practices)

r2_1_p <- performance::r2(m1_p)$R2_marginal
r2_2_p <- performance::r2(m2_p)$R2_marginal
r2_3_p <- performance::r2(m3_p)$R2_marginal

cat("\nPredicting TRIMP:\n")
cat("  Duration alone:        R² =", round(r2_1_p, 4), "\n")
cat("  sRPE-TL (original):    R² =", round(r2_2_p, 4), " (Δ =", round(r2_2_p - r2_1_p, 4), ")\n")
cat("  sRPE_intensity_TL:     R² =", round(r2_3_p, 4), " (Δ =", round(r2_3_p - r2_1_p, 4), ")\n")

# --- Step 5: Compare predictive validity for HR80 ---
cat("\n=== Step 4: Compare predictive validity for HR80 ===\n")

m1_p_hr <- lmer(HR80_mins ~ duration_within_z + duration_between_z + (1 | athlete_name),
                data = practices)

m2_p_hr <- lmer(HR80_mins ~ sRPE_TL_within_z + sRPE_TL_between_z + (1 | athlete_name),
                data = practices)

m3_p_hr <- lmer(HR80_mins ~ sRPE_int_TL_within_z + sRPE_int_TL_between_z + (1 | athlete_name),
                data = practices)

r2_1_p_hr <- performance::r2(m1_p_hr)$R2_marginal
r2_2_p_hr <- performance::r2(m2_p_hr)$R2_marginal
r2_3_p_hr <- performance::r2(m3_p_hr)$R2_marginal

cat("\nPredicting HR80:\n")
cat("  Duration alone:        R² =", round(r2_1_p_hr, 4), "\n")
cat("  sRPE-TL (original):    R² =", round(r2_2_p_hr, 4), " (Δ =", round(r2_2_p_hr - r2_1_p_hr, 4), ")\n")
cat("  sRPE_intensity_TL:     R² =", round(r2_3_p_hr, 4), " (Δ =", round(r2_3_p_hr - r2_1_p_hr, 4), ")\n")

# =============================================================================
# 3. GAMES
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# GAMES                                                                   #\n")
cat("###########################################################################\n")

games <- data %>%
  filter(context == "game", !is.na(toi_total_mins)) %>%
  filter(!is.na(sRPE), !is.na(TRIMP), !is.na(HR80_mins))

cat("\nN =", nrow(games), "observations\n")

# --- Step 1: Residualize sRPE for TOI ---
cat("\n=== Step 1: Extract intensity component from sRPE ===\n")

m_srpe_toi_g <- lmer(
  sRPE ~ toi_total_mins + (1 | athlete_name),
  data = games
)

games$sRPE_intensity <- residuals(m_srpe_toi_g)

cor_check_g <- cor(games$sRPE_intensity, games$toi_total_mins)
cat("Correlation(sRPE_intensity, TOI) =", round(cor_check_g, 4), "(should be ~0)\n")

# --- Step 2: Create load metrics ---
games <- games %>%
  mutate(
    sRPE_TL = sRPE * toi_total_mins,
    sRPE_intensity_TL = sRPE_intensity * toi_total_mins
  )

# --- Step 3: Within/between decomposition ---
games <- games %>%
  group_by(athlete_name) %>%
  mutate(
    toi_within = toi_total_mins - mean(toi_total_mins),
    toi_between = mean(toi_total_mins),
    sRPE_TL_within = sRPE_TL - mean(sRPE_TL),
    sRPE_TL_between = mean(sRPE_TL),
    sRPE_int_TL_within = sRPE_intensity_TL - mean(sRPE_intensity_TL),
    sRPE_int_TL_between = mean(sRPE_intensity_TL)
  ) %>%
  ungroup() %>%
  mutate(
    toi_within_z = scale(toi_within)[,1],
    toi_between_z = scale(toi_between)[,1],
    sRPE_TL_within_z = scale(sRPE_TL_within)[,1],
    sRPE_TL_between_z = scale(sRPE_TL_between)[,1],
    sRPE_int_TL_within_z = scale(sRPE_int_TL_within)[,1],
    sRPE_int_TL_between_z = scale(sRPE_int_TL_between)[,1]
  )

# --- Compare predictive validity for TRIMP ---
cat("\n=== Step 3: Compare predictive validity for TRIMP ===\n")

m1_g <- lmer(TRIMP ~ toi_within_z + toi_between_z + (1 | athlete_name),
             data = games)

m2_g <- lmer(TRIMP ~ sRPE_TL_within_z + sRPE_TL_between_z + (1 | athlete_name),
             data = games)

m3_g <- lmer(TRIMP ~ sRPE_int_TL_within_z + sRPE_int_TL_between_z + (1 | athlete_name),
             data = games)

r2_1_g <- performance::r2(m1_g)$R2_marginal
r2_2_g <- performance::r2(m2_g)$R2_marginal
r2_3_g <- performance::r2(m3_g)$R2_marginal

cat("\nPredicting TRIMP:\n")
cat("  TOI alone:             R² =", round(r2_1_g, 4), "\n")
cat("  sRPE-TL (original):    R² =", round(r2_2_g, 4), " (Δ =", round(r2_2_g - r2_1_g, 4), ")\n")
cat("  sRPE_intensity_TL:     R² =", round(r2_3_g, 4), " (Δ =", round(r2_3_g - r2_1_g, 4), ")\n")

# --- Compare predictive validity for HR80 ---
cat("\n=== Step 4: Compare predictive validity for HR80 ===\n")

m1_g_hr <- lmer(HR80_mins ~ toi_within_z + toi_between_z + (1 | athlete_name),
                data = games)

m2_g_hr <- lmer(HR80_mins ~ sRPE_TL_within_z + sRPE_TL_between_z + (1 | athlete_name),
                data = games)

m3_g_hr <- lmer(HR80_mins ~ sRPE_int_TL_within_z + sRPE_int_TL_between_z + (1 | athlete_name),
                data = games)

r2_1_g_hr <- performance::r2(m1_g_hr)$R2_marginal
r2_2_g_hr <- performance::r2(m2_g_hr)$R2_marginal
r2_3_g_hr <- performance::r2(m3_g_hr)$R2_marginal

cat("\nPredicting HR80:\n")
cat("  TOI alone:             R² =", round(r2_1_g_hr, 4), "\n")
cat("  sRPE-TL (original):    R² =", round(r2_2_g_hr, 4), " (Δ =", round(r2_2_g_hr - r2_1_g_hr, 4), ")\n")
cat("  sRPE_intensity_TL:     R² =", round(r2_3_g_hr, 4), " (Δ =", round(r2_3_g_hr - r2_1_g_hr, 4), ")\n")

# =============================================================================
# 4. SUMMARY TABLE
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# SUMMARY: Decomposing sRPE-TL Improvement                                #\n")
cat("###########################################################################\n")

results_decomp <- data.frame(
  context = c("Practice", "Practice", "Game", "Game"),
  outcome = c("TRIMP", "HR80", "TRIMP", "HR80"),
  duration_only_R2 = c(r2_1_p, r2_1_p_hr, r2_1_g, r2_1_g_hr),
  sRPE_TL_R2 = c(r2_2_p, r2_2_p_hr, r2_2_g, r2_2_g_hr),
  sRPE_intensity_TL_R2 = c(r2_3_p, r2_3_p_hr, r2_3_g, r2_3_g_hr),
  
  # Improvements
  sRPE_TL_improvement = c(r2_2_p - r2_1_p, r2_2_p_hr - r2_1_p_hr, 
                          r2_2_g - r2_1_g, r2_2_g_hr - r2_1_g_hr),
  intensity_only_improvement = c(r2_3_p - r2_1_p, r2_3_p_hr - r2_1_p_hr,
                                 r2_3_g - r2_1_g, r2_3_g_hr - r2_1_g_hr)
)

# Calculate what % of sRPE-TL's improvement comes from intensity vs duration
results_decomp <- results_decomp %>%
  mutate(
    pct_from_intensity = (intensity_only_improvement / sRPE_TL_improvement) * 100,
    pct_from_duration = 100 - pct_from_intensity
  )

cat("\n")
print(results_decomp %>% mutate(across(where(is.numeric), ~round(., 4))))

# =============================================================================
# 5. INTERPRETATION
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# INTERPRETATION                                                          #\n")
cat("###########################################################################\n")

cat("\nKey Question: When sRPE-TL improves prediction over duration alone,\n")
cat("              how much is from INTENSITY vs DURATION contamination?\n\n")

for (i in 1:nrow(results_decomp)) {
  row <- results_decomp[i, ]
  cat("---", row$context, "-", row$outcome, "---\n")
  cat("  Duration alone:        R² =", round(row$duration_only_R2, 3), "\n")
  cat("  sRPE-TL (original):    R² =", round(row$sRPE_TL_R2, 3), 
      " (+", round(row$sRPE_TL_improvement, 3), ")\n")
  cat("  sRPE_intensity_TL:     R² =", round(row$sRPE_intensity_TL_R2, 3),
      " (+", round(row$intensity_only_improvement, 3), ")\n")
  cat("  --> ", round(row$pct_from_intensity, 1), "% of improvement from INTENSITY\n")
  cat("  --> ", round(row$pct_from_duration, 1), "% of improvement from DURATION (redundant)\n\n")
}

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================
cat("\n=== SAVING RESULTS ===\n")

write.csv(results_decomp, "results_srpe_decomposition.csv", row.names = FALSE)

save(m1_p, m2_p, m3_p, m1_p_hr, m2_p_hr, m3_p_hr,
     m1_g, m2_g, m3_g, m1_g_hr, m2_g_hr, m3_g_hr,
     file = "aim7_decomposition_models.RData")

cat("Files saved:\n")
cat("  - results_srpe_decomposition.csv\n")
cat("  - aim7_decomposition_models.RData\n")

# =============================================================================
# 7. FINAL CONCLUSION
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# FINAL CONCLUSION                                                        #\n")
cat("###########################################################################\n")

avg_pct_intensity <- mean(results_decomp$pct_from_intensity, na.rm = TRUE)

cat("\nAverage across all comparisons:\n")
cat("  ", round(avg_pct_intensity, 1), "% of sRPE-TL's predictive value comes from INTENSITY\n")
cat("  ", round(100 - avg_pct_intensity, 1), "% comes from DURATION (double-counting)\n\n")

if (avg_pct_intensity > 75) {
  cat("CONCLUSION: sRPE-TL's value is primarily from capturing INTENSITY.\n")
  cat("  The improvement over duration alone reflects real perceptual information,\n")
  cat("  not just duration being counted twice.\n")
} else if (avg_pct_intensity > 50) {
  cat("CONCLUSION: sRPE-TL's value is MOSTLY from intensity, but duration\n")
  cat("  contamination contributes meaningfully to the apparent improvement.\n")
} else {
  cat("CONCLUSION: sRPE-TL's apparent improvement is MOSTLY from duration\n")
  cat("  being counted twice. The actual intensity component of sRPE\n")
  cat("  adds limited value beyond duration alone.\n")
}

cat("\n=============================================================================\n")