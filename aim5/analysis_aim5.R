# =============================================================================
# AIM 5 STATISTICAL ANALYSIS: Acute Fatigue Effects on sRPE-HR Relationships
# =============================================================================
# Location: C:\Users\avirgile\Dropbox\Adam Only\School\PhD\UVM Coursework\Dissertation\data\aim5
#
# Research Question:
#   Does acute fatigue from back-to-back competitions alter the coupling 
#   between perceived exertion (sRPE) and HR-derived metrics?
#
# Key Analyses:
#   1. Primary model: sRPE ~ HR_within × GamePosition + HR_between + GamePosition + (1|Athlete) + (1|Sequence)
#   2. Test multiple HR metrics: TRIMPmin, HRmean, HR80%, HRmax
#   3. GamePosition main effect: Does sRPE differ at same HR level in G2 vs G1?
#   4. Interaction: Does sRPE-HR coupling change from G1 to G2?
#   5. Ratio metrics: sRPE/TRIMP, sRPE/TRIMPmin as practitioner summaries
#   6. Random slopes sensitivity
#
# Inputs:
#   - ../merged_data_full.csv
#
# Outputs:
#   - b2b_sequences.csv (identified back-to-back sequences)
#   - b2b_paired_data.csv (athlete-level paired data for B2B games)
#   - descriptives_by_game_position.csv
#   - results_primary_models.csv (HR × GamePosition interaction models)
#   - results_game_position_effects.csv (main effects of G1 vs G2)
#   - results_ratio_metrics.csv (sRPE/TRIMP analysis)
#   - results_random_slopes.csv (sensitivity analysis)
#   - summary_aim5.csv
#   - aim5_models.RData
# =============================================================================

setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim5")

# =============================================================================
# 1. LOAD PACKAGES
# =============================================================================
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("=============================================================================\n")
cat("AIM 5: Acute Fatigue Effects on sRPE-HR Relationships\n")
cat("=============================================================================\n\n")

data <- read.csv("../merged_data_full.csv")
data$session_date <- as.Date(data$session_date)

# Filter to games with full participation and valid data
games <- data %>%
  filter(context == "game",
         full_participation == 1,
         !is.na(sRPE),
         !is.na(TRIMP_min))

cat("Data loaded:\n")
cat("  Total game observations:", nrow(games), "\n")
cat("  Athletes:", n_distinct(games$athlete_name), "\n")
cat("  Game dates:", n_distinct(games$session_date), "\n")

# =============================================================================
# 3. IDENTIFY BACK-TO-BACK SEQUENCES
# =============================================================================
cat("\n=== IDENTIFYING BACK-TO-BACK SEQUENCES ===\n")

# Get unique game dates
game_dates <- sort(unique(games$session_date))

# Find consecutive game dates (back-to-back)
b2b_sequences <- data.frame()
seq_id <- 1

i <- 1
while (i < length(game_dates)) {
  date1 <- game_dates[i]
  date2 <- game_dates[i + 1]
  
  if (as.numeric(date2 - date1) == 1) {
    b2b_sequences <- rbind(b2b_sequences, data.frame(
      sequence_id = seq_id,
      game1_date = date1,
      game2_date = date2
    ))
    seq_id <- seq_id + 1
    i <- i + 2  # Skip past both games
  } else {
    i <- i + 1
  }
}

cat("\nBack-to-back sequences identified:", nrow(b2b_sequences), "\n")
print(b2b_sequences)

# =============================================================================
# 4. CREATE PAIRED DATASET
# =============================================================================
cat("\n=== CREATING PAIRED DATASET ===\n")

# Tag each game observation with sequence info
games$sequence_id <- NA
games$game_position <- NA

for (i in 1:nrow(b2b_sequences)) {
  seq <- b2b_sequences[i, ]
  
  # Game 1
  games$sequence_id[games$session_date == seq$game1_date] <- seq$sequence_id
  games$game_position[games$session_date == seq$game1_date] <- 1
  
  # Game 2
  games$sequence_id[games$session_date == seq$game2_date] <- seq$sequence_id
  games$game_position[games$session_date == seq$game2_date] <- 2
}

# Filter to only B2B games
b2b_games <- games %>%
  filter(!is.na(sequence_id))

cat("Observations in B2B games:", nrow(b2b_games), "\n")

# Keep only athletes with BOTH games in a sequence
b2b_paired <- b2b_games %>%
  group_by(athlete_name, sequence_id) %>%
  filter(n() == 2) %>%  # Must have both G1 and G2
  ungroup()

cat("Paired observations (athletes with both G1 and G2):", nrow(b2b_paired), "\n")
cat("Unique athletes:", n_distinct(b2b_paired$athlete_name), "\n")
cat("Unique sequences with paired data:", n_distinct(b2b_paired$sequence_id), "\n")

# Summary by sequence
seq_summary <- b2b_paired %>%
  group_by(sequence_id) %>%
  summarise(
    game1_date = min(session_date),
    game2_date = max(session_date),
    n_athletes = n_distinct(athlete_name),
    .groups = "drop"
  )

cat("\nAthletes per sequence:\n")
print(seq_summary)

# =============================================================================
# 5. PREPARE VARIABLES
# =============================================================================
cat("\n=== PREPARING VARIABLES ===\n")

# GamePosition as factor (G1 = reference)
b2b_paired$game_position_f <- factor(b2b_paired$game_position, 
                                     levels = c(1, 2), 
                                     labels = c("Game1", "Game2"))

# Calculate within/between decomposition for HR metrics
# "Within" = deviation from athlete's mean across ALL their B2B games
# "Between" = athlete's overall mean

hr_vars <- c("TRIMP_min", "HR_avg_pct_max", "HR80_pct_fb", "HR_peak_pct_max", "TRIMP")

b2b_paired <- b2b_paired %>%
  group_by(athlete_name) %>%
  mutate(
    # TRIMPmin
    TRIMP_min_mean = mean(TRIMP_min, na.rm = TRUE),
    TRIMP_min_within = TRIMP_min - TRIMP_min_mean,
    TRIMP_min_between = TRIMP_min_mean,
    
    # HR average (% max)
    HR_avg_mean = mean(HR_avg_pct_max, na.rm = TRUE),
    HR_avg_within = HR_avg_pct_max - HR_avg_mean,
    HR_avg_between = HR_avg_mean,
    
    # HR80%
    HR80_mean = mean(HR80_pct_fb, na.rm = TRUE),
    HR80_within = HR80_pct_fb - HR80_mean,
    HR80_between = HR80_mean,
    
    # HR peak (% max)
    HR_peak_mean = mean(HR_peak_pct_max, na.rm = TRUE),
    HR_peak_within = HR_peak_pct_max - HR_peak_mean,
    HR_peak_between = HR_peak_mean,
    
    # TRIMP (for ratio metrics)
    TRIMP_mean = mean(TRIMP, na.rm = TRUE),
    TRIMP_within = TRIMP - TRIMP_mean,
    TRIMP_between = TRIMP_mean,
    
    # Ratio metrics
    sRPE_TRIMP_ratio = sRPE / TRIMP,
    sRPE_TRIMPmin_ratio = sRPE / TRIMP_min
  ) %>%
  ungroup()

# Z-score within components for standardized interpretation
b2b_paired <- b2b_paired %>%
  mutate(
    TRIMP_min_within_z = scale(TRIMP_min_within)[,1],
    TRIMP_min_between_z = scale(TRIMP_min_between)[,1],
    HR_avg_within_z = scale(HR_avg_within)[,1],
    HR_avg_between_z = scale(HR_avg_between)[,1],
    HR80_within_z = scale(HR80_within)[,1],
    HR80_between_z = scale(HR80_between)[,1],
    HR_peak_within_z = scale(HR_peak_within)[,1],
    HR_peak_between_z = scale(HR_peak_between)[,1]
  )

cat("Variables prepared.\n")

# =============================================================================
# 6. DESCRIPTIVES BY GAME POSITION
# =============================================================================
cat("\n=== DESCRIPTIVES BY GAME POSITION ===\n")

descriptives <- b2b_paired %>%
  group_by(game_position_f) %>%
  summarise(
    n = n(),
    sRPE_mean = mean(sRPE, na.rm = TRUE),
    sRPE_sd = sd(sRPE, na.rm = TRUE),
    TRIMP_mean = mean(TRIMP, na.rm = TRUE),
    TRIMP_sd = sd(TRIMP, na.rm = TRUE),
    TRIMP_min_mean = mean(TRIMP_min, na.rm = TRUE),
    TRIMP_min_sd = sd(TRIMP_min, na.rm = TRUE),
    HR_avg_mean = mean(HR_avg_pct_max, na.rm = TRUE),
    HR_avg_sd = sd(HR_avg_pct_max, na.rm = TRUE),
    HR80_mean = mean(HR80_pct_fb, na.rm = TRUE),
    HR80_sd = sd(HR80_pct_fb, na.rm = TRUE),
    HR_peak_mean = mean(HR_peak_pct_max, na.rm = TRUE),
    HR_peak_sd = sd(HR_peak_pct_max, na.rm = TRUE),
    sRPE_TRIMP_ratio_mean = mean(sRPE_TRIMP_ratio, na.rm = TRUE),
    sRPE_TRIMP_ratio_sd = sd(sRPE_TRIMP_ratio, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nDescriptives:\n")
print(descriptives %>% mutate(across(where(is.numeric), ~round(., 3))))

# =============================================================================
# 7. HELPER FUNCTIONS
# =============================================================================
safe_lmer <- function(formula, data) {
  tryCatch(
    lmer(formula, data = data, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))),
    error = function(e) NULL
  )
}

extract_model_results <- function(model, model_name) {
  if (is.null(model)) {
    return(data.frame(model = model_name, term = NA, estimate = NA, SE = NA, t = NA, p = NA))
  }
  
  fe <- summary(model)$coefficients
  ci <- tryCatch(confint(model, method = "Wald", parm = "beta_"), error = function(e) NULL)
  
  results <- data.frame(
    model = model_name,
    term = rownames(fe),
    estimate = fe[, "Estimate"],
    SE = fe[, "Std. Error"],
    t = fe[, "t value"],
    p = fe[, "Pr(>|t|)"],
    ci_lower = NA,
    ci_upper = NA,
    row.names = NULL
  )
  
  if (!is.null(ci)) {
    for (i in 1:nrow(results)) {
      if (results$term[i] %in% rownames(ci)) {
        results$ci_lower[i] <- ci[results$term[i], 1]
        results$ci_upper[i] <- ci[results$term[i], 2]
      }
    }
  }
  
  r2 <- tryCatch(performance::r2(model), error = function(e) NULL)
  results$R2_marginal <- if (!is.null(r2)) r2$R2_marginal else NA
  results$R2_conditional <- if (!is.null(r2)) r2$R2_conditional else NA
  
  return(results)
}

# =============================================================================
# 8. PRIMARY MODELS: HR × GamePosition Interaction
# =============================================================================
cat("\n=== PRIMARY MODELS: HR × GamePosition Interaction ===\n")
cat("Testing: sRPE ~ HR_within × GamePosition + HR_between + GamePosition + (1|Athlete) + (1|Sequence)\n\n")

# Model 1: TRIMPmin
m_trimp_min <- safe_lmer(
  
  sRPE ~ TRIMP_min_within_z * game_position_f + TRIMP_min_between_z + 
    (1 | athlete_name) + (1 | sequence_id),
  data = b2b_paired
)

# Model 2: HR average
m_hr_avg <- safe_lmer(
  sRPE ~ HR_avg_within_z * game_position_f + HR_avg_between_z + 
    (1 | athlete_name) + (1 | sequence_id),
  data = b2b_paired
)

# Model 3: HR80%
m_hr80 <- safe_lmer(
  sRPE ~ HR80_within_z * game_position_f + HR80_between_z + 
    (1 | athlete_name) + (1 | sequence_id),
  data = b2b_paired
)

# Model 4: HR peak
m_hr_peak <- safe_lmer(
  sRPE ~ HR_peak_within_z * game_position_f + HR_peak_between_z + 
    (1 | athlete_name) + (1 | sequence_id),
  data = b2b_paired
)

# Print summaries
cat("--- TRIMPmin Model ---\n")
if (!is.null(m_trimp_min)) print(summary(m_trimp_min))

cat("\n--- HR Average Model ---\n")
if (!is.null(m_hr_avg)) print(summary(m_hr_avg))

cat("\n--- HR80% Model ---\n")
if (!is.null(m_hr80)) print(summary(m_hr80))

cat("\n--- HR Peak Model ---\n")
if (!is.null(m_hr_peak)) print(summary(m_hr_peak))

# Combine results
results_primary <- bind_rows(
  extract_model_results(m_trimp_min, "TRIMPmin"),
  extract_model_results(m_hr_avg, "HR_avg"),
  extract_model_results(m_hr80, "HR80"),
  extract_model_results(m_hr_peak, "HR_peak")
)

# =============================================================================
# 9. EXTRACT KEY EFFECTS
# =============================================================================
cat("\n=== KEY EFFECTS SUMMARY ===\n")

extract_key_effects <- function(model, hr_name) {
  if (is.null(model)) return(NULL)
  
  fe <- summary(model)$coefficients
  
  # Find the interaction term (contains ":")
  interaction_term <- rownames(fe)[grep(":", rownames(fe))]
  hr_within_term <- rownames(fe)[grep("within_z$", rownames(fe))]
  game_term <- rownames(fe)[grep("game_position", rownames(fe))][1]
  
  data.frame(
    HR_metric = hr_name,
    
    # GamePosition main effect (is sRPE higher in G2 at same HR?)
    GamePosition_beta = fe[game_term, "Estimate"],
    GamePosition_SE = fe[game_term, "Std. Error"],
    GamePosition_p = fe[game_term, "Pr(>|t|)"],
    
    # HR within effect (baseline coupling in G1)
    HR_within_beta = fe[hr_within_term, "Estimate"],
    HR_within_SE = fe[hr_within_term, "Std. Error"],
    HR_within_p = fe[hr_within_term, "Pr(>|t|)"],
    
    # Interaction (does coupling change in G2?)
    Interaction_beta = fe[interaction_term, "Estimate"],
    Interaction_SE = fe[interaction_term, "Std. Error"],
    Interaction_p = fe[interaction_term, "Pr(>|t|)"]
  )
}

key_effects <- bind_rows(
  extract_key_effects(m_trimp_min, "TRIMPmin"),
  extract_key_effects(m_hr_avg, "HR_avg"),
  extract_key_effects(m_hr80, "HR80"),
  extract_key_effects(m_hr_peak, "HR_peak")
)

cat("\nKey Effects:\n")
print(key_effects %>% mutate(across(where(is.numeric), ~round(., 4))))

# =============================================================================
# 10. RATIO METRICS
# =============================================================================
cat("\n=== RATIO METRICS: sRPE/TRIMP ===\n")

# Model: Does sRPE/TRIMP ratio differ between G1 and G2?
m_ratio_trimp <- safe_lmer(
  sRPE_TRIMP_ratio ~ game_position_f + (1 | athlete_name) + (1 | sequence_id),
  data = b2b_paired %>% filter(is.finite(sRPE_TRIMP_ratio))
)

m_ratio_trimp_min <- safe_lmer(
  sRPE_TRIMPmin_ratio ~ game_position_f + (1 | athlete_name) + (1 | sequence_id),
  data = b2b_paired %>% filter(is.finite(sRPE_TRIMPmin_ratio))
)

if (!is.null(m_ratio_trimp)) {
  cat("\nsRPE/TRIMP Ratio Model:\n")
  print(summary(m_ratio_trimp))
}

if (!is.null(m_ratio_trimp_min)) {
  cat("\nsRPE/TRIMPmin Ratio Model:\n")
  print(summary(m_ratio_trimp_min))
}

# Extract ratio results
results_ratio <- bind_rows(
  extract_model_results(m_ratio_trimp, "sRPE_TRIMP_ratio"),
  extract_model_results(m_ratio_trimp_min, "sRPE_TRIMPmin_ratio")
)

# Ratio descriptives
ratio_desc <- b2b_paired %>%
  filter(is.finite(sRPE_TRIMP_ratio), is.finite(sRPE_TRIMPmin_ratio)) %>%
  group_by(game_position_f) %>%
  summarise(
    sRPE_TRIMP_mean = mean(sRPE_TRIMP_ratio, na.rm = TRUE),
    sRPE_TRIMP_sd = sd(sRPE_TRIMP_ratio, na.rm = TRUE),
    sRPE_TRIMPmin_mean = mean(sRPE_TRIMPmin_ratio, na.rm = TRUE),
    sRPE_TRIMPmin_sd = sd(sRPE_TRIMPmin_ratio, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nRatio Descriptives by Game Position:\n")
print(ratio_desc %>% mutate(across(where(is.numeric), ~round(., 4))))

# =============================================================================
# 11. RANDOM SLOPES SENSITIVITY
# =============================================================================
cat("\n=== RANDOM SLOPES SENSITIVITY ===\n")

# Try adding random slopes for HR_within by athlete
m_rs_trimp_min <- safe_lmer(
  sRPE ~ TRIMP_min_within_z * game_position_f + TRIMP_min_between_z + 
    (1 + TRIMP_min_within_z | athlete_name) + (1 | sequence_id),
  data = b2b_paired
)

rs_results <- data.frame(
  model = "TRIMPmin_RandomSlopes",
  converged = !is.null(m_rs_trimp_min),
  singular = if (!is.null(m_rs_trimp_min)) isSingular(m_rs_trimp_min) else NA
)

if (!is.null(m_rs_trimp_min) && !isSingular(m_rs_trimp_min)) {
  cat("Random slopes model converged successfully.\n")
  cat("AIC comparison: RI =", round(AIC(m_trimp_min), 1), 
      ", RS =", round(AIC(m_rs_trimp_min), 1), "\n")
  rs_results$AIC_RI <- AIC(m_trimp_min)
  rs_results$AIC_RS <- AIC(m_rs_trimp_min)
} else {
  cat("Random slopes model failed or singular. Using random intercepts only.\n")
  rs_results$AIC_RI <- if (!is.null(m_trimp_min)) AIC(m_trimp_min) else NA
  rs_results$AIC_RS <- NA
}

# =============================================================================
# 12. CREATE SUMMARY TABLE
# =============================================================================
cat("\n=== CREATING SUMMARY ===\n")

# Extract key findings
get_effect <- function(df, metric, effect_col) {
  val <- df[[effect_col]][df$HR_metric == metric]
  if (length(val) == 0 || is.na(val)) return(NA)
  val
}

summary_aim5 <- data.frame(
  finding = c(
    "N_sequences",
    "N_paired_observations",
    "N_athletes",
    "sRPE_Game1_mean",
    "sRPE_Game2_mean",
    "sRPE_diff_G2_minus_G1",
    "GamePosition_effect_TRIMPmin (β)",
    "GamePosition_effect_TRIMPmin (p)",
    "Interaction_TRIMPmin (β)",
    "Interaction_TRIMPmin (p)",
    "GamePosition_effect_HRavg (β)",
    "GamePosition_effect_HRavg (p)",
    "Interaction_HRavg (β)",
    "Interaction_HRavg (p)",
    "sRPE_TRIMP_ratio_G1",
    "sRPE_TRIMP_ratio_G2"
  ),
  value = c(
    nrow(b2b_sequences),
    nrow(b2b_paired),
    n_distinct(b2b_paired$athlete_name),
    descriptives$sRPE_mean[1],
    descriptives$sRPE_mean[2],
    descriptives$sRPE_mean[2] - descriptives$sRPE_mean[1],
    get_effect(key_effects, "TRIMPmin", "GamePosition_beta"),
    get_effect(key_effects, "TRIMPmin", "GamePosition_p"),
    get_effect(key_effects, "TRIMPmin", "Interaction_beta"),
    get_effect(key_effects, "TRIMPmin", "Interaction_p"),
    get_effect(key_effects, "HR_avg", "GamePosition_beta"),
    get_effect(key_effects, "HR_avg", "GamePosition_p"),
    get_effect(key_effects, "HR_avg", "Interaction_beta"),
    get_effect(key_effects, "HR_avg", "Interaction_p"),
    ratio_desc$sRPE_TRIMP_mean[1],
    ratio_desc$sRPE_TRIMP_mean[2]
  )
)

# =============================================================================
# 13. SAVE OUTPUTS
# =============================================================================
cat("\n=== SAVING OUTPUTS ===\n")

write.csv(b2b_sequences, "b2b_sequences.csv", row.names = FALSE)
write.csv(b2b_paired, "b2b_paired_data.csv", row.names = FALSE)
write.csv(descriptives, "descriptives_by_game_position.csv", row.names = FALSE)
write.csv(results_primary, "results_primary_models.csv", row.names = FALSE)
write.csv(key_effects, "results_game_position_effects.csv", row.names = FALSE)
write.csv(results_ratio, "results_ratio_metrics.csv", row.names = FALSE)
write.csv(rs_results, "results_random_slopes.csv", row.names = FALSE)
write.csv(summary_aim5, "summary_aim5.csv", row.names = FALSE)

save(m_trimp_min, m_hr_avg, m_hr80, m_hr_peak,
     m_ratio_trimp, m_ratio_trimp_min, m_rs_trimp_min,
     file = "aim5_models.RData")

cat("\nFiles saved to:", getwd(), "\n")
cat("  - b2b_sequences.csv\n")
cat("  - b2b_paired_data.csv\n")
cat("  - descriptives_by_game_position.csv\n")
cat("  - results_primary_models.csv\n")
cat("  - results_game_position_effects.csv\n")
cat("  - results_ratio_metrics.csv\n")
cat("  - results_random_slopes.csv\n")
cat("  - summary_aim5.csv\n")
cat("  - aim5_models.RData\n")

# =============================================================================
# 14. FINAL SUMMARY
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# AIM 5 RESULTS SUMMARY                                                   #\n")
cat("###########################################################################\n")

cat("\n=== DATA ===\n")
cat("Back-to-back sequences:", nrow(b2b_sequences), "\n")
cat("Paired observations:", nrow(b2b_paired), "\n")
cat("Athletes:", n_distinct(b2b_paired$athlete_name), "\n")

cat("\n=== DESCRIPTIVES ===\n")
cat("Game 1 sRPE:", round(descriptives$sRPE_mean[1], 2), "±", round(descriptives$sRPE_sd[1], 2), "\n")
cat("Game 2 sRPE:", round(descriptives$sRPE_mean[2], 2), "±", round(descriptives$sRPE_sd[2], 2), "\n")
cat("Difference:", round(descriptives$sRPE_mean[2] - descriptives$sRPE_mean[1], 2), "units\n")

cat("\n=== KEY QUESTION 1: Is sRPE higher in Game 2 at the same HR level? ===\n")
if (!is.null(key_effects)) {
  gp_beta <- key_effects$GamePosition_beta[key_effects$HR_metric == "TRIMPmin"]
  gp_p <- key_effects$GamePosition_p[key_effects$HR_metric == "TRIMPmin"]
  cat("GamePosition effect (TRIMPmin model):\n")
  cat("  β =", round(gp_beta, 3), ", p =", 
      ifelse(gp_p < 0.001, "< .001", round(gp_p, 3)), "\n")
  if (gp_p < 0.05 && gp_beta > 0) {
    cat("  --> YES: Athletes rate G2 higher than G1 at the same physiological load.\n")
    cat("      This suggests fatigue inflation of sRPE.\n")
  } else if (gp_p >= 0.05) {
    cat("  --> NO: sRPE does not systematically differ between G1 and G2.\n")
  }
}

cat("\n=== KEY QUESTION 2: Does sRPE-HR coupling change from G1 to G2? ===\n")
if (!is.null(key_effects)) {
  int_beta <- key_effects$Interaction_beta[key_effects$HR_metric == "TRIMPmin"]
  int_p <- key_effects$Interaction_p[key_effects$HR_metric == "TRIMPmin"]
  cat("HR × GamePosition interaction (TRIMPmin model):\n")
  cat("  β =", round(int_beta, 3), ", p =", 
      ifelse(int_p < 0.001, "< .001", round(int_p, 3)), "\n")
  if (int_p < 0.05) {
    cat("  --> YES: The sRPE-HR relationship differs between G1 and G2.\n")
  } else {
    cat("  --> NO: The sRPE-HR coupling is similar in G1 and G2.\n")
  }
}

cat("\n=== RATIO METRICS ===\n")
cat("sRPE/TRIMP ratio:\n")
cat("  Game 1:", round(ratio_desc$sRPE_TRIMP_mean[1], 4), "\n")
cat("  Game 2:", round(ratio_desc$sRPE_TRIMP_mean[2], 4), "\n")
cat("  Change:", round(ratio_desc$sRPE_TRIMP_mean[2] - ratio_desc$sRPE_TRIMP_mean[1], 4), "\n")

cat("\n=============================================================================\n")
cat("AIM 5 ANALYSIS COMPLETE\n")
cat("=============================================================================\n")