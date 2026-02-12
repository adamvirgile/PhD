# =============================================================================
# AIM 2 STATISTICAL ANALYSIS: Stability of sRPE-M Under Matched Demands
# =============================================================================
# Location: C:\Users\avirgile\Dropbox\Adam Only\School\PhD\UVM Coursework\Dissertation\data\aim2
#
# This script:
#   1. Algorithmically identifies matched practice sessions
#      - First practice of week (GAME-GAME-OFF-PRACTICE pattern)
#      - Duration-matched within 10 minutes
#      - Physiologically similar (TRIMP, HR80, HRavg)
#   2. Reports matching quality
#   3. Calculates stability metrics (ICC, Bland-Altman, SEM, MDC, Spearman)
#
# Inputs (relative to aim2 folder):
#   - ../aim1/aim1_analysis_data.csv (preprocessed data from Aim 1)
#   - ../team_schedule.xlsx (for identifying session patterns)
#
# Outputs (saved to aim2 folder):
#   - session_pair_selection.csv (all candidate pairs ranked)
#   - selected_pair_descriptives.csv (matching quality report)
#   - athlete_paired_data.csv (paired sRPE data for selected sessions)
#   - stability_results.csv (ICC, SEM, MDC, Spearman)
#   - bland_altman_results.csv (bias and LOA)
#   - selected_pair_info.csv (selected session details)
#   - bland_altman_plot.png
# =============================================================================

# Set working directory to aim2 folder
setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim2")

# =============================================================================
# 1. LOAD PACKAGES
# =============================================================================
library(dplyr)
library(tidyr)
library(lubridate)
library(readxl)
library(irr)       # For ICC
library(ggplot2)

# Check for blandr package
if (!requireNamespace("blandr", quietly = TRUE)) {
  message("Package 'blandr' not installed. Install with: install.packages('blandr')")
  message("Proceeding with manual Bland-Altman calculations...")
  use_blandr <- FALSE
} else {
  library(blandr)
  use_blandr <- TRUE
}

# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("=============================================================================\n")
cat("AIM 2: Stability of sRPE-M Under Matched Physiological Demands\n")
cat("=============================================================================\n\n")

# Load the preprocessed Aim 1 data (from ../aim1/)
aim1_data <- read.csv("../aim1/aim1_analysis_data.csv")

# Load team schedule (from parent /data/ folder)
schedule <- read_excel("../team_schedule.xlsx")

cat("Data loaded:\n")
cat("  Aim 1 observations:", nrow(aim1_data), "\n")
cat("  Schedule entries:", nrow(schedule), "\n")

# =============================================================================
# 3. PREPARE SCHEDULE DATA
# =============================================================================
cat("\n=== PREPARING SCHEDULE DATA ===\n")

# Clean schedule
schedule_clean <- schedule %>%
  mutate(
    session_date = as.Date(Date),
    event_type = tolower(trimws(Event))
  ) %>%
  arrange(session_date) %>%
  select(session_date, event_type)

# Add day of week
schedule_clean <- schedule_clean %>%
  mutate(
    day_of_week = weekdays(session_date),
    day_num = wday(session_date)  # 1 = Sunday, 7 = Saturday
  )

# Create lagged variables to identify preceding pattern
schedule_clean <- schedule_clean %>%
  mutate(
    lag1_event = lag(event_type, 1),  # Yesterday
    lag2_event = lag(event_type, 2),  # 2 days ago
    lag3_event = lag(event_type, 3)   # 3 days ago
  )

# Identify "first practice of week" pattern: GAME-GAME-OFF-PRACTICE
# This means: lag3=game, lag2=game, lag1=off, current=practice
schedule_clean <- schedule_clean %>%
  mutate(
    is_first_practice_of_week = (
      event_type == "practice" &
        lag1_event == "off" &
        lag2_event == "game" &
        lag3_event == "game"
    )
  )

# Get candidate practice dates
candidate_dates <- schedule_clean %>%
  filter(is_first_practice_of_week == TRUE) %>%
  pull(session_date)

cat("Candidate 'first practice of week' sessions found:", length(candidate_dates), "\n")
cat("Dates:", paste(candidate_dates, collapse = ", "), "\n")

# =============================================================================
# 4. FILTER AIM1 DATA TO PRACTICES ONLY
# =============================================================================

practices <- aim1_data %>%
  filter(context == "practice") %>%
  mutate(session_date = as.Date(session_date))

cat("\nPractice sessions in data:", n_distinct(practices$session_date), "\n")
cat("Practice observations:", nrow(practices), "\n")

# Filter to candidate dates that exist in our data
candidate_dates_in_data <- candidate_dates[candidate_dates %in% unique(practices$session_date)]
cat("Candidate dates with data:", length(candidate_dates_in_data), "\n")

if (length(candidate_dates_in_data) < 2) {
  cat("\nWARNING: Not enough candidate sessions found with GAME-GAME-OFF-PRACTICE pattern.\n")
  cat("Relaxing criteria to all practice sessions for matching...\n")
  candidate_dates_in_data <- unique(practices$session_date)
}

# =============================================================================
# 5. CALCULATE SESSION-LEVEL SUMMARIES
# =============================================================================
cat("\n=== CALCULATING SESSION-LEVEL SUMMARIES ===\n")

session_summaries <- practices %>%
  filter(session_date %in% candidate_dates_in_data) %>%
  group_by(session_date) %>%
  summarise(
    n_athletes = n(),
    duration_mean = mean(duration_mins, na.rm = TRUE),
    duration_sd = sd(duration_mins, na.rm = TRUE),
    TRIMP_mean = mean(TRIMP, na.rm = TRUE),
    TRIMP_sd = sd(TRIMP, na.rm = TRUE),
    HR80_mins_mean = mean(HR80_mins, na.rm = TRUE),
    HR80_mins_sd = sd(HR80_mins, na.rm = TRUE),
    HR_avg_pct_max_mean = mean(HR_avg_pct_max, na.rm = TRUE),
    HR_avg_pct_max_sd = sd(HR_avg_pct_max, na.rm = TRUE),
    sRPE_mean = mean(sRPE, na.rm = TRUE),
    sRPE_sd = sd(sRPE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    day_of_week = weekdays(session_date)
  )

cat("Session summaries calculated for", nrow(session_summaries), "candidate sessions\n")

# =============================================================================
# 6. FIND BEST-MATCHED SESSION PAIR
# =============================================================================
cat("\n=== FINDING BEST-MATCHED SESSION PAIR ===\n")

# Create all possible pairs
session_dates <- session_summaries$session_date
n_sessions <- length(session_dates)

if (n_sessions < 2) {
  stop("Not enough sessions to create pairs!")
}

# Generate all pairs
pairs_list <- combn(session_dates, 2, simplify = FALSE)
cat("Total possible pairs:", length(pairs_list), "\n")

# Function to evaluate a session pair
evaluate_pair <- function(date1, date2, practices_df, session_summaries_df) {
  
  # Get session-level info
  s1 <- session_summaries_df %>% filter(session_date == date1)
  s2 <- session_summaries_df %>% filter(session_date == date2)
  
  # Duration difference (hard constraint)
  duration_diff <- abs(s1$duration_mean - s2$duration_mean)
  
  # Get athlete-level data for both sessions
  d1 <- practices_df %>% filter(session_date == date1) %>%
    select(athlete_name, sRPE, duration_mins, TRIMP, HR80_mins, HR_avg_pct_max)
  d2 <- practices_df %>% filter(session_date == date2) %>%
    select(athlete_name, sRPE, duration_mins, TRIMP, HR80_mins, HR_avg_pct_max)
  
  # Find athletes present in both sessions
  common_athletes <- intersect(d1$athlete_name, d2$athlete_name)
  n_common <- length(common_athletes)
  
  if (n_common < 5) {
    return(data.frame(
      date1 = date1, date2 = date2,
      days_apart = as.numeric(date2 - date1),
      same_day_of_week = weekdays(date1) == weekdays(date2),
      duration_diff_mins = duration_diff,
      n_athletes_common = n_common,
      n_athletes_matched = 0,
      pct_matched = 0,
      mean_euclidean_dist = NA,
      TRIMP_diff_pct = NA,
      HR80_diff_pct = NA,
      HRavg_diff_pct = NA,
      feasible = FALSE
    ))
  }
  
  # Merge athlete data
  paired <- d1 %>%
    inner_join(d2, by = "athlete_name", suffix = c("_s1", "_s2"))
  
  # Calculate % differences for each athlete
  paired <- paired %>%
    mutate(
      TRIMP_diff_pct = abs(TRIMP_s1 - TRIMP_s2) / ((TRIMP_s1 + TRIMP_s2) / 2) * 100,
      HR80_diff_pct = abs(HR80_mins_s1 - HR80_mins_s2) / ((HR80_mins_s1 + HR80_mins_s2) / 2 + 0.01) * 100,
      HRavg_diff_pct = abs(HR_avg_pct_max_s1 - HR_avg_pct_max_s2) / ((HR_avg_pct_max_s1 + HR_avg_pct_max_s2) / 2) * 100,
      duration_diff_pct = abs(duration_mins_s1 - duration_mins_s2) / ((duration_mins_s1 + duration_mins_s2) / 2) * 100
    )
  
  # Count athletes with TRIMP < 25% difference (primary matching criterion)
  # HR80 and HRAvg are reported but NOT used for matching threshold
  paired <- paired %>%
    mutate(
      is_matched = (TRIMP_diff_pct < 25)  # Simplified: only TRIMP matching
    )
  
  n_matched <- sum(paired$is_matched, na.rm = TRUE)
  pct_matched <- n_matched / n_common * 100
  
  # Calculate Euclidean distance on z-scored metrics
  # Z-score within the pair
  paired <- paired %>%
    mutate(
      TRIMP_z1 = scale(TRIMP_s1)[,1],
      TRIMP_z2 = scale(TRIMP_s2)[,1],
      HR80_z1 = scale(HR80_mins_s1)[,1],
      HR80_z2 = scale(HR80_mins_s2)[,1],
      HRavg_z1 = scale(HR_avg_pct_max_s1)[,1],
      HRavg_z2 = scale(HR_avg_pct_max_s2)[,1]
    )
  
  # Euclidean distance per athlete
  paired <- paired %>%
    mutate(
      euclidean_dist = sqrt(
        (TRIMP_s1 - TRIMP_s2)^2 + 
          (HR80_mins_s1 - HR80_mins_s2)^2 + 
          (HR_avg_pct_max_s1 - HR_avg_pct_max_s2)^2
      )
    )
  
  mean_euc_dist <- mean(paired$euclidean_dist, na.rm = TRUE)
  
  # Mean % differences
  mean_TRIMP_diff <- mean(paired$TRIMP_diff_pct, na.rm = TRUE)
  mean_HR80_diff <- mean(paired$HR80_diff_pct, na.rm = TRUE)
  mean_HRavg_diff <- mean(paired$HRavg_diff_pct, na.rm = TRUE)
  
  data.frame(
    date1 = date1,
    date2 = date2,
    days_apart = as.numeric(date2 - date1),
    same_day_of_week = weekdays(date1) == weekdays(date2),
    duration_diff_mins = duration_diff,
    n_athletes_common = n_common,
    n_athletes_matched = n_matched,
    pct_matched = pct_matched,
    mean_euclidean_dist = mean_euc_dist,
    TRIMP_diff_pct = mean_TRIMP_diff,
    HR80_diff_pct = mean_HR80_diff,
    HRavg_diff_pct = mean_HRavg_diff,
    feasible = TRUE
  )
}

# Evaluate all pairs
cat("Evaluating all session pairs...\n")
pair_evaluations <- bind_rows(lapply(pairs_list, function(p) {
  evaluate_pair(p[1], p[2], practices, session_summaries)
}))

# =============================================================================
# MATCHING CONSTRAINTS (per proposal requirements)
# =============================================================================
MAX_DAYS_APART <- 14        # Sessions must be within 2 weeks
MAX_DURATION_DIFF <- 10     # Duration difference in minutes
TRIMP_THRESHOLD <- 25       # % difference threshold for TRIMP matching

# Filter to feasible pairs
# SORT BY: 1) Lowest TRIMP diff, 2) Highest % matched, 3) Most athletes
feasible_pairs <- pair_evaluations %>%
  filter(feasible == TRUE & 
           duration_diff_mins <= MAX_DURATION_DIFF &
           days_apart <= MAX_DAYS_APART) %>%
  arrange(TRIMP_diff_pct, desc(pct_matched), desc(n_athletes_common))

cat("\nFeasible pairs (duration diff <=", MAX_DURATION_DIFF, "mins AND days apart <=", MAX_DAYS_APART, "):", 
    nrow(feasible_pairs), "\n")

if (nrow(feasible_pairs) == 0) {
  cat("WARNING: No pairs meet strict criteria. Relaxing temporal constraint to 21 days...\n")
  feasible_pairs <- pair_evaluations %>%
    filter(feasible == TRUE & duration_diff_mins <= MAX_DURATION_DIFF & days_apart <= 21) %>%
    arrange(TRIMP_diff_pct, desc(pct_matched), desc(n_athletes_common))
}

if (nrow(feasible_pairs) == 0) {
  cat("WARNING: Still no pairs. Using duration constraint only...\n")
  feasible_pairs <- pair_evaluations %>%
    filter(feasible == TRUE & duration_diff_mins <= MAX_DURATION_DIFF) %>%
    arrange(TRIMP_diff_pct, desc(pct_matched), desc(n_athletes_common))
}

# Select best pair
best_pair <- feasible_pairs[1, ]

cat("\n=== SELECTED SESSION PAIR ===\n")
cat("Session 1:", as.character(best_pair$date1), "(", weekdays(best_pair$date1), ")\n")
cat("Session 2:", as.character(best_pair$date2), "(", weekdays(best_pair$date2), ")\n")
cat("Days apart:", best_pair$days_apart, "\n")
cat("Same day of week:", best_pair$same_day_of_week, "\n")
cat("Duration difference:", round(best_pair$duration_diff_mins, 1), "minutes\n")
cat("Athletes in both sessions:", best_pair$n_athletes_common, "\n")
cat("Athletes with TRIMP <25% diff:", best_pair$n_athletes_matched, 
    "(", round(best_pair$pct_matched, 1), "%)\n")
cat("\nMatching criteria used:\n")
cat("  - Duration: within", MAX_DURATION_DIFF, "minutes\n")
cat("  - Temporal: within", MAX_DAYS_APART, "days\n")
cat("  - TRIMP: <25% difference (for athlete-level matching)\n")
cat("  - HR80 and HRAvg: reported for description only\n")

# =============================================================================
# 7. CREATE PAIRED DATASET FOR SELECTED SESSIONS
# =============================================================================
cat("\n=== CREATING PAIRED DATASET ===\n")

selected_date1 <- best_pair$date1
selected_date2 <- best_pair$date2

# Get data for selected sessions
session1_data <- practices %>%
  filter(session_date == selected_date1) %>%
  select(athlete_name, sRPE_s1 = sRPE, duration_s1 = duration_mins, 
         TRIMP_s1 = TRIMP, HR80_mins_s1 = HR80_mins, 
         HR_avg_pct_max_s1 = HR_avg_pct_max, TRIMP_min_s1 = TRIMP_min)

session2_data <- practices %>%
  filter(session_date == selected_date2) %>%
  select(athlete_name, sRPE_s2 = sRPE, duration_s2 = duration_mins, 
         TRIMP_s2 = TRIMP, HR80_mins_s2 = HR80_mins, 
         HR_avg_pct_max_s2 = HR_avg_pct_max, TRIMP_min_s2 = TRIMP_min)

# Merge
paired_data <- session1_data %>%
  inner_join(session2_data, by = "athlete_name")

cat("Paired observations:", nrow(paired_data), "athletes\n")

# Calculate differences
paired_data <- paired_data %>%
  mutate(
    sRPE_diff = sRPE_s2 - sRPE_s1,
    sRPE_mean = (sRPE_s1 + sRPE_s2) / 2,
    TRIMP_diff_pct = abs(TRIMP_s1 - TRIMP_s2) / ((TRIMP_s1 + TRIMP_s2) / 2) * 100,
    HR80_diff_pct = abs(HR80_mins_s1 - HR80_mins_s2) / ((HR80_mins_s1 + HR80_mins_s2) / 2 + 0.01) * 100,
    HRavg_diff_pct = abs(HR_avg_pct_max_s1 - HR_avg_pct_max_s2) / ((HR_avg_pct_max_s1 + HR_avg_pct_max_s2) / 2) * 100
  )

# =============================================================================
# 8. REPORT MATCHING QUALITY
# =============================================================================
cat("\n=== MATCHING QUALITY REPORT ===\n")

# Session-level descriptives
s1_summary <- session_summaries %>% filter(session_date == selected_date1)
s2_summary <- session_summaries %>% filter(session_date == selected_date2)

matching_quality <- data.frame(
  metric = c("Duration (mins)", "TRIMP", "HR80 (mins)", "HR_avg_%max", "sRPE", "N Athletes"),
  session1_mean = c(s1_summary$duration_mean, s1_summary$TRIMP_mean, s1_summary$HR80_mins_mean,
                    s1_summary$HR_avg_pct_max_mean, s1_summary$sRPE_mean, s1_summary$n_athletes),
  session1_sd = c(s1_summary$duration_sd, s1_summary$TRIMP_sd, s1_summary$HR80_mins_sd,
                  s1_summary$HR_avg_pct_max_sd, s1_summary$sRPE_sd, NA),
  session2_mean = c(s2_summary$duration_mean, s2_summary$TRIMP_mean, s2_summary$HR80_mins_mean,
                    s2_summary$HR_avg_pct_max_mean, s2_summary$sRPE_mean, s2_summary$n_athletes),
  session2_sd = c(s2_summary$duration_sd, s2_summary$TRIMP_sd, s2_summary$HR80_mins_sd,
                  s2_summary$HR_avg_pct_max_sd, s2_summary$sRPE_sd, NA)
) %>%
  mutate(
    difference = session2_mean - session1_mean,
    pct_difference = abs(difference) / ((session1_mean + session2_mean) / 2) * 100
  )

cat("\nSession Comparison:\n")
print(matching_quality %>% mutate(across(where(is.numeric), ~round(., 2))))

# Athlete-level matching summary
cat("\nAthlete-Level Physiological Matching:\n")
cat("  Mean TRIMP % diff:", round(mean(paired_data$TRIMP_diff_pct, na.rm = TRUE), 1), "%\n")
cat("  Mean HR80 % diff:", round(mean(paired_data$HR80_diff_pct, na.rm = TRUE), 1), "%\n")
cat("  Mean HRavg % diff:", round(mean(paired_data$HRavg_diff_pct, na.rm = TRUE), 1), "%\n")
cat("  Athletes with all metrics <25% diff:", 
    sum(paired_data$TRIMP_diff_pct < 25 & paired_data$HR80_diff_pct < 25 & 
          paired_data$HRavg_diff_pct < 25, na.rm = TRUE), 
    "of", nrow(paired_data), "\n")

# =============================================================================
# 9. CALCULATE STABILITY METRICS
# =============================================================================
cat("\n=== CALCULATING STABILITY METRICS ===\n")

# --- ICC(3,1): Two-way mixed, consistency, single measure ---
icc_data <- paired_data %>%
  select(sRPE_s1, sRPE_s2) %>%
  as.matrix()

icc_result <- icc(icc_data, model = "twoway", type = "consistency", unit = "single")

cat("\nICC(3,1) - Two-way mixed, consistency, single measure:\n")
cat("  ICC:", round(icc_result$value, 3), "\n")
cat("  95% CI: [", round(icc_result$lbound, 3), ",", round(icc_result$ubound, 3), "]\n")
cat("  F-value:", round(icc_result$Fvalue, 2), "\n")
cat("  p-value:", ifelse(icc_result$p.value < 0.001, "< .001", round(icc_result$p.value, 4)), "\n")

# --- SEM and MDC ---
sd_pooled <- sd(c(paired_data$sRPE_s1, paired_data$sRPE_s2), na.rm = TRUE)
sem <- sd_pooled * sqrt(1 - icc_result$value)
mdc95 <- sem * 1.96 * sqrt(2)

cat("\nMeasurement Error:\n")
cat("  Pooled SD:", round(sd_pooled, 3), "\n")
cat("  SEM:", round(sem, 3), "\n")
cat("  MDC95:", round(mdc95, 3), "\n")

# --- Bland-Altman Analysis ---
cat("\nBland-Altman Analysis:\n")

mean_diff <- mean(paired_data$sRPE_diff, na.rm = TRUE)
sd_diff <- sd(paired_data$sRPE_diff, na.rm = TRUE)
loa_lower <- mean_diff - 1.96 * sd_diff
loa_upper <- mean_diff + 1.96 * sd_diff

cat("  Mean Bias:", round(mean_diff, 3), "\n")
cat("  SD of Differences:", round(sd_diff, 3), "\n")
cat("  95% Limits of Agreement: [", round(loa_lower, 3), ",", round(loa_upper, 3), "]\n")

# 95% CI for mean bias
n_pairs <- nrow(paired_data)
se_bias <- sd_diff / sqrt(n_pairs)
bias_ci_lower <- mean_diff - 1.96 * se_bias
bias_ci_upper <- mean_diff + 1.96 * se_bias
cat("  95% CI for Bias: [", round(bias_ci_lower, 3), ",", round(bias_ci_upper, 3), "]\n")

# --- Spearman Correlation (rank-order stability) ---
spearman_result <- cor.test(paired_data$sRPE_s1, paired_data$sRPE_s2, method = "spearman")

cat("\nSpearman Correlation (rank-order stability):\n")
cat("  rho:", round(spearman_result$estimate, 3), "\n")
cat("  p-value:", ifelse(spearman_result$p.value < 0.001, "< .001", round(spearman_result$p.value, 4)), "\n")

# --- Pearson Correlation (for comparison) ---
pearson_result <- cor.test(paired_data$sRPE_s1, paired_data$sRPE_s2, method = "pearson")

cat("\nPearson Correlation:\n")
cat("  r:", round(pearson_result$estimate, 3), "\n")
cat("  95% CI: [", round(pearson_result$conf.int[1], 3), ",", round(pearson_result$conf.int[2], 3), "]\n")
cat("  p-value:", ifelse(pearson_result$p.value < 0.001, "< .001", round(pearson_result$p.value, 4)), "\n")

# =============================================================================
# 10. CREATE OUTPUT TABLES
# =============================================================================
cat("\n=== CREATING OUTPUT TABLES ===\n")

# Stability results
stability_results <- data.frame(
  metric = c("ICC(3,1)", "ICC_CI_lower", "ICC_CI_upper", "ICC_F", "ICC_p",
             "SD_pooled", "SEM", "MDC95",
             "Spearman_rho", "Spearman_p",
             "Pearson_r", "Pearson_CI_lower", "Pearson_CI_upper", "Pearson_p",
             "N_athletes"),
  value = c(icc_result$value, icc_result$lbound, icc_result$ubound, icc_result$Fvalue, icc_result$p.value,
            sd_pooled, sem, mdc95,
            spearman_result$estimate, spearman_result$p.value,
            pearson_result$estimate, pearson_result$conf.int[1], pearson_result$conf.int[2], pearson_result$p.value,
            n_pairs)
)

# Bland-Altman results
bland_altman_results <- data.frame(
  metric = c("Mean_Bias", "SD_Diff", "LOA_Lower", "LOA_Upper", 
             "Bias_CI_Lower", "Bias_CI_Upper", "N_athletes"),
  value = c(mean_diff, sd_diff, loa_lower, loa_upper,
            bias_ci_lower, bias_ci_upper, n_pairs)
)

# Selected pair info
selected_pair_info <- data.frame(
  item = c("Session_1_Date", "Session_2_Date", 
           "Session_1_DayOfWeek", "Session_2_DayOfWeek",
           "Days_Apart", "Same_Day_of_Week",
           "Duration_Diff_Mins", "N_Athletes_Common", 
           "N_Athletes_Matched_25pct", "Pct_Matched"),
  value = c(as.character(selected_date1), as.character(selected_date2),
            weekdays(selected_date1), weekdays(selected_date2),
            best_pair$days_apart, best_pair$same_day_of_week,
            round(best_pair$duration_diff_mins, 1), best_pair$n_athletes_common,
            best_pair$n_athletes_matched, round(best_pair$pct_matched, 1))
)

# =============================================================================
# 11. SAVE OUTPUTS
# =============================================================================

# All outputs saved to current directory (aim2)

write.csv(feasible_pairs, "session_pair_selection.csv", row.names = FALSE)
write.csv(matching_quality, "selected_pair_descriptives.csv", row.names = FALSE)
write.csv(paired_data, "athlete_paired_data.csv", row.names = FALSE)
write.csv(stability_results, "stability_results.csv", row.names = FALSE)
write.csv(bland_altman_results, "bland_altman_results.csv", row.names = FALSE)
write.csv(selected_pair_info, "selected_pair_info.csv", row.names = FALSE)

cat("\nFiles saved to:", getwd(), "\n")
cat("  - session_pair_selection.csv (all candidate pairs ranked)\n")
cat("  - selected_pair_descriptives.csv (matching quality)\n")
cat("  - athlete_paired_data.csv (paired sRPE data)\n")
cat("  - stability_results.csv (ICC, SEM, MDC, correlations)\n")
cat("  - bland_altman_results.csv (bias and LOA)\n")
cat("  - selected_pair_info.csv (selected session info)\n")

# =============================================================================
# 12. CREATE BLAND-ALTMAN PLOT
# =============================================================================
cat("\n=== CREATING BLAND-ALTMAN PLOT ===\n")

ba_plot <- ggplot(paired_data, aes(x = sRPE_mean, y = sRPE_diff)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = mean_diff, color = "blue", linewidth = 1) +
  geom_hline(yintercept = loa_lower, color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = loa_upper, color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dotted") +
  annotate("text", x = max(paired_data$sRPE_mean, na.rm = TRUE), y = mean_diff, 
           label = paste("Bias =", round(mean_diff, 2)), hjust = 1, vjust = -0.5, color = "blue") +
  annotate("text", x = max(paired_data$sRPE_mean, na.rm = TRUE), y = loa_upper, 
           label = paste("+1.96 SD =", round(loa_upper, 2)), hjust = 1, vjust = -0.5, color = "red") +
  annotate("text", x = max(paired_data$sRPE_mean, na.rm = TRUE), y = loa_lower, 
           label = paste("-1.96 SD =", round(loa_lower, 2)), hjust = 1, vjust = 1.5, color = "red") +
  labs(
    title = "Bland-Altman Plot: sRPE-M Stability",
    subtitle = paste("Session 1:", selected_date1, "vs Session 2:", selected_date2),
    x = "Mean sRPE [(Session 1 + Session 2) / 2]",
    y = "Difference (Session 2 - Session 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("bland_altman_plot.png", ba_plot, width = 8, height = 6, dpi = 300)
cat("Bland-Altman plot saved: bland_altman_plot.png\n")

# =============================================================================
# 13. PRINT FINAL SUMMARY
# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("# AIM 2 RESULTS SUMMARY                                                   #\n")
cat("###########################################################################\n")

cat("\n=== SELECTED SESSION PAIR ===\n")
cat("Session 1:", as.character(selected_date1), "(", weekdays(selected_date1), ")\n")
cat("Session 2:", as.character(selected_date2), "(", weekdays(selected_date2), ")\n")
cat("Days apart:", best_pair$days_apart, "\n")
cat("Duration difference:", round(best_pair$duration_diff_mins, 1), "minutes\n")

cat("\n=== STABILITY METRICS ===\n")
cat("ICC(3,1):", round(icc_result$value, 3), 
    "[", round(icc_result$lbound, 3), ",", round(icc_result$ubound, 3), "]\n")
cat("SEM:", round(sem, 3), "sRPE units\n")
cat("MDC95:", round(mdc95, 3), "sRPE units\n")

cat("\n=== AGREEMENT ===\n")
cat("Mean Bias:", round(mean_diff, 3), "\n")
cat("95% LOA: [", round(loa_lower, 3), ",", round(loa_upper, 3), "]\n")

cat("\n=== CORRELATIONS ===\n")
cat("Spearman rho:", round(spearman_result$estimate, 3), "\n")
cat("Pearson r:", round(pearson_result$estimate, 3), "\n")

cat("\n=============================================================================\n")
cat("AIM 2 ANALYSIS COMPLETE\n")
cat("=============================================================================\n")