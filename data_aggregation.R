# =============================================================================
# DATA AGGREGATION SCRIPT FOR AIM 1
# =============================================================================
# Location: C:\Users\avirgile\Dropbox\Adam Only\School\PhD\UVM Coursework\Dissertation\data
# 
# This script reads raw data files from this folder and outputs:
#   - merged_data_full.csv (stays here in /data)
#   - aim1_analysis_data.csv (goes to /data/aim1)
#
# KEY DECISIONS:
#   - Practice = "Coaches Practice" only (not Captain's Practice, etc.)
#   - HR signal quality: <=30% error for practices, <=50% for games
#   - Game exposure = time-on-ice (TOI), not standardized game duration
#   - sRPE_TL_exposure = sRPE * exposure time (TOI for games, duration for practices)
#   - Athlete eligibility: >=55 total sessions with valid data + full participation
# =============================================================================

setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data")

library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
to_minutes <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (inherits(x, "difftime")) return(as.numeric(x, units = "mins"))
  if (inherits(x, "hms")) return(as.numeric(x) / 60)
  if (inherits(x, "POSIXt")) return(hour(x) * 60 + minute(x) + second(x) / 60)
  if (is.character(x)) {
    x_trim <- trimws(x)
    out <- rep(NA_real_, length(x_trim))
    is_blank <- is.na(x_trim) | x_trim == ""
    if (all(is_blank)) return(out)
    parts <- str_split(x_trim[!is_blank], ":", simplify = TRUE)
    if (ncol(parts) == 2) {
      mm <- suppressWarnings(as.numeric(parts[, 1]))
      ss <- suppressWarnings(as.numeric(parts[, 2]))
      mins <- mm + ss / 60
    } else if (ncol(parts) == 3) {
      hh <- suppressWarnings(as.numeric(parts[, 1]))
      mm <- suppressWarnings(as.numeric(parts[, 2]))
      ss <- suppressWarnings(as.numeric(parts[, 3]))
      mins <- hh * 60 + mm + ss / 60
    } else {
      mins <- rep(NA_real_, nrow(parts))
    }
    out[!is_blank] <- mins
    return(out)
  }
  if (is.numeric(x)) {
    out <- as.numeric(x)
    mx <- suppressWarnings(max(out, na.rm = TRUE))
    mn <- suppressWarnings(min(out, na.rm = TRUE))
    if (is.finite(mx) && mx <= 1.5 && is.finite(mn) && mn >= 0) return(out * 1440)
    if (is.finite(mx) && mx > 1000 && mx < 100000) return(out / 60)
    return(out)
  }
  suppressWarnings(as.numeric(x))
}

rescale_to_pct_0_100 <- function(x) {
  if (!is.numeric(x)) x <- suppressWarnings(as.numeric(x))
  mx <- suppressWarnings(max(x, na.rm = TRUE))
  if (is.finite(mx) && mx <= 1.5) return(x * 100)
  x
}

clip <- function(x, lo = -Inf, hi = Inf) { x[x < lo] <- lo; x[x > hi] <- hi; x }

# =============================================================================
# READ DATA
# =============================================================================
cat("=== READING DATA FILES ===\n")
cat("Working directory:", getwd(), "\n\n")

rd_data <- read_excel("uvm_data.xlsx", sheet = "Sheet1")
fb_data <- read_excel("firstbeat_data.xlsx", sheet = "Measurements")
athlete_profiles <- read_excel("athlete_profiles.xlsx")
team_schedule <- read_excel("team_schedule.xlsx")

# =============================================================================
# ATHLETE PROFILES
# =============================================================================
profiles_clean <- athlete_profiles %>%
  mutate(athlete_name = paste(FirstName, LastName)) %>%
  select(athlete_name, position = Position, class = Class) %>%
  mutate(athlete_name = str_trim(athlete_name))

cat("=== ATHLETE PROFILES ===\n")
cat("Athletes:", nrow(profiles_clean), "\n")
print(table(profiles_clean$position))

# =============================================================================
# TEAM SCHEDULE
# =============================================================================
schedule_clean <- team_schedule %>%
  mutate(session_date = as.Date(Date), schedule_event = Event) %>%
  select(session_date, schedule_event, Opponent, `Home/Road`, `Season Phase`) %>%
  rename(opponent = Opponent, home_road = `Home/Road`, season_phase = `Season Phase`)

# =============================================================================
# FIRSTBEAT DATA (HR)
# =============================================================================
fb_clean <- fb_data %>%
  rename(
    athlete_name = `Athlete name`,
    session_datetime = `Start date (dd.mm.yyyy)`,
    duration_raw = `Duration (hh:mm:ss)`,
    session_type = Sport,
    measurement_error = `Measurement error (%)`,
    artifact_pct = `Artifact percentage (%)`,
    HR_avg = `Average HR (bpm)`,
    HR_avg_pct_max = `Average %HRmax (%)`,
    HR_peak = `Peak HR (bpm)`,
    HR_peak_pct_max = `Peak %HRmax (%)`,
    HR_min = `Minimum HR (bpm)`,
    TRIMP = `TRIMP (Index)`,
    TRIMP_min = `TRIMP/min (Index)`,
    training_effect = `Training Effect (0-5)`,
    aerobic_TE = `Aerobic TE (0.0 - 5.0)`,
    anaerobic_TE = `Anaerobic TE (0.0 - 5.0)`,
    recovery_time = `Recovery training (hh:mm:ss)`,
    aerobic_z1_time = `Aerobic zone 1 (hh:mm:ss)`,
    aerobic_z2_time = `Aerobic zone 2 (hh:mm:ss)`,
    anaerobic_threshold_time = `Anaerobic threshold zone (hh:mm:ss)`,
    high_intensity_time = `High intensity training (hh:mm:ss)`
  ) %>%
  mutate(
    athlete_name = str_trim(athlete_name),
    session_date = as.Date(session_datetime),
    session_category = case_when(
      session_type %in% c("Game", "Exh Game", "Alumni Scrimmage") ~ "game",
      session_type == "Coaches Practice" ~ "practice",
      session_type == "Pre-Game Skate" ~ "pregame_skate",
      TRUE ~ "other"
    )
  ) %>%
  filter(session_category %in% c("practice", "game")) %>%
  mutate(
    context = session_category,
    context_binary = ifelse(session_category == "game", 1, 0),
    duration_mins_fb = to_minutes(duration_raw),
    anaerobic_threshold_mins = to_minutes(anaerobic_threshold_time),
    high_intensity_mins = to_minutes(high_intensity_time),
    HR_avg_pct_max = rescale_to_pct_0_100(HR_avg_pct_max),
    HR_peak_pct_max = rescale_to_pct_0_100(HR_peak_pct_max)
  ) %>%
  mutate(
    duration_mins_fb = ifelse(is.na(duration_mins_fb) | duration_mins_fb <= 0 |
                                duration_mins_fb > 300, NA_real_, duration_mins_fb),
    anaerobic_threshold_mins = ifelse(is.na(anaerobic_threshold_mins) |
                                        anaerobic_threshold_mins < 0, NA_real_,
                                      anaerobic_threshold_mins),
    high_intensity_mins = ifelse(is.na(high_intensity_mins) |
                                   high_intensity_mins < 0, NA_real_,
                                 high_intensity_mins),
    HR80_mins = ifelse(is.na(anaerobic_threshold_mins) | is.na(high_intensity_mins),
                       NA_real_, anaerobic_threshold_mins + high_intensity_mins),
    HR90_mins = high_intensity_mins,
    HR80_pct_fb = ifelse(!is.na(HR80_mins) & !is.na(duration_mins_fb) &
                           duration_mins_fb > 0,
                         (HR80_mins / duration_mins_fb) * 100, NA_real_),
    HR90_pct_fb = ifelse(!is.na(HR90_mins) & !is.na(duration_mins_fb) &
                           duration_mins_fb > 0,
                         (HR90_mins / duration_mins_fb) * 100, NA_real_)
  ) %>%
  mutate(
    HR80_mins = ifelse(!is.na(HR80_mins) & !is.na(duration_mins_fb) &
                         HR80_mins > (duration_mins_fb + 0.25), NA_real_, HR80_mins),
    HR90_mins = ifelse(!is.na(HR90_mins) & !is.na(duration_mins_fb) &
                         HR90_mins > (duration_mins_fb + 0.25), NA_real_, HR90_mins),
    HR80_pct_fb = ifelse(!is.na(HR80_pct_fb), clip(HR80_pct_fb, 0, 100), NA_real_),
    HR90_pct_fb = ifelse(!is.na(HR90_pct_fb), clip(HR90_pct_fb, 0, 100), NA_real_)
  ) %>%
  select(
    athlete_name, session_date, session_type, session_category, context, context_binary,
    measurement_error, artifact_pct,
    HR_avg, HR_avg_pct_max, HR_peak, HR_peak_pct_max, HR_min,
    TRIMP, TRIMP_min, training_effect, aerobic_TE, anaerobic_TE,
    duration_mins_fb, HR80_mins, HR80_pct_fb, HR90_mins, HR90_pct_fb
  )

cat("\n=== FIRSTBEAT DATA ===\n")
cat("Total rows:", nrow(fb_clean), "\n")
cat("Games:", sum(fb_clean$session_category == "game"), "\n")
cat("Practices:", sum(fb_clean$session_category == "practice"), "\n")

# =============================================================================
# sRPE DATA
# =============================================================================
srpe_game <- rd_data %>%
  filter(MetricName == "sRPE Game") %>%
  select(athlete_name = PlayerName, session_date = Date, sRPE_TL_game = Value) %>%
  mutate(session_date = as.Date(session_date), athlete_name = str_trim(athlete_name))

game_duration <- rd_data %>%
  filter(MetricName == "GAME DURATION") %>%
  select(athlete_name = PlayerName, session_date = Date, game_duration = Value) %>%
  mutate(session_date = as.Date(session_date), athlete_name = str_trim(athlete_name))

game_srpe <- srpe_game %>%
  left_join(game_duration, by = c("athlete_name", "session_date")) %>%
  filter(!is.na(game_duration) & game_duration > 0) %>%
  mutate(sRPE = sRPE_TL_game / game_duration, session_category = "game") %>%
  select(athlete_name, session_date, session_category, sRPE,
         sRPE_TL = sRPE_TL_game, duration_mins = game_duration)

srpe_practice <- rd_data %>%
  filter(MetricName == "sRPE Practice") %>%
  select(athlete_name = PlayerName, session_date = Date, sRPE_TL_practice = Value) %>%
  mutate(session_date = as.Date(session_date), athlete_name = str_trim(athlete_name))

practice_duration <- rd_data %>%
  filter(MetricName == "PRACTICE DURATION") %>%
  select(athlete_name = PlayerName, session_date = Date, practice_duration = Value) %>%
  mutate(session_date = as.Date(session_date), athlete_name = str_trim(athlete_name))

practice_srpe <- srpe_practice %>%
  left_join(practice_duration, by = c("athlete_name", "session_date")) %>%
  filter(!is.na(practice_duration) & practice_duration > 0) %>%
  mutate(sRPE = sRPE_TL_practice / practice_duration, session_category = "practice") %>%
  select(athlete_name, session_date, session_category, sRPE,
         sRPE_TL = sRPE_TL_practice, duration_mins = practice_duration)

srpe_combined <- bind_rows(game_srpe, practice_srpe)

cat("\n=== sRPE DATA ===\n")
cat("Game sRPE rows:", nrow(game_srpe), "\n")
cat("Practice sRPE rows:", nrow(practice_srpe), "\n")

# =============================================================================
# WELLBEING DATA
# =============================================================================
wellbeing_metrics <- c("FATIGUE LEVEL", "SLEEP QUALITY", "STRESS LEVEL",
                       "EXCITEMENT FOR TRAINING",
                       "CONFIDENCE FOR SPORT-SPECIFIC PERFORMANCE",
                       "ARE YOU FEELING SICK?")

wellbeing_data <- rd_data %>%
  filter(MetricName %in% wellbeing_metrics) %>%
  select(athlete_name = PlayerName, session_date = Date,
         metric_name = MetricName, value = Value) %>%
  mutate(
    session_date = as.Date(session_date),
    athlete_name = str_trim(athlete_name),
    metric_clean = case_when(
      metric_name == "FATIGUE LEVEL" ~ "fatigue",
      metric_name == "SLEEP QUALITY" ~ "sleep_quality",
      metric_name == "STRESS LEVEL" ~ "stress",
      metric_name == "EXCITEMENT FOR TRAINING" ~ "excitement",
      metric_name == "CONFIDENCE FOR SPORT-SPECIFIC PERFORMANCE" ~ "confidence",
      metric_name == "ARE YOU FEELING SICK?" ~ "sick"
    )
  ) %>%
  select(-metric_name) %>%
  group_by(athlete_name, session_date, metric_clean) %>%
  summarise(value = first(value), .groups = "drop") %>%
  pivot_wider(names_from = metric_clean, values_from = value)

# =============================================================================
# PARTICIPATION STATUS
# =============================================================================
participation_data <- rd_data %>%
  filter(MetricName == "PARTICIPATION STATUS") %>%
  select(athlete_name = PlayerName, session_date = Date,
         participation_value = Value, participation_display = Display) %>%
  mutate(session_date = as.Date(session_date), athlete_name = str_trim(athlete_name)) %>%
  group_by(athlete_name, session_date) %>%
  summarise(participation_value = first(participation_value),
            participation_display = first(participation_display), .groups = "drop") %>%
  mutate(full_participation = (participation_value == 4))

# =============================================================================
# TIME ON ICE DATA (Games only)
# =============================================================================
toi_data <- rd_data %>%
  filter(MetricName %in% c("GAME TOTAL TOI", "GAME TOTAL SHIFTS",
                           "GAME P1 TOI", "GAME P2 TOI", "GAME P3 TOI",
                           "GAME TOTAL EV TOI", "GAME TOTAL PP TOI", "GAME TOTAL PK TOI")) %>%
  select(athlete_name = PlayerName, session_date = Date,
         metric_name = MetricName, value = Value) %>%
  mutate(
    session_date = as.Date(session_date),
    athlete_name = str_trim(athlete_name),
    metric_clean = case_when(
      metric_name == "GAME TOTAL TOI" ~ "toi_total",
      metric_name == "GAME TOTAL SHIFTS" ~ "shifts_total",
      metric_name == "GAME P1 TOI" ~ "toi_p1",
      metric_name == "GAME P2 TOI" ~ "toi_p2",
      metric_name == "GAME P3 TOI" ~ "toi_p3",
      metric_name == "GAME TOTAL EV TOI" ~ "toi_ev",
      metric_name == "GAME TOTAL PP TOI" ~ "toi_pp",
      metric_name == "GAME TOTAL PK TOI" ~ "toi_pk"
    )
  ) %>%
  select(-metric_name) %>%
  group_by(athlete_name, session_date, metric_clean) %>%
  summarise(value = first(value), .groups = "drop") %>%
  pivot_wider(names_from = metric_clean, values_from = value) %>%
  mutate(
    toi_total_mins = toi_total / 60000,
    toi_p1_mins = toi_p1 / 60000,
    toi_p2_mins = toi_p2 / 60000,
    toi_p3_mins = toi_p3 / 60000,
    toi_ev_mins = toi_ev / 60000,
    toi_pp_mins = toi_pp / 60000,
    toi_pk_mins = toi_pk / 60000
  ) %>%
  select(athlete_name, session_date,
         toi_total_mins, toi_p1_mins, toi_p2_mins, toi_p3_mins,
         toi_ev_mins, toi_pp_mins, toi_pk_mins, shifts_total)

# =============================================================================
# PERFORMANCE RATINGS
# =============================================================================
self_performance <- rd_data %>%
  filter(MetricName %in% c("Rate Your Performance", "GRADE YOUR PERFORMANCE")) %>%
  select(athlete_name = PlayerName, session_date = Date, self_performance = Value) %>%
  mutate(session_date = as.Date(session_date), athlete_name = str_trim(athlete_name)) %>%
  group_by(athlete_name, session_date) %>%
  summarise(self_performance = first(self_performance), .groups = "drop")

coach_performance <- rd_data %>%
  filter(MetricName == "Grade Player Performance") %>%
  select(athlete_name = PlayerName, session_date = Date, coach_performance = Value) %>%
  mutate(session_date = as.Date(session_date), athlete_name = str_trim(athlete_name)) %>%
  group_by(athlete_name, session_date) %>%
  summarise(coach_performance = first(coach_performance), .groups = "drop")

# =============================================================================
# MERGE ALL DATA SOURCES
# =============================================================================
merged_data <- fb_clean %>%
  left_join(srpe_combined, by = c("athlete_name", "session_date", "session_category")) %>%
  left_join(wellbeing_data, by = c("athlete_name", "session_date")) %>%
  left_join(participation_data, by = c("athlete_name", "session_date")) %>%
  left_join(self_performance, by = c("athlete_name", "session_date")) %>%
  left_join(coach_performance, by = c("athlete_name", "session_date")) %>%
  left_join(toi_data, by = c("athlete_name", "session_date")) %>%
  left_join(profiles_clean, by = "athlete_name") %>%
  left_join(schedule_clean, by = "session_date")

# =============================================================================
# DERIVED METRICS
# =============================================================================
# KEY: For games, exposure = time-on-ice (TOI). For practices, exposure = duration.
merged_data <- merged_data %>%
  mutate(
    HR80_pct_rd = ifelse(!is.na(HR80_mins) & !is.na(duration_mins) & duration_mins > 0,
                         (HR80_mins / duration_mins) * 100, NA_real_),
    HR90_pct_rd = ifelse(!is.na(HR90_mins) & !is.na(duration_mins) & duration_mins > 0,
                         (HR90_mins / duration_mins) * 100, NA_real_),
    # Exposure-corrected duration (PRIMARY)
    duration_mins_exposure = case_when(
      context == "game" & !is.na(toi_total_mins) & toi_total_mins > 0 ~ toi_total_mins,
      context == "practice" & !is.na(duration_mins) & duration_mins > 0 ~ duration_mins,
      TRUE ~ NA_real_
    ),
    # sRPE-TL variants
    sRPE_TL_exposure = ifelse(!is.na(sRPE) & !is.na(duration_mins_exposure),
                              sRPE * duration_mins_exposure, NA_real_),
    sRPE_TL_duration = ifelse(!is.na(sRPE) & !is.na(duration_mins) & duration_mins > 0,
                              sRPE * duration_mins, NA_real_),
    sRPE_TL_toi = case_when(
      context == "game" & !is.na(sRPE) & !is.na(toi_total_mins) & toi_total_mins > 0
      ~ sRPE * toi_total_mins,
      TRUE ~ NA_real_
    ),
    category_mismatch = case_when(
      session_category == "game" & schedule_event == "Game" ~ FALSE,
      session_category == "practice" & schedule_event == "Practice" ~ FALSE,
      is.na(schedule_event) ~ NA,
      TRUE ~ TRUE
    )
  )

cat("\n=== MERGED DATA ===\n")
cat("Total rows:", nrow(merged_data), "\n")
cat("Unique athletes:", n_distinct(merged_data$athlete_name), "\n")

# =============================================================================
# AIM 1 ANALYSIS-READY DATASET
# =============================================================================

# --- Session-level filters ---
n_before_hr_filter <- nrow(merged_data %>%
                             filter(!is.na(sRPE) & !is.na(TRIMP_min)) %>%
                             filter(full_participation == TRUE) %>%
                             filter(!is.na(position) & position != "Goalie"))

aim1_data <- merged_data %>%
  filter(!is.na(sRPE) & !is.na(TRIMP_min)) %>%
  filter(full_participation == TRUE) %>%
  filter(!is.na(position) & position != "Goalie") %>%
  filter(
    is.na(measurement_error) |
      (context == "practice" & measurement_error <= 0.30) |
      (context == "game"     & measurement_error <= 0.50)
  )

n_after_hr_filter <- nrow(aim1_data)
cat("\n=== HR SIGNAL QUALITY FILTER ===\n")
cat("Before:", n_before_hr_filter, "| Excluded:", n_before_hr_filter - n_after_hr_filter,
    "| After:", n_after_hr_filter, "\n")

# --- Athlete-level eligibility (Section 2.5: >=55 total sessions) ---
athlete_counts <- aim1_data %>%
  group_by(athlete_name) %>%
  summarise(n_total = n(), n_practice = sum(context == "practice"),
            n_game = sum(context == "game"), .groups = "drop")

MIN_SESSIONS <- 55
eligible <- athlete_counts %>% filter(n_total >= MIN_SESSIONS) %>% pull(athlete_name)

cat("\n=== ATHLETE ELIGIBILITY (>= ", MIN_SESSIONS, " sessions) ===\n")
cat("Eligible:", length(eligible), "of", nrow(athlete_counts), "athletes\n")
cat("Excluded athletes:\n")
print(athlete_counts %>% filter(n_total < MIN_SESSIONS))

aim1_data <- aim1_data %>% filter(athlete_name %in% eligible)

# --- Within/between decomposition ---
aim1_data <- aim1_data %>%
  group_by(athlete_name) %>%
  mutate(
    n_sessions = n(),
    TRIMP_min_between = mean(TRIMP_min, na.rm = TRUE),
    HR_avg_pct_max_between = mean(HR_avg_pct_max, na.rm = TRUE),
    HR_peak_between = mean(HR_peak, na.rm = TRUE),
    HR80_pct_rd_between = mean(HR80_pct_rd, na.rm = TRUE),
    HR80_mins_between = mean(HR80_mins, na.rm = TRUE),
    TRIMP_min_within = TRIMP_min - TRIMP_min_between,
    HR_avg_pct_max_within = HR_avg_pct_max - HR_avg_pct_max_between,
    HR_peak_within = HR_peak - HR_peak_between,
    HR80_pct_rd_within = HR80_pct_rd - HR80_pct_rd_between,
    HR80_mins_within = HR80_mins - HR80_mins_between
  ) %>%
  ungroup()

cat("\n=== AIM 1 FINAL DATASET ===\n")
cat("Observations:", nrow(aim1_data), "| Athletes:", n_distinct(aim1_data$athlete_name), "\n")
cat("Practices:", sum(aim1_data$context == "practice"),
    "| Games:", sum(aim1_data$context == "game"), "\n")
cat("Non-missing sRPE_TL_exposure:", sum(!is.na(aim1_data$sRPE_TL_exposure)), "\n")

# =============================================================================
# SAVE
# =============================================================================
write.csv(merged_data, "merged_data_full.csv", row.names = FALSE)
write.csv(aim1_data, "aim1/aim1_analysis_data.csv", row.names = FALSE)
cat("\nSaved: merged_data_full.csv, aim1/aim1_analysis_data.csv\n")
cat("=== DATA AGGREGATION COMPLETE ===\n")
