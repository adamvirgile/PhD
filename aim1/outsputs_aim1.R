# =============================================================================
# AIM 1: FOCUSED OUTPUTS (Scatterplots + Table 1 + rmcorr Tables)
# =============================================================================
# Run AFTER analysis_aim1.R (needs aim1_analysis_data.csv in working dir)
#
# Produces:
#   SCATTERPLOTS (8 PNGs): 4 intensity + 2 load + 2 volume-unadjusted
#   TABLES: table1, rmcorr_overall, rmcorr_practice, rmcorr_game
# =============================================================================

setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim1")

library(dplyr)
library(tidyr)
library(ggplot2)
library(rmcorr)

# =============================================================================
# LOAD + PREP DATA
# =============================================================================
aim1_data <- read.csv("aim1_analysis_data.csv")

aim1_data <- aim1_data %>%
  mutate(
    position_clean = tolower(trimws(as.character(position))),
    position_group = case_when(
      position_clean %in% c("forward", "forwards", "f", "fw", "fwd") ~ "Forward",
      position_clean %in% c("defense", "defenceman", "defenseman", "defensemen", "d") ~ "Defenseman",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(position_group))

# De-identified IDs
athlete_ids <- aim1_data %>%
  distinct(athlete_name) %>% arrange(athlete_name) %>%
  mutate(athlete_id = sprintf("P%02d", row_number()))
write.csv(athlete_ids, "PRIVATE_athlete_id_mapping.csv", row.names = FALSE)

aim1_data <- aim1_data %>% left_join(athlete_ids, by = "athlete_name")

aim1_data <- aim1_data %>%
  mutate(Context = factor(ifelse(context == "practice", "Practice", "Game"),
                          levels = c("Practice", "Game")))

df <- aim1_data %>% filter(position_group %in% c("Forward", "Defenseman"))

# --- Exposure-corrected sRPE-TL (TOI for games, duration for practices) ---
if (!"duration_mins_exposure" %in% names(df)) {
  df <- df %>% mutate(
    duration_mins_exposure = case_when(
      context == "game" & !is.na(toi_total_mins) & toi_total_mins > 0 ~ toi_total_mins,
      context == "practice" & !is.na(duration_mins) & duration_mins > 0 ~ duration_mins,
      TRUE ~ NA_real_))
}
if (!"sRPE_TL_exposure" %in% names(df)) {
  df <- df %>% mutate(
    sRPE_TL_exposure = ifelse(!is.na(sRPE) & !is.na(duration_mins_exposure),
                              sRPE * duration_mins_exposure, NA_real_))
}

cat("=== Data loaded ===\n")
cat("N obs:", nrow(df), "| Athletes:", n_distinct(df$athlete_id), "\n")
cat("Practices:", sum(df$context == "practice"), "| Games:", sum(df$context == "game"), "\n\n")

# =============================================================================
# STYLE
# =============================================================================
theme_scatter <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.title = element_text(size = 10), axis.text = element_text(size = 9),
        legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray92"),
        plot.margin = margin(12, 15, 12, 12))

ctx_colors <- c("Practice" = "#2166AC", "Game" = "#B2182B")

make_scatter <- function(data, x_var, y_var, x_lab, y_lab, title, subtitle,
                         filename, width = 7, height = 6) {
  d <- data %>% filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))
  if (nrow(d) == 0) { cat("  Skipped:", filename, "\n"); return(invisible(NULL)) }
  x_max <- max(d[[x_var]], na.rm = TRUE) * 1.05
  y_max <- max(d[[y_var]], na.rm = TRUE) * 1.05
  p <- ggplot(d, aes(x = .data[[x_var]], y = .data[[y_var]], color = Context)) +
    geom_point(alpha = 0.2, size = 1, stroke = 0) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
    scale_color_manual(values = ctx_colors) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab, color = "Context") +
    theme_scatter
  ggsave(filename, p, width = width, height = height, dpi = 300, bg = "white")
  cat("  Saved:", filename, "\n")
}

# =============================================================================
# 1. SCATTERPLOTS
# =============================================================================
cat("=== Generating scatterplots ===\n")
n_ath <- n_distinct(df$athlete_id)
n_obs <- nrow(df)
sub_base <- paste0("N = ", n_ath, " athletes, ", n_obs, " sessions. Lines = linear fit +/- 95% CI.")
sub_tl <- paste0(sub_base, "\nNote: sRPE-TL uses practice duration (practice) and TOI (game).")

# Intensity pairs (proposal paragraph 380: TRIMPmin, HRmean, HR80%, HRmax)
make_scatter(df, "TRIMP_min", "sRPE", "TRIMP/min (AU/min)", "sRPE-M (0-10)",
             "TRIMP/min vs sRPE-M", sub_base, "scatter_01_TRIMPmin_vs_sRPE.png")

make_scatter(df, "HR80_pct_rd", "sRPE", "Time >80% HRmax (% of session)", "sRPE-M (0-10)",
             "Time >80% HRmax (%) vs sRPE-M", sub_base, "scatter_02_pctHR80_vs_sRPE.png")

make_scatter(df, "HR_avg_pct_max", "sRPE", "Mean HR (% HRmax)", "sRPE-M (0-10)",
             "Mean HR vs sRPE-M", sub_base, "scatter_03_HRmean_vs_sRPE.png")

make_scatter(df, "HR_peak", "sRPE", "Peak HR (bpm)", "sRPE-M (0-10)",
             "Peak HR vs sRPE-M", sub_base, "scatter_04_HRpeak_vs_sRPE.png")

# Load pairs (proposal: sRPE-TL with TRIMP and HR80)
make_scatter(df, "TRIMP", "sRPE_TL_exposure", "TRIMP (AU)", "sRPE-TL (AU)",
             "TRIMP vs sRPE-TL", sub_tl, "scatter_05_TRIMP_vs_sRPETL.png")

make_scatter(df, "HR80_mins", "sRPE_TL_exposure", "Time >80% HRmax (min)", "sRPE-TL (AU)",
             "Time >80% HRmax (min) vs sRPE-TL", sub_tl, "scatter_06_HR80mins_vs_sRPETL.png")

# Volume-unadjusted load pairs (TRIMP/HR80 vs raw sRPE)
make_scatter(df, "TRIMP", "sRPE", "TRIMP (AU)", "sRPE-M (0-10)",
             "TRIMP vs sRPE-M", sub_base, "scatter_07_TRIMP_vs_sRPE.png")

make_scatter(df, "HR80_mins", "sRPE", "Time >80% HRmax (min)", "sRPE-M (0-10)",
             "Time >80% HRmax (min) vs sRPE-M", sub_base, "scatter_08_HR80mins_vs_sRPE.png")

# =============================================================================
# 2. TABLE 1: PARTICIPANT DESCRIPTIVES
# =============================================================================
cat("\n=== Generating Table 1 ===\n")

fmt_mean_nr <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  paste0(round(m, 1), " [", round(m - s, 1), ", ", round(m + s, 1), "]")
}

vars_to_summarize <- c("duration_mins_exposure", "sRPE", "sRPE_TL_exposure",
                       "TRIMP", "TRIMP_min", "HR_avg_pct_max", "HR_peak",
                       "HR80_pct_rd", "HR80_mins")

athlete_desc <- df %>%
  group_by(athlete_id, position_group, Context) %>%
  summarise(n_sessions = n(),
            across(all_of(vars_to_summarize), fmt_mean_nr, .names = "{.col}"),
            .groups = "drop")

t1_prac <- athlete_desc %>% filter(Context == "Practice") %>%
  rename_with(~ paste0(.x, "_Practice"), all_of(c("n_sessions", vars_to_summarize))) %>%
  select(-Context)
t1_game <- athlete_desc %>% filter(Context == "Game") %>%
  rename_with(~ paste0(.x, "_Game"), all_of(c("n_sessions", vars_to_summarize))) %>%
  select(-Context)
table1 <- t1_prac %>% left_join(t1_game, by = c("athlete_id", "position_group")) %>%
  arrange(position_group, athlete_id)

# Position totals
pos_tot <- df %>%
  group_by(position_group, Context) %>%
  summarise(n_sessions = n(),
            across(all_of(vars_to_summarize), fmt_mean_nr, .names = "{.col}"),
            .groups = "drop") %>%
  mutate(athlete_id = paste0(position_group, " Total"))

pt_prac <- pos_tot %>% filter(Context == "Practice") %>%
  rename_with(~ paste0(.x, "_Practice"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
pt_game <- pos_tot %>% filter(Context == "Game") %>%
  rename_with(~ paste0(.x, "_Game"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
pos_combined <- pt_prac %>% left_join(pt_game, by = c("athlete_id", "position_group"))

# Grand totals
grand_tot <- df %>%
  group_by(Context) %>%
  summarise(n_sessions = n(),
            across(all_of(vars_to_summarize), fmt_mean_nr, .names = "{.col}"),
            .groups = "drop") %>%
  mutate(athlete_id = "All Skaters", position_group = "All")

gt_prac <- grand_tot %>% filter(Context == "Practice") %>%
  rename_with(~ paste0(.x, "_Practice"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
gt_game <- grand_tot %>% filter(Context == "Game") %>%
  rename_with(~ paste0(.x, "_Game"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
grand_combined <- gt_prac %>% left_join(gt_game, by = c("athlete_id", "position_group"))

table1_full <- bind_rows(table1, pos_combined, grand_combined) %>%
  rename(Participant = athlete_id, Position = position_group)

write.csv(table1_full, "table1_participant_descriptives.csv", row.names = FALSE)
cat("  Saved: table1_participant_descriptives.csv\n")

# =============================================================================
# 3. RMCORR TABLES (Overall, Practice, Game)
# =============================================================================
cat("\n=== Generating rmcorr tables ===\n")

# ALL 8 pairs matching the scatterplots
rmcorr_pairs <- list(
  list(x = "TRIMP_min",       y = "sRPE",             label = "TRIMP/min vs sRPE"),
  list(x = "HR80_pct_rd",     y = "sRPE",             label = "%HR80 vs sRPE"),
  list(x = "HR_avg_pct_max",  y = "sRPE",             label = "HRmean vs sRPE"),
  list(x = "HR_peak",         y = "sRPE",             label = "HRpeak vs sRPE"),
  list(x = "TRIMP",           y = "sRPE",             label = "TRIMP vs sRPE"),
  list(x = "HR80_mins",       y = "sRPE",             label = "HR80 mins vs sRPE"),
  list(x = "TRIMP",           y = "sRPE_TL_exposure", label = "TRIMP vs sRPE-TL"),
  list(x = "HR80_mins",       y = "sRPE_TL_exposure", label = "HR80 mins vs sRPE-TL")
)

compute_athlete_r <- function(data, x_var, y_var, pair_label) {
  data %>%
    filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]])) %>%
    group_by(athlete_id, position_group) %>%
    summarise(
      n_pair = n(),
      r = tryCatch({ if (n() < 5) return(NA_real_)
        unname(cor.test(.data[[x_var]], .data[[y_var]])$estimate) },
        error = function(e) NA_real_),
      ci_lo = tryCatch({ if (n() < 5) return(NA_real_)
        cor.test(.data[[x_var]], .data[[y_var]])$conf.int[1] },
        error = function(e) NA_real_),
      ci_hi = tryCatch({ if (n() < 5) return(NA_real_)
        cor.test(.data[[x_var]], .data[[y_var]])$conf.int[2] },
        error = function(e) NA_real_),
      .groups = "drop") %>%
    mutate(pair = pair_label)
}

compute_group_rmcorr <- function(data, x_var, y_var, pair_label, group_label) {
  d <- data %>% filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]])) %>%
    mutate(participant = factor(athlete_id))
  n_subj <- n_distinct(d$participant)
  if (nrow(d) < 6 || n_subj < 3) {
    return(data.frame(group = group_label, pair = pair_label,
                      n_obs = nrow(d), n_athletes = n_subj,
                      r = NA, ci_lo = NA, ci_hi = NA, p = NA))
  }
  out <- tryCatch(
    rmcorr(participant = participant, measure1 = d[[x_var]],
           measure2 = d[[y_var]], dataset = d),
    error = function(e) NULL)
  if (is.null(out)) {
    return(data.frame(group = group_label, pair = pair_label,
                      n_obs = nrow(d), n_athletes = n_subj,
                      r = NA, ci_lo = NA, ci_hi = NA, p = NA))
  }
  data.frame(group = group_label, pair = pair_label,
             n_obs = nrow(d), n_athletes = n_subj,
             r = round(out$r, 3), ci_lo = round(out$CI[1], 3),
             ci_hi = round(out$CI[2], 3), p = out$p)
}

build_rmcorr_table <- function(data, context_label) {
  athlete_r <- bind_rows(lapply(rmcorr_pairs, function(p)
    compute_athlete_r(data, p$x, p$y, p$label)))
  
  athlete_r_fmt <- athlete_r %>%
    mutate(r_fmt = ifelse(is.na(r),
                          ifelse(n_pair < 5, paste0("n=", n_pair, "; insufficient"), "insufficient"),
                          paste0("n=", n_pair, "; ", round(r, 2), " [", round(ci_lo, 2), ", ", round(ci_hi, 2), "]")))
  
  athlete_wide <- athlete_r_fmt %>%
    select(athlete_id, position_group, pair, r_fmt) %>%
    pivot_wider(id_cols = c(athlete_id, position_group), names_from = pair, values_from = r_fmt) %>%
    arrange(position_group, athlete_id)
  
  # Helper for n on summary rows
  get_pair_n <- function(d, x, y) d %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>% nrow()
  
  # Position summary
  pos_summary <- athlete_r %>% filter(!is.na(r)) %>%
    group_by(position_group, pair) %>%
    summarise(
      k = n(),
      r_summary = {
        m <- mean(r, na.rm = TRUE); s <- sd(r, na.rm = TRUE)
        paste0("k=", n(), "; ", round(m, 2), " [", round(m - s, 2), ", ", round(m + s, 2), "]")
      }, .groups = "drop") %>%
    select(-k) %>%
    pivot_wider(id_cols = position_group, names_from = pair, values_from = r_summary) %>%
    mutate(athlete_id = paste0(position_group, " Mean [+/-1SD]"))
  pos_summary <- pos_summary %>%
    select(athlete_id, position_group, everything())
  
  # Overall summary
  overall_summary <- athlete_r %>% filter(!is.na(r)) %>%
    group_by(pair) %>%
    summarise(
      k = n(),
      r_summary = {
        m <- mean(r, na.rm = TRUE); s <- sd(r, na.rm = TRUE)
        paste0("k=", n(), "; ", round(m, 2), " [", round(m - s, 2), ", ", round(m + s, 2), "]")
      }, .groups = "drop") %>%
    select(-k) %>%
    pivot_wider(names_from = pair, values_from = r_summary) %>%
    mutate(athlete_id = "All Skaters Mean [+/-1SD]", position_group = "All")
  
  # Group-level rmcorr
  rmcorr_all <- bind_rows(lapply(rmcorr_pairs, function(p)
    compute_group_rmcorr(data, p$x, p$y, p$label, "All Skaters rmcorr")))
  rmcorr_fwd <- bind_rows(lapply(rmcorr_pairs, function(p)
    compute_group_rmcorr(data %>% filter(position_group == "Forward"),
                         p$x, p$y, p$label, "Forwards rmcorr")))
  rmcorr_def <- bind_rows(lapply(rmcorr_pairs, function(p)
    compute_group_rmcorr(data %>% filter(position_group == "Defenseman"),
                         p$x, p$y, p$label, "Defensemen rmcorr")))
  
  rmcorr_combined <- bind_rows(rmcorr_all, rmcorr_fwd, rmcorr_def) %>%
    mutate(display = ifelse(is.na(r), "insufficient",
                            paste0("n=", n_obs, "; ", r, " [", ci_lo, ", ", ci_hi, "]")))
  
  rmcorr_wide <- rmcorr_combined %>%
    select(group, pair, display) %>%
    pivot_wider(id_cols = group, names_from = pair, values_from = display) %>%
    rename(athlete_id = group) %>% mutate(position_group = athlete_id)
  
  final <- bind_rows(athlete_wide, pos_summary, overall_summary, rmcorr_wide) %>%
    rename(Participant = athlete_id, Position = position_group)
  return(final)
}

# Build and save 3 tables
table_overall <- build_rmcorr_table(df, "Overall")
write.csv(table_overall, "table_rmcorr_overall.csv", row.names = FALSE)
cat("  Saved: table_rmcorr_overall.csv (", nrow(table_overall), "rows)\n")

df_practice <- df %>% filter(context == "practice")
table_practice <- build_rmcorr_table(df_practice, "Practice")
write.csv(table_practice, "table_rmcorr_practice.csv", row.names = FALSE)
cat("  Saved: table_rmcorr_practice.csv (", nrow(table_practice), "rows)\n")

df_game <- df %>% filter(context == "game")
table_game <- build_rmcorr_table(df_game, "Game")
write.csv(table_game, "table_rmcorr_game.csv", row.names = FALSE)
cat("  Saved: table_rmcorr_game.csv (", nrow(table_game), "rows)\n")

# Raw values for further analysis
athlete_r_all <- bind_rows(
  bind_rows(lapply(rmcorr_pairs, function(p)
    compute_athlete_r(df, p$x, p$y, p$label))) %>% mutate(context = "overall"),
  bind_rows(lapply(rmcorr_pairs, function(p)
    compute_athlete_r(df_practice, p$x, p$y, p$label))) %>% mutate(context = "practice"),
  bind_rows(lapply(rmcorr_pairs, function(p)
    compute_athlete_r(df_game, p$x, p$y, p$label))) %>% mutate(context = "game")
)
write.csv(athlete_r_all, "table_athlete_r_raw_all_contexts.csv", row.names = FALSE)
cat("  Saved: table_athlete_r_raw_all_contexts.csv\n")

cat("\n=== ALL OUTPUTS COMPLETE ===\n")
cat("Scatterplots: scatter_01 through scatter_08\n")
cat("Table 1: table1_participant_descriptives.csv\n")
cat("rmcorr: table_rmcorr_overall/practice/game.csv\n")
