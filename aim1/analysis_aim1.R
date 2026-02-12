# =============================================================================
# AIM 1: COMPLETE ANALYSIS + OUTPUTS
# Construct Validity of sRPE-M in NCAA D1 Men's Ice Hockey
# =============================================================================
#
# This is a SINGLE script that runs the full Aim 1 pipeline:
#   PART A: Data prep + filtering (one truth)
#   PART B: Mixed-effects models (intensity, load, TOI sensitivity)
#   PART C: Result extraction + summary tables
#   PART D: Scatterplots
#   PART E: Table 1 + rmcorr tables
#
# KEY: Game sessions are restricted to those with valid TOI data throughout.
#
# Inputs:
#   aim1_analysis_data.csv (produced by data_aggregation.R)
#
# Outputs:
#   Model CSVs:  ALL_/FWD_/DEF_ prefixed results, summaries, diagnostics
#   Figures:     scatter_01 through scatter_08 (.png)
#   Tables:      table1_participant_descriptives.csv
#                table_rmcorr_overall/practice/game.csv
#                table_athlete_r_raw_all_contexts.csv
#   Models:      ALL_/FWD_/DEF_ aim1_models.RData
# =============================================================================

setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim1")

# =============================================================================
# PACKAGES
# =============================================================================
library(lme4)
library(lmerTest)
library(performance)
library(dplyr)
library(tidyr)
library(broom.mixed)
library(emmeans)
library(rmcorr)
library(ggplot2)

# =============================================================================
# PART A: DATA LOAD + PREP
# =============================================================================
aim1_data <- read.csv("aim1_analysis_data.csv")

cat("=== DATA LOADED ===\n")
cat("Working directory:", getwd(), "\n")
cat("Observations:", nrow(aim1_data), "\n")
cat("Athletes:", length(unique(aim1_data$athlete_name)), "\n")

# --- Position labels ---
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

# --- De-identified IDs ---
athlete_ids <- aim1_data %>%
  distinct(athlete_name) %>% arrange(athlete_name) %>%
  mutate(athlete_id = sprintf("P%02d", row_number()))
write.csv(athlete_ids, "PRIVATE_athlete_id_mapping.csv", row.names = FALSE)
aim1_data <- aim1_data %>% left_join(athlete_ids, by = "athlete_name")

# --- Context factors ---
aim1_data <- aim1_data %>%
  mutate(
    context_f = factor(context, levels = c("practice", "game")),
    Context = factor(ifelse(context == "practice", "Practice", "Game"),
                     levels = c("Practice", "Game"))
  )

# --- Exposure-corrected sRPE-TL (TOI for games, duration for practices) ---
if (!"duration_mins_exposure" %in% names(aim1_data)) {
  aim1_data <- aim1_data %>% mutate(
    duration_mins_exposure = case_when(
      context == "game" & !is.na(toi_total_mins) & toi_total_mins > 0 ~ toi_total_mins,
      context == "practice" & !is.na(duration_mins) & duration_mins > 0 ~ duration_mins,
      TRUE ~ NA_real_))
}
if (!"sRPE_TL_exposure" %in% names(aim1_data)) {
  aim1_data <- aim1_data %>% mutate(
    sRPE_TL_exposure = ifelse(!is.na(sRPE) & !is.na(duration_mins_exposure),
                              sRPE * duration_mins_exposure, NA_real_))
}

# --- Z-score intensity within/between ---
aim1_data <- aim1_data %>%
  mutate(
    TRIMP_min_within_z = scale(TRIMP_min_within)[, 1],
    HR_avg_pct_max_within_z = scale(HR_avg_pct_max_within)[, 1],
    HR80_pct_rd_within_z = scale(HR80_pct_rd_within)[, 1],
    HR_peak_within_z = scale(HR_peak_within)[, 1],
    TRIMP_min_between_z = scale(TRIMP_min_between)[, 1],
    HR_avg_pct_max_between_z = scale(HR_avg_pct_max_between)[, 1],
    HR80_pct_rd_between_z = scale(HR80_pct_rd_between)[, 1],
    HR_peak_between_z = scale(HR_peak_between)[, 1]
  )

# --- Z-score load within/between ---
aim1_data <- aim1_data %>%
  group_by(athlete_name) %>%
  mutate(
    TRIMP_within = TRIMP - mean(TRIMP, na.rm = TRUE),
    TRIMP_between = mean(TRIMP, na.rm = TRUE),
    HR80_mins_within_new = HR80_mins - mean(HR80_mins, na.rm = TRUE),
    HR80_mins_between_new = mean(HR80_mins, na.rm = TRUE)
  ) %>% ungroup() %>%
  mutate(
    TRIMP_within_z = scale(TRIMP_within)[, 1],
    TRIMP_between_z = scale(TRIMP_between)[, 1],
    HR80_mins_within_z = scale(HR80_mins_within_new)[, 1],
    HR80_mins_between_z = scale(HR80_mins_between_new)[, 1]
  )

# ===========================================================================
# CRITICAL FILTER: For games, only keep sessions with valid TOI
# This ensures intensity models and load models use the SAME game subset.
# Practices are unaffected (they don't use TOI).
# ===========================================================================
n_before_toi <- nrow(aim1_data)
aim1_data <- aim1_data %>%
  filter(context == "practice" | (!is.na(toi_total_mins) & toi_total_mins > 0))

cat("\n=== TOI FILTER (games only) ===\n")
cat("Before:", n_before_toi, "| After:", nrow(aim1_data),
    "| Dropped:", n_before_toi - nrow(aim1_data), "games without TOI\n")
cat("Practices:", sum(aim1_data$context == "practice"),
    "| Games:", sum(aim1_data$context == "game"), "\n")
cat("Athletes:", n_distinct(aim1_data$athlete_name), "\n")
cat("sRPE_TL_exposure non-missing:", sum(!is.na(aim1_data$sRPE_TL_exposure)), "\n\n")

# =============================================================================
# PART A2: HELPER FUNCTIONS
# =============================================================================

safe_lmer <- function(formula, data, REML = TRUE) {
  warnings_log <- character(0)
  out <- tryCatch(
    withCallingHandlers(
      lmer(formula, data = data, REML = REML),
      warning = function(w) {
        warnings_log <<- c(warnings_log, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      warnings_log <<- c(warnings_log, paste0("ERROR: ", conditionMessage(e)))
      return(NULL)
    }
  )
  attr(out, "fit_warnings") <- warnings_log
  return(out)
}

fit_rs_with_fallback <- function(formula_rs, formula_rs_uncorr, data, REML = TRUE) {
  m <- safe_lmer(formula_rs, data = data, REML = REML)
  warns <- attr(m, "fit_warnings")
  if (!is.null(m) && !isSingular(m, tol = 1e-4) && length(warns) == 0) {
    return(list(model = m, type = "RS_correlated", warnings = warns, fallback_used = "none"))
  }
  cat("    -> RS correlated had issues, trying uncorrelated...\n")
  m_uc <- safe_lmer(formula_rs_uncorr, data = data, REML = REML)
  warns_uc <- attr(m_uc, "fit_warnings")
  if (!is.null(m_uc) && !isSingular(m_uc, tol = 1e-4) && length(warns_uc) == 0) {
    return(list(model = m_uc, type = "RS_uncorrelated", warnings = warns_uc, fallback_used = "uncorrelated"))
  }
  if (!is.null(m)) return(list(model = m, type = "RS_correlated_with_issues", warnings = warns, fallback_used = "kept_with_issues"))
  if (!is.null(m_uc)) return(list(model = m_uc, type = "RS_uncorrelated_with_issues", warnings = warns_uc, fallback_used = "kept_uc_with_issues"))
  return(list(model = NULL, type = "RS_failed", warnings = c(warns, warns_uc), fallback_used = "failed"))
}

extract_model_results <- function(model, model_name, outcome_var, model_type) {
  if (is.null(model)) {
    return(data.frame(model = model_name, outcome = outcome_var, model_type = model_type,
                      term = NA_character_, estimate = NA_real_, std_error = NA_real_,
                      t_value = NA_real_, p_value = NA_real_, ci_lower = NA_real_,
                      ci_upper = NA_real_, R2_marginal = NA_real_, R2_conditional = NA_real_,
                      AIC = NA_real_, BIC = NA_real_, random_intercept_var = NA_real_,
                      random_slope_var = NA_real_, random_intercept_slope_cor = NA_real_,
                      residual_var = NA_real_))
  }
  fe <- summary(model)$coefficients
  ci <- tryCatch(confint(model, method = "Wald", parm = "beta_"), error = function(e) NULL)
  r2_vals <- tryCatch(performance::r2(model), error = function(e) NULL)
  vc <- as.data.frame(VarCorr(model))
  ri <- vc %>% filter(grp == "athlete_name" & var1 == "(Intercept)" & is.na(var2))
  rs <- vc %>% filter(grp == "athlete_name" & var1 != "(Intercept)" & is.na(var2))
  cor_is <- vc %>% filter(grp == "athlete_name" & !is.na(var2))
  resid <- vc %>% filter(grp == "Residual")
  results <- data.frame(model = model_name, outcome = outcome_var, model_type = model_type,
                        term = rownames(fe), estimate = fe[, "Estimate"],
                        std_error = fe[, "Std. Error"], t_value = fe[, "t value"],
                        p_value = fe[, "Pr(>|t|)"], ci_lower = NA_real_, ci_upper = NA_real_,
                        row.names = NULL)
  if (!is.null(ci)) { rn <- results$term; ok <- rn %in% rownames(ci)
  results$ci_lower[ok] <- ci[rn[ok], 1]; results$ci_upper[ok] <- ci[rn[ok], 2] }
  results$R2_marginal <- if (!is.null(r2_vals)) r2_vals$R2_marginal else NA_real_
  results$R2_conditional <- if (!is.null(r2_vals)) r2_vals$R2_conditional else NA_real_
  results$AIC <- AIC(model); results$BIC <- BIC(model)
  results$random_intercept_var <- if (nrow(ri) > 0) ri$vcov[1] else NA_real_
  results$random_slope_var <- if (nrow(rs) > 0) rs$vcov[1] else NA_real_
  results$random_intercept_slope_cor <- if (nrow(cor_is) > 0) cor_is$vcov[1] else NA_real_
  results$residual_var <- if (nrow(resid) > 0) resid$vcov[1] else NA_real_
  return(results)
}

run_pearson_cor <- function(df, x, y, label_context, label_group) {
  d <- df %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  if (nrow(d) < 3) return(data.frame(group = label_group, context = label_context,
                                     x = x, y = y, n = nrow(d), r = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "pearson"))
  data.frame(group = label_group, context = label_context, x = x, y = y,
             n = nrow(d), r = unname(ct$estimate), p = ct$p.value)
}

run_rmcorr <- function(df, x, y, label_context, label_group) {
  d <- df %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>%
    select(athlete_name, all_of(x), all_of(y))
  n_subj <- n_distinct(d$athlete_name)
  if (nrow(d) < 6 || n_subj < 3) {
    return(data.frame(group = label_group, context = label_context, x = x, y = y,
                      n_obs = nrow(d), n_subjects = n_subj,
                      r_rm = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p = NA_real_))
  }
  out <- tryCatch(rmcorr::rmcorr(participant = athlete_name, measure1 = d[[x]],
                                 measure2 = d[[y]], dataset = d), error = function(e) NULL)
  if (is.null(out)) {
    return(data.frame(group = label_group, context = label_context, x = x, y = y,
                      n_obs = nrow(d), n_subjects = n_subj,
                      r_rm = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p = NA_real_))
  }
  data.frame(group = label_group, context = label_context, x = x, y = y,
             n_obs = nrow(d), n_subjects = n_subj,
             r_rm = out$r, ci_lower = out$CI[1], ci_upper = out$CI[2], p = out$p)
}

extract_emmeans_interaction <- function(model, model_name, within_var) {
  if (is.null(model)) return(data.frame(model = model_name, context = NA, estimate = NA,
                                        SE = NA, lower.CL = NA, upper.CL = NA))
  emm <- tryCatch(emtrends(model, ~ context_f, var = within_var), error = function(e) NULL)
  if (is.null(emm)) return(data.frame(model = model_name, context = NA, estimate = NA,
                                      SE = NA, lower.CL = NA, upper.CL = NA))
  s <- as.data.frame(summary(emm))
  trend_col <- grep("\\.trend$", names(s), value = TRUE)
  if (length(trend_col) == 0) trend_col <- names(s)[2]
  data.frame(model = model_name, context = as.character(s$context_f),
             estimate = s[[trend_col]], SE = s$SE, lower.CL = s$lower.CL, upper.CL = s$upper.CL)
}

# =============================================================================
# PART B + C: CORE RUNNER (models + results + outputs for a subgroup)
# =============================================================================
run_all_aim1 <- function(df, prefix) {
  
  cat("\n\n=============================================================\n")
  cat("RUNNING AIM 1 FOR GROUP:", prefix, "\n")
  cat("Observations:", nrow(df), " | Athletes:", n_distinct(df$athlete_name), "\n")
  cat("Practices:", sum(df$context == "practice"), " | Games:", sum(df$context == "game"), "\n")
  cat("=============================================================\n")
  
  # ----------------------------
  # INTENSITY MODELS (sRPE ~ HR)
  # ----------------------------
  cat("\n=== RUNNING INTENSITY MODELS ===\n")
  
  model_trimp_ri <- safe_lmer(sRPE ~ TRIMP_min_within_z + TRIMP_min_between_z +
                                context_f + TRIMP_min_within_z:context_f + (1 | athlete_name), df, REML = TRUE)
  model_hravg_ri <- safe_lmer(sRPE ~ HR_avg_pct_max_within_z + HR_avg_pct_max_between_z +
                                context_f + HR_avg_pct_max_within_z:context_f + (1 | athlete_name), df, REML = TRUE)
  model_hr80_ri <- safe_lmer(sRPE ~ HR80_pct_rd_within_z + HR80_pct_rd_between_z +
                               context_f + HR80_pct_rd_within_z:context_f + (1 | athlete_name), df, REML = TRUE)
  model_hrpeak_ri <- safe_lmer(sRPE ~ HR_peak_within_z + HR_peak_between_z +
                                 context_f + HR_peak_within_z:context_f + (1 | athlete_name), df, REML = TRUE)
  
  cat("  Fitting RS sensitivity models with fallback...\n")
  rs_trimp <- fit_rs_with_fallback(
    sRPE ~ TRIMP_min_within_z + TRIMP_min_between_z + context_f + TRIMP_min_within_z:context_f + (1 + TRIMP_min_within_z | athlete_name),
    sRPE ~ TRIMP_min_within_z + TRIMP_min_between_z + context_f + TRIMP_min_within_z:context_f + (1 | athlete_name) + (0 + TRIMP_min_within_z | athlete_name),
    data = df, REML = TRUE); model_trimp_rs <- rs_trimp$model
  rs_hravg <- fit_rs_with_fallback(
    sRPE ~ HR_avg_pct_max_within_z + HR_avg_pct_max_between_z + context_f + HR_avg_pct_max_within_z:context_f + (1 + HR_avg_pct_max_within_z | athlete_name),
    sRPE ~ HR_avg_pct_max_within_z + HR_avg_pct_max_between_z + context_f + HR_avg_pct_max_within_z:context_f + (1 | athlete_name) + (0 + HR_avg_pct_max_within_z | athlete_name),
    data = df, REML = TRUE); model_hravg_rs <- rs_hravg$model
  rs_hr80 <- fit_rs_with_fallback(
    sRPE ~ HR80_pct_rd_within_z + HR80_pct_rd_between_z + context_f + HR80_pct_rd_within_z:context_f + (1 + HR80_pct_rd_within_z | athlete_name),
    sRPE ~ HR80_pct_rd_within_z + HR80_pct_rd_between_z + context_f + HR80_pct_rd_within_z:context_f + (1 | athlete_name) + (0 + HR80_pct_rd_within_z | athlete_name),
    data = df, REML = TRUE); model_hr80_rs <- rs_hr80$model
  rs_hrpeak <- fit_rs_with_fallback(
    sRPE ~ HR_peak_within_z + HR_peak_between_z + context_f + HR_peak_within_z:context_f + (1 + HR_peak_within_z | athlete_name),
    sRPE ~ HR_peak_within_z + HR_peak_between_z + context_f + HR_peak_within_z:context_f + (1 | athlete_name) + (0 + HR_peak_within_z | athlete_name),
    data = df, REML = TRUE); model_hrpeak_rs <- rs_hrpeak$model
  
  # ----------------------------
  # LOAD MODELS (sRPE-TL ~ load)
  # ----------------------------
  cat("\n=== RUNNING LOAD MODELS (exposure-corrected sRPE-TL) ===\n")
  model_load_trimp_ri <- safe_lmer(sRPE_TL_exposure ~ TRIMP_within_z + TRIMP_between_z +
                                     context_f + TRIMP_within_z:context_f + (1 | athlete_name), df, REML = TRUE)
  model_load_hr80_ri <- safe_lmer(sRPE_TL_exposure ~ HR80_mins_within_z + HR80_mins_between_z +
                                    context_f + HR80_mins_within_z:context_f + (1 | athlete_name), df, REML = TRUE)
  rs_load_trimp <- fit_rs_with_fallback(
    sRPE_TL_exposure ~ TRIMP_within_z + TRIMP_between_z + context_f + TRIMP_within_z:context_f + (1 + TRIMP_within_z | athlete_name),
    sRPE_TL_exposure ~ TRIMP_within_z + TRIMP_between_z + context_f + TRIMP_within_z:context_f + (1 | athlete_name) + (0 + TRIMP_within_z | athlete_name),
    data = df, REML = TRUE); model_load_trimp_rs <- rs_load_trimp$model
  rs_load_hr80 <- fit_rs_with_fallback(
    sRPE_TL_exposure ~ HR80_mins_within_z + HR80_mins_between_z + context_f + HR80_mins_within_z:context_f + (1 + HR80_mins_within_z | athlete_name),
    sRPE_TL_exposure ~ HR80_mins_within_z + HR80_mins_between_z + context_f + HR80_mins_within_z:context_f + (1 | athlete_name) + (0 + HR80_mins_within_z | athlete_name),
    data = df, REML = TRUE); model_load_hr80_rs <- rs_load_hr80$model
  
  # ----------------------------
  # TOI MODELS (Games only — sensitivity)
  # ----------------------------
  cat("\n=== RUNNING TOI SENSITIVITY MODELS (Games only) ===\n")
  games_toi <- df %>% filter(context == "game") %>%
    group_by(athlete_name) %>%
    mutate(TRIMP_within_game = TRIMP - mean(TRIMP, na.rm = TRUE),
           HR80_mins_within_game = HR80_mins - mean(HR80_mins, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(TRIMP_within_game_z = scale(TRIMP_within_game)[, 1],
           HR80_mins_within_game_z = scale(HR80_mins_within_game)[, 1])
  model_toi_trimp <- safe_lmer(sRPE_TL_toi ~ TRIMP_within_game_z + (1 | athlete_name), games_toi, REML = TRUE)
  model_toi_hr80 <- safe_lmer(sRPE_TL_toi ~ HR80_mins_within_game_z + (1 | athlete_name), games_toi, REML = TRUE)
  model_dur_trimp_games <- safe_lmer(sRPE_TL_duration ~ TRIMP_within_game_z + (1 | athlete_name), games_toi, REML = TRUE)
  model_dur_hr80_games <- safe_lmer(sRPE_TL_duration ~ HR80_mins_within_game_z + (1 | athlete_name), games_toi, REML = TRUE)
  
  # =============================================================================
  # EXTRACT RESULTS
  # =============================================================================
  cat("\n=== EXTRACTING RESULTS ===\n")
  results_intensity <- bind_rows(
    extract_model_results(model_trimp_ri, "TRIMP_min", "sRPE", "RI"),
    extract_model_results(model_hravg_ri, "HR_avg_pct_max", "sRPE", "RI"),
    extract_model_results(model_hr80_ri, "HR80_pct_rd", "sRPE", "RI"),
    extract_model_results(model_hrpeak_ri, "HR_peak", "sRPE", "RI"),
    extract_model_results(model_trimp_rs, "TRIMP_min", "sRPE", "RS"),
    extract_model_results(model_hravg_rs, "HR_avg_pct_max", "sRPE", "RS"),
    extract_model_results(model_hr80_rs, "HR80_pct_rd", "sRPE", "RS"),
    extract_model_results(model_hrpeak_rs, "HR_peak", "sRPE", "RS"))
  results_load <- bind_rows(
    extract_model_results(model_load_trimp_ri, "TRIMP", "sRPE_TL_exposure", "RI"),
    extract_model_results(model_load_hr80_ri, "HR80_mins", "sRPE_TL_exposure", "RI"),
    extract_model_results(model_load_trimp_rs, "TRIMP", "sRPE_TL_exposure", "RS"),
    extract_model_results(model_load_hr80_rs, "HR80_mins", "sRPE_TL_exposure", "RS"))
  results_toi <- bind_rows(
    extract_model_results(model_toi_trimp, "TRIMP_TOI", "sRPE_TL_toi", "RI"),
    extract_model_results(model_toi_hr80, "HR80_mins_TOI", "sRPE_TL_toi", "RI"),
    extract_model_results(model_dur_trimp_games, "TRIMP_Duration", "sRPE_TL_duration", "RI"),
    extract_model_results(model_dur_hr80_games, "HR80_mins_Duration", "sRPE_TL_duration", "RI"))
  
  # --- Summary tables ---
  summary_validity <- results_intensity %>% filter(model_type == "RI") %>%
    filter(grepl("_within_z$", term) & !grepl(":", term)) %>%
    mutate(HR_metric = model, beta_within = round(estimate, 3), SE = round(std_error, 3),
           CI_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
           p = ifelse(is.na(p_value), NA, ifelse(p_value < 0.001, "< .001", round(p_value, 3))),
           R2_m = round(R2_marginal, 3), AIC = round(AIC, 1)) %>%
    select(HR_metric, beta_within, SE, CI_95, p, R2_m, AIC) %>% arrange(AIC)
  
  summary_moderation <- results_intensity %>% filter(model_type == "RI") %>%
    filter(grepl(":context_fgame", term)) %>%
    mutate(HR_metric = model, beta_interaction = round(estimate, 3), SE = round(std_error, 3),
           CI_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
           p = ifelse(is.na(p_value), NA, ifelse(p_value < 0.001, "< .001", round(p_value, 3)))) %>%
    select(HR_metric, beta_interaction, SE, CI_95, p)
  
  summary_context <- results_intensity %>% filter(model_type == "RI") %>%
    filter(term == "context_fgame") %>%
    mutate(HR_metric = model, beta_game = round(estimate, 3), SE = round(std_error, 3),
           CI_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
           p = ifelse(is.na(p_value), NA, ifelse(p_value < 0.001, "< .001", round(p_value, 3)))) %>%
    select(HR_metric, beta_game, SE, CI_95, p)
  
  summary_comparison <- results_intensity %>% filter(term == "(Intercept)") %>%
    select(model, model_type, AIC, BIC, R2_marginal, R2_conditional) %>%
    group_by(model) %>%
    mutate(AIC = round(AIC, 1), BIC = round(BIC, 1), R2_marginal = round(R2_marginal, 3),
           R2_conditional = round(R2_conditional, 3), delta_AIC = round(AIC - min(AIC, na.rm = TRUE), 1)) %>%
    ungroup() %>% rename(HR_metric = model) %>% arrange(HR_metric, model_type)
  
  summary_toi <- results_toi %>% filter(grepl("_within", term)) %>%
    mutate(method = ifelse(grepl("TOI", model), "TOI", "Duration"),
           HR_metric = gsub("_TOI|_Duration", "", model),
           beta = round(estimate, 3), SE = round(std_error, 3),
           p = ifelse(is.na(p_value), NA, ifelse(p_value < 0.001, "< .001", round(p_value, 3))),
           R2 = round(R2_marginal, 3)) %>%
    select(HR_metric, method, outcome, beta, SE, p, R2) %>% arrange(HR_metric, method)
  
  # --- emmeans ---
  cat("\n=== PROBING INTERACTIONS WITH EMMEANS ===\n")
  emmeans_results <- bind_rows(
    extract_emmeans_interaction(model_trimp_ri, "TRIMP_min", "TRIMP_min_within_z"),
    extract_emmeans_interaction(model_hravg_ri, "HR_avg_pct_max", "HR_avg_pct_max_within_z"),
    extract_emmeans_interaction(model_hr80_ri, "HR80_pct_rd", "HR80_pct_rd_within_z"),
    extract_emmeans_interaction(model_hrpeak_ri, "HR_peak", "HR_peak_within_z")) %>%
    mutate(estimate = round(estimate, 3), SE = round(SE, 3),
           lower.CL = round(lower.CL, 3), upper.CL = round(upper.CL, 3),
           CI_95 = paste0("[", lower.CL, ", ", upper.CL, "]"))
  
  # --- Pattern consistency ---
  cat("\n=== BUILDING PATTERN-CONSISTENCY SUMMARY ===\n")
  pattern_consistency <- results_intensity %>% filter(model_type == "RI") %>%
    filter(grepl("_within_z$", term) & !grepl(":", term)) %>%
    mutate(HR_metric = model, direction = ifelse(estimate > 0, "positive", "negative"),
           ci_excludes_zero = ifelse(!is.na(ci_lower) & !is.na(ci_upper),
                                     ifelse(ci_lower > 0 | ci_upper < 0, "yes", "no"), NA_character_),
           p_below_05 = ifelse(!is.na(p_value), ifelse(p_value < 0.05, "yes", "no"), NA_character_)) %>%
    select(HR_metric, estimate, ci_lower, ci_upper, direction, ci_excludes_zero, p_below_05)
  all_same_dir <- length(unique(pattern_consistency$direction)) == 1
  all_ci_excl <- all(pattern_consistency$ci_excludes_zero == "yes", na.rm = TRUE)
  pattern_consistency <- bind_rows(pattern_consistency,
                                   data.frame(HR_metric = "OVERALL_PATTERN", estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                                              direction = ifelse(all_same_dir, paste0("consistent_", unique(pattern_consistency$direction)), "mixed"),
                                              ci_excludes_zero = ifelse(all_ci_excl, "all_yes", "not_all"), p_below_05 = NA_character_))
  
  moderation_pattern <- results_intensity %>% filter(model_type == "RI") %>%
    filter(grepl(":context_fgame", term)) %>%
    mutate(HR_metric = model, direction = ifelse(estimate > 0, "positive", "negative"),
           ci_excludes_zero = ifelse(!is.na(ci_lower) & !is.na(ci_upper),
                                     ifelse(ci_lower > 0 | ci_upper < 0, "yes", "no"), NA_character_),
           p_below_05 = ifelse(!is.na(p_value), ifelse(p_value < 0.05, "yes", "no"), NA_character_),
           term_type = "interaction") %>%
    select(HR_metric, estimate, direction, ci_excludes_zero, p_below_05, term_type)
  
  # --- Descriptives ---
  descriptives <- df %>% group_by(context) %>%
    summarise(n = n(), n_athletes = n_distinct(athlete_name),
              sRPE_mean = round(mean(sRPE, na.rm = TRUE), 2), sRPE_sd = round(sd(sRPE, na.rm = TRUE), 2),
              TRIMP_min_mean = round(mean(TRIMP_min, na.rm = TRUE), 2), TRIMP_min_sd = round(sd(TRIMP_min, na.rm = TRUE), 2),
              HR_avg_pct_max_mean = round(mean(HR_avg_pct_max, na.rm = TRUE), 2), HR_avg_pct_max_sd = round(sd(HR_avg_pct_max, na.rm = TRUE), 2),
              HR80_pct_rd_mean = round(mean(HR80_pct_rd, na.rm = TRUE), 2), HR80_pct_rd_sd = round(sd(HR80_pct_rd, na.rm = TRUE), 2),
              HR_peak_mean = round(mean(HR_peak, na.rm = TRUE), 2), HR_peak_sd = round(sd(HR_peak, na.rm = TRUE), 2),
              TRIMP_mean = round(mean(TRIMP, na.rm = TRUE), 2), TRIMP_sd = round(sd(TRIMP, na.rm = TRUE), 2),
              HR80_mins_mean = round(mean(HR80_mins, na.rm = TRUE), 2), HR80_mins_sd = round(sd(HR80_mins, na.rm = TRUE), 2),
              duration_mins_mean = round(mean(duration_mins, na.rm = TRUE), 2), duration_mins_sd = round(sd(duration_mins, na.rm = TRUE), 2),
              duration_exposure_mean = round(mean(duration_mins_exposure, na.rm = TRUE), 2), duration_exposure_sd = round(sd(duration_mins_exposure, na.rm = TRUE), 2),
              sRPE_TL_exposure_mean = round(mean(sRPE_TL_exposure, na.rm = TRUE), 2), sRPE_TL_exposure_sd = round(sd(sRPE_TL_exposure, na.rm = TRUE), 2),
              .groups = "drop")
  
  # --- Correlations ---
  cat("\n=== RUNNING CORRELATIONS ===\n")
  corr_pairs_intensity <- list(c("sRPE","TRIMP_min"), c("sRPE","HR_avg_pct_max"), c("sRPE","HR80_pct_rd"), c("sRPE","HR_peak"))
  corr_pairs_load <- list(c("sRPE_TL_exposure","TRIMP"), c("sRPE_TL_exposure","HR80_mins"))
  all_pairs <- c(corr_pairs_intensity, corr_pairs_load)
  correlations_pearson <- bind_rows(
    lapply(all_pairs, function(p) run_pearson_cor(df, p[1], p[2], "all", prefix)),
    lapply(all_pairs, function(p) run_pearson_cor(df %>% filter(context == "practice"), p[1], p[2], "practice", prefix)),
    lapply(all_pairs, function(p) run_pearson_cor(df %>% filter(context == "game"), p[1], p[2], "game", prefix))) %>%
    mutate(r = round(r, 3), p = ifelse(is.na(p), NA, ifelse(p < 0.001, "< .001", round(p, 3))))
  correlations_rmcorr <- bind_rows(
    lapply(all_pairs, function(p) run_rmcorr(df, p[1], p[2], "all", prefix)),
    lapply(all_pairs, function(p) run_rmcorr(df %>% filter(context == "practice"), p[1], p[2], "practice", prefix)),
    lapply(all_pairs, function(p) run_rmcorr(df %>% filter(context == "game"), p[1], p[2], "game", prefix))) %>%
    mutate(r_rm = round(r_rm, 3), ci_lower = round(ci_lower, 3), ci_upper = round(ci_upper, 3),
           p = ifelse(is.na(p), NA, ifelse(p < 0.001, "< .001", round(p, 3))))
  
  # --- RS diagnostics ---
  rs_info <- list(TRIMP_min = rs_trimp, HR_avg_pct_max = rs_hravg, HR80_pct_rd = rs_hr80,
                  HR_peak = rs_hrpeak, TRIMP_load = rs_load_trimp, HR80_mins_load = rs_load_hr80)
  random_slope_diagnostics <- bind_rows(lapply(names(rs_info), function(nm) {
    info <- rs_info[[nm]]; m <- info$model
    if (is.null(m)) return(data.frame(model = nm, fit_ok = FALSE, is_singular = NA, converged = NA,
                                      rs_type = info$type, fallback_used = info$fallback_used,
                                      warnings = paste(info$warnings, collapse = " | "), n_obs = nrow(df), n_athletes = n_distinct(df$athlete_name)))
    opt <- m@optinfo; conv <- TRUE
    if (!is.null(opt$conv$lme4$messages)) conv <- FALSE
    data.frame(model = nm, fit_ok = TRUE, is_singular = lme4::isSingular(m, tol = 1e-4), converged = conv,
               rs_type = info$type, fallback_used = info$fallback_used,
               warnings = paste(info$warnings, collapse = " | "), n_obs = nrow(df), n_athletes = n_distinct(df$athlete_name))
  }))
  
  # =============================================================================
  # SAVE MODEL OUTPUTS
  # =============================================================================
  cat("\n=== SAVING MODEL OUTPUTS FOR:", prefix, "===\n")
  write.csv(results_intensity, paste0(prefix, "results_intensity_full.csv"), row.names = FALSE)
  write.csv(results_load, paste0(prefix, "results_load_full.csv"), row.names = FALSE)
  write.csv(results_toi, paste0(prefix, "results_toi_full.csv"), row.names = FALSE)
  write.csv(summary_validity, paste0(prefix, "summary_validity.csv"), row.names = FALSE)
  write.csv(summary_moderation, paste0(prefix, "summary_moderation.csv"), row.names = FALSE)
  write.csv(summary_context, paste0(prefix, "summary_context_effect.csv"), row.names = FALSE)
  write.csv(summary_comparison, paste0(prefix, "summary_model_comparison.csv"), row.names = FALSE)
  write.csv(summary_toi, paste0(prefix, "summary_toi_comparison.csv"), row.names = FALSE)
  write.csv(descriptives, paste0(prefix, "descriptives.csv"), row.names = FALSE)
  write.csv(emmeans_results, paste0(prefix, "summary_emmeans_interactions.csv"), row.names = FALSE)
  write.csv(pattern_consistency, paste0(prefix, "summary_pattern_consistency.csv"), row.names = FALSE)
  write.csv(moderation_pattern, paste0(prefix, "summary_moderation_pattern.csv"), row.names = FALSE)
  write.csv(correlations_pearson, paste0(prefix, "correlations_pearson.csv"), row.names = FALSE)
  write.csv(correlations_rmcorr, paste0(prefix, "correlations_rmcorr.csv"), row.names = FALSE)
  write.csv(random_slope_diagnostics, paste0(prefix, "random_slope_diagnostics.csv"), row.names = FALSE)
  save(model_trimp_ri, model_hravg_ri, model_hr80_ri, model_hrpeak_ri,
       model_trimp_rs, model_hravg_rs, model_hr80_rs, model_hrpeak_rs,
       model_load_trimp_ri, model_load_hr80_ri, model_load_trimp_rs, model_load_hr80_rs,
       model_toi_trimp, model_toi_hr80, model_dur_trimp_games, model_dur_hr80_games,
       file = paste0(prefix, "aim1_models.RData"))
  cat("Saved all outputs with prefix:", prefix, "\n")
}

# =============================================================================
# RUN MODELS: ALL / FWD / DEF
# =============================================================================
df_all <- aim1_data %>% filter(position_group %in% c("Forward", "Defenseman"))
run_all_aim1(df_all, "ALL_")
df_fwd <- aim1_data %>% filter(position_group == "Forward")
run_all_aim1(df_fwd, "FWD_")
df_def <- aim1_data %>% filter(position_group == "Defenseman")
run_all_aim1(df_def, "DEF_")

cat("\n\n=== MODELS COMPLETE ===\n")

# =============================================================================
# PART D: SCATTERPLOTS
# =============================================================================
cat("\n=== GENERATING SCATTERPLOTS ===\n")

df <- df_all  # Use same data as ALL_ models

theme_scatter <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.title = element_text(size = 10), axis.text = element_text(size = 9),
        legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray92"),
        plot.margin = margin(12, 15, 12, 12))
ctx_colors <- c("Practice" = "#2166AC", "Game" = "#B2182B")

make_scatter <- function(data, x_var, y_var, x_lab, y_lab, title, subtitle, filename, width = 7, height = 6) {
  d <- data %>% filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]))
  if (nrow(d) == 0) { cat("  Skipped:", filename, "\n"); return(invisible(NULL)) }
  x_max <- max(d[[x_var]], na.rm = TRUE) * 1.05; y_max <- max(d[[y_var]], na.rm = TRUE) * 1.05
  p <- ggplot(d, aes(x = .data[[x_var]], y = .data[[y_var]], color = Context)) +
    geom_point(alpha = 0.2, size = 1, stroke = 0) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.15) +
    scale_color_manual(values = ctx_colors) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab, color = "Context") + theme_scatter
  ggsave(filename, p, width = width, height = height, dpi = 300, bg = "white")
  cat("  Saved:", filename, "\n")
}

n_ath <- n_distinct(df$athlete_id); n_obs <- nrow(df)
sub_base <- paste0("N = ", n_ath, " athletes, ", n_obs, " sessions. Lines = linear fit +/- 95% CI.")
sub_tl <- paste0(sub_base, "\nNote: sRPE-TL uses practice duration (practice) and TOI (game).")

make_scatter(df, "TRIMP_min", "sRPE", "TRIMP/min (AU/min)", "sRPE-M (0-10)", "TRIMP/min vs sRPE-M", sub_base, "scatter_01_TRIMPmin_vs_sRPE.png")
make_scatter(df, "HR80_pct_rd", "sRPE", "Time >80% HRmax (% of session)", "sRPE-M (0-10)", "Time >80% HRmax (%) vs sRPE-M", sub_base, "scatter_02_pctHR80_vs_sRPE.png")
make_scatter(df, "HR_avg_pct_max", "sRPE", "Mean HR (% HRmax)", "sRPE-M (0-10)", "Mean HR vs sRPE-M", sub_base, "scatter_03_HRmean_vs_sRPE.png")
make_scatter(df, "HR_peak", "sRPE", "Peak HR (bpm)", "sRPE-M (0-10)", "Peak HR vs sRPE-M", sub_base, "scatter_04_HRpeak_vs_sRPE.png")
make_scatter(df, "TRIMP", "sRPE_TL_exposure", "TRIMP (AU)", "sRPE-TL (AU)", "TRIMP vs sRPE-TL", sub_tl, "scatter_05_TRIMP_vs_sRPETL.png")
make_scatter(df, "HR80_mins", "sRPE_TL_exposure", "Time >80% HRmax (min)", "sRPE-TL (AU)", "Time >80% HRmax (min) vs sRPE-TL", sub_tl, "scatter_06_HR80mins_vs_sRPETL.png")
make_scatter(df, "TRIMP", "sRPE", "TRIMP (AU)", "sRPE-M (0-10)", "TRIMP vs sRPE-M", sub_base, "scatter_07_TRIMP_vs_sRPE.png")
make_scatter(df, "HR80_mins", "sRPE", "Time >80% HRmax (min)", "sRPE-M (0-10)", "Time >80% HRmax (min) vs sRPE-M", sub_base, "scatter_08_HR80mins_vs_sRPE.png")

# =============================================================================
# PART E: TABLE 1 + RMCORR TABLES
# =============================================================================
cat("\n=== GENERATING TABLE 1 ===\n")

fmt_mean_nr <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  paste0(round(m, 1), " [", round(m - s, 1), ", ", round(m + s, 1), "]")
}
vars_to_summarize <- c("duration_mins_exposure", "sRPE", "sRPE_TL_exposure",
                       "TRIMP", "TRIMP_min", "HR_avg_pct_max", "HR_peak", "HR80_pct_rd", "HR80_mins")

athlete_desc <- df %>% group_by(athlete_id, position_group, Context) %>%
  summarise(n_sessions = n(), across(all_of(vars_to_summarize), fmt_mean_nr, .names = "{.col}"), .groups = "drop")
t1_prac <- athlete_desc %>% filter(Context == "Practice") %>%
  rename_with(~ paste0(.x, "_Practice"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
t1_game <- athlete_desc %>% filter(Context == "Game") %>%
  rename_with(~ paste0(.x, "_Game"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
table1 <- t1_prac %>% left_join(t1_game, by = c("athlete_id", "position_group")) %>% arrange(position_group, athlete_id)

pos_tot <- df %>% group_by(position_group, Context) %>%
  summarise(n_sessions = n(), across(all_of(vars_to_summarize), fmt_mean_nr, .names = "{.col}"), .groups = "drop") %>%
  mutate(athlete_id = paste0(position_group, " Total"))
pt_prac <- pos_tot %>% filter(Context == "Practice") %>%
  rename_with(~ paste0(.x, "_Practice"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
pt_game <- pos_tot %>% filter(Context == "Game") %>%
  rename_with(~ paste0(.x, "_Game"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
pos_combined <- pt_prac %>% left_join(pt_game, by = c("athlete_id", "position_group"))
grand_tot <- df %>% group_by(Context) %>%
  summarise(n_sessions = n(), across(all_of(vars_to_summarize), fmt_mean_nr, .names = "{.col}"), .groups = "drop") %>%
  mutate(athlete_id = "All Skaters", position_group = "All")
gt_prac <- grand_tot %>% filter(Context == "Practice") %>%
  rename_with(~ paste0(.x, "_Practice"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
gt_game <- grand_tot %>% filter(Context == "Game") %>%
  rename_with(~ paste0(.x, "_Game"), all_of(c("n_sessions", vars_to_summarize))) %>% select(-Context)
grand_combined <- gt_prac %>% left_join(gt_game, by = c("athlete_id", "position_group"))
table1_full <- bind_rows(table1, pos_combined, grand_combined) %>% rename(Participant = athlete_id, Position = position_group)
write.csv(table1_full, "table1_participant_descriptives.csv", row.names = FALSE)
cat("  Saved: table1_participant_descriptives.csv\n")

# --- RMCORR TABLES ---
cat("\n=== GENERATING RMCORR TABLES ===\n")
rmcorr_pairs <- list(
  list(x = "TRIMP_min", y = "sRPE", label = "TRIMP/min vs sRPE"),
  list(x = "HR80_pct_rd", y = "sRPE", label = "%HR80 vs sRPE"),
  list(x = "HR_avg_pct_max", y = "sRPE", label = "HRmean vs sRPE"),
  list(x = "HR_peak", y = "sRPE", label = "HRpeak vs sRPE"),
  list(x = "TRIMP", y = "sRPE", label = "TRIMP vs sRPE"),
  list(x = "HR80_mins", y = "sRPE", label = "HR80 mins vs sRPE"),
  list(x = "TRIMP", y = "sRPE_TL_exposure", label = "TRIMP vs sRPE-TL"),
  list(x = "HR80_mins", y = "sRPE_TL_exposure", label = "HR80 mins vs sRPE-TL"))

compute_athlete_r <- function(data, x_var, y_var, pair_label) {
  data %>% filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]])) %>%
    group_by(athlete_id, position_group) %>%
    summarise(
      n_pair = n(),
      r = tryCatch({ if (n() < 5) return(NA_real_); unname(cor.test(.data[[x_var]], .data[[y_var]])$estimate) }, error = function(e) NA_real_),
      ci_lo = tryCatch({ if (n() < 5) return(NA_real_); cor.test(.data[[x_var]], .data[[y_var]])$conf.int[1] }, error = function(e) NA_real_),
      ci_hi = tryCatch({ if (n() < 5) return(NA_real_); cor.test(.data[[x_var]], .data[[y_var]])$conf.int[2] }, error = function(e) NA_real_),
      .groups = "drop") %>% mutate(pair = pair_label)
}

compute_group_rmcorr <- function(data, x_var, y_var, pair_label, group_label) {
  d <- data %>% filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]])) %>% mutate(participant = factor(athlete_id))
  n_subj <- n_distinct(d$participant)
  if (nrow(d) < 6 || n_subj < 3) return(data.frame(group = group_label, pair = pair_label, n_obs = nrow(d), n_athletes = n_subj, r = NA, ci_lo = NA, ci_hi = NA, p = NA))
  out <- tryCatch(rmcorr(participant = participant, measure1 = d[[x_var]], measure2 = d[[y_var]], dataset = d), error = function(e) NULL)
  if (is.null(out)) return(data.frame(group = group_label, pair = pair_label, n_obs = nrow(d), n_athletes = n_subj, r = NA, ci_lo = NA, ci_hi = NA, p = NA))
  data.frame(group = group_label, pair = pair_label, n_obs = nrow(d), n_athletes = n_subj,
             r = round(out$r, 3), ci_lo = round(out$CI[1], 3), ci_hi = round(out$CI[2], 3), p = out$p)
}

build_rmcorr_table <- function(data, context_label) {
  athlete_r <- bind_rows(lapply(rmcorr_pairs, function(p) compute_athlete_r(data, p$x, p$y, p$label)))
  athlete_r_fmt <- athlete_r %>%
    mutate(r_fmt = ifelse(is.na(r), "insufficient", paste0(round(r, 2), " [", round(ci_lo, 2), ", ", round(ci_hi, 2), "]")))
  athlete_n <- athlete_r %>% group_by(athlete_id, position_group) %>%
    summarise(n = suppressWarnings(min(n_pair, na.rm = TRUE)), .groups = "drop") %>%
    mutate(n = ifelse(is.infinite(n), NA_integer_, as.integer(n)))
  athlete_wide <- athlete_r_fmt %>%
    select(athlete_id, position_group, pair, r_fmt) %>%
    pivot_wider(id_cols = c(athlete_id, position_group), names_from = pair, values_from = r_fmt) %>%
    left_join(athlete_n, by = c("athlete_id", "position_group")) %>%
    relocate(n, .after = position_group) %>% arrange(position_group, athlete_id)
  get_pair_n <- function(d, x, y) d %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>% nrow()
  pos_summary <- athlete_r %>% filter(!is.na(r)) %>% group_by(position_group, pair) %>%
    summarise(r_summary = { m <- mean(r, na.rm = TRUE); s <- sd(r, na.rm = TRUE)
    paste0(round(m, 2), " [", round(m - s, 2), ", ", round(m + s, 2), "]") }, .groups = "drop") %>%
    pivot_wider(id_cols = position_group, names_from = pair, values_from = r_summary) %>%
    mutate(athlete_id = paste0(position_group, " Mean [+/-1SD]"))
  pos_n <- bind_rows(lapply(unique(data$position_group), function(pg) {
    dpg <- data %>% filter(position_group == pg)
    data.frame(position_group = pg, n = as.integer(min(sapply(rmcorr_pairs, function(p) get_pair_n(dpg, p$x, p$y))))) }))
  pos_summary <- pos_summary %>% left_join(pos_n, by = "position_group") %>% relocate(n, .after = position_group) %>%
    select(athlete_id, position_group, n, everything())
  overall_summary <- athlete_r %>% filter(!is.na(r)) %>% group_by(pair) %>%
    summarise(r_summary = { m <- mean(r, na.rm = TRUE); s <- sd(r, na.rm = TRUE)
    paste0(round(m, 2), " [", round(m - s, 2), ", ", round(m + s, 2), "]") }, .groups = "drop") %>%
    pivot_wider(names_from = pair, values_from = r_summary) %>%
    mutate(athlete_id = "All Skaters Mean [+/-1SD]", position_group = "All",
           n = as.integer(min(sapply(rmcorr_pairs, function(p) get_pair_n(data, p$x, p$y))))) %>%
    relocate(n, .after = position_group)
  rmcorr_all <- bind_rows(lapply(rmcorr_pairs, function(p) compute_group_rmcorr(data, p$x, p$y, p$label, "All Skaters rmcorr")))
  rmcorr_fwd <- bind_rows(lapply(rmcorr_pairs, function(p) compute_group_rmcorr(data %>% filter(position_group == "Forward"), p$x, p$y, p$label, "Forwards rmcorr")))
  rmcorr_def <- bind_rows(lapply(rmcorr_pairs, function(p) compute_group_rmcorr(data %>% filter(position_group == "Defenseman"), p$x, p$y, p$label, "Defensemen rmcorr")))
  rmcorr_combined <- bind_rows(rmcorr_all, rmcorr_fwd, rmcorr_def) %>%
    mutate(display = ifelse(is.na(r), "insufficient", paste0(r, " [", ci_lo, ", ", ci_hi, "]")))
  rmcorr_n <- rmcorr_combined %>% group_by(group) %>%
    summarise(n = as.integer(suppressWarnings(min(n_obs, na.rm = TRUE))), .groups = "drop")
  rmcorr_wide <- rmcorr_combined %>% select(group, pair, display) %>%
    pivot_wider(id_cols = group, names_from = pair, values_from = display) %>%
    left_join(rmcorr_n, by = "group") %>% rename(athlete_id = group) %>%
    mutate(position_group = athlete_id) %>% relocate(n, .after = position_group)
  final <- bind_rows(athlete_wide, pos_summary, overall_summary, rmcorr_wide) %>%
    rename(Participant = athlete_id, Position = position_group)
  return(final)
}

table_overall <- build_rmcorr_table(df, "Overall")
write.csv(table_overall, "table_rmcorr_overall.csv", row.names = FALSE)
cat("  Saved: table_rmcorr_overall.csv\n")

df_practice <- df %>% filter(context == "practice")
table_practice <- build_rmcorr_table(df_practice, "Practice")
write.csv(table_practice, "table_rmcorr_practice.csv", row.names = FALSE)
cat("  Saved: table_rmcorr_practice.csv\n")

df_game <- df %>% filter(context == "game")
table_game <- build_rmcorr_table(df_game, "Game")
write.csv(table_game, "table_rmcorr_game.csv", row.names = FALSE)
cat("  Saved: table_rmcorr_game.csv\n")

athlete_r_all <- bind_rows(
  bind_rows(lapply(rmcorr_pairs, function(p) compute_athlete_r(df, p$x, p$y, p$label))) %>% mutate(context = "overall"),
  bind_rows(lapply(rmcorr_pairs, function(p) compute_athlete_r(df_practice, p$x, p$y, p$label))) %>% mutate(context = "practice"),
  bind_rows(lapply(rmcorr_pairs, function(p) compute_athlete_r(df_game, p$x, p$y, p$label))) %>% mutate(context = "game"))
write.csv(athlete_r_all, "table_athlete_r_raw_all_contexts.csv", row.names = FALSE)
cat("  Saved: table_athlete_r_raw_all_contexts.csv\n")

cat("\n=== AIM 1 COMPLETE (ALL / FWD / DEF + Figures + Tables) ===\n")
