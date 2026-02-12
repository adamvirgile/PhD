# =============================================================================
# AIM 1 STATISTICAL ANALYSIS: Construct Validity of sRPE-M (Revised)
#
# Core questions:
#   Q1. Does sRPE-M demonstrate convergent validity with HR-derived
#       physiological indicators (within-athlete association)?
#   Q2. Does this relationship differ between training and competition?
#
# Analyses:
#   - Primary: Mixed-effects models with within/between decomposition
#   - Sensitivity: Random slopes (with pre-specified convergence fallback)
#   - Descriptive: rmcorr + Pearson correlations by context
#   - Secondary: Position subgroup (Forwards vs Defensemen)
#   - Supplementary: emmeans interaction probing, pattern-consistency summary
#   - Volume-inclusive: sRPE-TL ~ TRIMP / HR80_mins load models
#   - TOI comparison: TOI vs Duration for game sRPE-TL
# =============================================================================
# Location: C:\Users\avirgile\Dropbox\Adam Only\School\PhD\UVM Coursework\Dissertation\data\aim1
#
# Inputs:
#   - aim1_analysis_data.csv (in this folder)
#
# Outputs (in this folder; each prefixed with ALL_, FWD_, DEF_):
#   - results_intensity_full.csv
#   - results_load_full.csv
#   - results_toi_full.csv
#   - summary_validity.csv
#   - summary_moderation.csv
#   - summary_context_effect.csv
#   - summary_model_comparison.csv
#   - summary_toi_comparison.csv
#   - summary_pattern_consistency.csv    [NEW]
#   - summary_emmeans_interactions.csv   [NEW]
#   - descriptives.csv
#   - correlations_pearson.csv
#   - correlations_rmcorr.csv
#   - random_slope_diagnostics.csv
#
# Also saves:
#   - <PREFIX>aim1_models.RData
# =============================================================================

# Set working directory to aim1 folder
setwd("C:/Users/avirgile/Dropbox/Adam Only/School/PhD/UVM Coursework/Dissertation/data/aim1")

# =============================================================================
# 1. LOAD PACKAGES
# =============================================================================
library(lme4)
library(lmerTest)
library(performance)
library(dplyr)
library(tidyr)
library(broom.mixed)
library(emmeans)       # [CHANGE 1] Now loaded and used for interaction probing

# For repeated-measures correlations (nested validity-style correlations)
suppressWarnings({
  if (!requireNamespace("rmcorr", quietly = TRUE)) {
    message("Package 'rmcorr' not installed. Install with: install.packages('rmcorr')")
  }
})
library(rmcorr)

# =============================================================================
# 2. LOAD DATA
# =============================================================================
aim1_data <- read.csv("aim1_analysis_data.csv")

cat("=== DATA LOADED ===\n")
cat("Working directory:", getwd(), "\n")
cat("Observations:", nrow(aim1_data), "\n")
cat("Athletes:", length(unique(aim1_data$athlete_name)), "\n")
cat("Games:", sum(aim1_data$context == "game"), "\n")
cat("Practices:", sum(aim1_data$context == "practice"), "\n")

# =============================================================================
# 3. PREPARE VARIABLES
# =============================================================================

# Standardize position labels defensively (your profiles already have position)
aim1_data <- aim1_data %>%
  mutate(
    position_clean = tolower(trimws(as.character(position))),
    position_group = case_when(
      position_clean %in% c("forward", "forwards", "f", "fw", "fwd") ~ "Forward",
      position_clean %in% c("defense", "defenceman", "defenseman", "defensemen", "d") ~ "Defenseman",
      TRUE ~ NA_character_
    )
  )

# NOTE: Your aggregation already excluded goalies; we still protect against weird labels
aim1_data <- aim1_data %>%
  filter(is.na(position_group) == FALSE)

# Core factors + z-scoring already in your script; kept as-is
aim1_data <- aim1_data %>%
  mutate(
    context_f = factor(context, levels = c("practice", "game")),
    
    # Standardize within-athlete intensity metrics
    TRIMP_min_within_z = scale(TRIMP_min_within)[, 1],
    HR_avg_pct_max_within_z = scale(HR_avg_pct_max_within)[, 1],
    HR80_pct_rd_within_z = scale(HR80_pct_rd_within)[, 1],
    HR_peak_within_z = scale(HR_peak_within)[, 1],
    
    # Standardize between-athlete intensity metrics
    TRIMP_min_between_z = scale(TRIMP_min_between)[, 1],
    HR_avg_pct_max_between_z = scale(HR_avg_pct_max_between)[, 1],
    HR80_pct_rd_between_z = scale(HR80_pct_rd_between)[, 1],
    HR_peak_between_z = scale(HR_peak_between)[, 1]
  )

# Create within/between for load metrics (your original approach)
# --- Exposure-corrected sRPE-TL (PRIMARY for load analyses) ---
# Games: sRPE * TOI (actual playing time); Practices: sRPE * practice duration
# If sRPE_TL_exposure wasn't created in aggregation, compute it here
if (!"sRPE_TL_exposure" %in% names(aim1_data)) {
  cat("  Creating sRPE_TL_exposure (TOI for games, duration for practices)...\n")
  aim1_data <- aim1_data %>%
    mutate(
      duration_mins_exposure = case_when(
        context == "game" & !is.na(toi_total_mins) & toi_total_mins > 0 ~ toi_total_mins,
        context == "practice" & !is.na(duration_mins) & duration_mins > 0 ~ duration_mins,
        TRUE ~ NA_real_
      ),
      sRPE_TL_exposure = ifelse(!is.na(sRPE) & !is.na(duration_mins_exposure),
                                sRPE * duration_mins_exposure, NA_real_)
    )
} else {
  cat("  sRPE_TL_exposure already present in data.\n")
}

# Ensure duration_mins_exposure exists
if (!"duration_mins_exposure" %in% names(aim1_data)) {
  aim1_data <- aim1_data %>%
    mutate(
      duration_mins_exposure = case_when(
        context == "game" & !is.na(toi_total_mins) & toi_total_mins > 0 ~ toi_total_mins,
        context == "practice" & !is.na(duration_mins) & duration_mins > 0 ~ duration_mins,
        TRUE ~ NA_real_
      )
    )
}

cat("  sRPE_TL_exposure: non-missing =", sum(!is.na(aim1_data$sRPE_TL_exposure)),
    "of", nrow(aim1_data), "\n")
cat("  duration_mins_exposure: non-missing =", sum(!is.na(aim1_data$duration_mins_exposure)),
    "of", nrow(aim1_data), "\n")

aim1_data <- aim1_data %>%
  group_by(athlete_name) %>%
  mutate(
    TRIMP_within = TRIMP - mean(TRIMP, na.rm = TRUE),
    TRIMP_between = mean(TRIMP, na.rm = TRUE),
    HR80_mins_within_new = HR80_mins - mean(HR80_mins, na.rm = TRUE),
    HR80_mins_between_new = mean(HR80_mins, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    TRIMP_within_z = scale(TRIMP_within)[, 1],
    TRIMP_between_z = scale(TRIMP_between)[, 1],
    HR80_mins_within_z = scale(HR80_mins_within_new)[, 1],
    HR80_mins_between_z = scale(HR80_mins_between_new)[, 1]
  )

# =============================================================================
# 4. HELPERS
# =============================================================================

# [CHANGE 2] safe_lmer now captures and LOGS warnings instead of silently
# swallowing them. This is critical for the pre-specified convergence decision
# rule: we need to know *what* went wrong (convergence failure, singular fit)
# so we can apply the correct fallback (uncorrelated slopes -> RI-only).
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

# [CHANGE 3] New helper to implement the Section 2.4.1 pre-specified decision
# rule for random slopes:
#   1. Fit RI model (primary)
#   2. Attempt RS model (correlated slopes)
#   3. If singular/convergence failure -> try uncorrelated slopes
#   4. If still failing -> report RI and flag limitation
fit_rs_with_fallback <- function(formula_rs, formula_rs_uncorr, data, REML = TRUE) {
  
  # Step 2: Try correlated random slopes
  m <- safe_lmer(formula_rs, data = data, REML = REML)
  warns <- attr(m, "fit_warnings")
  
  if (!is.null(m) && !isSingular(m, tol = 1e-4) && length(warns) == 0) {
    return(list(model = m, type = "RS_correlated", warnings = warns,
                fallback_used = "none"))
  }
  
  # Step 3: Try uncorrelated random slopes
  cat("    -> RS correlated had issues, trying uncorrelated...\n")
  m_uc <- safe_lmer(formula_rs_uncorr, data = data, REML = REML)
  warns_uc <- attr(m_uc, "fit_warnings")
  
  if (!is.null(m_uc) && !isSingular(m_uc, tol = 1e-4) && length(warns_uc) == 0) {
    return(list(model = m_uc, type = "RS_uncorrelated", warnings = warns_uc,
                fallback_used = "uncorrelated"))
  }
  
  # Step 4: Return best available (prefer correlated if it at least converged)
  if (!is.null(m)) {
    return(list(model = m, type = "RS_correlated_with_issues",
                warnings = warns, fallback_used = "kept_with_issues"))
  }
  if (!is.null(m_uc)) {
    return(list(model = m_uc, type = "RS_uncorrelated_with_issues",
                warnings = warns_uc, fallback_used = "kept_uc_with_issues"))
  }
  
  # Total failure
  return(list(model = NULL, type = "RS_failed",
              warnings = c(warns, warns_uc), fallback_used = "failed"))
}

extract_model_results <- function(model, model_name, outcome_var, model_type) {
  if (is.null(model)) {
    return(data.frame(
      model = model_name,
      outcome = outcome_var,
      model_type = model_type,
      term = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      R2_marginal = NA_real_,
      R2_conditional = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      random_intercept_var = NA_real_,
      random_slope_var = NA_real_,
      random_intercept_slope_cor = NA_real_,
      residual_var = NA_real_
    ))
  }
  
  fe <- summary(model)$coefficients
  
  # Wald CI for fixed effects only
  ci <- tryCatch(
    confint(model, method = "Wald", parm = "beta_"),
    error = function(e) NULL
  )
  
  r2_vals <- tryCatch(performance::r2(model), error = function(e) NULL)
  vc <- as.data.frame(VarCorr(model))
  
  # Pull random components
  ri <- vc %>% filter(grp == "athlete_name" & var1 == "(Intercept)" & is.na(var2))
  rs <- vc %>% filter(grp == "athlete_name" & var1 != "(Intercept)" & is.na(var2))
  cor_is <- vc %>% filter(grp == "athlete_name" & !is.na(var2))
  resid <- vc %>% filter(grp == "Residual")
  
  # Build results row-per-fixed-term
  results <- data.frame(
    model = model_name,
    outcome = outcome_var,
    model_type = model_type,
    term = rownames(fe),
    estimate = fe[, "Estimate"],
    std_error = fe[, "Std. Error"],
    t_value = fe[, "t value"],
    p_value = fe[, "Pr(>|t|)"],
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    row.names = NULL
  )
  
  if (!is.null(ci)) {
    rn <- results$term
    ok <- rn %in% rownames(ci)
    results$ci_lower[ok] <- ci[rn[ok], 1]
    results$ci_upper[ok] <- ci[rn[ok], 2]
  }
  
  results$R2_marginal <- if (!is.null(r2_vals)) r2_vals$R2_marginal else NA_real_
  results$R2_conditional <- if (!is.null(r2_vals)) r2_vals$R2_conditional else NA_real_
  results$AIC <- AIC(model)
  results$BIC <- BIC(model)
  
  results$random_intercept_var <- if (nrow(ri) > 0) ri$vcov[1] else NA_real_
  results$random_slope_var <- if (nrow(rs) > 0) rs$vcov[1] else NA_real_
  results$random_intercept_slope_cor <- if (nrow(cor_is) > 0) cor_is$vcov[1] else NA_real_
  results$residual_var <- if (nrow(resid) > 0) resid$vcov[1] else NA_real_
  
  return(results)
}

# Pearson correlations (optionally within a filter)
run_pearson_cor <- function(df, x, y, label_context, label_group) {
  d <- df %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  if (nrow(d) < 3) {
    return(data.frame(
      group = label_group, context = label_context, x = x, y = y,
      n = nrow(d), r = NA_real_, p = NA_real_
    ))
  }
  ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "pearson"))
  data.frame(
    group = label_group, context = label_context, x = x, y = y,
    n = nrow(d),
    r = unname(ct$estimate),
    p = ct$p.value
  )
}

# Repeated-measures correlations (rmcorr) across athletes
run_rmcorr <- function(df, x, y, label_context, label_group) {
  d <- df %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>%
    select(athlete_name, all_of(x), all_of(y))
  
  n_subj <- n_distinct(d$athlete_name)
  if (nrow(d) < 6 || n_subj < 3) {
    return(data.frame(
      group = label_group, context = label_context, x = x, y = y,
      n_obs = nrow(d), n_subjects = n_subj,
      r_rm = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p = NA_real_
    ))
  }
  
  out <- tryCatch(
    rmcorr::rmcorr(participant = athlete_name, measure1 = d[[x]], measure2 = d[[y]], dataset = d),
    error = function(e) NULL
  )
  
  if (is.null(out)) {
    return(data.frame(
      group = label_group, context = label_context, x = x, y = y,
      n_obs = nrow(d), n_subjects = n_subj,
      r_rm = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p = NA_real_
    ))
  }
  
  # [CHANGE 4] Now extracting rmcorr CI for richer reporting
  data.frame(
    group = label_group, context = label_context, x = x, y = y,
    n_obs = nrow(d), n_subjects = n_subj,
    r_rm = out$r,
    ci_lower = out$CI[1],
    ci_upper = out$CI[2],
    p = out$p
  )
}

# [CHANGE 5] New helper: extract emmeans simple slopes for interaction probing.
# For each RI model, this estimates the within-athlete HR slope at each context
# level (practice, game), giving a direct answer to "how strong is the
# sRPE-HR coupling in practice vs games?"
extract_emmeans_interaction <- function(model, model_name, within_var) {
  if (is.null(model)) {
    return(data.frame(
      model = model_name, context = NA, estimate = NA, SE = NA,
      lower.CL = NA, upper.CL = NA, t.ratio = NA, p.value = NA
    ))
  }
  
  emm <- tryCatch({
    emtrends(model, ~ context_f, var = within_var)
  }, error = function(e) NULL)
  
  if (is.null(emm)) {
    return(data.frame(
      model = model_name, context = NA, estimate = NA, SE = NA,
      lower.CL = NA, upper.CL = NA, t.ratio = NA, p.value = NA
    ))
  }
  
  s <- as.data.frame(summary(emm))
  # emtrends returns columns: context_f, <within_var>.trend, SE, df, lower.CL, upper.CL
  trend_col <- grep("\\.trend$", names(s), value = TRUE)
  if (length(trend_col) == 0) trend_col <- names(s)[2]  # fallback
  
  data.frame(
    model = model_name,
    context = as.character(s$context_f),
    estimate = s[[trend_col]],
    SE = s$SE,
    lower.CL = s$lower.CL,
    upper.CL = s$upper.CL
  )
}


# =============================================================================
# 5. CORE RUNNER (runs everything for a subgroup)
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
  
  # ---- Random intercept (primary) ----
  model_trimp_ri <- safe_lmer(
    sRPE ~ TRIMP_min_within_z + TRIMP_min_between_z +
      context_f + TRIMP_min_within_z:context_f +
      (1 | athlete_name),
    df, REML = TRUE
  )
  
  model_hravg_ri <- safe_lmer(
    sRPE ~ HR_avg_pct_max_within_z + HR_avg_pct_max_between_z +
      context_f + HR_avg_pct_max_within_z:context_f +
      (1 | athlete_name),
    df, REML = TRUE
  )
  
  model_hr80_ri <- safe_lmer(
    sRPE ~ HR80_pct_rd_within_z + HR80_pct_rd_between_z +
      context_f + HR80_pct_rd_within_z:context_f +
      (1 | athlete_name),
    df, REML = TRUE
  )
  
  model_hrpeak_ri <- safe_lmer(
    sRPE ~ HR_peak_within_z + HR_peak_between_z +
      context_f + HR_peak_within_z:context_f +
      (1 | athlete_name),
    df, REML = TRUE
  )
  
  # ---- Random slope (sensitivity, with Section 2.4.1 fallback) ----
  # [CHANGE 6] Using fit_rs_with_fallback to implement the pre-specified
  # decision rule: correlated RS -> uncorrelated RS -> report RI + flag.
  
  cat("  Fitting RS sensitivity models with fallback...\n")
  
  rs_trimp <- fit_rs_with_fallback(
    formula_rs       = sRPE ~ TRIMP_min_within_z + TRIMP_min_between_z +
      context_f + TRIMP_min_within_z:context_f +
      (1 + TRIMP_min_within_z | athlete_name),
    formula_rs_uncorr = sRPE ~ TRIMP_min_within_z + TRIMP_min_between_z +
      context_f + TRIMP_min_within_z:context_f +
      (1 | athlete_name) + (0 + TRIMP_min_within_z | athlete_name),
    data = df, REML = TRUE
  )
  model_trimp_rs <- rs_trimp$model
  
  rs_hravg <- fit_rs_with_fallback(
    formula_rs       = sRPE ~ HR_avg_pct_max_within_z + HR_avg_pct_max_between_z +
      context_f + HR_avg_pct_max_within_z:context_f +
      (1 + HR_avg_pct_max_within_z | athlete_name),
    formula_rs_uncorr = sRPE ~ HR_avg_pct_max_within_z + HR_avg_pct_max_between_z +
      context_f + HR_avg_pct_max_within_z:context_f +
      (1 | athlete_name) + (0 + HR_avg_pct_max_within_z | athlete_name),
    data = df, REML = TRUE
  )
  model_hravg_rs <- rs_hravg$model
  
  rs_hr80 <- fit_rs_with_fallback(
    formula_rs       = sRPE ~ HR80_pct_rd_within_z + HR80_pct_rd_between_z +
      context_f + HR80_pct_rd_within_z:context_f +
      (1 + HR80_pct_rd_within_z | athlete_name),
    formula_rs_uncorr = sRPE ~ HR80_pct_rd_within_z + HR80_pct_rd_between_z +
      context_f + HR80_pct_rd_within_z:context_f +
      (1 | athlete_name) + (0 + HR80_pct_rd_within_z | athlete_name),
    data = df, REML = TRUE
  )
  model_hr80_rs <- rs_hr80$model
  
  rs_hrpeak <- fit_rs_with_fallback(
    formula_rs       = sRPE ~ HR_peak_within_z + HR_peak_between_z +
      context_f + HR_peak_within_z:context_f +
      (1 + HR_peak_within_z | athlete_name),
    formula_rs_uncorr = sRPE ~ HR_peak_within_z + HR_peak_between_z +
      context_f + HR_peak_within_z:context_f +
      (1 | athlete_name) + (0 + HR_peak_within_z | athlete_name),
    data = df, REML = TRUE
  )
  model_hrpeak_rs <- rs_hrpeak$model
  
  # ----------------------------
  # LOAD MODELS (sRPE-TL ~ load)
  # ----------------------------
  # NOTE: sRPE_TL_exposure uses TOI for games, practice duration for practices.
  # This is the PRIMARY load metric.
  cat("\n=== RUNNING LOAD MODELS (exposure-corrected sRPE-TL) ===\n")
  
  model_load_trimp_ri <- safe_lmer(
    sRPE_TL_exposure ~ TRIMP_within_z + TRIMP_between_z +
      context_f + TRIMP_within_z:context_f +
      (1 | athlete_name),
    df, REML = TRUE
  )
  
  model_load_hr80_ri <- safe_lmer(
    sRPE_TL_exposure ~ HR80_mins_within_z + HR80_mins_between_z +
      context_f + HR80_mins_within_z:context_f +
      (1 | athlete_name),
    df, REML = TRUE
  )
  
  # Load RS with fallback
  rs_load_trimp <- fit_rs_with_fallback(
    formula_rs       = sRPE_TL_exposure ~ TRIMP_within_z + TRIMP_between_z +
      context_f + TRIMP_within_z:context_f +
      (1 + TRIMP_within_z | athlete_name),
    formula_rs_uncorr = sRPE_TL_exposure ~ TRIMP_within_z + TRIMP_between_z +
      context_f + TRIMP_within_z:context_f +
      (1 | athlete_name) + (0 + TRIMP_within_z | athlete_name),
    data = df, REML = TRUE
  )
  model_load_trimp_rs <- rs_load_trimp$model
  
  rs_load_hr80 <- fit_rs_with_fallback(
    formula_rs       = sRPE_TL_exposure ~ HR80_mins_within_z + HR80_mins_between_z +
      context_f + HR80_mins_within_z:context_f +
      (1 + HR80_mins_within_z | athlete_name),
    formula_rs_uncorr = sRPE_TL_exposure ~ HR80_mins_within_z + HR80_mins_between_z +
      context_f + HR80_mins_within_z:context_f +
      (1 | athlete_name) + (0 + HR80_mins_within_z | athlete_name),
    data = df, REML = TRUE
  )
  model_load_hr80_rs <- rs_load_hr80$model
  
  # ----------------------------
  # TOI MODELS (Games only)
  # ----------------------------
  cat("\n=== RUNNING TOI MODELS (Games only) ===\n")
  
  games_toi <- df %>%
    filter(context == "game" & !is.na(sRPE_TL_toi)) %>%
    group_by(athlete_name) %>%
    mutate(
      TRIMP_within_game = TRIMP - mean(TRIMP, na.rm = TRUE),
      HR80_mins_within_game = HR80_mins - mean(HR80_mins, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      TRIMP_within_game_z = scale(TRIMP_within_game)[, 1],
      HR80_mins_within_game_z = scale(HR80_mins_within_game)[, 1]
    )
  
  model_toi_trimp <- safe_lmer(sRPE_TL_toi ~ TRIMP_within_game_z + (1 | athlete_name), games_toi, REML = TRUE)
  model_toi_hr80 <- safe_lmer(sRPE_TL_toi ~ HR80_mins_within_game_z + (1 | athlete_name), games_toi, REML = TRUE)
  
  model_dur_trimp_games <- safe_lmer(sRPE_TL_duration ~ TRIMP_within_game_z + (1 | athlete_name), games_toi, REML = TRUE)
  model_dur_hr80_games <- safe_lmer(sRPE_TL_duration ~ HR80_mins_within_game_z + (1 | athlete_name), games_toi, REML = TRUE)
  
  # =============================================================================
  # 6. EXTRACT ALL RESULTS
  # =============================================================================
  cat("\n=== EXTRACTING RESULTS ===\n")
  
  results_intensity <- bind_rows(
    extract_model_results(model_trimp_ri,  "TRIMP_min",       "sRPE", "RI"),
    extract_model_results(model_hravg_ri,  "HR_avg_pct_max",  "sRPE", "RI"),
    extract_model_results(model_hr80_ri,   "HR80_pct_rd",     "sRPE", "RI"),
    extract_model_results(model_hrpeak_ri, "HR_peak",         "sRPE", "RI"),
    
    extract_model_results(model_trimp_rs,  "TRIMP_min",       "sRPE", "RS"),
    extract_model_results(model_hravg_rs,  "HR_avg_pct_max",  "sRPE", "RS"),
    extract_model_results(model_hr80_rs,   "HR80_pct_rd",     "sRPE", "RS"),
    extract_model_results(model_hrpeak_rs, "HR_peak",         "sRPE", "RS")
  )
  
  results_load <- bind_rows(
    extract_model_results(model_load_trimp_ri, "TRIMP",     "sRPE_TL_exposure", "RI"),
    extract_model_results(model_load_hr80_ri,  "HR80_mins", "sRPE_TL_exposure", "RI"),
    extract_model_results(model_load_trimp_rs, "TRIMP",     "sRPE_TL_exposure", "RS"),
    extract_model_results(model_load_hr80_rs,  "HR80_mins", "sRPE_TL_exposure", "RS")
  )
  
  results_toi <- bind_rows(
    extract_model_results(model_toi_trimp,       "TRIMP_TOI",        "sRPE_TL_toi",      "RI"),
    extract_model_results(model_toi_hr80,        "HR80_mins_TOI",     "sRPE_TL_toi",      "RI"),
    extract_model_results(model_dur_trimp_games, "TRIMP_Duration",    "sRPE_TL_duration", "RI"),
    extract_model_results(model_dur_hr80_games,  "HR80_mins_Duration","sRPE_TL_duration", "RI")
  )
  
  # =============================================================================
  # 7. SUMMARY TABLES
  # =============================================================================
  cat("\n=== CREATING SUMMARY TABLES ===\n")
  
  # Within-athlete validity (RI primary)
  summary_validity <- results_intensity %>%
    filter(model_type == "RI") %>%
    filter(grepl("_within_z$", term) & !grepl(":", term)) %>%
    select(model, estimate, std_error, ci_lower, ci_upper, p_value, R2_marginal, AIC) %>%
    mutate(
      HR_metric = model,
      beta_within = round(estimate, 3),
      SE = round(std_error, 3),
      CI_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
      p = ifelse(is.na(p_value), NA, ifelse(p_value < 0.001, "< .001", round(p_value, 3))),
      R2_m = round(R2_marginal, 3),
      AIC = round(AIC, 1)
    ) %>%
    select(HR_metric, beta_within, SE, CI_95, p, R2_m, AIC) %>%
    arrange(AIC)
  
  # Context moderation (RI primary)
  summary_moderation <- results_intensity %>%
    filter(model_type == "RI") %>%
    filter(grepl(":context_fgame", term)) %>%
    select(model, estimate, std_error, ci_lower, ci_upper, p_value) %>%
    mutate(
      HR_metric = model,
      beta_interaction = round(estimate, 3),
      SE = round(std_error, 3),
      CI_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
      p = ifelse(is.na(p_value), NA, ifelse(p_value < 0.001, "< .001", round(p_value, 3)))
    ) %>%
    select(HR_metric, beta_interaction, SE, CI_95, p)
  
  # Context main effect (RI primary)
  summary_context <- results_intensity %>%
    filter(model_type == "RI") %>%
    filter(term == "context_fgame") %>%
    select(model, estimate, std_error, ci_lower, ci_upper, p_value) %>%
    mutate(
      HR_metric = model,
      beta_game = round(estimate, 3),
      SE = round(std_error, 3),
      CI_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
      p = ifelse(is.na(p_value), NA, ifelse(p_value < 0.001, "< .001", round(p_value, 3)))
    ) %>%
    select(HR_metric, beta_game, SE, CI_95, p)
  
  # Model comparison for intensity (RI vs RS)
  summary_comparison <- results_intensity %>%
    filter(term == "(Intercept)") %>%
    select(model, model_type, AIC, BIC, R2_marginal, R2_conditional) %>%
    group_by(model) %>%
    mutate(
      AIC = round(AIC, 1),
      BIC = round(BIC, 1),
      R2_marginal = round(R2_marginal, 3),
      R2_conditional = round(R2_conditional, 3),
      delta_AIC = round(AIC - min(AIC, na.rm = TRUE), 1)
    ) %>%
    ungroup() %>%
    rename(HR_metric = model) %>%
    arrange(HR_metric, model_type)
  
  # TOI vs Duration summary
  summary_toi <- results_toi %>%
    filter(grepl("_within", term)) %>%
    select(model, outcome, estimate, std_error, p_value, R2_marginal) %>%
    mutate(
      method = ifelse(grepl("TOI", model), "TOI", "Duration"),
      HR_metric = gsub("_TOI|_Duration", "", model),
      beta = round(estimate, 3),
      SE = round(std_error, 3),
      p = ifelse(is.na(p_value), NA, ifelse(p_value < 0.001, "< .001", round(p_value, 3))),
      R2 = round(R2_marginal, 3)
    ) %>%
    select(HR_metric, method, outcome, beta, SE, p, R2) %>%
    arrange(HR_metric, method)
  
  # =============================================================================
  # 7b. [CHANGE 7] EMMEANS INTERACTION PROBING
  # =============================================================================
  # The dissertation says to use emmeans to probe the within x context
  # interaction. This gives the simple slope of each HR metric at each
  # context level -- directly answering "how strong is the sRPE-physiology
  # coupling in practice vs games?"
  cat("\n=== PROBING INTERACTIONS WITH EMMEANS ===\n")
  
  emmeans_results <- bind_rows(
    extract_emmeans_interaction(model_trimp_ri,  "TRIMP_min",       "TRIMP_min_within_z"),
    extract_emmeans_interaction(model_hravg_ri,  "HR_avg_pct_max",  "HR_avg_pct_max_within_z"),
    extract_emmeans_interaction(model_hr80_ri,   "HR80_pct_rd",     "HR80_pct_rd_within_z"),
    extract_emmeans_interaction(model_hrpeak_ri, "HR_peak",         "HR_peak_within_z")
  )
  
  emmeans_results <- emmeans_results %>%
    mutate(
      estimate = round(estimate, 3),
      SE = round(SE, 3),
      lower.CL = round(lower.CL, 3),
      upper.CL = round(upper.CL, 3),
      CI_95 = paste0("[", lower.CL, ", ", upper.CL, "]")
    )
  
  # =============================================================================
  # 7c. [CHANGE 8] PATTERN-CONSISTENCY SUMMARY
  # =============================================================================
  # The dissertation emphasizes: "conclusions about construct validity will be
  # based on pattern consistency across indicators and precision of within-
  # athlete estimates, rather than isolated significance testing of any single
  # metric." This table synthesizes across all 4 HR metrics to support that
  # interpretive framework.
  cat("\n=== BUILDING PATTERN-CONSISTENCY SUMMARY ===\n")
  
  pattern_consistency <- results_intensity %>%
    filter(model_type == "RI") %>%
    filter(grepl("_within_z$", term) & !grepl(":", term)) %>%
    select(model, estimate, ci_lower, ci_upper, p_value) %>%
    mutate(
      HR_metric = model,
      direction = ifelse(estimate > 0, "positive", "negative"),
      ci_excludes_zero = ifelse(!is.na(ci_lower) & !is.na(ci_upper),
                                ifelse(ci_lower > 0 | ci_upper < 0, "yes", "no"),
                                NA_character_),
      p_below_05 = ifelse(!is.na(p_value), ifelse(p_value < 0.05, "yes", "no"), NA_character_)
    ) %>%
    select(HR_metric, estimate, ci_lower, ci_upper, direction, ci_excludes_zero, p_below_05)
  
  # Overall consistency flag
  all_same_direction <- length(unique(pattern_consistency$direction)) == 1
  all_ci_exclude_zero <- all(pattern_consistency$ci_excludes_zero == "yes", na.rm = TRUE)
  
  pattern_summary_row <- data.frame(
    HR_metric = "OVERALL_PATTERN",
    estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    direction = ifelse(all_same_direction,
                       paste0("consistent_", unique(pattern_consistency$direction)),
                       "mixed"),
    ci_excludes_zero = ifelse(all_ci_exclude_zero, "all_yes", "not_all"),
    p_below_05 = NA_character_
  )
  
  pattern_consistency <- bind_rows(pattern_consistency, pattern_summary_row)
  
  # Also build pattern summary for the interaction terms (context moderation)
  moderation_pattern <- results_intensity %>%
    filter(model_type == "RI") %>%
    filter(grepl(":context_fgame", term)) %>%
    select(model, estimate, ci_lower, ci_upper, p_value) %>%
    mutate(
      HR_metric = model,
      direction = ifelse(estimate > 0, "positive", "negative"),
      ci_excludes_zero = ifelse(!is.na(ci_lower) & !is.na(ci_upper),
                                ifelse(ci_lower > 0 | ci_upper < 0, "yes", "no"),
                                NA_character_),
      p_below_05 = ifelse(!is.na(p_value), ifelse(p_value < 0.05, "yes", "no"), NA_character_),
      term_type = "interaction"
    ) %>%
    select(HR_metric, estimate, direction, ci_excludes_zero, p_below_05, term_type)
  
  # =============================================================================
  # 8. DESCRIPTIVES
  # =============================================================================
  descriptives <- df %>%
    group_by(context) %>%
    summarise(
      n = n(),
      n_athletes = n_distinct(athlete_name),
      sRPE_mean = round(mean(sRPE, na.rm = TRUE), 2),
      sRPE_sd = round(sd(sRPE, na.rm = TRUE), 2),
      TRIMP_min_mean = round(mean(TRIMP_min, na.rm = TRUE), 2),
      TRIMP_min_sd = round(sd(TRIMP_min, na.rm = TRUE), 2),
      HR_avg_pct_max_mean = round(mean(HR_avg_pct_max, na.rm = TRUE), 2),
      HR_avg_pct_max_sd = round(sd(HR_avg_pct_max, na.rm = TRUE), 2),
      HR80_pct_rd_mean = round(mean(HR80_pct_rd, na.rm = TRUE), 2),
      HR80_pct_rd_sd = round(sd(HR80_pct_rd, na.rm = TRUE), 2),
      HR_peak_mean = round(mean(HR_peak, na.rm = TRUE), 2),
      HR_peak_sd = round(sd(HR_peak, na.rm = TRUE), 2),
      TRIMP_mean = round(mean(TRIMP, na.rm = TRUE), 2),
      TRIMP_sd = round(sd(TRIMP, na.rm = TRUE), 2),
      HR80_mins_mean = round(mean(HR80_mins, na.rm = TRUE), 2),
      HR80_mins_sd = round(sd(HR80_mins, na.rm = TRUE), 2),
      duration_mins_mean = round(mean(duration_mins, na.rm = TRUE), 2),
      duration_mins_sd = round(sd(duration_mins, na.rm = TRUE), 2),
      duration_exposure_mean = round(mean(duration_mins_exposure, na.rm = TRUE), 2),
      duration_exposure_sd = round(sd(duration_mins_exposure, na.rm = TRUE), 2),
      sRPE_TL_exposure_mean = round(mean(sRPE_TL_exposure, na.rm = TRUE), 2),
      sRPE_TL_exposure_sd = round(sd(sRPE_TL_exposure, na.rm = TRUE), 2),
      .groups = "drop"
    )
  
  # =============================================================================
  # 9. CORRELATIONS (Pearson + rmcorr)
  # =============================================================================
  cat("\n=== RUNNING CORRELATIONS ===\n")
  
  corr_pairs_intensity <- list(
    c("sRPE", "TRIMP_min"),
    c("sRPE", "HR_avg_pct_max"),
    c("sRPE", "HR80_pct_rd"),
    c("sRPE", "HR_peak")
  )
  
  corr_pairs_load <- list(
    c("sRPE_TL_exposure", "TRIMP"),
    c("sRPE_TL_exposure", "HR80_mins")
  )
  
  # Pearson: overall + by context
  pearson_all <- bind_rows(
    lapply(c(corr_pairs_intensity, corr_pairs_load), function(p) run_pearson_cor(df, p[1], p[2], "all", prefix))
  )
  
  pearson_practice <- bind_rows(
    lapply(c(corr_pairs_intensity, corr_pairs_load), function(p) run_pearson_cor(df %>% filter(context == "practice"), p[1], p[2], "practice", prefix))
  )
  
  pearson_game <- bind_rows(
    lapply(c(corr_pairs_intensity, corr_pairs_load), function(p) run_pearson_cor(df %>% filter(context == "game"), p[1], p[2], "game", prefix))
  )
  
  correlations_pearson <- bind_rows(pearson_all, pearson_practice, pearson_game) %>%
    mutate(
      r = round(r, 3),
      p = ifelse(is.na(p), NA, ifelse(p < 0.001, "< .001", round(p, 3)))
    )
  
  # rmcorr: overall + by context (nested)
  rmc_all <- bind_rows(
    lapply(c(corr_pairs_intensity, corr_pairs_load), function(p) run_rmcorr(df, p[1], p[2], "all", prefix))
  )
  
  rmc_practice <- bind_rows(
    lapply(c(corr_pairs_intensity, corr_pairs_load), function(p) run_rmcorr(df %>% filter(context == "practice"), p[1], p[2], "practice", prefix))
  )
  
  rmc_game <- bind_rows(
    lapply(c(corr_pairs_intensity, corr_pairs_load), function(p) run_rmcorr(df %>% filter(context == "game"), p[1], p[2], "game", prefix))
  )
  
  correlations_rmcorr <- bind_rows(rmc_all, rmc_practice, rmc_game) %>%
    mutate(
      r_rm = round(r_rm, 3),
      ci_lower = round(ci_lower, 3),
      ci_upper = round(ci_upper, 3),
      p = ifelse(is.na(p), NA, ifelse(p < 0.001, "< .001", round(p, 3)))
    )
  
  # =============================================================================
  # 10. RANDOM SLOPE DIAGNOSTICS (enhanced)
  # =============================================================================
  cat("\n=== RANDOM SLOPE DIAGNOSTICS ===\n")
  
  # [CHANGE 9] Diagnostics now include the fallback type and captured warnings,
  # so you have a transparent audit trail for which decision-rule branch was taken.
  rs_info <- list(
    TRIMP_min       = rs_trimp,
    HR_avg_pct_max  = rs_hravg,
    HR80_pct_rd     = rs_hr80,
    HR_peak         = rs_hrpeak,
    TRIMP_load      = rs_load_trimp,
    HR80_mins_load  = rs_load_hr80
  )
  
  random_slope_diagnostics <- bind_rows(lapply(names(rs_info), function(nm) {
    info <- rs_info[[nm]]
    m <- info$model
    if (is.null(m)) {
      return(data.frame(
        model = nm, fit_ok = FALSE, is_singular = NA,
        converged = NA, rs_type = info$type,
        fallback_used = info$fallback_used,
        warnings = paste(info$warnings, collapse = " | "),
        n_obs = nrow(df), n_athletes = n_distinct(df$athlete_name)
      ))
    }
    opt <- m@optinfo
    conv <- TRUE
    if (!is.null(opt$conv$lme4$messages)) conv <- FALSE
    data.frame(
      model = nm,
      fit_ok = TRUE,
      is_singular = lme4::isSingular(m, tol = 1e-4),
      converged = conv,
      rs_type = info$type,
      fallback_used = info$fallback_used,
      warnings = paste(info$warnings, collapse = " | "),
      n_obs = nrow(df),
      n_athletes = n_distinct(df$athlete_name)
    )
  }))
  
  # =============================================================================
  # 11. SAVE OUTPUTS
  # =============================================================================
  cat("\n=== SAVING OUTPUTS FOR:", prefix, "===\n")
  
  write.csv(results_intensity, paste0(prefix, "results_intensity_full.csv"), row.names = FALSE)
  write.csv(results_load,      paste0(prefix, "results_load_full.csv"), row.names = FALSE)
  write.csv(results_toi,       paste0(prefix, "results_toi_full.csv"), row.names = FALSE)
  
  write.csv(summary_validity,    paste0(prefix, "summary_validity.csv"), row.names = FALSE)
  write.csv(summary_moderation,  paste0(prefix, "summary_moderation.csv"), row.names = FALSE)
  write.csv(summary_context,     paste0(prefix, "summary_context_effect.csv"), row.names = FALSE)
  write.csv(summary_comparison,  paste0(prefix, "summary_model_comparison.csv"), row.names = FALSE)
  write.csv(summary_toi,         paste0(prefix, "summary_toi_comparison.csv"), row.names = FALSE)
  write.csv(descriptives,        paste0(prefix, "descriptives.csv"), row.names = FALSE)
  
  # New outputs
  write.csv(emmeans_results,       paste0(prefix, "summary_emmeans_interactions.csv"), row.names = FALSE)
  write.csv(pattern_consistency,   paste0(prefix, "summary_pattern_consistency.csv"), row.names = FALSE)
  write.csv(moderation_pattern,    paste0(prefix, "summary_moderation_pattern.csv"), row.names = FALSE)
  
  write.csv(correlations_pearson, paste0(prefix, "correlations_pearson.csv"), row.names = FALSE)
  write.csv(correlations_rmcorr,  paste0(prefix, "correlations_rmcorr.csv"), row.names = FALSE)
  
  write.csv(random_slope_diagnostics, paste0(prefix, "random_slope_diagnostics.csv"), row.names = FALSE)
  
  # Save models
  save(
    model_trimp_ri, model_hravg_ri, model_hr80_ri, model_hrpeak_ri,
    model_trimp_rs, model_hravg_rs, model_hr80_rs, model_hrpeak_rs,
    model_load_trimp_ri, model_load_hr80_ri,
    model_load_trimp_rs, model_load_hr80_rs,
    model_toi_trimp, model_toi_hr80, model_dur_trimp_games, model_dur_hr80_games,
    file = paste0(prefix, "aim1_models.RData")
  )
  
  cat("Saved outputs with prefix:", prefix, "\n")
}

# =============================================================================
# 6. RUN: ALL SKATERS TOGETHER + POSITION SUBGROUPS
# =============================================================================

# Together (Forwards + Defensemen)
df_all <- aim1_data %>% filter(position_group %in% c("Forward", "Defenseman"))
run_all_aim1(df_all, "ALL_")

# Forwards only
df_fwd <- aim1_data %>% filter(position_group == "Forward")
run_all_aim1(df_fwd, "FWD_")

# Defensemen only
df_def <- aim1_data %>% filter(position_group == "Defenseman")
run_all_aim1(df_def, "DEF_")

cat("\n\n=== AIM 1 ANALYSIS COMPLETE (ALL / FWD / DEF) ===\n")