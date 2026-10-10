## ============================================================================
## Combined multilevel logistic-regression analysis
##
## Part 1 (first): Factors associated with IPV prevalence among women in the
##                 "child age [1,3]" group
## Part 2 (second): Moderation analysis comparing "child age [1,3]" with
##                  "no child" using group_bin interactions
##
## ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr)
  library(readxl); library(writexl)
  library(lme4); library(rlang)
})
## -------------------------
## 1) Settings
## -------------------------

data_path <- "..."
out_dir   <- "..."

P13_LABEL      <- "child age [1,3]"
NO_CHILD_LABEL <- "no child"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

prev_dir <- file.path(out_dir, "01_prevalence_models")
pd_dir   <- file.path(out_dir, "02_pd_models")
prev_rds_dir <- file.path(prev_dir, "models")
pd_rds_dir   <- file.path(pd_dir, "models")

for (x in c(prev_dir, pd_dir, prev_rds_dir, pd_rds_dir)) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
}

## -------------------------
## 2) Save/load helpers
## -------------------------

safe_save_rds <- function(obj, path) {
  tryCatch({
    saveRDS(obj, path, compress = "xz")
    message("Saved: ", path)
    TRUE
  }, error = function(e) {
    message("ERROR saving ", path, ": ", e$message)
    FALSE
  })
}

safe_load_rds <- function(path) {
  tryCatch({
    readRDS(path)
  }, error = function(e) {
    message("ERROR loading ", path, ": ", e$message)
    NULL
  })
}

safe_write_xlsx <- function(list_of_dfs, path) {
  tryCatch({
    writexl::write_xlsx(list_of_dfs, path = path)
    message("Excel written: ", path)
    TRUE
  }, error = function(e) {
    message("ERROR writing Excel ", path, ": ", e$message)
    FALSE
  })
}

safe_write_csv <- function(df, path) {
  tryCatch({
    write.csv(df, path, row.names = FALSE)
    message("CSV written: ", path)
    TRUE
  }, error = function(e) {
    message("ERROR writing CSV ", path, ": ", e$message)
    FALSE
  })
}

log_error <- function(msg) {
  cat(
    sprintf("[%s] %s\n", Sys.time(), msg),
    file = file.path(out_dir, "model_errors.log"),
    append = TRUE
  )
}

safe_rename <- function(df, map_old_to_new) {
  old <- intersect(names(df), names(map_old_to_new))
  if (length(old) == 0) return(df)
  new <- unname(map_old_to_new[old])
  dplyr::rename(df, !!!setNames(old, new))
}

## -------------------------
## 3) Read and prepare data
## -------------------------

message("Reading data from: ", data_path)
data <- read_excel(data_path)
names(data) <- tolower(names(data))

data$country     <- substr(data$v000, 1, 2)
data$DHS_version <- suppressWarnings(as.integer(substr(data$v000, 3, 3)))

if ("mv714" %in% names(data)) data$mv714[data$mv714 == 9] <- NA
if ("v714"  %in% names(data)) data$v714[data$v714 == 9] <- NA

rename_map <- c(
  v012  = "age",
  v730  = "partner_age",
  v511  = "age_first_marriage",
  v025  = "residence",
  v133  = "edu_years",
  v218  = "n_children",
  v212  = "age_first_birth",
  v714  = "job",
  mv714 = "partner_job",
  v702  = "partner_edu_yrs",
  v007  = "interview_year"
)

data <- safe_rename(data, rename_map)

data <- data %>%
  mutate(
    age = ifelse(!is.na(age) & age < 15, NA, age),
    age_first_birth = ifelse(
      !is.na(age_first_birth) & age_first_birth < 15,
      NA,
      age_first_birth
    ),
    
    ## Divide selected age variables by 5 so that their ORs are per 5 years.
    partner_age        = suppressWarnings(as.numeric(partner_age)) / 5,
    age_first_marriage = suppressWarnings(as.numeric(age_first_marriage)) / 5,
    age_first_birth    = suppressWarnings(as.numeric(age_first_birth)) / 5,
    age                = suppressWarnings(as.numeric(age)) / 5,
    
    first_child = case_when(
      n_children == 1 ~ 1,
      n_children > 1  ~ 0,
      TRUE            ~ NA_real_
    )
  )

data <- data %>%
  mutate(
    partner_edu_yrs = suppressWarnings(as.numeric(partner_edu_yrs)),
    edu_years       = suppressWarnings(as.numeric(edu_years)),
    age_diff        = suppressWarnings(as.numeric(partner_age)) -
      suppressWarnings(as.numeric(age)),
    edu_year_diff   = edu_years - partner_edu_yrs,
    interview_year  = suppressWarnings(as.numeric(interview_year))
  )

data <- data %>%
  mutate(
    gii_std = ifelse(!is.na(gii) & gii <= 1, gii * 100, gii),
    gdi_std = ifelse(!is.na(gdi) & gdi <= 1, gdi * 100, gdi),
    lfpr_f_std = ifelse(
      !is.na(lfpr_f) & lfpr_f <= 1,
      lfpr_f * 100,
      lfpr_f
    ),
    lfpr_m_std = ifelse(
      !is.na(lfpr_m) & lfpr_m <= 1,
      lfpr_m * 100,
      lfpr_m
    ),
    lfpr_diff_std = ifelse(
      !is.na(lfpr_diff) & abs(lfpr_diff) <= 1,
      lfpr_diff * 100,
      lfpr_diff
    ),
    paid_leave_f_weeks = ifelse(
      !is.na(paid_leave_f_days),
      paid_leave_f_days / 7,
      NA_real_
    ),
    paid_leave_m_weeks = ifelse(
      !is.na(paid_leave_m_days),
      paid_leave_m_days / 7,
      NA_real_
    )
  )

data <- data %>%
  mutate(
    group_code = as.character(group_code),
    group_bin = case_when(
      group_code == NO_CHILD_LABEL ~ 0,
      group_code == P13_LABEL ~ 1,
      TRUE ~ NA_real_
    )
  )

## -------------------------
## 4) Model variables and analysis samples
## -------------------------

outcomes <- c("any_ipv", "phy_ipv", "emo_ipv", "sex_ipv")

## Part 1: prevalence-factor models in the child age [1,3] subgroup.
prev_individual_predictors <- c(
  "residence", "job", "partner_job", "first_child",
  "age", "age_first_birth", "age_first_marriage", "interview_year",
  "edu_years", "age_diff", "edu_year_diff"
)

prev_country_predictors <- c(
  "dv_law_explicit", "gii_std", "gdi_std", "female_gni_ppp_k",
  "eys_f", "eys_diff", "lfpr_f_std", "lfpr_diff_std",
  "paid_leave_f_weeks", "paid_leave_m_weeks"
)

## Part 2: PD moderation models comparing child age [1,3] with no child.
pd_individual_predictors <- c(
  "residence", "job", "partner_job",
  "age", "age_first_marriage", "interview_year",
  "edu_years", "age_diff", "edu_year_diff"
)

pd_country_predictors <- c(
  "dv_law_explicit", "gii_std", "gdi_std", "female_gni_ppp_k",
  "eys_f", "eys_diff", "lfpr_f_std", "lfpr_diff_std"
)

categorical_variables <- c(
  "residence", "job", "partner_job", "first_child", "dv_law_explicit"
)

for (v in categorical_variables) {
  if (v %in% names(data) && !is.factor(data[[v]])) data[[v]] <- factor(data[[v]])
}

numeric_variables <- unique(c(
  prev_individual_predictors,
  prev_country_predictors,
  pd_individual_predictors,
  pd_country_predictors
))
numeric_variables <- setdiff(numeric_variables, categorical_variables)

for (v in numeric_variables) {
  if (v %in% names(data)) data[[v]] <- suppressWarnings(as.numeric(data[[v]]))
}

data_prevalence <- subset(data, group_code == P13_LABEL)
if (nrow(data_prevalence) == 0) {
  stop("The P13 prevalence-analysis sample is empty. Check group_code.")
}

data_pd <- subset(data, group_bin %in% c(0, 1))
if (nrow(data_pd) == 0) {
  stop("The PD-analysis sample is empty. Check group_code/group_bin.")
}

## -------------------------
## 5) GLMM helper functions
## -------------------------

ctrl <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

tidy_glmer_all <- function(fit, analysis, model_set, outcome) {
  if (is.null(fit)) return(data.frame())
  
  cf <- summary(fit)$coefficients
  rn <- rownames(cf)
  rn <- rn[rn != "(Intercept)"]
  
  out <- lapply(rn, function(r) {
    est <- cf[r, "Estimate"]
    se  <- cf[r, "Std. Error"]
    p   <- cf[r, "Pr(>|z|)"]
    OR  <- exp(est)
    CI  <- exp(est + c(-1, 1) * 1.96 * se)
    
    data.frame(
      Analysis = analysis,
      Model_Set = model_set,
      Outcome = outcome,
      Term = r,
      Beta = round(est, 3),
      OR = round(OR, 3),
      CI_low = round(CI[1], 3),
      CI_high = round(CI[2], 3),
      p_value = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
      N = nrow(model.frame(fit)),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

fit_and_save_big <- function(
    formula,
    model_data,
    analysis,
    model_set,
    outcome,
    name_tag,
    model_dir,
    result_dir) {
  
  rds_path  <- file.path(model_dir, paste0(name_tag, ".rds"))
  csv_path  <- file.path(result_dir, paste0(name_tag, "_results.csv"))
  xlsx_path <- file.path(result_dir, paste0(name_tag, "_results.xlsx"))
  
  if (file.exists(rds_path)) {
    message("Loading existing model: ", rds_path)
    fit <- safe_load_rds(rds_path)
    
    if (!is.null(fit)) {
      result <- tidy_glmer_all(fit, analysis, model_set, outcome)
      if (nrow(result) > 0) {
        safe_write_csv(result, csv_path)
        safe_write_xlsx(list(Result = result), xlsx_path)
      }
      return(list(model = fit, result = result))
    }
  }
  
  message("Fitting BIG model: ", name_tag)
  
  fit <- tryCatch({
    glmer(
      formula,
      data = model_data,
      family = binomial,
      control = ctrl
    )
  }, error = function(e) {
    log_error(paste("BIG fit error:", name_tag, e$message))
    message("BIG model fit error for ", name_tag, ": ", e$message)
    NULL
  })
  
  if (!is.null(fit)) {
    safe_save_rds(fit, rds_path)
    result <- tidy_glmer_all(fit, analysis, model_set, outcome)
    
    if (nrow(result) > 0) {
      safe_write_csv(result, csv_path)
      safe_write_xlsx(list(Result = result), xlsx_path)
    }
    
    return(list(model = fit, result = result))
  }
  
  list(model = NULL, result = data.frame())
}

## -------------------------
## 6) Part 1: prevalence-factor BIG models
## -------------------------

prev_formula <- function(outcome, predictors) {
  as.formula(
    paste(
      outcome,
      "~",
      paste(predictors, collapse = " + "),
      "+ (1 | country/v000)"
    )
  )
}

prev_specs <- list(
  BIG_Individual = prev_individual_predictors,
  BIG_Country = prev_country_predictors
)

prev_runs <- list()

for (model_set in names(prev_specs)) {
  predictors <- prev_specs[[model_set]]
  
  for (outcome in outcomes) {
    name_tag <- paste0("Prevalence_", model_set, "_", outcome)
    
    prev_runs[[name_tag]] <- fit_and_save_big(
      formula = prev_formula(outcome, predictors),
      model_data = data_prevalence,
      analysis = "Prevalence factors",
      model_set = model_set,
      outcome = outcome,
      name_tag = name_tag,
      model_dir = prev_rds_dir,
      result_dir = prev_dir
    )
  }
}

prevalence_results <- bind_rows(lapply(prev_runs, `[[`, "result"))

if (nrow(prevalence_results) > 0) {
  safe_write_csv(
    prevalence_results,
    file.path(prev_dir, "Prevalence_BIG_All_OR.csv")
  )
  safe_write_xlsx(
    list(Prevalence_BIG = prevalence_results),
    file.path(prev_dir, "Prevalence_BIG_All_OR.xlsx")
  )
}

## -------------------------
## 7) Part 2: PD moderation BIG models
## -------------------------

pd_formula <- function(outcome, predictors) {
  as.formula(
    paste0(
      outcome,
      " ~ group_bin * (",
      paste(predictors, collapse = " + "),
      ") + (1 | country/v000)"
    )
  )
}

pd_specs <- list(
  BIG_Individual = pd_individual_predictors,
  BIG_Country = pd_country_predictors
)

pd_runs <- list()

for (model_set in names(pd_specs)) {
  predictors <- pd_specs[[model_set]]
  
  for (outcome in outcomes) {
    name_tag <- paste0("PD_", model_set, "_", outcome)
    
    pd_runs[[name_tag]] <- fit_and_save_big(
      formula = pd_formula(outcome, predictors),
      model_data = data_pd,
      analysis = "PD moderation",
      model_set = model_set,
      outcome = outcome,
      name_tag = name_tag,
      model_dir = pd_rds_dir,
      result_dir = pd_dir
    )
  }
}

pd_results <- bind_rows(lapply(pd_runs, `[[`, "result"))

if (nrow(pd_results) > 0) {
  safe_write_csv(
    pd_results,
    file.path(pd_dir, "PD_BIG_All_OR.csv")
  )
  safe_write_xlsx(
    list(PD_BIG = pd_results),
    file.path(pd_dir, "PD_BIG_All_OR.xlsx")
  )
}

## -------------------------
## 8) Combined outputs
## -------------------------

all_results <- bind_rows(prevalence_results, pd_results)

if (nrow(all_results) > 0) {
  safe_write_csv(
    all_results,
    file.path(out_dir, "BIG_All_Models_OR.csv")
  )
  
  combined_sheets <- list()
  if (nrow(prevalence_results) > 0) {
    combined_sheets$Prevalence_OR <- prevalence_results
  }
  if (nrow(pd_results) > 0) {
    combined_sheets$PD_OR <- pd_results
  }
  
  safe_write_xlsx(
    combined_sheets,
    file.path(out_dir, "BIG_All_Models_OR.xlsx")
  )
}

safe_write_csv(
  data_prevalence,
  file.path(out_dir, "data_for_prevalence_analysis.csv")
)

safe_write_csv(
  data_pd,
  file.path(out_dir, "data_for_PD_analysis.csv")
)

## -------------------------
## 9) Session information
## -------------------------

session_info_path <- file.path(
  out_dir,
  paste0("session_info_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

tryCatch({
  sink(session_info_path)
  cat("Timestamp: ", as.character(Sys.time()), "\n\n")
  print(sessionInfo())
  sink()
  message("Saved session info: ", session_info_path)
}, error = function(e) {
  message("ERROR saving session info: ", e$message)
})

message("All done. Outputs under: ", out_dir)
