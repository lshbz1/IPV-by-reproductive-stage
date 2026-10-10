# ============================================================
# Mutually exclusive reproductive-stage IPV analysis by parity
# ============================================================


# ============================================================
# 0) Packages
# ============================================================

packages <- c(
  "dplyr", "tidyr", "stringr", "haven",
  "lme4", "emmeans", "writexl", "purrr", "tibble"
)

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) install.packages(pkg)
}

library(dplyr)
library(tidyr)
library(stringr)
library(haven)
library(lme4)
library(emmeans)
library(writexl)
library(purrr)
library(tibble)


# ============================================================
# 1) Paths and settings
# ============================================================

data_directory <- "/Users/apple/Desktop/STATA data"

out_file <- "/Users/apple/Desktop/IPV_mutually_exclusive_stage_results_FIXED.xlsx"

save_excel <- TRUE
save_intermediate_csv <- FALSE

combined_csv_out <- "/Users/apple/Desktop/combined_data_DHS_mutually_exclusive_FIXED.csv"
clean_csv_out    <- "/Users/apple/Desktop/cleaned_data_DHS_mutually_exclusive_FIXED.csv"

# Surveys to remove
remove_codes <- c("BO4", "BO5", "ZM4", "MZ6")

# Model combinations
run_tags <- c(
  "noAge_noWt",
  "noAge_wt",
  "ageAdj_noWt"
)

# Main model for publication table
primary_model_tag <- "noAge_noWt"

# Reference stage
reference_stage <- "No child, not pregnant"

# IPV outcomes
ipvs <- c(
  "any_IPV",
  "physex_IPV",
  "phy_IPV",
  "emo_IPV",
  "sex_IPV"
)

ipv_labels <- c(
  any_IPV    = "Any IPV",
  physex_IPV = "Physical or sexual IPV",
  phy_IPV    = "Physical IPV",
  emo_IPV    = "Emotional IPV",
  sex_IPV    = "Sexual IPV"
)


# ============================================================
# 2) Mutually exclusive stage levels
# ============================================================

stage_levels <- c(
  "No child, not pregnant",
  "First pregnancy",
  "First child, 0-3 years",
  "First child, 3-6 years",
  "First child, >6 years",
  "Second pregnancy",
  "Second child, 0-3 years",
  "Second child, 3-6 years",
  "Second child, >6 years",
  "Third pregnancy",
  "Third child, 0-3 years",
  "Third child, 3-6 years",
  "Third child, >6 years",
  "Fourth-or-higher pregnancy",
  "Fourth-or-higher child, 0-3 years",
  "Fourth-or-higher child, 3-6 years",
  "Fourth-or-higher child, >6 years"
)


# ============================================================
# 3) Helper functions
# ============================================================

to_num <- function(x) {
  x <- haven::zap_labels(x)
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

to_chr <- function(x) {
  x <- haven::zap_labels(x)
  as.character(x)
}

norm_group <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[\u2012\u2013\u2014\u2212]", "-")
  x <- str_replace_all(x, "[()]", "")
  str_squish(x)
}

prepare_weights <- function(df) {
  if (!"d005" %in% names(df)) return(df)
  
  df %>%
    mutate(
      w_raw = to_num(d005) / 1e6,
      country = substr(v000, 1, 2)
    ) %>%
    group_by(v000) %>%
    mutate(
      w_mean = mean(w_raw[w_raw > 0], na.rm = TRUE),
      w_norm = ifelse(is.finite(w_mean) & w_mean > 0, w_raw / w_mean, 1)
    ) %>%
    ungroup() %>%
    select(-w_mean)
}

pick_age_var <- function(df, age_var = "v012") {
  if (age_var %in% names(df)) age_var else NULL
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

get_ci_cols <- function(df) {
  lower_candidates <- c("asymp.LCL", "lower.CL", "LCL")
  upper_candidates <- c("asymp.UCL", "upper.CL", "UCL")
  
  lower_col <- lower_candidates[lower_candidates %in% names(df)][1]
  upper_col <- upper_candidates[upper_candidates %in% names(df)][1]
  
  if (is.na(lower_col) || is.na(upper_col)) {
    stop(
      "Cannot find CI columns. Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  list(lower = lower_col, upper = upper_col)
}

get_response_col <- function(df) {
  response_candidates <- c("prob", "response", "emmean")
  response_col <- response_candidates[response_candidates %in% names(df)][1]
  
  if (is.na(response_col)) {
    stop(
      "Cannot find response column. Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  response_col
}

make_contrast_list <- function(grid_stages, ref_idx) {
  target_idx <- seq_along(grid_stages)[-ref_idx]
  target_stages <- grid_stages[target_idx]
  
  contrast_list <- lapply(target_idx, function(i) {
    v <- rep(0, length(grid_stages))
    v[i] <- 1
    v[ref_idx] <- -1
    v
  })
  
  names(contrast_list) <- target_stages
  
  list(
    contrast_list = contrast_list,
    target_stages = target_stages
  )
}


# ============================================================
# 4) Read and combine raw DHS .DTA files
# ============================================================

core_vars <- c(
  "v044", "v222", "v213", "v201", "v007", "v000", "d005",
  "d105a", "d105b", "d105c", "d105d", "d105e",
  "d105f", "d105j", "d103a", "d103b", "d103c",
  "d105h", "d105i", "d105k", "d130a", "d130b", "d130c"
)

table_vars <- c(
  "d118a", "d118l", "v214", "v012", "v730",
  paste0("b8_", sprintf("%02d", 1:20)),
  paste0("b4_", sprintf("%02d", 1:20)),
  "v511", "v501", "v025", "v632",
  "v744a", "v744b", "v744c", "v744d", "v744e",
  "v106", "v133", "v155", "v218", "v212",
  "v217", "v702", "v701", "v714", "mv714"
)

all_vars <- unique(c(core_vars, table_vars))

data_files <- list.files(
  path = data_directory,
  pattern = "\\.DTA$|\\.dta$",
  full.names = TRUE
)

if (length(data_files) == 0) {
  stop("No .DTA files found in: ", data_directory)
}

selected_data_list <- list()

for (file in data_files) {
  message("Processing: ", basename(file))
  
  data_i <- read_dta(file, encoding = "latin1")
  names(data_i) <- tolower(names(data_i))
  
  selected <- data_i[, intersect(names(data_i), all_vars), drop = FALSE]
  
  missing_vars <- setdiff(all_vars, names(selected))
  for (var in missing_vars) {
    selected[[var]] <- NA
  }
  
  selected <- selected[, all_vars]
  selected$source_file <- basename(file)
  
  selected_data_list[[length(selected_data_list) + 1]] <- selected
}

combined_data <- bind_rows(selected_data_list) %>%
  mutate(v000 = to_chr(v000))

cat("Combined data rows:", nrow(combined_data), "\n")

if (save_intermediate_csv) {
  write.csv(combined_data, combined_csv_out, row.names = FALSE)
}


# ============================================================
# 5) Keep Domestic Violence Module sample
# ============================================================

data_clean <- combined_data %>%
  mutate(
    v000 = to_chr(v000),
    across(
      any_of(c(
        "v044", "v222", "v213", "v201", "v007", "d005",
        "v214", "v012", "v218"
      )),
      to_num
    )
  ) %>%
  filter(v044 == 1)

cat("Domestic Violence Module sample rows:", nrow(data_clean), "\n")

# Fix selected survey years
fix_map <- tibble::tibble(
  v000 = c("AF7", "NP6", "NP8"),
  v007_fix = c(2015, 2011, 2022)
)

data_clean <- data_clean %>%
  left_join(fix_map, by = "v000") %>%
  mutate(v007 = dplyr::coalesce(v007_fix, v007)) %>%
  select(-v007_fix)


# ============================================================
# 6) Clean IPV raw variables
# ============================================================

freq_vars <- c(
  "d105a", "d105b", "d105c", "d105d", "d105e", "d105f", "d105j",
  "d103a", "d103b", "d103c",
  "d105h", "d105i", "d105k"
)

d130_vars <- c("d130a", "d130b", "d130c")

data_clean <- data_clean %>%
  mutate(
    across(
      all_of(c(freq_vars, d130_vars)),
      to_num
    )
  ) %>%
  mutate(
    across(
      all_of(freq_vars),
      ~ ifelse(.x %in% 0:4, .x, 999)
    )
  ) %>%
  mutate(
    across(
      all_of(d130_vars),
      ~ ifelse(is.na(.x), 999, .x)
    )
  )


# ============================================================
# 7) Create IPV outcomes
# ============================================================

data_clean <- data_clean %>%
  mutate(
    phy_IPV = case_when(
      d105a %in% 1:2 | d105b %in% 1:2 | d105c %in% 1:2 |
        d105d %in% 1:2 | d105e %in% 1:2 | d105f %in% 1:2 |
        d105j %in% 1:2 | d130a == 1 ~ 1,
      
      d105a == 999 & d105b == 999 & d105c == 999 &
        d105d == 999 & d105e == 999 & d105f == 999 &
        d105j == 999 & d130a == 999 ~ NA_real_,
      
      TRUE ~ 0
    ),
    
    emo_IPV = case_when(
      d103a %in% 1:2 | d103b %in% 1:2 | d103c %in% 1:2 |
        d130c == 1 ~ 1,
      
      d103a == 999 & d103b == 999 & d103c == 999 &
        d130c == 999 ~ NA_real_,
      
      TRUE ~ 0
    ),
    
    sex_IPV = case_when(
      d105h %in% 1:2 | d105i %in% 1:2 | d105k %in% 1:2 |
        d130b == 1 ~ 1,
      
      d105h == 999 & d105i == 999 & d105k == 999 &
        d130b == 999 ~ NA_real_,
      
      TRUE ~ 0
    ),
    
    physex_IPV = case_when(
      phy_IPV == 1 | sex_IPV == 1 ~ 1,
      is.na(phy_IPV) & is.na(sex_IPV) ~ NA_real_,
      TRUE ~ 0
    ),
    
    any_IPV = case_when(
      phy_IPV == 1 | emo_IPV == 1 | sex_IPV == 1 ~ 1,
      is.na(phy_IPV) & is.na(emo_IPV) & is.na(sex_IPV) ~ NA_real_,
      TRUE ~ 0
    )
  )

cat("\nIPV outcome distributions:\n")
for (x in ipvs) {
  cat("\n", x, ":\n", sep = "")
  print(table(data_clean[[x]], useNA = "ifany"))
}


# ============================================================
# 8) Remove problematic surveys
# ============================================================

summary_anyipv_na_before <- data_clean %>%
  group_by(v000) %>%
  summarise(
    n_total = n(),
    n_na = sum(is.na(any_IPV)),
    pct_na = round(100 * n_na / n_total, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_na))

data_clean <- data_clean %>%
  filter(!v000 %in% remove_codes)

summary_anyipv_na_after <- data_clean %>%
  group_by(v000) %>%
  summarise(
    n_total = n(),
    n_na = sum(is.na(any_IPV)),
    pct_na = round(100 * n_na / n_total, 2),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_na))

cat("\nRows after removing problematic surveys:", nrow(data_clean), "\n")


# ============================================================
# 9) Create mutually exclusive reproductive-stage groups
# ============================================================
# Logic:
#
# v201 = number of children ever born
# v213 = currently pregnant
# v214 = duration of current pregnancy, months
# v222 = months since last birth
#
# Pregnancy clean window:
#   currently pregnant and v214 >= 8
#   OR 0-1 months after childbirth
#
# Post-childbirth clean windows:
#   0-3 years: 12-36 months after most recent birth
#   3-6 years: 48-72 months after most recent birth
#   >6 years:  >=84 months after most recent birth
#
# Ambiguous windows excluded:
#   2-11 months after birth
#   37-47 months after birth
#   73-83 months after birth
# ============================================================

data_clean <- data_clean %>%
  mutate(
    country = substr(v000, 1, 2),
    DHS_version = suppressWarnings(as.integer(substr(v000, 3, 3))),
    
    parity_cat = case_when(
      v201 == 0 ~ "0 children ever born",
      v201 == 1 ~ "1 child ever born",
      v201 == 2 ~ "2 children ever born",
      v201 == 3 ~ "3 children ever born",
      v201 >= 4 ~ "4+ children ever born",
      TRUE ~ NA_character_
    ),
    
    # Clean-window flags
    is_no_child_clean = v201 == 0 & v213 == 0,
    
    is_pregnancy_clean = (
      (v213 == 1 & !is.na(v214) & v214 >= 8) |
        (v213 == 0 & !is.na(v222) & v222 >= 0 & v222 <= 1)
    ),
    
    is_pp_0_3_clean = (
      v213 == 0 & !is.na(v222) &
        v222 >= 12 & v222 <= 36
    ),
    
    is_pp_3_6_clean = (
      v213 == 0 & !is.na(v222) &
        v222 >= 48 & v222 <= 72
    ),
    
    is_pp_gt6_clean = (
      v213 == 0 & !is.na(v222) &
        v222 >= 84
    ),
    
    # Mutually exclusive reproductive stages
    reproductive_stage = case_when(
      # No child
      v201 == 0 & v213 == 0 ~
        "No child, not pregnant",
      
      # First pregnancy:
      # currently pregnant with first child OR 0-1 months after first birth
      (
        v201 == 0 & v213 == 1 &
          !is.na(v214) & v214 >= 8
      ) |
        (
          v201 == 1 & v213 == 0 &
            !is.na(v222) & v222 >= 0 & v222 <= 1
        ) ~
        "First pregnancy",
      
      # First child after birth
      v201 == 1 & v213 == 0 &
        !is.na(v222) & v222 >= 12 & v222 <= 36 ~
        "First child, 0-3 years",
      
      v201 == 1 & v213 == 0 &
        !is.na(v222) & v222 >= 48 & v222 <= 72 ~
        "First child, 3-6 years",
      
      v201 == 1 & v213 == 0 &
        !is.na(v222) & v222 >= 84 ~
        "First child, >6 years",
      
      # Second pregnancy:
      # currently pregnant with second child OR 0-1 months after second birth
      (
        v201 == 1 & v213 == 1 &
          !is.na(v214) & v214 >= 8
      ) |
        (
          v201 == 2 & v213 == 0 &
            !is.na(v222) & v222 >= 0 & v222 <= 1
        ) ~
        "Second pregnancy",
      
      # Second child after birth
      v201 == 2 & v213 == 0 &
        !is.na(v222) & v222 >= 12 & v222 <= 36 ~
        "Second child, 0-3 years",
      
      v201 == 2 & v213 == 0 &
        !is.na(v222) & v222 >= 48 & v222 <= 72 ~
        "Second child, 3-6 years",
      
      v201 == 2 & v213 == 0 &
        !is.na(v222) & v222 >= 84 ~
        "Second child, >6 years",
      
      # Third pregnancy:
      # currently pregnant with third child OR 0-1 months after third birth
      (
        v201 == 2 & v213 == 1 &
          !is.na(v214) & v214 >= 8
      ) |
        (
          v201 == 3 & v213 == 0 &
            !is.na(v222) & v222 >= 0 & v222 <= 1
        ) ~
        "Third pregnancy",
      
      # Third child after birth
      v201 == 3 & v213 == 0 &
        !is.na(v222) & v222 >= 12 & v222 <= 36 ~
        "Third child, 0-3 years",
      
      v201 == 3 & v213 == 0 &
        !is.na(v222) & v222 >= 48 & v222 <= 72 ~
        "Third child, 3-6 years",
      
      v201 == 3 & v213 == 0 &
        !is.na(v222) & v222 >= 84 ~
        "Third child, >6 years",
      
      # Fourth-or-higher pregnancy:
      # currently pregnant with fourth-or-higher child
      # OR 0-1 months after fourth-or-higher birth
      (
        v201 >= 3 & v213 == 1 &
          !is.na(v214) & v214 >= 8
      ) |
        (
          v201 >= 4 & v213 == 0 &
            !is.na(v222) & v222 >= 0 & v222 <= 1
        ) ~
        "Fourth-or-higher pregnancy",
      
      # Fourth-or-higher child after birth
      v201 >= 4 & v213 == 0 &
        !is.na(v222) & v222 >= 12 & v222 <= 36 ~
        "Fourth-or-higher child, 0-3 years",
      
      v201 >= 4 & v213 == 0 &
        !is.na(v222) & v222 >= 48 & v222 <= 72 ~
        "Fourth-or-higher child, 3-6 years",
      
      v201 >= 4 & v213 == 0 &
        !is.na(v222) & v222 >= 84 ~
        "Fourth-or-higher child, >6 years",
      
      TRUE ~ NA_character_
    ),
    
    reproductive_stage = factor(
      reproductive_stage,
      levels = stage_levels
    ),
    
    mapping_exclusion_reason = case_when(
      !is.na(reproductive_stage) ~
        "included in mutually exclusive clean mapping",
      
      v213 == 1 & (is.na(v214) | v214 < 8) ~
        "excluded: early/mid pregnancy; recall window mostly before pregnancy",
      
      v213 == 0 & !is.na(v222) & v222 >= 2 & v222 < 12 ~
        "excluded: early postpartum; recall window overlaps pregnancy/pre-childbirth",
      
      v213 == 0 & !is.na(v222) & v222 > 36 & v222 < 48 ~
        "excluded: recall window overlaps 0-3y and 3-6y postpartum",
      
      v213 == 0 & !is.na(v222) & v222 > 72 & v222 < 84 ~
        "excluded: recall window overlaps 3-6y and >6y postpartum",
      
      TRUE ~
        "excluded: unclassified or missing reproductive timing/parity"
    )
  )

cat("\nMutually exclusive reproductive-stage counts:\n")
print(table(data_clean$reproductive_stage, useNA = "ifany"))

cat("\nStage by parity:\n")
print(table(data_clean$reproductive_stage, data_clean$parity_cat, useNA = "ifany"))

cat("\nMapping exclusion reasons:\n")
print(table(data_clean$mapping_exclusion_reason, useNA = "ifany"))

if (save_intermediate_csv) {
  write.csv(data_clean, clean_csv_out, row.names = FALSE)
}


# ============================================================
# 10) Descriptive summary
# ============================================================

summarise_by_stage <- function(df,
                               stage_var = "reproductive_stage",
                               stage_levels = stage_levels,
                               ipvs = ipvs,
                               ipv_labels = ipv_labels,
                               include_total = TRUE) {
  
  df_included <- df %>%
    filter(!is.na(.data[[stage_var]]))
  
  long <- df_included %>%
    mutate(stage_code = as.character(.data[[stage_var]])) %>%
    select(stage_code, all_of(ipvs)) %>%
    pivot_longer(
      cols = all_of(ipvs),
      names_to = "IPV",
      values_to = "val"
    ) %>%
    group_by(stage_code, IPV) %>%
    summarise(
      sample_size = sum(!is.na(val)),
      cases = sum(val == 1, na.rm = TRUE),
      prob = cases / sample_size,
      .groups = "drop"
    ) %>%
    mutate(
      prob = ifelse(is.nan(prob), NA_real_, prob),
      IPV_label = ipv_labels[IPV]
    )
  
  if (include_total) {
    total_long <- df_included %>%
      select(all_of(ipvs)) %>%
      pivot_longer(
        cols = all_of(ipvs),
        names_to = "IPV",
        values_to = "val"
      ) %>%
      group_by(IPV) %>%
      summarise(
        sample_size = sum(!is.na(val)),
        cases = sum(val == 1, na.rm = TRUE),
        prob = cases / sample_size,
        .groups = "drop"
      ) %>%
      mutate(
        prob = ifelse(is.nan(prob), NA_real_, prob),
        stage_code = "Total",
        IPV_label = ipv_labels[IPV]
      ) %>%
      select(stage_code, IPV, IPV_label, sample_size, cases, prob)
    
    long <- bind_rows(long, total_long)
  }
  
  long <- long %>%
    mutate(
      stage_code = factor(stage_code, levels = c(stage_levels, "Total")),
      IPV = factor(IPV, levels = ipvs),
      IPV_label = factor(IPV_label, levels = ipv_labels[ipvs]),
      prevalence_fmt = sprintf("%.2f%%", prob * 100)
    ) %>%
    arrange(IPV, stage_code)
  
  wide <- long %>%
    select(stage_code, IPV_label, sample_size, cases, prob, prevalence_fmt) %>%
    pivot_wider(
      id_cols = stage_code,
      names_from = IPV_label,
      values_from = c(sample_size, cases, prob, prevalence_fmt),
      names_sep = "__"
    ) %>%
    arrange(factor(stage_code, levels = c(stage_levels, "Total"), ordered = TRUE))
  
  list(long = long, wide = wide)
}

desc <- summarise_by_stage(
  data_clean,
  stage_var = "reproductive_stage",
  stage_levels = stage_levels,
  ipvs = ipvs,
  ipv_labels = ipv_labels,
  include_total = TRUE
)

desc_long <- desc$long
desc_wide <- desc$wide

cat("\nDescriptive summary:\n")
print(desc_wide)


# ============================================================
# 11) GLMM model functions
# ============================================================

make_stage_for_model <- function(df,
                                 stage_var = "reproductive_stage",
                                 stage_levels = stage_levels) {
  df %>%
    mutate(
      stage_code = factor(
        as.character(.data[[stage_var]]),
        levels = stage_levels
      ),
      country = substr(v000, 1, 2),
      DHS_version = suppressWarnings(as.integer(substr(v000, 3, 3)))
    )
}

fit_ipv_once <- function(df,
                         ipv,
                         control_age = FALSE,
                         weighted = FALSE,
                         ref_stage = reference_stage,
                         stage_var = "reproductive_stage",
                         stage_levels = stage_levels,
                         age_var = "v012") {
  
  dat <- make_stage_for_model(
    df = df,
    stage_var = stage_var,
    stage_levels = stage_levels
  ) %>%
    prepare_weights() %>%
    filter(!is.na(stage_code), !is.na(.data[[ipv]])) %>%
    droplevels()
  
  if (nrow(dat) == 0) {
    stop("No data available after filtering for stage and outcome: ", ipv)
  }
  
  if (!ref_stage %in% levels(dat$stage_code)) {
    stop("Reference stage is not available in model data: ", ref_stage)
  }
  
  if (n_distinct(dat$stage_code) < 2) {
    stop("Fewer than two stages available for model: ", ipv)
  }
  
  if (n_distinct(dat[[ipv]], na.rm = TRUE) < 2) {
    stop("Outcome has no variation for model: ", ipv)
  }
  
  use_age <- if (control_age) pick_age_var(dat, age_var) else NULL
  
  rhs <- "stage_code"
  if (!is.null(use_age)) rhs <- paste(rhs, "+", use_age)
  
  fmla <- as.formula(paste(ipv, "~", rhs, "+ (1 | country/v000)"))
  
  ctrl <- glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
  
  fit <- if (weighted && "w_norm" %in% names(dat)) {
    glmer(
      fmla,
      data = dat,
      family = binomial,
      weights = w_norm,
      control = ctrl
    )
  } else {
    glmer(
      fmla,
      data = dat,
      family = binomial,
      control = ctrl
    )
  }
  
  # ------------------------------------------------------------
  # Adjusted prevalence
  # ------------------------------------------------------------
  
  emm_resp <- emmeans(fit, ~ stage_code, type = "response")
  
  prev_raw <- as.data.frame(summary(emm_resp))
  prev_ci <- get_ci_cols(prev_raw)
  prev_response_col <- get_response_col(prev_raw)
  
  prev_df <- prev_raw %>%
    transmute(
      IPV = ipv,
      IPV_Type = ipv_labels[ipv],
      stage_code = as.character(stage_code),
      mean = .data[[prev_response_col]],
      lower = .data[[prev_ci$lower]],
      upper = .data[[prev_ci$upper]]
    )
  
  grid_stages <- as.character(emm_resp@grid$stage_code)
  ref_idx <- which(grid_stages == ref_stage)
  
  if (length(ref_idx) != 1) {
    stop(
      "Reference stage not found uniquely in emmeans grid: ",
      ref_stage,
      ". Grid stages: ",
      paste(grid_stages, collapse = " | ")
    )
  }
  
  contrast_info <- make_contrast_list(
    grid_stages = grid_stages,
    ref_idx = ref_idx
  )
  
  contrast_list <- contrast_info$contrast_list
  target_stages <- contrast_info$target_stages
  
  # ------------------------------------------------------------
  # Prevalence difference on response scale
  # ------------------------------------------------------------
  
  emm_resp2 <- regrid(emm_resp, transform = "response")
  
  pd_contr <- contrast(
    emm_resp2,
    method = contrast_list
  )
  
  pd_raw <- as.data.frame(summary(pd_contr, infer = c(TRUE, TRUE)))
  pd_ci <- get_ci_cols(pd_raw)
  
  if (nrow(pd_raw) != length(target_stages)) {
    stop(
      "PD contrast rows do not match target stages for ",
      ipv,
      ". nrow(pd_raw) = ", nrow(pd_raw),
      "; target stages = ", length(target_stages)
    )
  }
  
  pd_df <- pd_raw %>%
    mutate(
      stage = target_stages,
      ref_stage = ref_stage
    ) %>%
    transmute(
      IPV = ipv,
      IPV_Type = ipv_labels[ipv],
      stage = as.character(stage),
      ref_stage = as.character(ref_stage),
      pd = estimate,
      pd_lower = .data[[pd_ci$lower]],
      pd_upper = .data[[pd_ci$upper]],
      p_value = p.value
    )
  
  # ------------------------------------------------------------
  # Odds ratio on logit scale
  # ------------------------------------------------------------
  
  emm_link <- emmeans(fit, ~ stage_code)
  
  or_contr <- contrast(
    emm_link,
    method = contrast_list
  )
  
  or_raw <- as.data.frame(summary(or_contr, infer = c(TRUE, TRUE)))
  or_ci <- get_ci_cols(or_raw)
  
  if (nrow(or_raw) != length(target_stages)) {
    stop(
      "OR contrast rows do not match target stages for ",
      ipv,
      ". nrow(or_raw) = ", nrow(or_raw),
      "; target stages = ", length(target_stages)
    )
  }
  
  or_df <- or_raw %>%
    mutate(
      stage = target_stages,
      ref_stage = ref_stage
    ) %>%
    transmute(
      IPV = ipv,
      IPV_Type = ipv_labels[ipv],
      stage = as.character(stage),
      ref_stage = as.character(ref_stage),
      OR = exp(estimate),
      OR_lower = exp(.data[[or_ci$lower]]),
      OR_upper = exp(.data[[or_ci$upper]]),
      OR_p = p.value
    )
  
  pd_or <- pd_df %>%
    left_join(
      or_df,
      by = c("IPV", "IPV_Type", "stage", "ref_stage")
    )
  
  list(
    prevalence = prev_df,
    pd_or = pd_or,
    model = fit
  )
}


run_one_combo <- function(df,
                          ipvs,
                          control_age,
                          weighted,
                          ref_stage = reference_stage,
                          stage_var = "reproductive_stage",
                          stage_levels = stage_levels,
                          age_var = "v012",
                          print_n = 20) {
  
  tag <- paste0(
    ifelse(control_age, "ageAdj", "noAge"),
    "_",
    ifelse(weighted, "wt", "noWt")
  )
  
  cat("\n==============================\n")
  cat("Running: ", tag, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  
  model_failures <- list()
  
  res_list <- lapply(
    ipvs,
    function(x) {
      tryCatch(
        fit_ipv_once(
          df = df,
          ipv = x,
          control_age = control_age,
          weighted = weighted,
          ref_stage = ref_stage,
          stage_var = stage_var,
          stage_levels = stage_levels,
          age_var = age_var
        ),
        error = function(e) {
          msg <- conditionMessage(e)
          message("Model failed: ", tag, " | ", x, " | ", msg)
          
          model_failures[[length(model_failures) + 1]] <<- tibble(
            model_tag = tag,
            IPV = x,
            IPV_Type = ipv_labels[x],
            error_message = msg
          )
          
          list(
            prevalence = tibble(
              IPV = character(),
              IPV_Type = character(),
              stage_code = character(),
              mean = numeric(),
              lower = numeric(),
              upper = numeric()
            ),
            pd_or = tibble(
              IPV = character(),
              IPV_Type = character(),
              stage = character(),
              ref_stage = character(),
              pd = numeric(),
              pd_lower = numeric(),
              pd_upper = numeric(),
              p_value = numeric(),
              OR = numeric(),
              OR_lower = numeric(),
              OR_upper = numeric(),
              OR_p = numeric()
            ),
            model = NULL
          )
        }
      )
    }
  )
  
  names(res_list) <- ipvs
  
  prevalence_table <- bind_rows(lapply(res_list, `[[`, "prevalence"))
  pd_or_table <- bind_rows(lapply(res_list, `[[`, "pd_or"))
  model_failures_table <- bind_rows(model_failures)
  
  if (nrow(prevalence_table) == 0) {
    final_table <- tibble()
  } else {
    
    prevalence_table2 <- prevalence_table %>%
      mutate(
        IPV = as.character(IPV),
        IPV_Type = as.character(IPV_Type),
        stage_code_n = norm_group(stage_code),
        is_ref = stage_code_n == norm_group(ref_stage)
      )
    
    pd_or_table2 <- pd_or_table %>%
      mutate(
        IPV = as.character(IPV),
        IPV_Type = as.character(IPV_Type),
        stage_n = norm_group(stage),
        ref_stage_n = norm_group(ref_stage)
      )
    
    pd_nums <- pd_or_table2 %>%
      select(
        IPV, IPV_Type, stage_n, ref_stage_n,
        pd, pd_lower, pd_upper, p_value,
        OR, OR_lower, OR_upper, OR_p
      )
    
    final_table <- prevalence_table2 %>%
      left_join(
        pd_nums,
        by = c(
          "IPV" = "IPV",
          "IPV_Type" = "IPV_Type",
          "stage_code_n" = "stage_n"
        )
      ) %>%
      mutate(
        model_tag = tag,
        reference_stage = ref_stage,
        
        prevalence_fmt = sprintf(
          "%.2f%% (%.2f–%.2f)",
          mean * 100, lower * 100, upper * 100
        ),
        
        pd_fmt = case_when(
          is_ref ~ "Reference",
          !is.na(pd) ~ sprintf(
            "%.2f%% (%.2f – %.2f)",
            pd * 100, pd_lower * 100, pd_upper * 100
          ),
          TRUE ~ "Contrast unavailable"
        ),
        
        or_fmt = case_when(
          is_ref ~ "Reference",
          !is.na(OR) ~ sprintf(
            "%.2f (%.2f–%.2f)",
            OR, OR_lower, OR_upper
          ),
          TRUE ~ "Contrast unavailable"
        ),
        
        p_fmt = case_when(
          is_ref ~ "",
          TRUE ~ format_p(OR_p)
        )
      ) %>%
      select(
        model_tag,
        IPV,
        `IPV Type` = IPV_Type,
        `Reproductive Stage` = stage_code,
        `Reference Stage` = reference_stage,
        `Prevalence (with 95% CI)` = prevalence_fmt,
        `Prevalence Difference vs Reference (with 95% CI)` = pd_fmt,
        `Odds Ratio vs Reference (with 95% CI)` = or_fmt,
        `P value` = p_fmt,
        mean, lower, upper,
        pd, pd_lower, pd_upper,
        OR, OR_lower, OR_upper,
        OR_p
      ) %>%
      mutate(
        IPV = factor(IPV, levels = ipvs),
        `IPV Type` = factor(`IPV Type`, levels = ipv_labels[ipvs]),
        `Reproductive Stage` = factor(`Reproductive Stage`, levels = stage_levels)
      ) %>%
      arrange(IPV, `Reproductive Stage`)
  }
  
  cat("Finished: ", tag, "\n", sep = "")
  print(dplyr::slice_head(final_table, n = print_n))
  
  list(
    tag = tag,
    prevalence = prevalence_table,
    PD_OR = pd_or_table,
    final = final_table,
    model_failures = model_failures_table
  )
}


# ============================================================
# 12) Run GLMM models
# ============================================================

all_results <- list()
final_all_models_long <- list()
model_failures_all <- list()

for (tg in run_tags) {
  
  parts <- strsplit(tg, "_", fixed = TRUE)[[1]]
  
  control_age <- parts[1] == "ageAdj"
  weighted <- parts[2] == "wt"
  
  res <- run_one_combo(
    df = data_clean,
    ipvs = ipvs,
    control_age = control_age,
    weighted = weighted,
    ref_stage = reference_stage,
    stage_var = "reproductive_stage",
    stage_levels = stage_levels,
    age_var = "v012",
    print_n = 20
  )
  
  all_results[[paste0("prev_", res$tag)]] <- res$prevalence
  all_results[[paste0("PDOR_", res$tag)]] <- res$PD_OR
  all_results[[paste0("final_", res$tag)]] <- res$final
  
  final_all_models_long[[res$tag]] <- res$final
  model_failures_all[[res$tag]] <- res$model_failures
}

final_all_models_long <- bind_rows(final_all_models_long)
model_failures_all <- bind_rows(model_failures_all)


# ============================================================
# 13) Publication-style summary tables
# ============================================================

final_primary_long <- final_all_models_long %>%
  filter(model_tag == primary_model_tag) %>%
  select(
    `IPV Type`,
    `Reproductive Stage`,
    `Reference Stage`,
    `Prevalence (with 95% CI)`,
    `Prevalence Difference vs Reference (with 95% CI)`,
    `Odds Ratio vs Reference (with 95% CI)`,
    `P value`
  ) %>%
  arrange(`IPV Type`, `Reproductive Stage`)

final_primary_wide <- final_primary_long %>%
  pivot_wider(
    id_cols = c(`Reproductive Stage`, `Reference Stage`),
    names_from = `IPV Type`,
    values_from = c(
      `Prevalence (with 95% CI)`,
      `Prevalence Difference vs Reference (with 95% CI)`,
      `Odds Ratio vs Reference (with 95% CI)`,
      `P value`
    ),
    names_glue = "{`IPV Type`}__{.value}"
  ) %>%
  arrange(factor(`Reproductive Stage`, levels = stage_levels, ordered = TRUE))

final_all_models_compact <- final_all_models_long %>%
  select(
    model_tag,
    `IPV Type`,
    `Reproductive Stage`,
    `Reference Stage`,
    `Prevalence (with 95% CI)`,
    `Prevalence Difference vs Reference (with 95% CI)`,
    `Odds Ratio vs Reference (with 95% CI)`,
    `P value`
  ) %>%
  arrange(model_tag, `IPV Type`, `Reproductive Stage`)


# ============================================================
# 14) Diagnostics
# ============================================================

mapping_counts <- data_clean %>%
  summarise(
    n_total = n(),
    n_included = sum(!is.na(reproductive_stage)),
    n_excluded = sum(is.na(reproductive_stage)),
    pct_included = round(100 * n_included / n_total, 2),
    pct_excluded = round(100 * n_excluded / n_total, 2)
  )

stage_counts <- data_clean %>%
  count(reproductive_stage, name = "n") %>%
  mutate(
    percent = round(100 * n / sum(n, na.rm = TRUE), 2)
  ) %>%
  arrange(factor(reproductive_stage, levels = stage_levels, ordered = TRUE))

stage_by_parity <- data_clean %>%
  count(reproductive_stage, parity_cat, name = "n") %>%
  group_by(reproductive_stage) %>%
  mutate(
    percent_within_stage = round(100 * n / sum(n, na.rm = TRUE), 2)
  ) %>%
  ungroup()

exclusion_summary <- data_clean %>%
  count(mapping_exclusion_reason, sort = TRUE) %>%
  mutate(
    percent = round(100 * n / sum(n), 2)
  )

missing_summary <- data.frame(
  variable = c(
    "reproductive_stage",
    "v012", "v000", "d005", "country",
    "v222", "v213", "v201", "v214", "v218",
    "phy_IPV", "sex_IPV", "physex_IPV", "emo_IPV", "any_IPV"
  ),
  missing_count = c(
    sum(is.na(data_clean$reproductive_stage)),
    sum(is.na(data_clean$v012)),
    sum(is.na(data_clean$v000)),
    sum(is.na(data_clean$d005)),
    sum(is.na(data_clean$country)),
    sum(is.na(data_clean$v222)),
    sum(is.na(data_clean$v213)),
    sum(is.na(data_clean$v201)),
    sum(is.na(data_clean$v214)),
    sum(is.na(data_clean$v218)),
    sum(is.na(data_clean$phy_IPV)),
    sum(is.na(data_clean$sex_IPV)),
    sum(is.na(data_clean$physex_IPV)),
    sum(is.na(data_clean$emo_IPV)),
    sum(is.na(data_clean$any_IPV))
  ),
  total_count = nrow(data_clean)
) %>%
  mutate(
    missing_percent = round(100 * missing_count / total_count, 2)
  )

# Survey-level sample table
analysis_df_for_sample_table <- data_clean %>%
  filter(!is.na(reproductive_stage), !is.na(any_IPV))

survey_sample_long <- analysis_df_for_sample_table %>%
  mutate(
    survey_code = as.character(v000),
    country_code = substr(survey_code, 1, 2),
    survey_year = as.numeric(v007)
  ) %>%
  group_by(survey_code, country_code) %>%
  summarise(
    sample_size = n(),
    year_min = min(survey_year, na.rm = TRUE),
    year_max = max(survey_year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    year = ifelse(
      is.finite(year_min) & is.finite(year_max) & year_min != year_max,
      paste0(year_min, "–", year_max),
      as.character(year_min)
    )
  ) %>%
  arrange(country_code, survey_code) %>%
  mutate(No = row_number()) %>%
  select(
    No,
    `survey code` = survey_code,
    country_code,
    year,
    `sample size` = sample_size
  )


# ============================================================
# 15) Export
# ============================================================

export_list <- c(
  list(
    mapping_counts = mapping_counts,
    stage_counts = stage_counts,
    stage_by_parity = stage_by_parity,
    exclusion_summary = exclusion_summary,
    missing_summary = missing_summary,
    model_failures = model_failures_all,
    anyipv_na_before = summary_anyipv_na_before,
    anyipv_na_after = summary_anyipv_na_after,
    desc_long = desc_long,
    desc_wide = desc_wide,
    final_primary_long = final_primary_long,
    final_primary_wide = final_primary_wide,
    final_all_models_compact = final_all_models_compact,
    final_all_models_long = final_all_models_long,
    survey_sample_long = survey_sample_long
  ),
  all_results
)

if (save_excel) {
  write_xlsx(export_list, out_file)
  cat("\nAll results exported to: ", out_file, "\n", sep = "")
}

cat("\nDone.\n")