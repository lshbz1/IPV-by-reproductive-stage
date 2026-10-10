# ============================================================
# Clean-window mapping IPV analysis: 5 reproductive-stage groups
# Groups:
#   no child
#   pregnancy
#   0-3 years postpartum
#   3-6 years after childbirth
#   >6 years after childbirth
#
# Outcomes:
#   Any IPV
#   Physical or sexual IPV
#   Physical IPV
#   Emotional IPV
#   Sexual IPV
#
# ============================================================

# ============================================================
# 0) Packages
# ============================================================

packages <- c(
  "dplyr", "tidyr", "stringr", "haven",
  "lme4", "emmeans", "writexl", "purrr"
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

# ============================================================
# 1) Paths and settings
# ============================================================

data_directory <- "/Users/apple/Desktop/STATA data"

out_file <- "/Users/apple/Desktop/IPV_clean_window_5groups_physex_results.xlsx"

save_excel <- TRUE
save_intermediate_csv <- FALSE

combined_csv_out <- "/Users/apple/Desktop/combined_data_DHS_clean_window_5groups_physex.csv"
clean_csv_out    <- "/Users/apple/Desktop/cleaned_data_DHS_clean_window_5groups_physex.csv"


# Model combinations
run_tags <- c(
  "noAge_noWt",
  "noAge_wt",
  "ageAdj_noWt"
)

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

# Five clean-window reproductive-stage groups
group5_levels <- c(
  "no child",
  "pregnancy",
  "0-3 years postpartum",
  "3-6 years after childbirth",
  ">6 years after childbirth"
)

primary_model_tag <- "noAge_noWt"

# ============================================================
# 2) Helper functions
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

# ============================================================
# 3) Read and combine raw DHS .DTA files
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
# 4) Keep Domestic Violence Module sample
# ============================================================

data_clean <- combined_data %>%
  mutate(
    v000 = to_chr(v000),
    across(
      any_of(c("v044", "v222", "v213", "v201", "v007", "d005", "v214", "v012")),
      to_num
    )
  ) %>%
  filter(v044 == 1)

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
# 5) Clean IPV raw variables
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
# 6) Create IPV outcomes
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
    
    # Physical or sexual IPV
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
# 7) Remove problematic surveys
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

# Surveys to remove
remove_codes <- c("BO4", "BO5", "ZM4", "MZ6")

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

# ============================================================
# 8) Five-group clean-window mapping
# ============================================================
# Mapping logic:
#
# no child:
#   v201 == 0 & v213 == 0
#
# pregnancy approximation:
#   v213 == 1 & v214 >= 8
#   OR v213 == 0 & v222 0-1 months postpartum
#
# 0-3 years postpartum occurrence:
#   interviewed at 12-36 months postpartum
#
# 3-6 years after childbirth occurrence:
#   interviewed at 48-72 months postpartum
#
# >6 years after childbirth occurrence:
#   interviewed at >=84 months postpartum
#
# Ambiguous windows are excluded:
#   2-11 months postpartum
#   37-47 months postpartum
#   73-83 months postpartum
# ============================================================

data_clean <- data_clean %>%
  mutate(
    country = substr(v000, 1, 2),
    DHS_version = suppressWarnings(as.integer(substr(v000, 3, 3))),
    
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
    
    group5_clean = case_when(
      is_no_child_clean ~ "no child",
      is_pregnancy_clean ~ "pregnancy",
      is_pp_0_3_clean ~ "0-3 years postpartum",
      is_pp_3_6_clean ~ "3-6 years after childbirth",
      is_pp_gt6_clean ~ ">6 years after childbirth",
      TRUE ~ NA_character_
    ),
    
    group5_clean = factor(
      group5_clean,
      levels = group5_levels
    ),
    
    mapping_exclusion_reason = case_when(
      !is.na(group5_clean) ~ "included in 5-group clean mapping",
      
      v213 == 1 & (is.na(v214) | v214 < 8) ~
        "excluded: early/mid pregnancy; recall window mostly before pregnancy",
      
      v213 == 0 & !is.na(v222) & v222 >= 2 & v222 < 12 ~
        "excluded: early postpartum; recall window overlaps pregnancy/pre-childbirth",
      
      v213 == 0 & !is.na(v222) & v222 > 36 & v222 < 48 ~
        "excluded: recall window overlaps 0-3y and 3-6y postpartum",
      
      v213 == 0 & !is.na(v222) & v222 > 72 & v222 < 84 ~
        "excluded: recall window overlaps 3-6y and >6y postpartum",
      
      TRUE ~ "excluded: unclassified or missing reproductive timing"
    )
  )

cat("\n5-group clean mapping counts:\n")
print(table(data_clean$group5_clean, useNA = "ifany"))

cat("\nMapping exclusion reasons:\n")
print(table(data_clean$mapping_exclusion_reason, useNA = "ifany"))

if (save_intermediate_csv) {
  write.csv(data_clean, clean_csv_out, row.names = FALSE)
}

# ============================================================
# 9) Descriptive summary by 5 groups
# ============================================================

summarise_by_group <- function(df,
                               group_var = "group5_clean",
                               ipvs = ipvs,
                               ipv_labels = ipv_labels,
                               include_total = TRUE) {
  
  group_sym <- rlang::sym(group_var)
  
  long <- df %>%
    filter(!is.na(!!group_sym)) %>%
    select(!!group_sym, all_of(ipvs)) %>%
    pivot_longer(
      cols = all_of(ipvs),
      names_to = "IPV",
      values_to = "val"
    ) %>%
    group_by(!!group_sym, IPV) %>%
    summarise(
      sample_size = sum(!is.na(val)),
      cases = sum(val == 1, na.rm = TRUE),
      prob = cases / sample_size,
      .groups = "drop"
    ) %>%
    mutate(
      prob = ifelse(is.nan(prob), NA_real_, prob),
      IPV_label = ipv_labels[IPV]
    ) %>%
    rename(group_code = !!group_sym)
  
  if (include_total) {
    total_long <- df %>%
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
        group_code = "Total",
        IPV_label = ipv_labels[IPV]
      ) %>%
      select(group_code, IPV, IPV_label, sample_size, cases, prob)
    
    long <- bind_rows(long, total_long)
  }
  
  long <- long %>%
    mutate(
      group_code = as.character(group_code),
      group_code = factor(group_code, levels = c(group5_levels, "Total")),
      IPV = factor(IPV, levels = ipvs),
      IPV_label = factor(IPV_label, levels = ipv_labels[ipvs]),
      prevalence_fmt = sprintf("%.2f%%", prob * 100)
    ) %>%
    arrange(IPV, group_code)
  
  wide <- long %>%
    select(group_code, IPV_label, sample_size, cases, prob, prevalence_fmt) %>%
    pivot_wider(
      id_cols = group_code,
      names_from = IPV_label,
      values_from = c(sample_size, cases, prob, prevalence_fmt),
      names_sep = "__"
    ) %>%
    arrange(factor(group_code, levels = c(group5_levels, "Total"), ordered = TRUE))
  
  list(long = long, wide = wide)
}

sum5 <- summarise_by_group(
  data_clean,
  group_var = "group5_clean",
  ipvs = ipvs,
  ipv_labels = ipv_labels,
  include_total = TRUE
)

summary_5_long <- sum5$long
summary_5_wide <- sum5$wide

cat("\nDescriptive summary: 5-group clean mapping\n")
print(summary_5_wide)

# ============================================================
# 10) GLMM model functions
# ============================================================

make_groups_for_model <- function(df) {
  df %>%
    mutate(
      group_code = factor(
        as.character(group5_clean),
        levels = group5_levels
      ),
      country = substr(v000, 1, 2),
      DHS_version = suppressWarnings(as.integer(substr(v000, 3, 3)))
    )
}

fit_ipv_once <- function(df,
                         ipv,
                         control_age = FALSE,
                         weighted = FALSE,
                         ref_group = "no child",
                         age_var = "v012") {
  
  dat <- make_groups_for_model(df) %>%
    prepare_weights() %>%
    filter(!is.na(group_code), !is.na(.data[[ipv]])) %>%
    droplevels()
  
  if (!ref_group %in% levels(dat$group_code)) {
    stop("Reference group is not available: ", ref_group)
  }
  
  use_age <- if (control_age) pick_age_var(dat, age_var) else NULL
  
  rhs <- "group_code"
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
  
  # Prevalence on response scale
  emm_resp <- emmeans(fit, ~ group_code, type = "response")
  
  prev_df <- as.data.frame(summary(emm_resp)) %>%
    transmute(
      IPV = ipv,
      IPV_Type = ipv_labels[ipv],
      group_code = as.character(group_code),
      mean = prob,
      lower = asymp.LCL,
      upper = asymp.UCL
    )
  
  ref_idx <- which(levels(emm_resp@grid$group_code) == ref_group)
  if (length(ref_idx) == 0) stop("Reference group not found in emmeans grid.")
  
  # Prevalence difference on response scale
  emm_resp2 <- regrid(emm_resp, transform = "response")
  pd_contr <- contrast(emm_resp2, method = "trt.vs.ctrl", ref = ref_idx)
  
  pd_df <- as.data.frame(summary(pd_contr, infer = c(TRUE, TRUE))) %>%
    tidyr::separate(
      contrast,
      into = c("group", "ref_group"),
      sep = " - ",
      remove = TRUE,
      extra = "merge",
      fill = "right"
    ) %>%
    transmute(
      IPV = ipv,
      IPV_Type = ipv_labels[ipv],
      group = as.character(group),
      ref_group = as.character(ref_group),
      pd = estimate,
      pd_lower = asymp.LCL,
      pd_upper = asymp.UCL,
      p_value = p.value
    ) %>%
    filter(group != ref_group)
  
  # OR on logit scale
  emm_link <- emmeans(fit, ~ group_code)
  or_contr <- contrast(emm_link, method = "trt.vs.ctrl", ref = ref_idx)
  
  or_df <- as.data.frame(summary(or_contr, infer = c(TRUE, TRUE))) %>%
    tidyr::separate(
      contrast,
      into = c("group", "ref_group"),
      sep = " - ",
      remove = TRUE,
      extra = "merge",
      fill = "right"
    ) %>%
    transmute(
      IPV = ipv,
      IPV_Type = ipv_labels[ipv],
      group = as.character(group),
      ref_group = as.character(ref_group),
      OR = exp(estimate),
      OR_lower = exp(asymp.LCL),
      OR_upper = exp(asymp.UCL),
      OR_p = p.value
    ) %>%
    filter(group != ref_group)
  
  pd_or <- pd_df %>%
    left_join(or_df, by = c("IPV", "IPV_Type", "group", "ref_group"))
  
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
                          ref_group = "no child",
                          age_var = "v012",
                          print_n = 12) {
  
  tag <- paste0(
    ifelse(control_age, "ageAdj", "noAge"),
    "_",
    ifelse(weighted, "wt", "noWt")
  )
  
  cat("\n==============================\n")
  cat("Running: ", tag, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  
  res_list <- lapply(
    ipvs,
    function(x) {
      fit_ipv_once(
        df = df,
        ipv = x,
        control_age = control_age,
        weighted = weighted,
        ref_group = ref_group,
        age_var = age_var
      )
    }
  )
  
  names(res_list) <- ipvs
  
  prevalence_table <- bind_rows(lapply(res_list, `[[`, "prevalence"))
  pd_or_table <- bind_rows(lapply(res_list, `[[`, "pd_or"))
  
  prevalence_table2 <- prevalence_table %>%
    mutate(
      IPV = as.character(IPV),
      IPV_Type = as.character(IPV_Type),
      group_code_n = norm_group(group_code),
      is_ref = group_code_n == norm_group(ref_group)
    )
  
  pd_or_table2 <- pd_or_table %>%
    mutate(
      IPV = as.character(IPV),
      IPV_Type = as.character(IPV_Type),
      group_n = norm_group(group),
      ref_group_n = norm_group(ref_group)
    )
  
  pd_nums <- pd_or_table2 %>%
    select(
      IPV, IPV_Type, group_n, ref_group_n,
      pd, pd_lower, pd_upper, p_value,
      OR, OR_lower, OR_upper, OR_p
    )
  
  final_table <- prevalence_table2 %>%
    left_join(
      pd_nums,
      by = c(
        "IPV" = "IPV",
        "IPV_Type" = "IPV_Type",
        "group_code_n" = "group_n"
      )
    ) %>%
    mutate(
      model_tag = tag,
      
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
        TRUE ~ "Not estimable"
      ),
      
      or_fmt = case_when(
        is_ref ~ "Reference",
        !is.na(OR) ~ sprintf(
          "%.2f (%.2f–%.2f)",
          OR, OR_lower, OR_upper
        ),
        TRUE ~ ""
      ),
      
      p_fmt = case_when(
        is_ref ~ "",
        TRUE ~ format_p(p_value)
      )
    ) %>%
    select(
      model_tag,
      IPV,
      `IPV Type` = IPV_Type,
      `Reproductive Stage` = group_code,
      `Prevalence (with 95% CI)` = prevalence_fmt,
      `Prevalence Difference vs no child (with 95% CI)` = pd_fmt,
      `Odds Ratio vs no child (with 95% CI)` = or_fmt,
      `P value` = p_fmt,
      mean, lower, upper,
      pd, pd_lower, pd_upper,
      OR, OR_lower, OR_upper,
      p_value
    ) %>%
    mutate(
      IPV = factor(IPV, levels = ipvs),
      `IPV Type` = factor(`IPV Type`, levels = ipv_labels[ipvs]),
      `Reproductive Stage` = factor(`Reproductive Stage`, levels = group5_levels)
    ) %>%
    arrange(IPV, `Reproductive Stage`)
  
  cat("Finished: ", tag, "\n", sep = "")
  print(dplyr::slice_head(final_table, n = print_n))
  
  list(
    tag = tag,
    prevalence = prevalence_table,
    PD_OR = pd_or_table,
    final = final_table
  )
}

# ============================================================
# 11) Run GLMM models
# ============================================================

all_results <- list()
final_all_models_long <- list()

for (tg in run_tags) {
  
  parts <- strsplit(tg, "_", fixed = TRUE)[[1]]
  
  control_age <- parts[1] == "ageAdj"
  weighted <- parts[2] == "wt"
  
  res <- run_one_combo(
    df = data_clean,
    ipvs = ipvs,
    control_age = control_age,
    weighted = weighted,
    ref_group = "no child",
    age_var = "v012",
    print_n = 12
  )
  
  all_results[[paste0("prev_", res$tag)]] <- res$prevalence
  all_results[[paste0("PDOR_", res$tag)]] <- res$PD_OR
  all_results[[paste0("final_", res$tag)]] <- res$final
  
  final_all_models_long[[res$tag]] <- res$final
}

final_all_models_long <- bind_rows(final_all_models_long)

# ============================================================
# 12) Publication-style summary tables
# ============================================================

final_primary_long <- final_all_models_long %>%
  filter(model_tag == primary_model_tag) %>%
  select(
    `IPV Type`,
    `Reproductive Stage`,
    `Prevalence (with 95% CI)`,
    `Prevalence Difference vs no child (with 95% CI)`,
    `Odds Ratio vs no child (with 95% CI)`,
    `P value`
  ) %>%
  arrange(`IPV Type`, `Reproductive Stage`)

final_primary_wide <- final_primary_long %>%
  pivot_wider(
    id_cols = `Reproductive Stage`,
    names_from = `IPV Type`,
    values_from = c(
      `Prevalence (with 95% CI)`,
      `Prevalence Difference vs no child (with 95% CI)`,
      `Odds Ratio vs no child (with 95% CI)`,
      `P value`
    ),
    names_glue = "{`IPV Type`}__{.value}"
  ) %>%
  arrange(factor(`Reproductive Stage`, levels = group5_levels, ordered = TRUE))

final_all_models_compact <- final_all_models_long %>%
  select(
    model_tag,
    `IPV Type`,
    `Reproductive Stage`,
    `Prevalence (with 95% CI)`,
    `Prevalence Difference vs no child (with 95% CI)`,
    `Odds Ratio vs no child (with 95% CI)`,
    `P value`
  ) %>%
  arrange(model_tag, `IPV Type`, `Reproductive Stage`)

# ============================================================
# 13) Diagnostics and export
# ============================================================

mapping_counts <- data_clean %>%
  summarise(
    n_total = n(),
    n_group5_included = sum(!is.na(group5_clean)),
    n_group5_excluded = sum(is.na(group5_clean)),
    n_no_child = sum(group5_clean == "no child", na.rm = TRUE),
    n_pregnancy = sum(group5_clean == "pregnancy", na.rm = TRUE),
    n_0_3_postpartum = sum(group5_clean == "0-3 years postpartum", na.rm = TRUE),
    n_3_6_after_childbirth = sum(group5_clean == "3-6 years after childbirth", na.rm = TRUE),
    n_gt6_after_childbirth = sum(group5_clean == ">6 years after childbirth", na.rm = TRUE)
  )

exclusion_summary <- data_clean %>%
  count(mapping_exclusion_reason, sort = TRUE)

missing_summary <- data.frame(
  variable = c(
    "group5_clean",
    "v012", "v000", "d005", "country",
    "v222", "v213", "v201", "v214",
    "phy_IPV", "sex_IPV", "physex_IPV", "emo_IPV", "any_IPV"
  ),
  missing_count = c(
    sum(is.na(data_clean$group5_clean)),
    sum(is.na(data_clean$v012)),
    sum(is.na(data_clean$v000)),
    sum(is.na(data_clean$d005)),
    sum(is.na(data_clean$country)),
    sum(is.na(data_clean$v222)),
    sum(is.na(data_clean$v213)),
    sum(is.na(data_clean$v201)),
    sum(is.na(data_clean$v214)),
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

export_list <- c(
  list(
    mapping_counts = mapping_counts,
    exclusion_summary = exclusion_summary,
    missing_summary = missing_summary,
    anyipv_na_before = summary_anyipv_na_before,
    anyipv_na_after = summary_anyipv_na_after,
    desc5_long = summary_5_long,
    desc5_wide = summary_5_wide,
    final_primary_long = final_primary_long,
    final_primary_wide = final_primary_wide,
    final_all_models_compact = final_all_models_compact,
    final_all_models_long = final_all_models_long
  ),
  all_results
)

if (save_excel) {
  write_xlsx(export_list, out_file)
  cat("\nAll results exported to: ", out_file, "\n", sep = "")
}

cat("\nDone.\n")

