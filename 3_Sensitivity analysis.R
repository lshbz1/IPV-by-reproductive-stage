# ============================================================
# Sensitivity analysis with 6-month reproductive stage intervals
# ============================================================

# ============================================================
# 0) Packages
# ============================================================

packages <- c(
  "dplyr", "tidyr", "stringr", "haven",
  "lme4", "emmeans", "writexl", "purrr",
  "ggplot2", "tibble", "scales"
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
library(ggplot2)
library(tibble)
library(scales)

# ============================================================
# 1) Paths and settings
# ============================================================

data_directory <- "/Users/apple/Desktop/STATA data"

save_intermediate_csv <- FALSE

combined_csv_out <- "/Users/apple/Desktop/combined_data_DHS_clean_window_5groups_physex.csv"
clean_csv_out    <- "/Users/apple/Desktop/cleaned_data_DHS_clean_window_5groups_physex.csv"

# Surveys to remove
remove_codes <- c("BO4", "BO5", "ZM4", "MZ6")


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
# 8) Reproductive-stage mapping:
#    no child + pregnancy unchanged;
#    postpartum every 6 months until 18Y;
#    >18Y combined
# ============================================================

library(ggplot2)

# 0 to 18 years = 0 to 216 months
pp_breaks_m <- seq(0, 216, by = 6)

pp_labels_m <- paste0(
  pp_breaks_m[-length(pp_breaks_m)],
  "-",
  pp_breaks_m[-1],
  "m postpartum"
)

format_year_num <- function(x) {
  x <- x / 12
  ifelse(x %% 1 == 0, as.character(x), sprintf("%.1f", x))
}

pp_labels_y <- paste0(
  format_year_num(pp_breaks_m[-length(pp_breaks_m)]),
  "-",
  format_year_num(pp_breaks_m[-1]),
  "Y"
)

data_clean <- data_clean %>%
  mutate(
    country = substr(v000, 1, 2),
    DHS_version = suppressWarnings(as.integer(substr(v000, 3, 3))),
    
    # keep no child unchanged
    is_no_child_clean = v201 == 0 & v213 == 0,
    
    # keep pregnancy / pre unchanged
    is_pregnancy_clean = (
      (v213 == 1 & !is.na(v214) & v214 >= 8) |
        (v213 == 0 & !is.na(v222) & v222 >= 0 & v222 <= 1)
    ),
    
    # postpartum grouped by interview timing, every 6 months
    pp6_group = case_when(
      v213 == 0 & !is.na(v222) & v222 >= 0 & v222 < 216 ~
        as.character(
          cut(
            v222,
            breaks = pp_breaks_m,
            labels = pp_labels_m,
            right = FALSE,
            include.lowest = TRUE
          )
        ),
      
      v213 == 0 & !is.na(v222) & v222 >= 216 ~
        ">18Y postpartum",
      
      TRUE ~ NA_character_
    ),
    
    group6 = case_when(
      is_no_child_clean ~ "no child",
      is_pregnancy_clean ~ "pregnancy",
      !is.na(pp6_group) ~ pp6_group,
      TRUE ~ NA_character_
    )
  )

group6_levels <- c(
  "no child",
  "pregnancy",
  pp_labels_m,
  ">18Y postpartum"
)

group6_axis_labels <- c(
  "No child",
  "Currently pregnant",
  pp_labels_y,
  ">18Y"
)

data_clean <- data_clean %>%
  mutate(
    group6 = factor(group6, levels = group6_levels),
    
    group6_axis_label = case_when(
      as.character(group6) == "no child" ~ "No child",
      as.character(group6) == "pregnancy" ~ "Currently pregnant",
      as.character(group6) == ">18Y postpartum" ~ ">18Y",
      as.character(group6) %in% pp_labels_m ~
        group6_axis_labels[match(as.character(group6), group6_levels)],
      TRUE ~ NA_character_
    ),
    
    group6_axis_label = factor(
      group6_axis_label,
      levels = group6_axis_labels
    ),
    
    mapping_exclusion_reason = case_when(
      !is.na(group6) ~ "included in 6-month mapping",
      v213 == 1 & (is.na(v214) | v214 < 8) ~ "excluded: early/mid pregnancy",
      TRUE ~ "excluded: unclassified or missing reproductive timing"
    )
  )

cat("\n6-month group counts:\n")
print(table(data_clean$group6, useNA = "ifany"))

cat("\nMapping exclusion reasons:\n")
print(table(data_clean$mapping_exclusion_reason, useNA = "ifany"))


# ============================================================
# 9) IPV outcomes for reference-style figure
# ============================================================

ipvs_for_plot <- c(
  "any_IPV",
  "phy_IPV",
  "emo_IPV",
  "sex_IPV"
)

ipv_labels_plot <- c(
  any_IPV = "Any IPV",
  phy_IPV = "Physical IPV",
  emo_IPV = "Emotional IPV",
  sex_IPV = "Sexual IPV"
)


# ============================================================
# 10) Model functions
# ============================================================

make_groups_for_model <- function(df) {
  df %>%
    mutate(
      group_code = factor(
        as.character(group6),
        levels = group6_levels
      ),
      country = substr(v000, 1, 2),
      DHS_version = suppressWarnings(as.integer(substr(v000, 3, 3)))
    )
}

fit_ipv_prevalence <- function(df,
                               ipv,
                               control_age = FALSE,
                               weighted = FALSE,
                               age_var = "v012") {
  
  dat <- make_groups_for_model(df) %>%
    prepare_weights() %>%
    filter(!is.na(group_code), !is.na(.data[[ipv]])) %>%
    droplevels()
  
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
  
  emm_resp <- emmeans(fit, ~ group_code, type = "response")
  
  prev_df <- as.data.frame(summary(emm_resp)) %>%
    transmute(
      IPV = ipv,
      IPV_Type = ipv_labels_plot[ipv],
      group_code = as.character(group_code),
      mean = prob,
      lower = asymp.LCL,
      upper = asymp.UCL
    )
  
  list(
    prevalence = prev_df,
    model = fit
  )
}

run_model_set <- function(df,
                          ipvs_for_plot,
                          control_age = FALSE,
                          weighted = FALSE,
                          age_var = "v012") {
  
  tag <- paste0(
    ifelse(control_age, "ageAdj", "noAge"),
    "_",
    ifelse(weighted, "wt", "noWt")
  )
  
  cat("\n==============================\n")
  cat("Running model: ", tag, "\n", sep = "")
  
  res_list <- lapply(
    ipvs_for_plot,
    function(x) {
      fit_ipv_prevalence(
        df = df,
        ipv = x,
        control_age = control_age,
        weighted = weighted,
        age_var = age_var
      )
    }
  )
  
  names(res_list) <- ipvs_for_plot
  
  prevalence_table <- bind_rows(lapply(res_list, `[[`, "prevalence")) %>%
    mutate(model_tag = tag)
  
  list(
    tag = tag,
    prevalence = prevalence_table,
    models = res_list
  )
}


# ============================================================
# 11) Run two models:
#     1. no age adjustment, no weight
#     2. age adjusted, no weight
# ============================================================

res_noAge_noWt <- run_model_set(
  df = data_clean,
  ipvs_for_plot = ipvs_for_plot,
  control_age = FALSE,
  weighted = FALSE,
  age_var = "v012"
)

res_ageAdj_noWt <- run_model_set(
  df = data_clean,
  ipvs_for_plot = ipvs_for_plot,
  control_age = TRUE,
  weighted = FALSE,
  age_var = "v012"
)


# ============================================================
# 12) Prepare plot data
# ============================================================

make_plot_data <- function(prevalence_df) {
  
  prevalence_df %>%
    mutate(
      group_code = factor(group_code, levels = group6_levels),
      
      x_label = case_when(
        as.character(group_code) == "no child" ~ "No child",
        as.character(group_code) == "pregnancy" ~ "Currently pregnant",
        as.character(group_code) == ">18Y postpartum" ~ ">18Y",
        as.character(group_code) %in% pp_labels_m ~
          group6_axis_labels[match(as.character(group_code), group6_levels)],
        TRUE ~ as.character(group_code)
      ),
      
      x_label = factor(x_label, levels = group6_axis_labels),
      x_num = as.numeric(x_label),
      
      prevalence_pct = mean * 100,
      lower_pct = lower * 100,
      upper_pct = upper * 100,
      
      IPV_Type = factor(
        IPV_Type,
        levels = ipv_labels_plot[ipvs_for_plot]
      )
    ) %>%
    arrange(IPV_Type, x_num)
}

plot_data_noAge_noWt <- make_plot_data(res_noAge_noWt$prevalence)
plot_data_ageAdj_noWt <- make_plot_data(res_ageAdj_noWt$prevalence)


# ============================================================
# 13) Updated clean reference-style plotting section
#     - keep original bottom x-axis labels
#     - add 0Y, 3Y, 6Y... above x-axis
#     - remove >18Y from plot
#     - black border on every point
#     - 6-month points
#     - 3-year average bars
#     - Any IPV average labels
# ============================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(scales)

# ------------------------------------------------------------
# 13.0 Check required objects
# ------------------------------------------------------------

required_objects <- c(
  "plot_data_noAge_noWt",
  "plot_data_ageAdj_noWt",
  "pp_labels_y",
  "group6_axis_labels"
)

missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    "\nPlease run the data-preparation/modeling sections before this plotting section."
  )
}

# ------------------------------------------------------------
# 13.1 Remove >18Y
# ------------------------------------------------------------

plot_data_noAge_noWt <- plot_data_noAge_noWt %>%
  filter(as.character(x_label) != ">18Y")

plot_data_ageAdj_noWt <- plot_data_ageAdj_noWt %>%
  filter(as.character(x_label) != ">18Y")

# Keep bottom x-axis labels, but remove >18Y label
group6_axis_labels_plot <- group6_axis_labels[
  as.character(group6_axis_labels) != ">18Y"
]

# ------------------------------------------------------------
# 13.2 IPV order and visual settings
# ------------------------------------------------------------

ipv_order <- c(
  "Any IPV",
  "Physical IPV",
  "Emotional IPV",
  "Sexual IPV"
)

ipv_color_values <- c(
  "Any IPV" = "grey15",
  "Physical IPV" = "grey45",
  "Emotional IPV" = "grey60",
  "Sexual IPV" = "grey35"
)

ipv_fill_values <- c(
  "Any IPV" = "grey15",
  "Physical IPV" = "grey65",
  "Emotional IPV" = "grey85",
  "Sexual IPV" = "white"
)

ipv_shape_values <- c(
  "Any IPV" = 21,
  "Physical IPV" = 22,
  "Emotional IPV" = 24,
  "Sexual IPV" = 23
)

plot_data_noAge_noWt <- plot_data_noAge_noWt %>%
  mutate(
    IPV_Type = factor(IPV_Type, levels = ipv_order)
  )

plot_data_ageAdj_noWt <- plot_data_ageAdj_noWt %>%
  mutate(
    IPV_Type = factor(IPV_Type, levels = ipv_order)
  )


# ------------------------------------------------------------
# 13.3 Calculate 3-year average bars
# ------------------------------------------------------------

make_bar_data_3y <- function(plot_data) {
  
  plot_data %>%
    mutate(
      pp_index = case_when(
        as.character(x_label) %in% pp_labels_y ~
          match(as.character(x_label), pp_labels_y),
        TRUE ~ NA_integer_
      ),
      
      bar_group = case_when(
        !is.na(pp_index) & pp_index >= 1  & pp_index <= 6  ~ "0-3Y",
        !is.na(pp_index) & pp_index >= 7  & pp_index <= 12 ~ "3-6Y",
        !is.na(pp_index) & pp_index >= 13 & pp_index <= 18 ~ "6-9Y",
        !is.na(pp_index) & pp_index >= 19 & pp_index <= 24 ~ "9-12Y",
        !is.na(pp_index) & pp_index >= 25 & pp_index <= 30 ~ "12-15Y",
        !is.na(pp_index) & pp_index >= 31 & pp_index <= 36 ~ "15-18Y",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(bar_group)) %>%
    group_by(IPV_Type, bar_group) %>%
    summarise(
      y_bar = mean(prevalence_pct, na.rm = TRUE),
      x_start = min(x_num, na.rm = TRUE) - 0.35,
      x_end = max(x_num, na.rm = TRUE) + 0.35,
      x_mid = mean(c(x_start, x_end)),
      .groups = "drop"
    ) %>%
    mutate(
      IPV_Type = factor(IPV_Type, levels = ipv_order),
      bar_group = factor(
        bar_group,
        levels = c(
          "0-3Y",
          "3-6Y",
          "6-9Y",
          "9-12Y",
          "12-15Y",
          "15-18Y"
        )
      )
    ) %>%
    arrange(IPV_Type, bar_group)
}

bar_data_noAge_noWt <- make_bar_data_3y(plot_data_noAge_noWt)
bar_data_ageAdj_noWt <- make_bar_data_3y(plot_data_ageAdj_noWt)


# ------------------------------------------------------------
# 13.4 Vertical dashed lines and year labels above x-axis
# ------------------------------------------------------------
# Keep 18Y as the final boundary, but do not plot >18Y.

vline_data <- tibble(
  xintercept = c(
    2.5,   # 0Y
    8.5,   # 3Y
    14.5,  # 6Y
    20.5,  # 9Y
    26.5,  # 12Y
    32.5,  # 15Y
    38.5   # 18Y
  ),
  year_label = c(
    "0Y",
    "3Y",
    "6Y",
    "9Y",
    "12Y",
    "15Y",
    "18Y"
  )
)

x_axis_marker_data <- tibble(
  x = c(
    2.5,
    8.5,
    14.5,
    20.5,
    26.5,
    32.5,
    38.5
  ),
  label = c(
    "0Y",
    "3Y",
    "6Y",
    "9Y",
    "12Y",
    "15Y",
    "18Y"
  )
)


# ------------------------------------------------------------
# 13.5 Any IPV labels on 3-year average bars
# ------------------------------------------------------------

make_any_bar_label_data <- function(bar_data) {
  
  bar_data %>%
    filter(IPV_Type == "Any IPV") %>%
    mutate(
      label = sprintf("%.1f%%", y_bar),
      y_label = y_bar + 0.75
    )
}

any_bar_labels_noAge_noWt <- make_any_bar_label_data(bar_data_noAge_noWt)
any_bar_labels_ageAdj_noWt <- make_any_bar_label_data(bar_data_ageAdj_noWt)


# ============================================================
# 14) Updated plot function
# ============================================================

plot_ipv_reference_style_clean <- function(plot_data,
                                           bar_data,
                                           any_bar_labels,
                                           title,
                                           subtitle,
                                           out_file) {
  
  y_max <- 35
  marker_y <- 0.9
  
  # x-axis breaks and labels:
  # keep original bottom labels, but no >18Y
  x_breaks <- seq_along(group6_axis_labels_plot)
  x_labels <- group6_axis_labels_plot
  
  p <- ggplot() +
    
    geom_vline(
      data = vline_data,
      aes(xintercept = xintercept),
      linetype = "dashed",
      linewidth = 0.38,
      color = "grey55",
      alpha = 0.85
    ) +
    
    geom_errorbar(
      data = plot_data,
      aes(
        x = x_num,
        ymin = lower_pct,
        ymax = upper_pct
      ),
      width = 0,
      linewidth = 0.32,
      color = "grey72",
      alpha = 0.55,
      show.legend = FALSE
    ) +
    
    geom_line(
      data = plot_data,
      aes(
        x = x_num,
        y = prevalence_pct,
        group = IPV_Type,
        color = IPV_Type
      ),
      linewidth = 0.62,
      alpha = 0.95
    ) +
    
    geom_point(
      data = plot_data,
      aes(
        x = x_num,
        y = prevalence_pct,
        fill = IPV_Type,
        shape = IPV_Type
      ),
      color = "black",
      size = 2.35,
      stroke = 0.70
    ) +
    
    geom_segment(
      data = bar_data,
      aes(
        x = x_start,
        xend = x_end,
        y = y_bar,
        yend = y_bar,
        color = IPV_Type
      ),
      linewidth = 1.15,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    
    geom_text(
      data = any_bar_labels,
      aes(
        x = x_mid,
        y = y_label,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.3,
      fontface = "bold",
      color = "grey10"
    ) +
    
    geom_text(
      data = x_axis_marker_data,
      aes(
        x = x,
        y = marker_y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.4,
      fontface = "bold",
      color = "grey10",
      vjust = 0
    ) +
    
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(0.7, 38.8),
      expand = expansion(mult = c(0.005, 0.01))
    ) +
    
    scale_y_continuous(
      breaks = seq(0, y_max, by = 5),
      labels = function(x) paste0(x, "%"),
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    scale_color_manual(
      values = ipv_color_values
    ) +
    
    scale_fill_manual(
      values = ipv_fill_values
    ) +
    
    scale_shape_manual(
      values = ipv_shape_values
    ) +
    
    guides(
      color = "none",
      fill = guide_legend(
        title = "IPV type",
        override.aes = list(
          shape = c(21, 22, 24, 23),
          color = "black",
          size = 3.0,
          stroke = 0.8
        )
      ),
      shape = "none"
    ) +
    
    labs(
      x = "Reproductive stage / time since childbirth",
      y = "Prevalence",
      title = title,
      subtitle = subtitle,
      fill = "IPV type"
    ) +
    
    coord_cartesian(
      ylim = c(0, y_max),
      clip = "off"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 17,
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 10.5,
        hjust = 0,
        color = "grey25"
      ),
      
      legend.position = "top",
      legend.justification = "center",
      legend.title = element_text(size = 10.5),
      legend.text = element_text(size = 10),
      
      axis.title.x = element_text(size = 12, margin = margin(t = 24)),
      axis.title.y = element_text(size = 12),
      
      # Bottom labels retained
      axis.text.x = element_text(
        size = 7,
        color = "grey25",
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = 10, color = "grey25"),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(
        color = "grey90",
        linewidth = 0.38
      ),
      
      axis.line.x = element_line(color = "grey20", linewidth = 0.45),
      axis.line.y = element_blank(),
      
      plot.margin = margin(t = 10, r = 20, b = 48, l = 20)
    )
  
  print(p)
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 18,
    height = 7.2,
    dpi = 300
  )
  
  cat("\nPlot saved to: ", out_file, "\n", sep = "")
  
  p
}


# ============================================================
# 15) Draw two separate figures
# ============================================================

p_noAge_noWt <- plot_ipv_reference_style_clean(
  plot_data = plot_data_noAge_noWt,
  bar_data = bar_data_noAge_noWt,
  any_bar_labels = any_bar_labels_noAge_noWt,
  title = "IPV Prevalence by Reproductive Stage and Time Since Childbirth",
  subtitle = "Points show 6-month intervals; dashed vertical lines and horizontal bars show 3-year averages; model: no age adjustment, no weights",
  out_file = "/Users/apple/Desktop/IPV_6month_clean_noAge_noWt.png"
)

p_ageAdj_noWt <- plot_ipv_reference_style_clean(
  plot_data = plot_data_ageAdj_noWt,
  bar_data = bar_data_ageAdj_noWt,
  any_bar_labels = any_bar_labels_ageAdj_noWt,
  title = "IPV Prevalence by Reproductive Stage and Time Since Childbirth",
  subtitle = "Points show 6-month intervals; dashed vertical lines and horizontal bars show 3-year averages; model: age-adjusted, no weights",
  out_file = "/Users/apple/Desktop/IPV_6month_clean_ageAdj_noWt.png"
)