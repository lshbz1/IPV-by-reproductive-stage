# ============================================================
# Four-axis plot with "No child, not pregnant" added as baseline point to axes 2, 3, and 4
# ============================================================

packages <- c(
  "readxl", "dplyr", "stringr", "ggplot2", "scales", "showtext", "sysfonts", "ggh4x"
)
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) install.packages(pkg)
}

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(showtext)
library(sysfonts)
library(ggh4x)

# -------------------------
# 0) Font settings
# -------------------------
title_font_ok <- tryCatch({
  font_add_google("Playfair Display", "playfair")
  TRUE
}, error = function(e) FALSE)

body_font_ok <- tryCatch({
  font_add_google("Inter", "inter")
  TRUE
}, error = function(e) FALSE)

showtext_auto()
output_dpi <- 320
showtext_opts(dpi = output_dpi)
title_font <- if (title_font_ok) "playfair" else "serif"
base_font  <- if (body_font_ok) "inter" else "sans"

# -------------------------
# 1) Read model results
# -------------------------
# Update this path to match the actual location of the results file.
result_file <- "/Users/apple/Desktop/IPV_mutually_exclusive_stage_results_FIXED.xlsx"

df <- read_excel(
  path = result_file,
  sheet = "final_all_models_long"
)

# -------------------------
# 2) Select the model and IPV outcome
# -------------------------
model_to_plot <- "noAge_noWt"
ipv_to_plot   <- "Any IPV"
plot_out_file <- "/Users/apple/Desktop/IPV_four_axes_noAge_noWeight_any_IPV.png"

# -------------------------
# 3) Prepare the main dataset
# -------------------------
main_df <- df %>%
  filter(
    model_tag == model_to_plot,
    `IPV Type` == ipv_to_plot
  ) %>%
  mutate(
    mean_pct  = mean  * 100,
    lower_pct = lower * 100,
    upper_pct = upper * 100
  )

if (nrow(main_df) == 0) {
  stop(
    "main_df is empty: model_to_plot = \"", model_to_plot,
    "\" or ipv_to_plot = \"", ipv_to_plot,
    "\" has no matching rows in the data.\n",
    "Please check the available values below:\n",
    "  Available model_tag values: ", paste(unique(df$model_tag), collapse = ", "), "\n",
    "  Available IPV Type values: ", paste(unique(df$`IPV Type`), collapse = ", ")
  )
}

no_child_row <- main_df %>%
  filter(`Reproductive Stage` == "No child, not pregnant") %>%
  slice(1)

baseline_value <- no_child_row$mean_pct[1]

# -------------------------
# 4) Observed stage data for the four axes
# -------------------------
real_plot_df <- main_df %>%
  mutate(
    transition_axis = case_when(
      `Reproductive Stage` %in% c(
        "No child, not pregnant",
        "First pregnancy",
        "First child, 0-3 years",
        "First child, 3-6 years",
        "First child, >6 years"
      ) ~ "1. First-child transition",
      
      `Reproductive Stage` %in% c(
        "Second pregnancy",
        "Second child, 0-3 years",
        "Second child, 3-6 years",
        "Second child, >6 years"
      ) ~ "2. Second-child transition",
      
      `Reproductive Stage` %in% c(
        "Third pregnancy",
        "Third child, 0-3 years",
        "Third child, 3-6 years",
        "Third child, >6 years"
      ) ~ "3. Third-child transition",
      
      `Reproductive Stage` %in% c(
        "Fourth-or-higher pregnancy",
        "Fourth-or-higher child, 0-3 years",
        "Fourth-or-higher child, 3-6 years",
        "Fourth-or-higher child, >6 years"
      ) ~ "4. Fourth-or-higher transition",
      
      TRUE ~ NA_character_
    ),
    
    stage_on_axis = case_when(
      `Reproductive Stage` == "No child, not pregnant" ~ "No child\nnot pregnant",
      str_detect(`Reproductive Stage`, "pregnancy")     ~ "Pregnancy",
      str_detect(`Reproductive Stage`, "0-3 years")     ~ "0-3 years",
      str_detect(`Reproductive Stage`, "3-6 years")     ~ "3-6 years",
      str_detect(`Reproductive Stage`, ">6 years")      ~ ">6 years",
      TRUE ~ as.character(`Reproductive Stage`)
    ),
    
    stage_order = case_when(
      stage_on_axis == "No child\nnot pregnant" ~ 0,
      stage_on_axis == "Pregnancy"               ~ 1,
      stage_on_axis == "0-3 years"               ~ 2,
      stage_on_axis == "3-6 years"               ~ 3,
      stage_on_axis == ">6 years"                ~ 4,
      TRUE ~ NA_real_
    ),
    
    point_type = ifelse(
      `Reproductive Stage` == "No child, not pregnant",
      "No-child baseline",
      "Observed stage"
    )
  ) %>%
  filter(!is.na(transition_axis))

# -------------------------
# 5) Duplicate the no-child baseline on axes 2, 3, and 4
# -------------------------
no_child_added_df <- bind_rows(
  no_child_row %>%
    mutate(transition_axis = "2. Second-child transition",
           stage_on_axis = "No child\nnot pregnant",
           stage_order = 0,
           point_type = "Added no-child baseline"),
  no_child_row %>%
    mutate(transition_axis = "3. Third-child transition",
           stage_on_axis = "No child\nnot pregnant",
           stage_order = 0,
           point_type = "Added no-child baseline"),
  no_child_row %>%
    mutate(transition_axis = "4. Fourth-or-higher transition",
           stage_on_axis = "No child\nnot pregnant",
           stage_order = 0,
           point_type = "Added no-child baseline")
)

axis_levels <- c(
  "1. First-child transition",
  "2. Second-child transition",
  "3. Third-child transition",
  "4. Fourth-or-higher transition"
)

plot_df <- bind_rows(real_plot_df, no_child_added_df) %>%
  mutate(transition_axis = factor(transition_axis, levels = axis_levels)) %>%
  arrange(transition_axis, stage_order)

# -------------------------
# 6) Solid-line data (observed stages)
# The first axis includes the no-child point; the other axes use a dashed segment instead
# -------------------------
solid_line_df <- real_plot_df %>%
  filter(
    # For the first axis, retain all observations, including the no-child group
    (transition_axis == "1. First-child transition") |
      # For the other axes, exclude the no-child group because it is a duplicated baseline
      (transition_axis != "1. First-child transition" & `Reproductive Stage` != "No child, not pregnant")
  ) %>%
  mutate(transition_axis = factor(transition_axis, levels = axis_levels)) %>%
  arrange(transition_axis, stage_order)

# -------------------------
# 7) Dashed-line data: baseline to pregnancy (axes 2, 3, and 4 only)
# -------------------------
dashed_line_df <- plot_df %>%
  filter(
    transition_axis %in% c(
      "2. Second-child transition",
      "3. Third-child transition",
      "4. Fourth-or-higher transition"
    ),
    stage_order %in% c(0, 1)
  ) %>%
  arrange(transition_axis, stage_order)

# -------------------------
# 8) Color palette
# -------------------------
stage_colors <- c(
  "No child\nnot pregnant" = "#9080B8",
  "Pregnancy"               = "#E8922E",
  "0-3 years"               = "#E0568F",
  "3-6 years"               = "#2E9E82",
  ">6 years"                = "#E8B93A"
)

stage_labels <- c(
  "No child\nnot pregnant" = "No children",
  "Pregnancy"               = "Pregnancy",
  "0-3 years"               = "0-3y after childbirth",
  "3-6 years"               = "3-6y after childbirth",
  ">6 years"                = ">6y after childbirth"
)

col_line      <- "#2C2A28"   # Dark solid line
col_dashed    <- "#B05A3A"   # Prominent dashed line
col_errorbar  <- "grey55"
col_grid      <- "#ECE6DC"
col_bg        <- "#FFFFFF"
col_text_dark <- "#33302B"
col_text_mid  <- "#6F6A61"

# Y-axis range: start at 5% to reduce empty space
y_min <- 5
y_max <- max(plot_df$upper_pct, na.rm = TRUE) * 1.05

plot_df <- plot_df %>%
  mutate(stage_on_axis = factor(stage_on_axis, levels = names(stage_colors)))

# -------------------------
# 8b) X-axis labels
# -------------------------
x_axis_labels <- list(
  "1. First-child transition" = c(
    "No child\nnot pregnant", "First\npregnancy",
    "First child\n0-3 years", "First child\n3-6 years", "First child\n>6 years"
  ),
  "2. Second-child transition" = c(
    "No child\nnot pregnant", "Second\npregnancy",
    "Second child\n0-3 years", "Second child\n3-6 years", "Second child\n>6 years"
  ),
  "3. Third-child transition" = c(
    "No child\nnot pregnant", "Third\npregnancy",
    "Third child\n0-3 years", "Third child\n3-6 years", "Third child\n>6 years"
  ),
  "4. Fourth-or-higher transition" = c(
    "No child\nnot pregnant", "Fourth-or-higher\npregnancy",
    "Fourth-or-higher child\n0-3 years", "Fourth-or-higher child\n3-6 years", "Fourth-or-higher child\n>6 years"
  )
)

x_scales <- lapply(x_axis_labels, function(labs) {
  scale_x_continuous(
    breaks = 0:4,
    labels = labs,
    expand = expansion(mult = c(0.09, 0.09))
  )
})

# -------------------------
# 9) Plotting
# -------------------------
p <- ggplot() +
  
  # Baseline reference line across panels
  geom_hline(
    yintercept = baseline_value,
    linetype = "dotted", linewidth = 0.6,
    color = "grey55", alpha = 0.6
  ) +
  
  # Dashed line: baseline to pregnancy (axes 2, 3, and 4)
  geom_line(
    data = dashed_line_df,
    aes(x = stage_order, y = mean_pct, group = transition_axis),
    linewidth = 1.0, linetype = "dashed", color = col_dashed, alpha = 0.9,
    lineend = "round"
  ) +
  
  # Solid line: observed trajectory (the first axis includes the full path from the no-child group)
  geom_line(
    data = solid_line_df,
    aes(x = stage_order, y = mean_pct, group = transition_axis),
    linewidth = 1.5, color = col_line, lineend = "round"
  ) +
  
  # Error bars
  geom_errorbar(
    data = plot_df,
    aes(x = stage_order, ymin = lower_pct, ymax = upper_pct,
        alpha = point_type),
    width = 0.10, linewidth = 0.55, color = col_errorbar
  ) +
  
  # Data points
  geom_point(
    data = plot_df,
    aes(x = stage_order, y = mean_pct,
        fill = stage_on_axis, alpha = point_type),
    size = 4.4, shape = 21, color = "white", stroke = 0.7
  ) +
  
  # Value labels (not bold)
  geom_text(
    data = plot_df,
    aes(x = stage_order, y = upper_pct, label = sprintf("%.1f%%", mean_pct)),
    vjust = -1.1, size = 4.0, family = base_font, color = col_text_dark
  ) +
  
  scale_alpha_manual(
    values = c("Observed stage" = 1, "Added no-child baseline" = 0.55),
    guide = "none"
  ) +
  
  scale_fill_manual(
    values = stage_colors,
    labels = stage_labels,
    name = "Stage",
    breaks = names(stage_colors)
  ) +
  
  facet_wrap(~ transition_axis, ncol = 1, scales = "free_x") +
  
  ggh4x::facetted_pos_scales(x = x_scales) +
  
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(y_min, y_max),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  
  labs(
    subtitle = paste0(
      "Model: ", model_to_plot,
      "  ·  point color = stage (see legend)  ·  faded point = repeated no-child baseline"
    ),
    x = NULL,
    y = "Adjusted Prevalence of IPV",
    caption = "Error bars represent 95% confidence intervals. Lines connect observed stage means; dashed lines indicate the hypothetical baseline-to-pregnancy segment for higher parities."
  ) +
  
  theme_minimal(base_size = 16, base_family = base_font) +
  theme(
    plot.subtitle = element_text(size = 12.5, color = col_text_mid, face = "italic",
                                 margin = margin(b = 10)),
    plot.caption = element_text(size = 10, color = col_text_mid, hjust = 0, margin = margin(t = 14)),
    
    legend.position = "top",
    legend.title = element_text(size = 12, color = col_text_dark),
    legend.text = element_text(size = 11, color = col_text_mid),
    legend.margin = margin(b = 10),
    legend.key = element_blank(),
    
    strip.text = element_text(family = base_font, size = 14,
                              color = col_text_dark, hjust = 0, margin = margin(t = 4, b = 10)),
    strip.background = element_blank(),
    
    panel.spacing = unit(2.0, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = col_grid, linewidth = 0.4),
    panel.background = element_rect(fill = col_bg, color = NA),
    panel.border = element_blank(),
    
    # Move y-axis labels closer to the axis by reducing the right margin; leave the x-axis unchanged
    axis.text.y = element_text(size = 11, color = col_text_mid, margin = margin(r = 2)),
    axis.title.y = element_text(size = 13, margin = margin(r = 6), color = col_text_mid),
    axis.text.x = element_text(size = 12, color = col_text_mid, lineheight = 0.9),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    # Slightly reduce the overall margins to give the panels more space
    plot.background = element_rect(fill = col_bg, color = NA),
    plot.margin = margin(20, 20, 16, 20)
  )

print(p)

# Keep the saved figure height at 18 inches so that each panel has sufficient vertical space
ggsave(
  filename = plot_out_file,
  plot = p,
  width = 10,
  height = 18,
  dpi = output_dpi,
  bg = col_bg
)