# ============================================================
# Figure: prevalence (Panel A) + OR forest plot (Panel B)
# Reference group for OR: No children
# ============================================================

# =========================
# 0) Packages
# =========================
library(readxl)
library(janitor)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(grid)

# =========================
# 1) Read results
# =========================
path_tbl <- "/Users/apple/Desktop/IPV_clean_window_5groups_no_child_ref_pairwise_OR_results.xlsx"

# Numeric columns are available in sheets such as final_noAge_noWt,
# final_ageAdj_noWt, final_noAge_wt, and final_ageAdj_wt.
model_to_plot <- "noAge_noWt"
sheet_name <- paste0("final_", model_to_plot)

df_raw <- read_excel(path_tbl, sheet = sheet_name) %>%
  clean_names()

# =========================
# 2) Standardise variables
# =========================
stage_levels <- c(
  "No children",
  "Pregnancy",
  "0–3y after childbirth",
  "3–6y after childbirth",
  ">6y after childbirth"
)

ipv_levels <- c(
  "Any IPV",
  "Emotional IPV",
  "Physical IPV",
  "Sexual IPV"
)

clean_stage <- function(x) {
  key <- x %>%
    as.character() %>%
    str_squish() %>%
    str_to_lower() %>%
    str_replace_all("[–—−]", "-") %>%
    str_remove_all("[()]" )

  case_when(
    key %in% c("no child", "no children") ~ "No children",
    key %in% c("pregnancy", "currently pregnant") ~ "Pregnancy",
    key %in% c("0-3 years postpartum", "0-3 years after childbirth") ~
      "0–3y after childbirth",
    key %in% c("3-6 years postpartum", "3-6 years after childbirth") ~
      "3–6y after childbirth",
    key %in% c(">6 years postpartum", ">6 years after childbirth") ~
      ">6y after childbirth",
    TRUE ~ as.character(x)
  )
}

dfp <- df_raw %>%
  transmute(
    model_tag,
    ipv_type = ipv_type,
    repro_stage = clean_stage(reproductive_stage),
    prev = mean,
    prev_lower = lower,
    prev_upper = upper,
    or = or,
    or_lower = or_lower,
    or_upper = or_upper,
    or_p = or_p
  ) %>%
  mutate(
    repro_stage = factor(repro_stage, levels = stage_levels),
    ipv_type = factor(ipv_type, levels = ipv_levels),

    # Explicitly retain the reference group in Panel B.
    or = if_else(repro_stage == "No children", 1, or),
    or_lower = if_else(repro_stage == "No children", 1, or_lower),
    or_upper = if_else(repro_stage == "No children", 1, or_upper),

    or_label = if_else(
      repro_stage == "No children",
      "Reference",
      sprintf("%.2f (%.2f–%.2f)", or, or_lower, or_upper)
    )
  ) %>%
  filter(
    !is.na(repro_stage),
    !is.na(ipv_type),
    ipv_type != "Physical or sexual IPV"
  ) %>%
  droplevels()

stopifnot(
  all(is.na(dfp$prev_lower) | is.na(dfp$prev_upper) |
        dfp$prev_upper >= dfp$prev_lower),
  all(is.na(dfp$or_lower) | is.na(dfp$or_upper) |
        dfp$or_upper >= dfp$or_lower)
)

# =========================
# 3) Shared style
# =========================
pal_stage <- c(
  "No children"             = "#7570B3",
  "Pregnancy"               = "#D95F02",
  "0–3y after childbirth"   = "#D84A8A",
  "3–6y after childbirth"   = "#1B9E77",
  ">6y after childbirth"    = "#E6AB02"
)

grey_line_col <- "grey55"
line_w <- 0.8
err_w <- 0.7
pt_size <- 3.0
pt_stroke <- 0.7
pd <- position_dodge(width = 0.38)

shape_map <- c(
  "Any IPV"       = 21,
  "Emotional IPV" = 24,
  "Physical IPV"  = 22,
  "Sexual IPV"    = 23
)

shared_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey85"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.tag = element_text(face = "bold", size = 13)
  )

# =========================
# 4) Panel A: prevalence
# =========================
p_a <- ggplot(
  dfp,
  aes(x = repro_stage, y = prev, group = ipv_type)
) +
  geom_line(
    color = grey_line_col,
    linewidth = line_w,
    alpha = 0.85,
    position = pd
  ) +
  geom_linerange(
    data = dfp %>% filter(!is.na(prev_lower), !is.na(prev_upper)),
    aes(ymin = prev_lower, ymax = prev_upper),
    color = grey_line_col,
    linewidth = err_w,
    alpha = 0.9,
    position = pd
  ) +
  geom_point(
    aes(fill = repro_stage, shape = ipv_type),
    size = pt_size,
    color = grey_line_col,
    stroke = pt_stroke,
    position = pd
  ) +
  scale_fill_manual(
    values = pal_stage,
    name = "Reproductive Stage",
    breaks = stage_levels,
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      order = 1,
      override.aes = list(
        shape = 21,
        size = 4,
        color = grey_line_col,
        stroke = pt_stroke
      )
    )
  ) +
  scale_shape_manual(
    values = shape_map,
    name = "IPV Type",
    breaks = ipv_levels,
    guide = guide_legend(nrow = 1, byrow = TRUE, order = 2)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 0.35, 0.05),
    limits = c(0, 0.35),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    title = "Model-estimated IPV prevalence",
    subtitle = paste0("Clean-window analysis; model: ", model_to_plot),
    x = NULL,
    y = "Prevalence",
    tag = "A"
  ) +
  shared_theme +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = unit(4, "pt"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(margin = margin(t = 6), size = 9),
    plot.margin = margin(t = 8, r = 70, b = 8, l = 12)
  )

# Right-side direct labels for Panel A
right_level <- tail(stage_levels, 1)
lab_data <- dfp %>%
  filter(repro_stage == right_level) %>%
  group_by(ipv_type) %>%
  summarise(
    x = as.numeric(repro_stage) + 0.45,
    y = mean(prev, na.rm = TRUE),
    label = first(ipv_type),
    .groups = "drop"
  )

p_a <- p_a +
  coord_cartesian(clip = "off") +
  ggrepel::geom_text_repel(
    data = lab_data,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    direction = "y",
    nudge_x = 0.08,
    min.segment.length = Inf,
    segment.color = NA,
    box.padding = 0.25,
    point.padding = 0.15,
    size = 3.0,
    color = "grey25"
  )

# =========================
# 5) Panel B: OR forest plot
#    Reference = No children
# =========================
# Panel order from top to bottom:
# Any IPV, Emotional IPV, Physical IPV, Sexual IPV
dfp_or <- dfp %>%
  mutate(
    ipv_type_plot = factor(ipv_type, levels = ipv_levels),
    stage_plot = factor(repro_stage, levels = rev(stage_levels))
  )

p_b <- ggplot(
  dfp_or,
  aes(x = or, y = stage_plot, fill = repro_stage)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey45"
  ) +
  geom_linerange(
    data = dfp_or %>% filter(repro_stage != "No children"),
    aes(xmin = or_lower, xmax = or_upper),
    color = grey_line_col,
    linewidth = err_w
  ) +
  geom_point(
    aes(shape = ipv_type),
    size = pt_size,
    color = grey_line_col,
    stroke = pt_stroke
  ) +
  geom_text(
    aes(x = 2.02, label = or_label),
    hjust = 0,
    size = 2.8,
    color = "grey25",
    show.legend = FALSE
  ) +
  facet_grid(
    rows = vars(ipv_type_plot),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(values = pal_stage, breaks = stage_levels) +
  scale_shape_manual(values = shape_map, breaks = ipv_levels) +
  scale_x_log10(
    breaks = c(0.75, 1, 1.25, 1.5, 1.75, 2),
    limits = c(0.72, 2.55),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Odds ratios for IPV by reproductive stage",
    subtitle = "Reference group: No children",
    x = "Odds ratio (log scale)",
    y = NULL,
    tag = "B"
  ) +
  shared_theme +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 9,
      hjust = 1
    ),
    axis.text.y = element_text(size = 8.8),
    plot.margin = margin(t = 8, r = 92, b = 8, l = 8)
  )

# =========================
# 6) Combine and save
# =========================
final_plot <- p_a + p_b +
  plot_layout(widths = c(1.08, 1.12))

print(final_plot)

ggsave(
  "/Users/apple/Desktop/ipv_prevalence_and_or_no_child_ref_4types.png",
  final_plot,
  width = 20,
  height = 7.8,
  dpi = 320,
  bg = "white"
)

ggsave(
  "/Users/apple/Desktop/ipv_prevalence_and_or_no_child_ref_4types.pdf",
  final_plot,
  width = 16,
  height = 7.8,
  device = cairo_pdf,
  bg = "white"
)
