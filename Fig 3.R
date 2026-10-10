# =========================================================
# Panel A: Any IPV (child 1–3y) — left map (bottom colorbar) + right histogram
# Panel B: PD = prev(child 1–3y) − prev(no child) — left map (bottom colorbar) + right histogram
# Fully transparent background; panels stacked vertically; histogram width ≈ 25% of map width (4:1)
# Legends are placed below each map and use continuous color scales
# =========================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)
library(stringr)
library(ggrepel)
library(scales)
library(patchwork)
library(grid)  # unit()

# -------------------------------
# 0) Parameters
# -------------------------------
path1 <- "/Users/apple/Desktop/data_for_analysis/cleaned_data_DHS_with_mod.csv"

# Panel A: Any IPV —  (yellow → orange → reddish brown)
prev_upper <- 0.60
prev_palette_steps <- c("#FFB347","#FF6F61","#C94CBE","#6A00FF","#3A00A1")
binwidth_prev <- 0.025  # 2.5%

# Panel B: PD —  (negative blue → positive red; coarse steps for histogram consistency)
pd_breaks     <- c(-0.10,  -0.05, 0, 0.05, 0.10,0.20,  0.30)
pd_palette_steps <- c("#2C5AA0","#6FA4D6","#BFD7EE","#E9E9E9","#F5B1B1","#E57373","#C73B3B")
pd_palette_steps <- c("#2C7FB8","#BFD7EE","#F7FCF5","#FDD0A2","#F5B1B1","#FB6A4A","#CB181D")

show_country_labels_on_map <- TRUE

# -------------------------------
# 1) Read data and construct the two latest-survey indicators
# -------------------------------
# Select data file
df <- readr::read_csv(path1, show_col_types = FALSE)

# ------------------------------------------------------------------
# Preprocessing: group by the alphabetic prefix of v000 (e.g., "NG") and retain rows with the largest numeric part (e.g., "7")
# If all v000 values within a prefix group have no numeric part (all NA), retain all rows in that group (usually a single survey round)
# ------------------------------------------------------------------
library(stringr)
df_latest_v000 <- df %>%
  mutate(
    v000 = as.character(v000),
    v000_letters = str_extract(v000, "^[A-Za-z]+"),
    v000_digits  = as.integer(str_extract(v000, "[0-9]+"))
  ) %>%
  group_by(v000_letters) %>%
  mutate(
    # If the group contains at least one non-missing numeric part, take the maximum; otherwise set to NA
    max_digit = if (all(is.na(v000_digits))) NA_integer_ else max(v000_digits, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(
    # Retain rows whose numeric part equals the group maximum; if the group has no numeric part (max_digit is NA), retain all rows
    ( !is.na(max_digit) & !is.na(v000_digits) & v000_digits == max_digit ) |
      is.na(max_digit)
  ) %>%
  select(-v000_letters, -v000_digits, -max_digit)

# -----------------------
# Any IPV (child 1–3y)
# -----------------------
any_prev_latest <- df_latest_v000 %>%
  filter(group_code == "child age [1,3]") %>%
  transmute(
    iso2 = str_sub(v000, 1, 2),
    any  = as.numeric(any_IPV)
  ) %>%
  mutate(
    country = countrycode::countrycode(iso2, "iso2c", "country.name"),
    iso3c   = countrycode::countrycode(iso2, "iso2c", "iso3c")
  ) %>%
  filter(!is.na(country)) %>%
  group_by(iso3c, country) %>%
  summarise(prev = mean(any, na.rm = TRUE), .groups = "drop") %>%
  arrange(country) %>%
  distinct(iso3c, .keep_all = TRUE) %>%   # Retain one row per country after applying the latest-v000 filtering above
  select(iso3c, country, prev)

# -----------------------
# PD = child 1–3y − no child
# -----------------------
pd_latest <- df_latest_v000 %>%
  filter(group_code %in% c("child age [1,3]", "no child")) %>%
  transmute(
    iso2 = str_sub(v000, 1, 2),
    grp  = group_code,
    any  = as.numeric(any_IPV)
  ) %>%
  mutate(
    country = countrycode::countrycode(iso2, "iso2c", "country.name"),
    iso3c   = countrycode::countrycode(iso2, "iso2c", "iso3c")
  ) %>%
  filter(!is.na(country)) %>%
  group_by(iso3c, country,grp) %>%
  summarise(prev = mean(any, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = grp, values_from = prev) %>%
  mutate(pd = `child age [1,3]` - `no child`) %>%
  arrange(country) %>%
  filter(!is.na(pd)) %>%
  distinct(iso3c, .keep_all = TRUE) %>%
  select(iso3c, country, pd)

# -------------------------------
# 2) Base map (Robinson projection)
# -------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% mutate(iso3c = iso_a3)
world_rob <- st_transform(world, "ESRI:54030")

map_prev_rob <- st_transform(world %>% left_join(any_prev_latest, by = "iso3c"), "ESRI:54030")
map_pd_rob   <- st_transform(world %>% left_join(pd_latest,      by = "iso3c"), "ESRI:54030")

lab_prev_df <- map_prev_rob %>% filter(!is.na(prev)) %>%
  st_point_on_surface() %>% mutate(country_short = str_trunc(country, 16)) %>%
  { cbind(st_drop_geometry(.), st_coordinates(.)) }

lab_pd_df <- map_pd_rob %>% filter(!is.na(pd)) %>%
  st_point_on_surface() %>% mutate(country_short = str_trunc(country, 16)) %>%
  { cbind(st_drop_geometry(.), st_coordinates(.)) }

grat_rob <- st_transform(st_graticule(lat = seq(-60, 80, 20), lon = seq(-160, 160, 20)), "ESRI:54030")

# -------------------------------
# 3) Themes and legend styling
# -------------------------------
theme_transparent_void <- function(base_size = 22){
  theme_void(base_size = base_size) %+replace%
    theme(plot.background = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA),
          legend.background = element_rect(fill = NA, color = NA),
          legend.key = element_rect(fill = NA, color = NA),
          plot.margin = margin(0,0,0,0))
}
theme_transparent_min <- function(base_size = 21){
  theme_minimal(base_size = base_size) %+replace%
    theme(plot.background = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA),
          legend.background = element_rect(fill = NA, color = NA),
          legend.key = element_rect(fill = NA, color = NA),
          panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0,0,0,0))
}
legend_bottom_small <- guides(
  fill = guide_colorbar(direction = "horizontal", title.position = "top",
                        barheight = unit(3, "mm"), barwidth = unit(70, "mm"))
)

# -------------------------------
# 4) Panel A: Any IPV (map with bottom colorbar + histogram on the right)
# -------------------------------
g_map_prev <- ggplot() +
  geom_sf(data = world_rob, fill = NA, color = "#E6E0D8", linewidth = 0.15) +
  geom_sf(data = map_prev_rob %>% filter(!is.na(prev)),
          aes(fill = prev),  # continuous color scale
          color = "#E6E0D8", linewidth = 0.15) +
  geom_sf(data = grat_rob, color = "#E6E0D8", linewidth = 0.15) +
  scale_fill_gradientn(
    colours = prev_palette_steps,
    limits  = c(0, prev_upper),
    oob     = scales::squish,
    name    = "Postpartum IPV prevalence",
    labels  = label_percent(accuracy = 1)
  ) +
  { if (show_country_labels_on_map)
    geom_text_repel(data = lab_prev_df, aes(X, Y, label = country_short),
                    size = 8, color = "black",
                    min.segment.length = 0, box.padding = 0.15, point.padding = 0.1,
                    max.overlaps = Inf) else NULL } +
  coord_sf(crs = "ESRI:54030", expand = FALSE) +
  theme_transparent_void(12) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 28),
        legend.text  = element_text(size = 28)) +
  legend_bottom_small

# Histogram (use the same continuous color scale, without a legend)
# === Histogram (use the same continuous color scale, without a legend) ===
dist_prev <- any_prev_latest %>% filter(!is.na(prev))
q_prev    <- quantile(dist_prev$prev, probs = c(0.25,0.5,0.75), na.rm = TRUE)
dens_prev <- density(dist_prev$prev, from = 0, to = prev_upper, na.rm = TRUE)
h_prev    <- hist(dist_prev$prev,
                  breaks = seq(0, prev_upper, by = binwidth_prev),
                  plot = FALSE, include.lowest = TRUE, right = FALSE)
dens_prev_scaled <- tibble(
  x = dens_prev$x,
  y = dens_prev$y * (max(h_prev$counts) / max(dens_prev$y))
)

g_hist_prev <- ggplot(dist_prev, aes(x = prev)) +
  annotate("rect", xmin = q_prev[1], xmax = q_prev[3], ymin = 0, ymax = Inf,
           alpha = 0.08, fill = "#888888") +
  geom_histogram(
    aes(fill = after_stat(x)),  # continuous
    binwidth = binwidth_prev, boundary = 0, closed = "left",
    color = alpha("white", 0.6), linewidth = 0.2
  ) +
  scale_fill_gradientn(
    colours = prev_palette_steps,
    limits = c(0, prev_upper),
    oob = scales::squish,
    guide = "none"
  ) +
  geom_line(data = dens_prev_scaled, aes(x = x, y = y),
            linewidth = 1.0, color = "#5A3A3A") +
  geom_vline(xintercept = q_prev[2], linetype = "dashed",
             linewidth = 0.6, color = "#5A3A3A") +
  annotate("text", x = q_prev[2], y = max(h_prev$counts) * 0.95,
           label = paste0("Median ", percent(q_prev[2], 1)),
           hjust = -0.05, vjust = 1, size = 8, color = "#5A3A3A") +  # enlarged to approximately 28-point text
  scale_x_continuous(
    limits = c(0, prev_upper),
    breaks = seq(0, prev_upper, by = 0.10),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = NULL,
    x = "Postpartum IPV prevalence",
    y = "Countries",
    subtitle = paste0("Bin width = ", percent(binwidth_prev, 0.1))
  ) +
  theme_transparent_min(11) +
  theme(
    text = element_text(size = 28),                 # global text size
    axis.title = element_text(size = 28),
    axis.text  = element_text(size = 28),
    plot.subtitle = element_text(size = 28),
    panel.grid.major = element_line(color = "#E2E2E2", linewidth = 0.25),
    legend.position = "none"
  )

# --------------------------------
# 5) Panel B: PD (map + histogram on the right)
# --------------------------------

# === Map ===
g_map_pd <- ggplot() +
  geom_sf(data = world_rob, fill = NA, color = "#E6E0D8", linewidth = 0.15) +
  geom_sf(data = map_pd_rob %>% filter(!is.na(pd)),
          aes(fill = pd),
          color = "#E6E0D8", linewidth = 0.15) +
  geom_sf(data = grat_rob, color = "#E6E0D8", linewidth = 0.15) +
  scale_fill_stepsn(
    colors = pd_palette_steps,          # symmetric blue–red palette
    breaks = pd_breaks,
    limits = range(pd_breaks),
    values = scales::rescale(pd_breaks),
    oob = scales::squish,
    labels = label_percent(accuracy = 1),
    name = "IPV PD (postpartum − no children)"
  ) +
  { if (show_country_labels_on_map)
    geom_text_repel(data = lab_pd_df, aes(X, Y, label = country_short),
                    size = 8, color = "black",
                    min.segment.length = 0, box.padding = 0.15, point.padding = 0.1,
                    max.overlaps = Inf) else NULL } +
  coord_sf(crs = "ESRI:54030", expand = FALSE) +
  theme_transparent_void(12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 28),
    legend.text  = element_text(size = 28)
  ) +
  legend_bottom_small


# === Histogram ===
dist_pd <- pd_latest %>% filter(!is.na(pd))
q_pd    <- quantile(dist_pd$pd, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
dens_pd <- density(dist_pd$pd, from = min(pd_breaks), to = max(pd_breaks), na.rm = TRUE)
h_pd    <- hist(dist_pd$pd,
                breaks = seq(min(pd_breaks), max(pd_breaks), by = 0.02),
                plot = FALSE, include.lowest = TRUE, right = FALSE)
dens_pd_scaled <- tibble(
  x = dens_pd$x,
  y = dens_pd$y * (max(h_pd$counts) / max(dens_pd$y))
)

g_hist_pd <- ggplot(dist_pd, aes(x = pd)) +
  annotate("rect",
           xmin = q_pd[1], xmax = q_pd[3],
           ymin = 0, ymax = Inf,
           alpha = 0.08, fill = "#888888") +
  geom_histogram(
    aes(fill = after_stat(x)),
    binwidth = 0.02, boundary = 0, closed = "left",
    color = alpha("white", 0.6), linewidth = 0.2
  ) +
  scale_fill_stepsn(
    colors = pd_palette_steps,
    breaks = pd_breaks,
    limits = range(pd_breaks),
    values = scales::rescale(pd_breaks),
    oob = scales::squish,
    guide = "none"
  ) +
  geom_line(data = dens_pd_scaled, aes(x = x, y = y),
            linewidth = 1.0, color = "#4A4A4A") +
  geom_vline(xintercept = q_pd[2], linetype = "dashed",
             linewidth = 0.6, color = "#4A4A4A") +
  annotate("text", x = q_pd[2], y = max(h_pd$counts) * 0.95,
           label = paste0("Median ", percent(q_pd[2], 1)),
           hjust = -0.05, vjust = 1, size = 8, color = "#4A4A4A") +  # enlarged annotation text
  scale_x_continuous(
    limits = range(pd_breaks),
    breaks = pd_breaks,
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = NULL,
    x = "IPV PD (postpartum − no children)",
    y = "Countries",
    subtitle = "Bin width = 2 pp"
  ) +
  theme_transparent_min(11) +
  theme(
    text = element_text(size = 28),                # global text size
    axis.title = element_text(size = 28),
    axis.text  = element_text(size = 28),
    plot.subtitle = element_text(size = 28),
    legend.position = "none",
    panel.grid.major = element_line(color = "#E2E2E2", linewidth = 0.25)
  )


# -------------------------------
# 6) Assemble panels (narrower histograms; legends already placed below each map)
# -------------------------------
make_panel <- function(g_map, g_hist){
  wrap_plots(
    A = g_map,
    B = g_hist,
    design = "
AB
AB
",
    widths  = c(4, 1),   # Map : histogram = 4 : 1 (histogram ≈ 25% of map width)
    heights = c(1, 1)
  )
}

panel_any <- make_panel(g_map_prev, g_hist_prev)
panel_pd  <- make_panel(g_map_pd,  g_hist_pd)

# Final layout: Any IPV on top, PD below (no timeline)
p_all <- panel_any / panel_pd + plot_layout(heights = c(1, 1))
print(p_all)

# Export as a transparent, high-resolution PNG
library(ragg)
ggsave("/Users/apple/Desktop/ipv_any_pd_continuous_legends.png",
       p_all, width = 10, height = 10, units = "in", dpi = 300,
       bg = "transparent", device = ragg::agg_png)


