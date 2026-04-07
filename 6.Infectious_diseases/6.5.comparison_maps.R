# Publication-style comparison maps:
# 1) side-by-side comparison
# 2) single-map overlay

setwd("/Users/egepehlivanoglu/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Anthropocene-maps")

if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian", repos = "https://cloud.r-project.org")
}

library(librarian)
shelf(terra, sf, ggplot2, dplyr, rnaturalearth, patchwork, viridis)

# ---- 1) source the summed binary raster workflow ----
source_raster_objects <- function(script_path, keep_names) {
  env <- new.env(parent = globalenv())
  sys.source(script_path, envir = env)
  rm(list = setdiff(ls(envir = env), keep_names), envir = env)
  env
}

binary_to_zero <- function(r) {
  terra::ifel(is.na(r), 0, r)
}

deforestation_env <- source_raster_objects(
  "1.Deforestation/1.4.Deforestation_bin.R",
  c("deforestation_bin_SR")
)
conflict_env <- source_raster_objects(
  "2.Conflict/2.5.Conflict_bin.R",
  c("conflict_bin_SR")
)
climate_env <- source_raster_objects(
  "4.Climate_change/4.4.Climate_change_bin_SR.R",
  c("clim_change_binSR")
)
pesticide_env <- source_raster_objects(
  "5.Pesticide_use/5.3.Pesticide_RiskScores_SRbin.R",
  c("pesticide_SR_bin")
)

r_list <- list(
  deforestation = get("deforestation_bin_SR", envir = deforestation_env),
  conflict = get("conflict_bin_SR", envir = conflict_env),
  climate = get("clim_change_binSR", envir = climate_env),
  pesticide = get("pesticide_SR_bin", envir = pesticide_env)
)

target_crs <- terra::crs(r_list$deforestation)
target_res <- c(
  max(vapply(r_list, function(r) terra::res(r)[1], numeric(1))),
  max(vapply(r_list, function(r) terra::res(r)[2], numeric(1)))
)

target_extent <- terra::ext(
  min(vapply(r_list, function(r) terra::ext(r)[1], numeric(1))),
  max(vapply(r_list, function(r) terra::ext(r)[2], numeric(1))),
  min(vapply(r_list, function(r) terra::ext(r)[3], numeric(1))),
  max(vapply(r_list, function(r) terra::ext(r)[4], numeric(1)))
)

common_template_raster <- terra::rast(
  ext = target_extent,
  resolution = target_res,
  crs = target_crs
)

r_aligned <- lapply(r_list, function(r) {
  terra::resample(r, common_template_raster, method = "near")
})
r_aligned_zero <- lapply(r_aligned, binary_to_zero)

binary_sum <- Reduce(`+`, r_aligned_zero)
names(binary_sum) <- "sum_binary"
binary_sum_df <- as.data.frame(binary_sum, xy = TRUE, na.rm = FALSE) %>%
  mutate(sum_binary_cat = factor(sum_binary, levels = c(0, 1, 2, 3, 4)))

# ---- 2) source the infectious disease threshold workflow ----
setwd("/Users/egepehlivanoglu/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Anthropocene-maps/6.Infectious_diseases")
source("6.2.Infectious_diseases.R")

hist_df <- fig3b_df %>%
  dplyr::filter(!is.na(risk) & risk > 0)

top25_cutoff <- stats::quantile(hist_df$risk, probs = 0.75, na.rm = TRUE)
top25_fig3b_df <- hist_df %>%
  dplyr::filter(risk >= top25_cutoff)

# Reproject infectious threshold map to Mollweide so both figures match
top25_vect <- terra::vect(top25_fig3b_df, geom = c("x", "y"), crs = "EPSG:4326")
top25_vect_moll <- terra::project(top25_vect, target_crs)
top25_raster <- terra::rasterize(top25_vect_moll, common_template_raster, field = "risk", fun = "max")
names(top25_raster) <- "risk"
top25_df_moll <- as.data.frame(top25_raster, xy = TRUE, na.rm = TRUE)

# Binary presence layer for overlay
top25_presence <- terra::ifel(is.na(top25_raster), NA, 1)
names(top25_presence) <- "presence"
top25_presence_df <- as.data.frame(top25_presence, xy = TRUE, na.rm = TRUE)
top25_presence_zero <- terra::ifel(is.na(top25_presence), 0, top25_presence)
names(top25_presence_zero) <- "presence"

# ---- 3) shared basemap ----
setwd("/Users/egepehlivanoglu/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Anthropocene-maps")
world_moll <- sf::st_transform(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"),
  crs = target_crs
)

# ---- 4) publication-style side-by-side figure ----
sum_map_plot <- ggplot() +
  geom_tile(
    data = binary_sum_df,
    aes(x = x, y = y, fill = sum_binary_cat)
  ) +
  geom_sf(data = world_moll, fill = NA, color = "grey15", linewidth = 0.05) +
  coord_sf(crs = st_crs(target_crs), expand = FALSE) +
  scale_fill_manual(
    values = c(
      "0" = "#f7f7f7",
      "1" = "#d9e6f2",
      "2" = "#92b4d6",
      "3" = "#477db3",
      "4" = "#123b6d"
    ),
    na.value = "white",
    name = "Summed\nbinary layers"
  ) +
  labs(
    title = "A. Summed binary raster surface",
    subtitle = "Count of aligned binary layers",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

infectious_map_plot <- ggplot() +
  geom_tile(
    data = top25_df_moll,
    aes(x = x, y = y, fill = risk)
  ) +
  geom_sf(data = world_moll, fill = NA, color = "grey15", linewidth = 0.05) +
  coord_sf(crs = st_crs(target_crs), expand = FALSE) +
  scale_fill_gradientn(
    colours = c("#fff7bc", "#fec44f", "#fe9929", "#d95f0e", "#993404"),
    name = "Top 25%\ninfectious risk"
  ) +
  labs(
    title = "B. Infectious disease top 25%",
    subtitle = "Upper quartile of EID risk index",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

comparison_side_by_side <- sum_map_plot + infectious_map_plot +
  plot_annotation(
    title = "Comparison of summed binary pressures and infectious disease hotspots",
    subtitle = "Both layers shown in Mollweide projection with matched extent"
  )

# ---- 5) publication-style single-map overlay ----
overlay_plot <- ggplot() +
  geom_tile(
    data = binary_sum_df,
    aes(x = x, y = y, fill = sum_binary_cat)
  ) +
  geom_tile(
    data = top25_presence_df,
    aes(x = x, y = y),
    fill = "black",
    alpha = 0.28
  ) +
  geom_sf(data = world_moll, fill = NA, color = "grey10", linewidth = 0.05) +
  coord_sf(crs = st_crs(target_crs), expand = FALSE) +
  scale_fill_manual(
    values = c(
      "0" = "#f7f7f7",
      "1" = "#d9e6f2",
      "2" = "#92b4d6",
      "3" = "#477db3",
      "4" = "#123b6d"
    ),
    na.value = "white",
    name = "Summed\nbinary layers"
  ) +
  labs(
    title = "Overlay of summed binary pressures and infectious disease hotspots",
    subtitle = "Black overlay marks cells in the top 25% of infectious disease risk",
    caption = "Recommended for publication only when accompanied by the side-by-side comparison panel.",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# ---- 6) bivariate comparison map ----
bivariate_stack <- c(binary_sum, top25_presence_zero)
names(bivariate_stack) <- c("sum_binary", "presence")
bivariate_df <- as.data.frame(bivariate_stack, xy = TRUE, na.rm = FALSE) %>%
  mutate(
    sum_class = if_else(sum_binary >= 2, "High summed pressure", "Low summed pressure"),
    inf_class = if_else(presence >= 1, "Top 25% infectious risk", "Lower infectious risk"),
    bivar_class = dplyr::case_when(
      sum_class == "Low summed pressure" & inf_class == "Lower infectious risk" ~ "1",
      sum_class == "High summed pressure" & inf_class == "Lower infectious risk" ~ "2",
      sum_class == "Low summed pressure" & inf_class == "Top 25% infectious risk" ~ "3",
      sum_class == "High summed pressure" & inf_class == "Top 25% infectious risk" ~ "4",
      TRUE ~ NA_character_
    ),
    bivar_class = factor(bivar_class, levels = c("1", "2", "3", "4"))
  )

bivariate_palette <- c(
  "1" = "#e8e8e8",
  "2" = "#73ae80",
  "3" = "#be64ac",
  "4" = "#2a5a5b"
)

bivariate_legend_df <- data.frame(
  x = c(1, 2, 1, 2),
  y = c(1, 1, 2, 2),
  bivar_class = factor(c("1", "2", "3", "4"), levels = c("1", "2", "3", "4"))
)

bivariate_map <- ggplot() +
  geom_tile(
    data = bivariate_df,
    aes(x = x, y = y, fill = bivar_class)
  ) +
  geom_sf(data = world_moll, fill = NA, color = "grey10", linewidth = 0.05) +
  coord_sf(crs = st_crs(target_crs), expand = FALSE) +
  scale_fill_manual(
    values = bivariate_palette,
    na.value = "white",
    guide = "none"
  ) +
  labs(
    title = "C. Bivariate comparison of cumulative pressure and infectious disease hotspots",
    subtitle = "2x2 broad classes following bivariate mapping guidance",
    caption = "Summed binary pressure is grouped into low (0-1) and high (2-4). Infectious disease risk is grouped as outside versus inside the top 25% of the distribution.",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

bivariate_legend_plot <- ggplot(
  bivariate_legend_df,
  aes(x = x, y = y, fill = bivar_class)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = bivariate_palette, guide = "none") +
  coord_equal(expand = FALSE) +
  labs(
    x = "Summed binary pressure",
    y = "Infectious disease risk"
  ) +
  annotate("text", x = 1, y = 0.55, label = "Low", size = 3) +
  annotate("text", x = 2, y = 0.55, label = "High", size = 3) +
  annotate("text", x = 0.55, y = 1, label = "Low", angle = 90, size = 3) +
  annotate("text", x = 0.55, y = 2, label = "High", angle = 90, size = 3) +
  theme_void(base_size = 10) +
  theme(
    plot.margin = margin(20, 20, 20, 5),
    axis.title.x = element_text(size = 9, margin = margin(t = 10)),
    axis.title.y = element_text(size = 9, angle = 90, margin = margin(r = 10))
  )

bivariate_plot <- bivariate_map + bivariate_legend_plot +
  patchwork::plot_layout(widths = c(5, 1.3))

comparison_side_by_side
overlay_plot
bivariate_plot

# Optional exports
# ggsave("6.5.side_by_side_comparison.png", comparison_side_by_side, width = 14, height = 6.5, dpi = 300)
# ggsave("6.5.overlay_comparison.png", overlay_plot, width = 8, height = 6.5, dpi = 300)
# ggsave("6.5.bivariate_comparison.png", bivariate_plot, width = 8, height = 6.5, dpi = 300)
