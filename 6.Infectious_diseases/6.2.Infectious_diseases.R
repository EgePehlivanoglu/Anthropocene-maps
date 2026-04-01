# Reproduce Fig. 3B from Allen et al. (2017)
# Paper: https://www.nature.com/articles/s41467-017-00923-8
# Fig. 3B = weighted model output reweighted by population
# Variable used here: bsm_weight_pop

if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian", repos = "https://cloud.r-project.org")
}

library(librarian)
shelf(raster, terra, sf, sp, rnaturalearth, tidyverse, viridis, maps)

# ---- 1) load archived model predictions ----
predictions_path <- "~/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Base maps/6.Infectious_diseases/Allen 2017/predictions.RData"
stopifnot(file.exists(path.expand(predictions_path)))
load(path.expand(predictions_path))  # loads object: predictions

predictions <- predictions %>%
  dplyr::select(gridid, lon, lat, bsm_weight_pop) %>%
  dplyr::filter(is.finite(lon), is.finite(lat), is.finite(bsm_weight_pop))

# ---- 2) helper functions following the authors' workflow ----
clip_at_sd <- function(x, multiple = 1) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE) * multiple
  x %>%
    pmin(m + s) %>%
    pmax(m - s)
}

template_raster <- function(df) {
  template_df <- df %>%
    dplyr::select(lon, lat, gridid)
  sp::coordinates(template_df) <- c("lon", "lat")
  sp::gridded(template_df) <- TRUE
  sp::proj4string(template_df) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  raster::raster(template_df)
}

# ---- 3) clip values to +/- 2.5 SD, as in the figure caption ----
eid_mean <- mean(predictions$bsm_weight_pop, na.rm = TRUE)
eid_sd <- sd(predictions$bsm_weight_pop, na.rm = TRUE)
eid_lower <- eid_mean - 2.5 * eid_sd
eid_upper <- eid_mean + 2.5 * eid_sd

predictions$bsm_weight_pop_clipped <- clip_at_sd(predictions$bsm_weight_pop, 2.5)

# ---- 4) build the raster the way the authors' scripts do ----
country_outlines <- terra::vect(
  sf::st_transform(
    ne_countries(scale = "medium", returnclass = "sf"),
    crs = "EPSG:4326"
  )
)

fig3b_raster <- predictions %>%
  dplyr::select(x = lon, y = lat, z = bsm_weight_pop_clipped) %>%
  raster::rasterFromXYZ(crs = raster::crs(template_raster(predictions))) %>%
  raster::disaggregate(2, method = "bilinear")

fig3b_raster <- terra::rast(fig3b_raster)
terra::crs(fig3b_raster) <- "EPSG:4326"
fig3b_raster <- terra::mask(fig3b_raster, country_outlines)
names(fig3b_raster) <- "risk"

fig3b_df <- as.data.frame(fig3b_raster, xy = TRUE) %>%
  stats::na.omit()

numcolors <- length(unique(fig3b_df$risk))
world_df <- ggplot2::map_data("world")

# ---- 5) plot ----
fig_3b_plot <- ggplot() +
  geom_polygon(
    aes(x = long, y = lat, group = group),
    data = world_df,
    inherit.aes = FALSE,
    fill = viridis::viridis(1)
  ) +
  geom_raster(
    aes(x = x, y = y, fill = risk),
    data = fig3b_df
  ) +
  geom_path(
    aes(x = long, y = lat, group = group),
    data = world_df,
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.1
  ) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-65, 90), expand = FALSE) +
  scale_fill_gradientn(
    colours = viridis::viridis(numcolors),
    guide = guide_colorbar(
      label = TRUE,
      label.position = "right",
      title = "EID Risk Index"
    )
  ) +
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black"),
    line = element_blank(),
    legend.title = element_text(color = "white", size = 8),
    legend.text = element_text(color = "white", size = 8),
    legend.title.align = 0,
    legend.background = element_blank(),
    legend.position = c(0.11, 0.45)
  ) +
  labs(x = NULL, y = NULL)

fig_3b_plot

# ---- 6) alternative plot with explicit SD-scaled palette ----
fig_3b_sd_plot <- ggplot() +
  geom_polygon(
    aes(x = long, y = lat, group = group),
    data = world_df,
    inherit.aes = FALSE,
    fill = viridis::viridis(1)
  ) +
  geom_raster(
    aes(x = x, y = y, fill = risk),
    data = fig3b_df
  ) +
  geom_path(
    aes(x = long, y = lat, group = group),
    data = world_df,
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.15
  ) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-65, 90), expand = FALSE) +
  scale_fill_gradientn(
    colours = viridis::viridis(256),
    limits = c(eid_lower, eid_upper),
    oob = scales::squish,
    guide = guide_colorbar(
      label = TRUE,
      label.position = "right",
      title = "EID Risk Index\n(+/- 2.5 SD)"
    )
  ) +
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black"),
    line = element_blank(),
    legend.title = element_text(color = "white", size = 8),
    legend.text = element_text(color = "white", size = 8),
    legend.title.align = 0,
    legend.background = element_blank(),
    legend.position = c(0.11, 0.45)
  ) +
  labs(
    title = "Figure 3B reproduction with explicit SD scaling",
    subtitle = "Colour palette limited to 2.5 SD above and below the mean",
    x = NULL, y = NULL
  )

fig_3b_sd_plot

# ---- 7) optional save ----
# ggsave("fig_3b_infectious_diseases.png", fig_3b_plot, width = 11, height = 6.5, dpi = 300)
