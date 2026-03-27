# Reproduce Fig. 3B from Allen et al. (2017)
# Paper: https://www.nature.com/articles/s41467-017-00923-8
# Fig. 3B = weighted model output reweighted by population
# Variable used here: bsm_weight_pop

# install librarian package quietly if needed
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian", repos = "https://cloud.r-project.org")
}

library(librarian)
shelf(terra, sf, rnaturalearth, tidyverse, viridis)

# ---- 1) load archived model predictions ----
predictions_path <- "~/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Base maps/6.Infectious_diseases/Allen 2017/predictions.RData"
stopifnot(file.exists(path.expand(predictions_path)))

load(path.expand(predictions_path))  # loads object: predictions

eid_pred <- predictions %>%
  select(gridid, lon, lat, bsm_weight_pop) %>%
  filter(is.finite(lon), is.finite(lat), is.finite(bsm_weight_pop))

# ---- 2) build raster from the archived grid ----
r_eid <- terra::rast(
  eid_pred[, c("lon", "lat", "bsm_weight_pop")],
  type = "xyz",
  crs = "EPSG:4326"
)
names(r_eid) <- "bsm_weight_pop"

# ---- 3) standard deviation scaling used in the paper ----
# Fig. 3 caption: palette scaled to 2.5 SD above and below the mean
eid_vals <- values(r_eid, na.rm = TRUE)
eid_mean <- mean(eid_vals, na.rm = TRUE)
eid_sd <- sd(eid_vals, na.rm = TRUE)
eid_min <- eid_mean - 2.5 * eid_sd
eid_max <- eid_mean + 2.5 * eid_sd

r_eid_scaled <- clamp(r_eid, lower = eid_min, upper = eid_max, values = TRUE)

# ---- 4) convert to plotting data ----
eid_df <- as.data.frame(r_eid_scaled, xy = TRUE, na.rm = TRUE)

# ---- 5) country boundaries ----
world <- ne_countries(scale = "medium", returnclass = "sf")

# ---- 6) plot ----
# Land mask (so country borders are visible)
land_sf     <- ne_countries(scale = "medium", returnclass = "sf") 
fig_3b_plot <- ggplot() +
  geom_raster(data = eid_df, aes(x = x, y = y, fill = bsm_weight_pop)) +
  geom_sf(data =land_sf, fill = NA, color = "white", size = 0.15) +
  #  coord_fixed() +
  ylim(-65, 90) +
  scale_fill_gradientn(colours = viridis(100),
                       guide = guide_colorbar(label = TRUE,
                                              label.position = "right",
                                              title = "Event probability\n(relative to\nreporting effort)")) +
  #  theme_black_legend() +
  theme(panel.background = element_rect(fill = "black"),
        plot.background  = element_rect(fill = "black"),
        line = element_blank(),
        legend.title = element_text(color = "white", size = 8),
        legend.text = element_text(color = "white", size = 8),
        legend.title.align = 0,
        legend.background = element_blank(),
        legend.position = c(0.11, 0.45)) +
  labs(x = NULL, y = NULL)

  # geom_sf(data = world, fill = NA, color = "white", linewidth = 0.15) +
  # coord_sf(
  #   xlim = c(-180, 180),
  #   ylim = c(-65, 90),
  #   expand = FALSE
  # ) +
  # scale_fill_gradientn(
  #   colours = viridis::viridis(256, option = "D", direction = 1),
  #   limits = c(eid_min, eid_max),
  #   oob = scales::squish,
  #   guide = guide_colorbar(
  #     title = "High",
  #     label.position = "right",
  #     barheight = unit(6, "cm")
  #   )
  # ) +
  # labs(
  #   title = "Predicted relative risk distribution of zoonotic emerging infectious diseases",
  #   subtitle = "Figure 3B reproduction after factoring out reporting bias",
  #   caption = "Based on Allen et al. (2017); weighted model output reweighted by population, scaled to +/- 2.5 SD around the mean.",
  #   x = NULL,
  #   y = NULL
  # ) +
  # theme_minimal(base_size = 11) +
  # theme(
  #   panel.background = element_rect(fill = "black", color = NA),
  #   plot.background = element_rect(fill = "black", color = NA),
  #   panel.grid = element_blank(),
  #   axis.text = element_blank(),
  #   axis.ticks = element_blank(),
  #   legend.position = c(0.11, 0.45),
  #   legend.title = element_text(color = "white", size = 8),
  #   legend.text = element_text(color = "white", size = 8),
  #   legend.background = element_blank(),
  #   plot.title = element_text(color = "white", face = "bold"),
  #   plot.subtitle = element_text(color = "white"),
  #   plot.caption = element_text(color = "white")
  # )

fig_3b_plot

# ---- 7) optional save ----
# ggsave("fig_3b_infectious_diseases.png", fig_3b_plot, width = 11, height = 6.5, dpi = 300)
