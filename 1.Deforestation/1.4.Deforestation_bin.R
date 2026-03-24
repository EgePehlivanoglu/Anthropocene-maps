
# install librarian package
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian", repos = "https://cloud.r-project.org")
}
library(librarian)
shelf(httr, terra, tmap, tidyverse, rnaturalearth, sf) # install and load packages

# 1) Download the public GeoTIFF (v1.2, 2001–2024) from Zenodo #######
url <- "https://zenodo.org/records/15366671/files/drivers_forest_loss_1km_2001_2024_v1_2.tif?download=1"  # Zenodo direct
tif_path <- tempfile(fileext = ".tif")
res <- GET(url, write_disk(tif_path, overwrite = TRUE))
stop_for_status(res)

# 2) Read the raster ########
# Band 1 = dominant class (uint8, 1..7); Bands 2..8 = class probabilities (0..250)
r_all <- rast(tif_path)
r <- r_all[[1]]          # dominant class only
r <- clamp(r, 1, 7,  values = FALSE)      # be safe: anything outside 1..7 to NA
# 3) Reproject to Mollweide (categorical -> nearest neighbor)
#r_moll <- project(r, "ESRI:54009", method = "near")
#use new projection since many points failed to transform (specify the resolution)
moll <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
r_moll <- terra::project(r, moll, method="near", res=10000)

# metadata: -1 = non-ag land, -2 = water/no data  -> set to NA
r_moll[r_moll %in% c(-1, -2)] <- NA              # keep raw counts for everything else
names(r_moll) <- "value"
# 3) Reproject ######## 
# target CRS: Mollweide, Using a standard proj string for world Mollweide
crs_moll <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"

# 4) Add country boundaries -> sf -> Mollweide ######
world <- ne_countries(scale = "medium", returnclass = "sf")
world_moll <- st_transform(world, crs = crs_moll)
gc() #important so the memory limit is OK
# 5) Create data frame for ggplot ----
# Attach human-readable labels from the dataset docs (v1.2)
# Codes: 1 Perm. agriculture, 2 Hard commodities, 3 Shifting cultivation,
# 4 Logging, 5 Wildfire, 6 Settlements & infrastructure, 7 Other natural disturbances
lev <- data.frame(
  ID   = 1:7,
  name = c("Permanent agriculture",
           "Hard commodities",
           "Shifting cultivation",
           "Logging",
           "Wildfire",
           "Settlements & infrastructure",
           "Other natural disturbances"),
  binary_codes =c(1,1,1,0,0,0,0)
)

# 6) Binary Transformation #######
# Reclassify the dominant-class raster directly to a binary raster:
# 1 = permanent disturbances, 0 = temporary disturbances
r_bin <- classify(
  r_moll,
  rcl = cbind(lev$ID, lev$binary_codes)
)
names(r_bin) <- "binary_codes"

# Convert the binary raster to a plotting data frame only after classification
df_bin <- as.data.frame(r_bin, xy = TRUE, na.rm = TRUE) %>%
  mutate(binary_codes = factor(binary_codes, levels = c(0, 1)))

deforestation_bin <- r_bin
world_moll_deforestation <- world_moll
# 7) Plotting with ggplot2 ######
bin_plot_hr <- ggplot() +
  geom_tile(data = df_bin, aes(x = x, y = y, fill = binary_codes)) +
  geom_sf(data = world_moll, fill = NA, color = "black", linewidth = 0.1) +
  coord_sf() +
  scale_fill_manual(
    name = "Disturbance Types",
    values = c("0" = "seagreen4",  # low risk (example yellow)
               "1" = "salmon4",  # high risk (example red)
               "NA" = "white"),  # explicitly map NA to white
    na.value = "white",           # ensures NAs are white
    labels = c("0" = "Temporary", "1" = "Permanent")
  ) +
  labs(
    title = "Deforestation Plot",
    subtitle = "Binary Transformed Data (10000 m)",
    caption = "Permanent Disturbances:Permanent agriculture, Hard commodities, Shifting cultivation. Temporary Disturbances: Logging, Wildfire, Settlements & infrastructure, and Other natural disturbances",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    # panel.background = element_rect(fill = "gray85", color = NA), # gray background (ocean)
    # plot.background = element_rect(fill = "gray85", color = NA),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.6, "cm"),
    legend.text = element_text(size = 9)
  )
bin_plot_hr
