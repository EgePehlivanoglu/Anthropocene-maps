# Workflow
# 1. Load data (tiff file)
# 2. Reproject to Mollweide
# 3. Binary transformation (low risk score vs high risk score)


# Load packages
# install librarian package (working better than pacman)
# From CRAN:
install.packages("librarian")
librarian::shelf(raster,       # for raster data handling
                 sf,           # for vector data (shapefiles)
                 rnaturalearth,# to get country boundaries
                 tidyverse,      # for plotting
              #  viridis, # for nice color scales (optional)
                 terra) # for rast function

# Binary Global pesticide risk scores ########
# *Assumptions*
#   + Value for non-agricultural land: -1, Value for water/no data: -2 given NA
# + Give NA to non-agricultural land (previously 0 on the original figure)
# + Give 0 for values low to medium (until 3)
# + Give 1 for values high (3 or above)

# ---- 1. load data ----
tif_path <- "~/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Base maps/5.Pesticide_use/Tang 2021/data/Global_pesticide_risk_scores.tif"  # <- your file
# ---- read raster in WGS84 and clean special values ----
r <- rast(tif_path)                     # should be WGS84 per metadata
# metadata: -1 = non-ag land, -2 = water/no data  -> set to NA
r[r %in% c(-1, -2)] <- NA              # keep raw counts for everything else
names(r) <- "value"

# ----- investigate the dataset -----
# r is your SpatRaster
# Resolution (pixel size in map units)
res(r)
# CRS / projection (WKT string)
crs(r)
# A friendlier summary: EPSG if recognized + proj string
terra::crs(r, describe=TRUE)
# Extent (bounding box) in CRS units
ext(r)
# Dimensions: nrows, ncols, nlyr
dim(r)
# Number of cells
ncell(r)

# ---- 2. reproject CRS: Mollweide ----
# Using a standard proj string for world Mollweide
crs_moll <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"

# ---- reproject raster (use nearest neighbor to preserve raw counts) ----
r_moll <- project(r, crs_moll, method = "near")

# ----- investigate the reporjected dataset -----
# r_moll is your SpatRaster
r_moll
# Resolution (pixel size in map units)
res(r_moll)
# CRS / projection (WKT string)
crs(r_moll)
# A friendlier summary: EPSG if recognized + proj string
terra::crs(r_moll, describe=TRUE)
# Extent (bounding box) in CRS units
ext(r_moll)
# Dimensions: nrows, ncols, nlyr
dim(r_moll)
# Number of cells
ncell(r_moll)


# ---- add country boundaries -> sf -> Mollweide ----
world <- ne_countries(scale = "medium", returnclass = "sf")
world_moll <- st_transform(world, crs = crs_moll)

# ---- 3. binary transformation ----
df <- as.data.frame(r_moll, xy = TRUE, na.rm = TRUE)  # columns: x, y, value

df_bin <- df %>% 
  mutate(RS_bin= case_when(value == 0 ~ NA,
                           value >0 & value <3  ~ "0",  
                           value >=3 ~ "1"), factor(RS_bin))

# ----- 4. plotting ----
bin_plot <- ggplot() +
  geom_tile(data = df_bin, aes(x = x, y = y, fill = RS_bin)) +
  geom_sf(data = world_moll, fill = NA, color = "black", linewidth = 0.2) +
  coord_sf() +
  scale_fill_manual(
    name = "Binary Risk Score",
    values = c("0" = "#f4d03f",  # low risk (example yellow)
               "1" = "#c0392b",  # high risk (example red)
               "NA" = "white"),  # explicitly map NA to white
    na.value = "white",           # ensures NAs are white
    labels = c("0" = "<3", "1" = "≥3")
  ) +
  labs(
    title = "Global Pesticide Risk Score",
    subtitle = "Figure 1 from Tang et al. (2021) reproduced with Mollweide projection and binary transformation",
    caption = "Non-agricultural lands and places with no data shown in white.",
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


# ---- save (optional) ----
# ggsave("binary_rs_mollweide.png", bin_plot, width = 11, height = 6.5, dpi = 300)



