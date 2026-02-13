# 2.4.Conflict after threshold--------
## Creating the binary maps
# *Assumptions*
# + "Data: starting raster object is r_sum_land.", 
# + "1st- it is clamped and all NAs are excluded.", 
# + "2nd- all 0s are excluded from df_bin and df_bin_no0 created.",
# + "3rd- new columns are created according to the percentiles (0.5, 0.75, 0.9, 0.95) and binary values are given (< threshold 0, >= threshold 1)")
# 
# We have decided to continue with 50th percentile dataset. 
# So anything below 50th percentile is 0, anything above is 1. 

# Load packages
# install.packages("pacman")
library(pacman)
p_load(tidyverse)
p_load(terra)
p_load(rnaturalearth)
p_load(sf)
p_load(maps)

## Load Data

GEDEvent_v25_1 <- GEDEvent_v25_1 <- readRDS("~/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Base maps/2.Conflict/GEDEvent_v25_1.rds")

# Manipulate Data  -----------------------------
### Parameters # -----------------------------
sum_field    <- "best"            # "best", "high", "low", etc.
cell_km      <- 50                # grid resolution (km): try 25/50/100
crs_moll     <- "ESRI:54009"      # Mollweide (equal-area)

## Data → points (WGS84 → Mollweide)# -----------------------------

dat_df <- GEDEvent_v25_1 %>%
  filter(year <= 2024, 
         year >= 2000,
         is.finite(longitude), is.finite(latitude)) %>%
  as.data.frame()

stopifnot(sum_field %in% names(dat_df))

pts_ll   <- vect(dat_df, geom = c("longitude","latitude"), crs = "EPSG:4326")
pts_moll <- project(pts_ll, crs_moll)

## Equal-area raster grid in Mollweide # -----------------------------

cell_m  <- cell_km * 1000
r_moll  <- rast(xmin = -18000000, xmax = 18000000,
                ymin =  -9000000, ymax =  9000000,
                resolution = cell_m, crs = crs_moll)

## Rasterize: SUM within each cell# -----------------------------

r_sum <- rasterize(pts_moll, r_moll, field = sum_field, fun = "sum", background = 0)

# Land mask (so oceans are blank)
land_sf     <- ne_countries(scale = "medium", returnclass = "sf") |> st_transform(crs_moll)
# create main data------
r_sum_land  <- mask(r_sum, vect(land_sf)) 

# this brings back all the 0s
bin_sum_land <- clamp(r_sum_land, values = TRUE)
# 1. All NAs excluded
df_bin <- as.data.frame(bin_sum_land, xy = TRUE, na.rm = TRUE)
# 2. All 0s excluded
df_bin_no0 <- df_bin %>%
  filter(sum != 0) # remove rows where sum == 0

# 3. create new columns with below threshold 0, above threshold 1
# thresholds are decided according to the quantiles (0.5, 0.75, 0.9, 0.95)
# make it factorial so it is binary not continuous
df_bin_no0 <- df_bin_no0 %>%             
  mutate(q50_bin = factor(ifelse(sum < quantile(sum, probs= 0.5), 0, 1)),
         q75_bin = factor(ifelse(sum < quantile(sum, probs= 0.75), 0, 1)),
         q90_bin = factor(ifelse(sum < quantile(sum, probs= 0.9), 0, 1)),
         q95_bin = factor(ifelse(sum < quantile(sum, probs= 0.95), 0, 1))
  )  
# Plot ############

p_50 <- ggplot() +
  geom_tile(data = df_bin_no0, aes(x = x, y = y, fill = q50_bin)) +
  geom_sf(data = land_sf, fill = NA, color = "lightgray", linewidth = 0.1) +
  coord_sf(crs = st_crs(crs_moll), expand = FALSE) +
  scale_fill_manual(values = c("0" = "#F8766D", "1" = "#00BFC4"),  
                    limits = c("0","1"),        # keep order
                    drop   = FALSE,             # show both even if one is absent
                    name = "Sum of Deaths",
                    labels = c("0" = "< threshold", "1" = "≥ threshold")) +
  labs(
    title = "Organized Violence (sum of estimated deaths)",
    subtitle = sprintf("2000-2024 • %dkm equal-area grid", cell_km),
    caption = "Projection: Mollweide (equal-area). Values included after 50th percentile."
  ) + 
  xlab("longitude") +
  ylab("latitude") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(size = 9, margin = margin(b = 6)),
    plot.caption    = element_text(size = 9, color = "#666666"),
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10, face = "bold")
  )
p_50
# pdf(file = "q50_bin.pdf")
# p_50
# dev.off()
# Save df
conflict_bin_df <- df_bin_no0 %>% select(c("x", "y", "q50_bin"))

# rm(list=setdiff(ls(), c("conflict_bin_df")))


# 4th of February
# MAKE NAs 0 also ######

### this is only for visualisation dataset hasn't changed
df_bin_NAsconverted0 <- df_bin_no0 %>% select(c(x,y, q50_bin)) %>% filter(q50_bin==1)

conflict_binary_single <- ggplot() +
  geom_tile(data = df_bin_NAsconverted0, aes(x = x, y = y, fill = q50_bin)) +
  geom_sf(data = land_sf, fill = NA, color = "lightgray", linewidth = 0.1) +
  coord_sf(crs = st_crs(crs_moll), expand = FALSE) +
  scale_fill_manual(values = c("1" = "#00BFC4"),  
                    limits = c("1"),        # keep order
                    drop   = FALSE,             # show both even if one is absent
                    name = "Sum of Deaths",
                    labels = c("1" = "≥ threshold")) +
  labs(
    title = "Organized Violence (sum of estimated deaths)",
    subtitle = sprintf("2000-2024 • %dkm equal-area grid", cell_km),
    caption = "Projection: Mollweide (equal-area). Values included after 50th percentile.NAs treated as 0"
  ) + 
  xlab("longitude") +
  ylab("latitude") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(size = 9, margin = margin(b = 6)),
    plot.caption    = element_text(size = 9, color = "#666666"),
    legend.position = "right",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10, face = "bold")
  )

# pdf(file = "q50_bin_NAsconverted0.pdf")
# conflict_binary_single
# dev.off()
