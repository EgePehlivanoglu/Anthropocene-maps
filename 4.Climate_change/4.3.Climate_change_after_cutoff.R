#Climate Change with bigger than 2 degreee cutoff


library(pacman)
p_load(terra)
p_load(R.utils)
p_load(zoo)
p_load(tidyverse)
p_load(scales)
p_load(rnaturalearth)
p_load(rnaturalearthdata)



# ---- 1) get GISTEMP monthly LOTI grid ----
url <- "https://data.giss.nasa.gov/pub/gistemp/gistemp1200_GHCNv4_ERSSTv5.nc.gz"
gzf <- file.path(tempdir(), basename(url))
download.file(url, gzf, mode = "wb")
ncfile <- R.utils::gunzip(gzf, overwrite = TRUE, remove = FALSE)

# ---- 2) load monthly raster ----
r_mon <- rast(ncfile)   # monthly 2°×2° anomaly grid (°C)
n <- nlyr(r_mon)

# ---- 3) build dates WITHOUT using time<- (GISTEMP monthly starts 1880-01) ----
# If you want to be extra safe, you can check the first valid year in the dataset,
# but for GISTEMP LOTI monthly grids 1880-01 is correct.
dates <- seq(as.Date("1880-01-01"), by = "1 month", length.out = n)
years <- format(dates, "%Y")

# ---- 4) annual means per cell (group by year labels) ----
r_ann <- tapp(r_mon, years, mean, na.rm = TRUE)
names(r_ann) <- sort(unique(years))

# ---- 5) 10-year running means ----
k <- 10
idx <- embed(1:nlyr(r_ann), k)[, k:1]
r_10yr <- rast(lapply(seq_len(nrow(idx)), function(i) {
  mean(r_ann[[idx[i,]]], na.rm = TRUE)
}))
names(r_10yr) <- apply(
  cbind(names(r_ann)[idx[,1]], names(r_ann)[idx[,k]]), 1,
  function(x) paste0(x[1], "-", x[2])
)

# Different Cut offs
## Exceeding 2 -keep this one

# Reproject the raster to Mollweide projection
mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
r_moll <- project(r_10yr[[nlyr(r_10yr)]], mollweide_crs)

# r_10yr[[nlyr(r_10yr)]] 
# note that it is updated now it spans 2017-2026 and hence the middle of russia is blank.

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world)
world_moll <- project(world_vect, mollweide_crs)


# 6) convert spatraster to a dataframe ---------

# Create a mask for values greater than 2
r_anomaly <- clamp(r_moll, lower = 2, value = FALSE)  # Sets values < 2 to NA
# The value column name is the layer name: "2017-2026"
val_col <- names(r_anomaly)[1]
df <- as.data.frame(r_anomaly, xy = TRUE, na.rm = FALSE)  # columns: x, y, value #### important that na.rm is false so it keeps NAs too.
df_bin <- df %>%
  mutate(
    bin = if_else(!is.na(.data[[val_col]]), 1L, 0L, missing = NA)
  )
df_bin$bin <- as.factor(df_bin$bin)

# ---- country boundaries -> sf -> Mollweide ----
world <- ne_countries(scale = "medium", returnclass = "sf")
world_moll <- st_transform(world, crs = mollweide_crs)

# 7) plot----------
bin_climate_change_plot <- ggplot() +
  geom_tile(data = df_bin, aes(x = x, y = y, fill = factor(bin))) +
  # optional if you already have it:
  geom_sf(data = world_moll, fill = NA, color = "black", linewidth = 0.1) +
  coord_sf(crs = mollweide_crs, expand = TRUE) +
  scale_fill_manual(values = c(`0` = "white", `1` = "red"), name = "Data Present") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  labs(
    title = "Areas with Anomaly > 2 Degree Celcius (2017-2026)",
    subtitle = "Binary Transformed Data with Mollweide Projection",
    caption = "Years for 2017-2026, if the increase is more than 2 degree Celcius it is shown in red",
    x = NULL, y = NULL
  ) 

# save the plot in pdf
# pdf(file = "binary climate change plot2017-26.pdf")
# bin_climate_change_plot
# dev.off()

