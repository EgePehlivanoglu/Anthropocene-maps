#Climate Change with bigger than 2 degreee cutoff

# Workflow
# 1. Get GISTEMP monthly LOTI grid (1880-2026)
# 2. Load monthly raster
# 3. Build dates WITHOUT using time<- (GISTEMP monthly starts 1880  -01)  
# 4. Annual means per cell (group by year labels) 
# 5. 10-year running means
# 6. Convert spatraster to a dataframe
# 7. Plot

# install librarian package (working better than pacman)
# From CRAN:
install.packages("librarian")
librarian::shelf(R.utils, terra, zoo, tidyverse, scales, sf, rnaturalearth, rnaturalearthdata) # install and load packages


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
r_moll <- project(r_10yr[[nlyr(r_10yr)-1]], mollweide_crs)

# r_10yr[[nlyr(r_10yr)]] is the last layer, which is 2017-2026
# note that it is updated now it spans 2016-2025 (first version)

# 6) create binary spatraster ---------

# Create a mask for values greater than 2
r_anomaly <- clamp(r_moll, lower = 2, value = FALSE)  # Sets values < 2 to NA

# Create a binary SpatRaster:
# 1 = anomaly > 2C, 0 = anomaly <= 2C or NA in the masked raster
r_anomaly_bin <- ifel(is.na(r_anomaly), 0, 1)
names(r_anomaly_bin) <- "bin"

# Convert the binary raster to a data frame only for plotting
df_bin <- as.data.frame(r_anomaly_bin, xy = TRUE, na.rm = FALSE)
df_bin$bin <- factor(df_bin$bin, levels = c(0, 1))

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
    title = "Areas with Anomaly > 2 Degree Celcius (2016-2025)",
    subtitle = "Binary Transformed Data with Mollweide Projection",
    caption = "Years for 2016-2025, if the increase is more than 2 degree Celcius it is shown in red",
    x = NULL, y = NULL
  ) 

# # save the plot in pdf
# pdf(file = "binary climate change plot2016-25.pdf")
# bin_climate_change_plot
# dev.off()

# rename the data for next steps
clim_change_binSR <- r_anomaly_bin
world_moll_clim <- world_moll