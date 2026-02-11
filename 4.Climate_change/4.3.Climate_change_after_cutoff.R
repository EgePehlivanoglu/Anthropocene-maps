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

# Create a mask for values greater than 2
r_anomaly <- clamp(r_moll, lower = 2, value = FALSE)  # Sets values < 2 to NA

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world)
world_moll <- project(world_vect, mollweide_crs)
# Create and reverse the "Heat 2" color palette
heat2_pal <- hcl.colors(100, "Heat 2")
reversed_heat2 <- rev(heat2_pal)  # Reverse the color order
#pdf(file = "Areas_withanomalybiggerthan2.pdf")
# Plot only the areas with anomaly > 2
plot(r_anomaly, 
     col = reversed_heat2,  # Using red colors for emphasis
     main = "Areas with Anomaly > 2 (2016-2025)", 
     plg = list(title = "Anomaly Value"),
     axes = FALSE)

# Add land borders for reference
plot(world_moll, add = TRUE, col = NA, border = "black", lwd = 0.1)
#dev.off()
