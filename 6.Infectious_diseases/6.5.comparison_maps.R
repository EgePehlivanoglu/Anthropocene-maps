# Publication-style comparison maps:
# 1) side-by-side comparison
# 2) single-map overlay

setwd("/Users/egepehlivanoglu/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Anthropocene-maps")

if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian", repos = "https://cloud.r-project.org")
}

library(librarian)
shelf(terra, sf, ggplot2, dplyr, rnaturalearth, viridis, knitr)

# ---- 1) source the summed binary raster workflow ----
load_summed_binary_context <- function() {
  required_names <- c("binary_sum", "target_crs", "target_res")
  context_rds <- "summed_binary_context.rds"

  if (file.exists(context_rds)) {
    sum_context <- readRDS(context_rds)
    sum_env <- new.env(parent = globalenv())
    list2env(sum_context, envir = sum_env)
    return(sum_env)
  }

  if (all(vapply(required_names, exists, logical(1), envir = .GlobalEnv, inherits = FALSE))) {
    return(.GlobalEnv)
  }

  sum_env <- new.env(parent = globalenv())
  sum_r_script <- tempfile(fileext = ".R")
  knitr::purl(
    input = "raster_math_binary_sum.Rmd",
    output = sum_r_script,
    quiet = TRUE
  )

  tryCatch(
    sys.source(sum_r_script, envir = sum_env),
    error = function(e) {
      stop(
        paste(
          "Could not load summed binary raster context automatically.",
          "Prefer running raster_math_binary_sum.Rmd first to create summed_binary_context.rds.",
          "If the fallback notebook sourcing is used, its data sources must be reachable.",
          "Original error:", conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  missing_required <- required_names[!vapply(required_names, exists, logical(1), envir = sum_env, inherits = FALSE)]
  if (length(missing_required) > 0) {
    stop(
      paste(
        "The sourced raster_math_binary_sum.Rmd did not create required objects:",
        paste(missing_required, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  sum_env
}

valid_spatraster <- function(x) {
  inherits(x, "SpatRaster") &&
    !inherits(try(terra::nlyr(x), silent = TRUE), "try-error")
}

rebuild_extent_from_df <- function(df, res_xy) {
  terra::ext(
    min(df$x, na.rm = TRUE) - res_xy[1] / 2,
    max(df$x, na.rm = TRUE) + res_xy[1] / 2,
    min(df$y, na.rm = TRUE) - res_xy[2] / 2,
    max(df$y, na.rm = TRUE) + res_xy[2] / 2
  )
}

sum_env <- load_summed_binary_context()

binary_sum <- if (exists("binary_sum", envir = sum_env, inherits = FALSE)) {
  get("binary_sum", envir = sum_env)
} else {
  NULL
}
target_crs <- if (exists("target_crs", envir = sum_env, inherits = FALSE)) {
  get("target_crs", envir = sum_env)
} else {
  terra::crs(binary_sum)
}
target_res <- if (exists("target_res", envir = sum_env, inherits = FALSE)) {
  get("target_res", envir = sum_env)
} else {
  terra::res(binary_sum)
}

if (exists("binary_sum_df", envir = sum_env, inherits = FALSE)) {
  binary_sum_df <- get("binary_sum_df", envir = sum_env)
} else {
  binary_sum_df <- as.data.frame(binary_sum, xy = TRUE, na.rm = FALSE)
}

if (!valid_spatraster(binary_sum) && exists("binary_sum_df", inherits = FALSE)) {
  binary_sum <- terra::rast(
    stats::na.omit(binary_sum_df[, c("x", "y", "sum_binary")]),
    type = "xyz",
    crs = target_crs
  )
  names(binary_sum) <- "sum_binary"
}

target_extent <- if (exists("target_extent", envir = sum_env, inherits = FALSE)) {
  do.call(terra::ext, as.list(get("target_extent", envir = sum_env)))
} else if (valid_spatraster(binary_sum)) {
  terra::ext(binary_sum)
} else {
  rebuild_extent_from_df(binary_sum_df, target_res)
}

common_template_raster <- terra::rast(
  ext = target_extent,
  resolution = target_res,
  crs = target_crs
)

if (exists("world_moll", envir = sum_env, inherits = FALSE)) {
  world_moll <- get("world_moll", envir = sum_env)
} else {
  world_moll <- sf::st_transform(
    rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"),
    crs = target_crs
  )
}

binary_sum_df <- binary_sum_df %>%
  mutate(sum_binary_cat = factor(sum_binary, levels = c(0, 1, 2, 3, 4)))

# ---- 2) source the infectious disease threshold workflow ----
setwd("/Users/egepehlivanoglu/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Anthropocene-maps/6.Infectious_diseases")
source("6.2.Infectious_diseases.R")

hist_df <- fig3b_df %>%
  dplyr::filter(!is.na(risk) & risk > 0)

top5_cutoff <- stats::quantile(hist_df$risk, probs = 0.95, na.rm = TRUE)
top5_fig3b_df <- hist_df %>%
  dplyr::filter(risk >= top5_cutoff)

# Reproject infectious threshold map to Mollweide so both figures match
top5_vect <- terra::vect(top5_fig3b_df, geom = c("x", "y"), crs = "EPSG:4326")
top5_vect_moll <- terra::project(top5_vect, target_crs)
top5_raster <- terra::rasterize(top5_vect_moll, common_template_raster, field = "risk", fun = "max")
names(top5_raster) <- "risk"
top5_df_moll <- as.data.frame(top5_raster, xy = TRUE, na.rm = TRUE)

# Binary presence layer for overlay
top5_presence <- terra::ifel(is.na(top5_raster), NA, 1)
names(top5_presence) <- "presence"
top5_presence_df <- as.data.frame(top5_presence, xy = TRUE, na.rm = TRUE)
top5_presence_zero <- terra::ifel(is.na(top5_presence), 0, top5_presence)
names(top5_presence_zero) <- "presence"

# ---- 3) publication-style single-map overlay ----
overlay_plot <- ggplot() +
  geom_tile(
    data = binary_sum_df,
    aes(x = x, y = y, fill = sum_binary_cat)
  ) +
  geom_tile(
    data = top5_presence_df,
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
    subtitle = "Black overlay marks cells in the top 5% of infectious disease risk",
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
overlay_plot

