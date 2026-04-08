# HISTOGRAMS AND THRESHOLD MAPS
## for fig 3B reproduction
setwd(dir = "/Users/egepehlivanoglu/Library/CloudStorage/OneDrive-StockholmUniversity/KVA backup/KVAOneDrive_backup_28Jan2026/2. Projects/4.Cascades/Editorial Review 20250807/map/Anthropocene-maps/6.Infectious_diseases")
source("6.2.Infectious_diseases.R")

if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny", repos = "https://cloud.r-project.org")
}
library(shiny)

# ---- 1) histogram of the Fig. 3B raster values ----
# Exclude 0s and NAs before calculating upper-tail thresholds
hist_df <- fig3b_df %>%
  dplyr::filter(!is.na(risk) & risk > 0)

cutoff_5 <- stats::quantile(hist_df$risk, probs = 0.95, na.rm = TRUE)
cutoff_10 <- stats::quantile(hist_df$risk, probs = 0.90, na.rm = TRUE)
cutoff_25 <- stats::quantile(hist_df$risk, probs = 0.75, na.rm = TRUE)

hist(
  hist_df$risk,
  xlab = "EID Risk Index",
  main = "Histogram of EID Risk Index",
  sub = "Thresholds calculated after excluding 0s and NAs",
  col = "steelblue",
  border = "white"
)

abline(v = cutoff_5, col = "#b30000", lwd = 2, lty = 1)
abline(v = cutoff_10, col = "#e34a33", lwd = 2, lty = 2)
abline(v = cutoff_25, col = "#fdbb84", lwd = 2, lty = 3)

legend(c("bottomright"),
       legend = c("Top 5%", "Top 10%", "Top 25%"),
       col = c("#b30000", "#e34a33", "#fdbb84"),
       lty = c(1, 2, 3),
       lwd = 2,
       bty = "n"
)

# ---- 2) helper for upper-tail thresholds ----
top_threshold_data <- function(df, percent = 10) {
  filtered_input <- df %>%
    dplyr::filter(!is.na(risk) & risk > 0)

  cutoff <- stats::quantile(filtered_input$risk, probs = 1 - percent / 100, na.rm = TRUE)
  filtered <- filtered_input %>%
    dplyr::filter(risk >= cutoff)

  list(
    data = filtered,
    cutoff = unname(cutoff)
  )
}

# ---- 3) top 10% data from the histogram distribution ----
top10_result <- top_threshold_data(hist_df, percent = 10)
top10_fig3b_df <- top10_result$data
top10_cutoff <- top10_result$cutoff

# ---- 4) map using the same plotting approach as fig_3b_plot ----
top10_fig3b_plot <- ggplot() +
  geom_polygon(
    aes(x = long, y = lat, group = group),
    data = world_df,
    inherit.aes = FALSE,
    fill = viridis::viridis(1)
  ) +
  geom_raster(
    aes(x = x, y = y, fill = risk),
    data = top10_fig3b_df
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
    colours = viridis::viridis(max(2, length(unique(top10_fig3b_df$risk)))),
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
  labs(
    title = "Top 10% of EID Risk Index",
    subtitle = paste0(
      "Cells at or above ",
      round(top10_cutoff, 6),
      " after excluding 0s and NAs"
    ),
    caption = "Threshold calculated from non-zero, non-missing cells only.",
    x = NULL,
    y = NULL
  )

top10_fig3b_plot

# ---- 5) interactive threshold map ----
threshold_map_plot <- function(percent = 10) {
  threshold_result <- top_threshold_data(hist_df, percent = percent)
  threshold_df <- threshold_result$data
  threshold_cutoff <- threshold_result$cutoff

  ggplot() +
    geom_polygon(
      aes(x = long, y = lat, group = group),
      data = world_df,
      inherit.aes = FALSE,
      fill = viridis::viridis(1)
    ) +
    geom_raster(
      aes(x = x, y = y, fill = risk),
      data = threshold_df
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
      colours = viridis::viridis(max(2, length(unique(threshold_df$risk)))),
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
    labs(
      title = paste0("Top ", percent, "% of EID Risk Index"),
      subtitle = paste0(
        "Cutoff value: ",
        round(threshold_cutoff, 6),
        " after excluding 0s and NAs"
      ),
      caption = "Interactive threshold map based on non-zero, non-missing cells only.",
      x = NULL,
      y = NULL
    )
}

threshold_app <- shinyApp(
  ui = fluidPage(
    titlePanel("Interactive threshold map for infectious disease risk"),
    p("Thresholds are calculated after excluding 0s and NAs from the risk distribution."),
    sidebarLayout(
      sidebarPanel(
        sliderInput(
          inputId = "threshold_percent",
          label = "Top percentage to retain",
          min = 1,
          max = 100,
          value = 10,
          step = 1,
          post = "%"
        )
      ),
      mainPanel(
        plotOutput("threshold_map", height = "700px")
      )
    )
  ),
  server = function(input, output, session) {
    output$threshold_map <- renderPlot({
      threshold_map_plot(input$threshold_percent)
    })
  }
)

# Run interactively when desired:
shiny::runApp(threshold_app)
