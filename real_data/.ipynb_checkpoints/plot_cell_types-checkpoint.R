.libPaths(c("~/R/library", .libPaths()))
library(tidyverse)
library(glue)
library(SingleCellExperiment)
library(lemur)
source("~/conformeR-paper/util.R")

# Load data.
fit10 <- readRDS("~/fitseed10_standard.rds")
reducedDim(fit10, "fit_al_umap")[,1] <- -reducedDim(fit10, "fit_al_umap")[,1]

cell_type_colors <- structure(c("#FAC1BD", "#8dd3c7", "#80b1d3", "#b3de69"), names = c("Tumor cells", "T cells", "Oligodendrocytes", "Myeloid cells"))
umap_fit <- reducedDim(fit10, "fit_al_umap") |> as.data.frame()

ct <- colData(fit10) |> as.data.frame() |> dplyr::select(cell_type) |> dplyr::rename("Cell-type"="cell_type") 

ct_plot_data <- ct |>
  rownames_to_column("cell") |>
  left_join(
    umap_fit |> rownames_to_column("cell"),
    by = "cell"
  )

ggplot(ct_plot_data, aes(x = V1, y = V2)) +
  ggrastr::rasterise(
    geom_point(
      aes(color = `Cell-type`),
      alpha = 0.2,
      size = 0.8,
      stroke = 0
    )
  ) +
  scale_x_continuous(limits = range(umap_fit[, 1]) * 1.1) +
  scale_y_continuous(limits = range(umap_fit[, 2]) * 1.1) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed() +
  facet_wrap(~`Cell-type`) +
  small_axis(
    fontsize = 18,
    xlim = range(umap_fit[, 1]),
    ylim = range(umap_fit[, 2]),
    arrow_length = 6
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_void(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    strip.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )