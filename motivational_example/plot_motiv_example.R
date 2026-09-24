.libPaths(c("~/R/library", .libPaths()))
library(tidyverse)
library(glue)
library(SingleCellExperiment)
library(lemur)
source("~/conformeR-paper/util.R")

# Load data.
fit10 <- readRDS("~/fitseed10_standard.rds")
fit11 <- readRDS("~/fitseed11_standard.rds")
nei10 <- readRDS("~/nei_standard10.rds")
nei11 <- readRDS("~/nei_standard11.rds")
pred_set_10 <- readRDS("~/pred_set_pan10.rds")
pred_set_11 <- readRDS("~/pred_set_pan11.rds")

reducedDim(fit10, "fit_al_umap")[,1] <- -reducedDim(fit10, "fit_al_umap")[,1]
reducedDim(fit11, "fit_al_umap")[,1] <- -reducedDim(fit11, "fit_al_umap")[,1]

genes_of_interest <- rowData(fit10) %>%
  as_tibble() %>%
  dplyr::select(gid, gene) %>%
  filter(gene %in% c("HBEGF","MBP"))

plot_umap_simple <- function(fit, nei, pred_set) {
  
  fit_small <- fit[genes_of_interest$gid, ]
  umap_fit <- reducedDim(fit_small, "fit_al_umap") |> as.data.frame()
  
  nei <- nei |>
    dplyr::filter(name %in% genes_of_interest$gid) |>
    mutate(
      neighborhood = map(neighborhood, ~ colnames(fit) %in% .x)
    ) |>
    dplyr::select(name, neighborhood) |>
    unnest(neighborhood) |>
    mutate(cell = rep(colnames(fit), nrow(genes_of_interest))) |>
    left_join(genes_of_interest, by = c("name" = "gid"))
  
  selection <- pred_set |>
    dplyr::filter(gene %in% genes_of_interest$gid) |>
    left_join(nei, by = c("cell", "gene" = "name"))
  
  de_plot_data <- as_tibble(
    colData(fit_small),
    rownames = "cell"
  ) |>
    mutate(
      umap = umap_fit,
      de = as_tibble(t(assay(fit_small, "DE_panobinostat")))
    ) |>
    unnest(de, names_sep = "-") |>
    pivot_longer(
      starts_with("de-"),
      names_sep = "-",
      values_to = "de",
      names_to = c(NA, "gid")
    ) |>
    inner_join(
      as_tibble(rowData(fit_small)),
      by = "gid"
    ) |>
    mutate(gene = factor(gene)) |>
    left_join(
      selection,
      by = c("gid" = "gene", "cell")
    )
  
  de_plot_data <- list(
    lemur = de_plot_data |>
      mutate(
        inside = factor(
          ifelse(neighborhood, "in", "out"),
          levels = c("in", "out")
        )
      ),
    cs = de_plot_data |>
      mutate(
        inside = factor(
          ifelse(inside_cs, "in", "out"),
          levels = c("in", "out")
        )
      )
  )
  
  plots <- lapply(de_plot_data, \(data) {
    
    abs_max <- max(abs(quantile(data$de, c(0.95, 0.05))))
    
    ggplot(data, aes(x = umap[, 1], y = umap[, 2])) +
      ggrastr::rasterise(
        geom_point(
          aes(
            color = de,
            alpha = scales::rescale(abs(de), to = c(0.05, 1))
          ),
          size = 0.5
        )
      ) +
      scale_alpha_identity() +
      scale_color_de_gradient(
        abs_max,
        mid_width = 0.2
      ) +
      facet_grid(
        rows = vars(inside),
        cols = vars(gene),
        labeller = labeller(
          inside = c(
            "in" = "Inside",
            "out" = "Outside"
          ),
          gene = as_labeller(\(x) paste0(x))
        ),
        switch = "y"
      ) +
      scale_x_continuous(
        limits = range(umap_fit[, 1]) * 1.1
      ) +
      scale_y_continuous(
        limits = range(umap_fit[, 2]) * 1.1,
        position = "left"
      ) +
      small_axis(
        fontsize = 18,
        xlim = range(umap_fit[, 1]),
        ylim = range(umap_fit[, 2]),
        arrow_length = 6
      ) +
      guides(
        color = guide_colorbar(
          barheight = unit(4, "mm"),
          barwidth = unit(30, "mm"),
          title = ""
        )
      ) +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(fill = "white", colour = NA),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(
          size = 20,
          angle = 90
        ),
        strip.placement = "outside",
        axis.text = element_text(size = 18),
        panel.spacing.x = unit(12, "mm"),
        panel.spacing.y = unit(4, "mm"),
        legend.text = element_text(size = 18),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(3, "mm"),
        legend.position = "bottom",
        plot.margin = margin(2, 2, 2, 2)
      )
  })
  
  list(
    plot_lemur = plots$lemur,
    plot_cs = plots$cs
  )
}

plot_seed10 <- plot_umap_simple(fit10, nei10, pred_set_10) 
plot_seed11 <- plot_umap_simple(fit11, nei11, pred_set_11) 