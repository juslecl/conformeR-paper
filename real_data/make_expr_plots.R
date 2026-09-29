library(tidyverse)
library(tibble)
library(SingleCellExperiment)
library(patchwork)
source("~/conformeR-paper/util.R")
## --- Put outputs together --------------------------
outputs <- lapply(1:651, function(i) {
  f <- file.path("~", paste0("results", i, ".rds"))
  if (!file.exists(f)) return(NULL)
  readRDS(f) 
}) |>
  dplyr::bind_rows()

#saveRDS(outputs, "res_pb.rds")
condition_colors <- structure(c("#fc8d62", "#8da0cb"), names = c("ctrl", "panobinostat"))

#res_pb <- readRDS("~/res_pb.rds") |> dplyr::rename("name"="gene", "neighborhood"="cs_list")
res_pb <- outputs |> dplyr::rename("name"="gene", "neighborhood"="cs_list")
nei <- readRDS("~/nei.rds")
fit <- readRDS("~/fit.rds")

pan_genes_order <- c("HIST3H2A", "MBP","SPATA13")
sel_patient <- c("PW032","PW030")
pan_genes_ids <- deframe(rowData(fit)[, c("gene", "gid")])[c(pan_genes_order)]
fit_small <- fit[pan_genes_ids,fit$colData$patient_id %in% sel_patient]
full_row_data <- rowData(fit_small)

## --- Plot function --------------------------

make_pb_plot_same_axis <- function(fit_small, set, y_max = NULL, return_data = FALSE){
  
  mask <- matrix(
    NA,
    nrow = length(pan_genes_ids),
    ncol = ncol(fit_small), 
    dimnames = list(pan_genes_ids, colnames(fit_small))
  )
  
  mask2 <- matrix(
    1,
    nrow = length(pan_genes_ids),
    ncol = ncol(fit_small), 
    dimnames = list(pan_genes_ids, colnames(fit_small))
  )
  
  for(id in pan_genes_ids){
    inside_cells <- intersect(
      colnames(fit_small),
      filter(set, name == id)$neighborhood[[1]]
    )
    
    mask[id, inside_cells] <- 1
    mask2[id, inside_cells] <- NA
  }
  
  psce <- glmGamPoi::pseudobulk(
    SingleCellExperiment(
      list(
        inside = as.matrix(
          logcounts(fit_small)[pan_genes_ids,] * mask
        ),
        outside = as.matrix(
          logcounts(fit_small)[pan_genes_ids,] * mask2
        )
      )
    ),
    group_by = vars(patient_id, condition, pat_cond),
    n = n(),
    aggregation_functions = list(
      .default = \(...) matrixStats::rowMeans2(..., na.rm = TRUE)
    ),
    col_data = as.data.frame(colData(fit_small))
  )
  
  comparison_data <- as_tibble(colData(psce)) %>%
    mutate(expr_in = as_tibble(t(assay(psce, "inside")))) %>%
    mutate(expr_out = as_tibble(t(assay(psce, "outside")))) %>%
    unpack(starts_with("expr"), names_sep = "-") %>%
    pivot_longer(
      starts_with("expr"),
      names_sep = "[-_]",
      names_to = c(".value", "inside", "gid")
    ) %>%
    left_join(as_tibble(full_row_data)) %>%
    filter(condition != "etoposide") %>%
    mutate(
      condition_short = case_when(
        condition == "ctrl" ~ "ctrl",
        condition == "panobinostat" ~ "pano."
      ),
      inside = factor(inside, levels = c("in", "out")),
      gene = factor(
        gene,
        levels = c(pan_genes_order[1], "LMO2", pan_genes_order[-1])
      )
    )
  
  if (exists("return_data") && return_data) {
    return(comparison_data)
  }
  
  comparison_plots <- comparison_data %>%
    group_by(gene) %>%
    group_map(\(data, key){
      
      data$gene <- key[[1]][1]
      
      # Get the y-axis maximum for this gene
      gene_y_max <- y_max[[as.character(key[[1]])]]
      
      ggplot(data, aes(x = condition_short, y = expr)) +
        geom_point(
          aes(
            group = condition,
            color = condition
          ),
          position = position_dodge(width = 0.7),
          size = 1.5,
          show.legend = FALSE
        ) +
        geom_line(
          aes(
            group = patient_id,
            x = condition_short,
            y = expr
          ),
          color = "lightgrey",
          linewidth = 0.7
        ) +
        scale_color_manual(values = condition_colors) +
        scale_y_continuous(
          limits = c(0, gene_y_max),
          expand = expansion(mult = c(0, 0.1)),
          n.breaks = 3
        ) +
        coord_cartesian(clip = "off") +
        ggh4x::facet_nested(
          ~ gene + inside,
          labeller = labeller(
            inside = as_labeller(
              c("in" = "Inside", "out" = "Outside")
            ),
            gene = as_labeller(\(x) paste0(x))
          ),
          strip = ggh4x::strip_nested(clip = "off")
        ) +
        guides(x = guide_axis(angle = 90)) +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.line = element_line(colour = "black"),
          strip.text.x = element_text(size = 12),
          strip.placement = "outside",
          panel.spacing.x = unit(2, "mm"),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          strip.background = element_rect(fill = "white", colour = NA),
          panel.grid = element_blank()
        )
    }) %>%
    cowplot::plot_grid(plotlist = ., nrow = 1)
  
  return(comparison_plots)
}

# Get data for both methods
data_lemur <- make_pb_plot_same_axis(fit_small, nei, return_data = TRUE)
data_cs <- make_pb_plot_same_axis(fit_small, res_pb, return_data = TRUE)

# Maximum expression per gene, taking the maximum across both datasets
y_max <- bind_rows(
  data_lemur %>% select(gene, expr),
  data_cs %>% select(gene, expr)
) %>%
  group_by(gene) %>%
  summarise(
    y_max = max(expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tibble::deframe()

# Make both plots with the same gene-specific limits
pb_lemur <- make_pb_plot_same_axis(
  fit_small,
  nei,
  y_max = y_max
)

pb_cs <- make_pb_plot_same_axis(
  fit_small,
  res_pb,
  y_max = y_max
)

plot_lemur <- pb_lemur +
  labs(title = "LEMUR")  +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5)
  )

plot_conf_selection <- pb_cs +
  labs(title = "Conformal selection on LEMUR")  +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5)
  )

plot_lemur / plot_conf_selection