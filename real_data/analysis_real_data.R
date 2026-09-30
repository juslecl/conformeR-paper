library(SingleCellExperiment)
set.seed(2002)

# Load data
sce <- readRDS("sce_5pat_2cond.RDS")

# Define the genes of interest
pan_affect_genes <- rowData(sce) %>%
  tibble::as_tibble() %>%
  dplyr::select(gid, gene) %>%
  dplyr::filter(gene %in% c("EPC1", "MBP", "HIST3H2A", "PTPRZ1", "HBEGF","SPATA13"))

rownames(sce) <- rowData(sce)$gene
# Run conformeR
res <- conformeR(
  sce,
  design_lemur = ~ patient_id + condition,
  contrast_column = "condition",
  genes_of_interest = pan_affect_genes$gene,
  what = c("conf_selection", "conf_clustering"),
  n_cores = 7,
  verbose = TRUE
)

# Plot
cp_res <- res$conf_results |>
  left_join(pan_affect_genes, by = c("gene" = "gid")) |>
  dplyr::select(-gene) |>
  dplyr::rename(gene = gene.y)

plotter_conformal_selection(res, pan_affect_genes$gene)
plotter_conformal_clustering(res, pan_affect_genes$gene)