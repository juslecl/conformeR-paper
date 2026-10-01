library(SingleCellExperiment)
library(dplyr)
set.seed(1)

# Load data
sce <- readRDS("sce_5pat_2cond.RDS")

# Define the genes of interest
pan_affect_genes <- c("EPC1", "MBP", "HIST3H2A", "PTPRZ1", "HBEGF","SPATA13")

rownames(sce) <- rowData(sce)$gene

# Run conformeR
res <- conformeR::conformeR(
  sce,
  design_lemur = ~ patient_id + condition,
  contrast_column = "condition",
  genes_of_interest = pan_affect_genes,
  what = c("conf_selection", "conf_clustering"),
  n_cores = 7,
  verbose = TRUE
)

# Plot
plot_cs <- conformeR::plotter_conformal_selection(res, pan_affect_genes)
# 2 plots for conformal clustering, so that each panels are not too small.
plot_cc1 <- conformeR::plotter_conformal_clustering(res, pan_affect_genes[1:3])
plot_cc2 <- conformeR::plotter_conformal_clustering(res, pan_affect_genes[4:6])
