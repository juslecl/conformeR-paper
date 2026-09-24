.libPaths(c("~/R/library", .libPaths()))
library(tidyverse)
library(glue)
library(SingleCellExperiment)
library(lemur)
library(rsample)
library(rlang)
library(simdata)
library(scater)
library(scuttle)
library(scran)
library(MASS)
library(tidyr)
library(tibble)
library(purrr)
library(scales)
library(BiocParallel)
library(conformeR)

# LOAD DATA
sce_full <- readRDS("~/sce_5pat_2cond.RDS")
pan_affect_genes <- c("HBEGF", "MBP")
pan_genes <- rowData(sce_full) |>
  as.data.frame() |>
  dplyr::select(gene, gid) |>
  dplyr::filter(gene %in% c(pan_affect_genes))

# RUN ONE SEED
run_replication <- function(seed, sce_full=sce_full, genes_of_interest=pan_genes$gid, obs_condition="condition",
                            replicate_id="patient_id", gene_batch_size = 7) {
  result <- tryCatch({
    colData_df <- as.data.frame(colData(sce_full)) |>
      dplyr::mutate(
        row = seq_len(nrow(colData(sce_full))),
        strata = interaction(colData(sce_full)[[obs_condition]], colData(sce_full)[[replicate_id]])
      )
    set.seed(seed)
    split1 <- initial_split(colData_df, prop = .8, strata = strata)
    select_idx <- training(split1)$row
    sce <- sce_full[, select_idx]
    set.seed(seed)
    conformeR::conformeR(sce,
                         design_lemur = patient_id+condition,
                         contrast_column = "condition",
                         genes_of_interest = genes_of_interest)

  }, error = function(e) {
    list(seed = seed, status = "error", message = conditionMessage(e))
  })

  result
}

seeds <- c(10,11)

results <- lapply(
  seeds,
  function(s) run_replication(seed=s)
)

saveRDS(results[[1]]$fit_lemur, "fitseed10_standard.rds")
saveRDS(results[[1]]$nei_lemur, "neiseed10_standard.rds")
saveRDS(results[[1]]$conf_results, "pred_set_pan10.rds")
saveRDS(results[[2]]$fit_lemur, "fitseed11_standard.rds")
saveRDS(results[[2]]$nei_lemur, "neiseed11_standard.rds")
saveRDS(results[[2]]$conf_results, "pred_set_pan11.rds")
