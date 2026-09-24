library(SingleCellExperiment)
library(tidyverse)
library(lemur)
library(SingleCellExperiment)
source("~/conformeR-paper/semi_synth_simulations/sce_simulation.R")
source("~/conformeR-paper/semi_synth_simulations/sim_helpers.R")

# ----------------------------------------------
# ----------------------------------------------
# LOAD DATA
# ----------------------------------------------
# ----------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
seed <- args[1]
sce <- readRDS("sce_5pat_2cond.RDS")
sce <- sce[, sce$condition == "ctrl"]

# ----------------------------------------------
# ----------------------------------------------
# APPLY CONFORMAL LAYERS ON THE SIMULATED CORRUPTED LABELS
# ----------------------------------------------
# ----------------------------------------------
sim_analysis <- function(sce, epsilon, gene_batch_size = 7,seed=seed) {
  result <- tryCatch({
    set.seed(seed)
    sce <- label_shifter(sce, epsilon)
    genes_of_interest <- rowData(sce)[rowData(sce)$is_simulated,]$name
    set.seed(seed)
    fit <- lemur(sce, design = ~ fake_condition, n_embedding = 60, 
                 use_assay = "logcounts", test_fraction = 0.5)
    fit <- align_harmony(fit)
    fit <- test_de(fit, contrast = cond(fake_condition = "fake_trt") - cond(fake_condition = "fake_ctrl"))
    rowData(fit)$is_de_cell[rowData(fit)$lfc==0] <- list(rep(FALSE,ncol(fit)))
    set.seed(seed)
    
    split <- conformeR::data_processing(fit, "sample", "fake_condition")
    
    pred_train <- split$train
    pred_cal <- split$cal
    pred_test <- split$test
    
    #Update rowData --> is_de_cell
    idx_test <- match(colnames(pred_test), colnames(fit))
    rowData(pred_test)$is_de_cell <-
      lapply(rowData(pred_test)$is_de_cell, `[`, idx_test)
    idx_train <- match(colnames(pred_train), colnames(fit))
    rowData(pred_train)$is_de_cell <-
      lapply(rowData(pred_train)$is_de_cell, `[`, idx_train)
    idx_cal <- match(colnames(pred_cal), colnames(fit))
    rowData(pred_cal)$is_de_cell <-
      lapply(rowData(pred_cal)$is_de_cell, `[`, idx_cal)
    
    #Update rowData --> corrup_labels
    idx_test <- match(colnames(pred_test), colnames(fit))
    rowData(pred_test)$corrup_labels <-
      lapply(rowData(pred_test)$corrup_labels, `[`, idx_test)
    idx_train <- match(colnames(pred_train), colnames(fit))
    rowData(pred_train)$corrup_labels <-
      lapply(rowData(pred_train)$corrup_labels, `[`, idx_train)
    idx_cal <- match(colnames(pred_cal), colnames(fit))
    rowData(pred_cal)$corrup_labels <-
      lapply(rowData(pred_cal)$corrup_labels, `[`, idx_cal)
    
    gene_chunks <- split(genes_of_interest, ceiling(seq_along(genes_of_interest) / gene_batch_size))
    
    set.seed(seed)
    param <- MulticoreParam(workers = 7)
    
    chunk_results <- bplapply(
      seq_along(gene_chunks),
      function(i) {
        
        chunk_genes <- gene_chunks[[i]]
        
        pred_train_chunk <- pred_train[chunk_genes, ]
        pred_cal_chunk   <- pred_cal[chunk_genes, ]
        pred_test_chunk  <- pred_test[chunk_genes, ]
        
        chunk_res <- conf_layer_shifted(
          pred_train_chunk,
          pred_cal_chunk,
          pred_test_chunk
        )
        
        rm(pred_train_chunk, pred_cal_chunk, pred_test_chunk)
        gc(FALSE)
        
        chunk_res
      },
      BPPARAM = param
    )
    
    pred_set_all <- dplyr::bind_rows(chunk_results) |> group_by(gene) |>
      dplyr::summarize(cc_list=list(cell[inside_cc]),
                       cc_list_out=list(cell[outside_cc]),
                       cs_list=list(cell[inside_cs]),
                       cs_list_adj=list(cell[inside_cs_adj]),
                       .groups="drop") 
    
    pred_test_small <- pred_test[genes_of_interest, ]
    res <- fdr_tpr_rate_cs(pred_set_all, pred_test_small)
    res2 <- fdr_tpr_rate_cc(pred_set_all, pred_test_small)
    saveRDS(res, paste0("results_cs", args[1] ,"eps",epsilon,".rds"))
    saveRDS(res2, paste0("results_cc", args[1] ,"eps",epsilon,".rds"))
    list(seed = seed, status = "ok")
    
  }, error = function(e) {message("ERROR: ", conditionMessage(e))
    list(status = "error", message = conditionMessage(e))
  })
  result
}

new_sce <- sce_simulation(sce, seed=seed)
for (eps in seq(0,0.15,0.01)){
  sim_analysis(new_sce, epsilon=eps, seed=seed)
}

