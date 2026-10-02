################################################################################
# data upload and preparation
################################################################################
sce <- readRDS("sce_5pat_2cond.RDS")
dec <- scran::modelGeneVar(sce)
hvgs <- scran::getTopHVGs(dec, n = 5000)
sce <- sce[hvgs, ]
all_genes <- c("HIST3H2A", "MBP","SPATA13")
gene_batch_size <- 7
rownames(sce) <- rowData(sce)$gene

################################################################################
# analysis
################################################################################
obs_condition <- "condition"
replicate_id  <- "patient_id"
set.seed(2002)
fit <- lemur::lemur(
  sce,
  design = ~ condition + patient_id,
  n_embedding = 60,
  test_fraction = 0.5
)
SingleCellExperiment::reducedDim(fit, "fit_umap") <-
  uwot::umap(t(fit$embedding))
fit <- lemur::align_harmony(fit)
SingleCellExperiment::reducedDim(fit, "fit_al_umap") <-
  uwot::umap(t(fit$embedding))
fit <- lemur::test_de(fit,
                      contrast = cond(condition = "panobinostat") -
                        cond(condition = "ctrl"))
set.seed(2002)
split <- conformeR::data_processing(
  fit,
  strat_by = c("patient_id", "condition"),
  size_train = 0.5,
  size_cal = 0.25
)

pred_train <- split$train
pred_cal   <- split$cal
pred_test  <- split$test

saveRDS(pred_test[all_genes,], "fit.rds")
rm(fit, sce)

nei_train <- lemur::find_de_neighborhoods(
  pred_train,
  group_by = dplyr::vars(patient_id, condition),
  de_mat = SummarizedExperiment::assay(pred_train, "DE"),
  test_method = "edgeR"
)

nei_cal <- lemur::find_de_neighborhoods(
  pred_cal,
  group_by = dplyr::vars(patient_id, condition),
  de_mat = SummarizedExperiment::assay(pred_cal, "DE"),
  test_method = "edgeR"
)

nei_test <- lemur::find_de_neighborhoods(
  pred_test,
  group_by = dplyr::vars(patient_id, condition),
  de_mat = SummarizedExperiment::assay(pred_test, "DE"),
  test_method = "edgeR"
)

saveRDS(nei_test |> dplyr::filter(name %in% all_genes), "nei.rds")

gene_chunks <- split(all_genes, ceiling(seq_along(all_genes) / gene_batch_size))

set.seed(2002)

param <- BiocParallel::MulticoreParam(workers = 7)

chunk_results <- BiocParallel::bplapply(seq_along(gene_chunks), function(i) {
  chunk_genes <- gene_chunks[[i]]
  
  pred_train_chunk <- pred_train[chunk_genes, ]
  pred_cal_chunk   <- pred_cal[chunk_genes, ]
  pred_test_chunk  <- pred_test[chunk_genes, ]
  
  nei_train_chunk <-
    nei_train |>
    dplyr::filter(name %in% chunk_genes) |>
    dplyr::mutate(neighborhood = purrr::map(neighborhood, ~ colnames(pred_train) %in% .x))
  
  nei_cal_chunk <-
    nei_cal |>
    dplyr::filter(name %in% chunk_genes) |>
    dplyr::mutate(neighborhood = purrr::map(neighborhood, ~ colnames(pred_cal) %in% .x))
  
  chunk_res <- conformeR::conformal_lemur(
    pred_train_chunk,
    pred_cal_chunk,
    pred_test_chunk,
    nei_train_chunk,
    nei_cal_chunk,
    what=c("conf_selection","conf_clustering"),
    alpha = 0.05
  )
  
  rm(
    pred_train_chunk,
    pred_cal_chunk,
    pred_test_chunk,
    nei_train_chunk,
    nei_cal_chunk
  )
  return(chunk_res)
}, BPPARAM = param)

pred_set_all <- dplyr::bind_rows(chunk_results)
saveRDS(pred_set_all, "pred_set_realdata.RDS")

reducedDim(pred_test, "fit_al_umap")[,2] <- -reducedDim(pred_test, "fit_al_umap")[,2]

plot_cs <- conformeR::plotter_conformal_selection(list(fit_lemur=pred_test[all_genes,],nei_lemur=nei_test |> dplyr::filter(name %in% all_genes), conf_results=pred_set_all), all_genes)
plot_cc <- conformeR::plotter_conformal_clustering(list(fit_lemur=pred_test[all_genes,],nei_lemur=nei_test |> dplyr::filter(name %in% all_genes), conf_results=pred_set_all), all_genes)
