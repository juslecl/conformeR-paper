.libPaths(c("~/R/library", .libPaths()))
library(SingleCellExperiment)
library(tidyverse)
library(glue)
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

#sce <- readRDS("~/glioblastoma_annotated_sce.rds")
#sce <- sce[,sce$condition!="etoposide"]
#sce <- sce[,sce$patient_id!="PW029"]

seed=2
home_dir <- Sys.getenv("HOME")
sce <- readRDS(paste0(home_dir,"/sce_5pat_2cond.RDS"))
#already <- readRDS(paste0(home_dir,"/res_pb0608.rds"))
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n = 5000)
sce <- sce[hvgs, ]

missing <- which(!(rownames(sce) %in% already$gene_id))
all_genes <- rownames(sce)[missing]

gene_batch_size=7

set.seed(2)

data_processing <- function(sce, replicate_id, obs_condition,
                            size_train = 0.3, size_cal = 0.25) {
  
  stopifnot(is(sce, "SingleCellExperiment"))
  colData(sce)[[replicate_id]] <- as.factor(colData(sce)[[replicate_id]])
  colData(sce)[[obs_condition]] <- as.factor(colData(sce)[[obs_condition]])
  colData_df <- as.data.frame(colData(sce)) |>
    dplyr::mutate(
      row = seq_len(nrow(colData(sce))),
      strata = interaction(colData(sce)[[obs_condition]], colData(sce)[[replicate_id]])
    )
  split1 <- initial_split(colData_df, prop = size_train, strata = strata)
  remaining <- training(split1)
  test_idx <- testing(split1)$row
  
  split2 <- initial_split(remaining, prop = 1 - size_cal, strata = strata)
  proper_idx <- training(split2)$row
  cal_idx <- testing(split2)$row
  
  train_set <- sce[, proper_idx]
  cal_set <- sce[, cal_idx]
  test_set <- sce[, test_idx]
  
  return(list(
    train = train_set,
    cal = cal_set,
    test = test_set))
}

make_cc <- function(pred_train, pred_cal, pred_test, nei_train, nei_cal) {
  genes <- nei_train$name
  softies <- lapply(
    genes,
    function(gene_name) {
      idx <- which(nei_train$name==gene_name)
      y <- as.factor(as.numeric(nei_train$neighborhood[[idx]]))
      
      expr_de <- assay(pred_train, "DE_panobinostat")[gene_name, ]
      
      embedding <- t(pred_train$embedding)
      colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
      
      X <- data.frame(
        expr_de = expr_de,
        embedding
      )
      
      pca <- prcomp(X, scale. = TRUE)
      
      var_expl <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
      k <- which(var_expl > 0.9)[1]
      
      dat <- data.frame(
        y = y,
        pca$x[, 1:k]
      )
      
      fit <- glm(y ~ ., data = dat, family = "binomial")
      
      list(
        fit = fit,
        pca = pca,
        k = k
      )
    }
  )
  
  pred_proba <- lapply(
    seq_along(genes),
    function(i) {
      
      gene_name <- genes[i]
      
      expr_de <- assay(pred_cal, "DE_panobinostat")[gene_name, ]
      
      embedding <- t(pred_cal$embedding)
      colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
      
      X_new <- data.frame(
        expr_de = expr_de,
        embedding
      )
      
      pca_scores <- predict(
        softies[[i]]$pca,
        newdata = X_new
      )
      
      dat_new <- data.frame(
        pca_scores[, seq_len(softies[[i]]$k), drop = FALSE]
      )
      
      pred <- predict(
        softies[[i]]$fit,
        newdata = dat_new,
        type = "response"
      )
      
      data.frame(
        pred = pred,
        cell_id = colnames(pred_cal),
        gene = gene_name
      )
    }
  )
  
  pred_proba_test <- lapply(
    seq_along(genes),
    function(i) {
      
      gene_name <- genes[i]
      
      expr_de <- assay(pred_test, "DE_panobinostat")[gene_name, ]
      
      embedding <- t(pred_test$embedding)
      colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
      
      X_new <- data.frame(
        expr_de = expr_de,
        embedding
      )
      
      pca_scores <- predict(
        softies[[i]]$pca,
        newdata = X_new
      )
      
      dat_new <- data.frame(
        pca_scores[, seq_len(softies[[i]]$k), drop = FALSE]
      )
      
      pred <- predict(
        softies[[i]]$fit,
        newdata = dat_new,
        type = "response"
      )
      
      data.frame(
        pred = pred,
        cell = colnames(pred_test),
        gene = gene_name
      )
    }
  )
  
  alpha <- 0.05
  pred_proba <- do.call(rbind.data.frame, pred_proba)
  
  gt_cal <- nei_cal |> 
    unnest(neighborhood) |> 
    mutate(cell_id=rep(colnames(pred_cal),nrow(nei_cal))) |>
    dplyr::rename("gene"="name")
  
  score_cal <- pred_proba |> left_join(gt_cal, by = c("cell_id", "gene"))
  
  score_cal <- score_cal |> mutate(scores_cal = 1000 *neighborhood - pred)
  
  scores_test <- lapply(
    seq_along(genes),
    function(row) {
      
      cbind.data.frame(
        scores_test = -pred_proba_test[[row]]$pred,
        gene = genes[row],
        cell = colnames(pred_test)
      )
    }
  )
  
  conf_pval <- lapply(
    seq_along(scores_test),
    function(gene) {
      
      current_gene <- unique(score_cal$gene)[gene]
      cal_scores <- score_cal$scores_cal[
        score_cal$gene == current_gene
      ]
      
      p <- vapply(
        scores_test[[gene]][, 1],
        function(score)
          (sum(cal_scores < score) + runif(1)*(1+sum(cal_scores = score))) /
          (ncol(pred_cal) + 1),
        numeric(1)
      )
      
      cbind.data.frame(
        conformal_p_val = p,
        cell = colnames(pred_test)
      )
    }
  )
  
  scores_test <- do.call(rbind.data.frame, scores_test)
  
  q <- 0.05
  m <- ncol(pred_test)
  k_val <- 1:m
  
  k_max <- lapply(conf_pval, function(mat_gene) {
    p_adj <- p.adjust(mat_gene[, 1], method = "BH")
    sum(p_adj <= q)
  }) |> unlist()
  
  label_set <- lapply(seq_along(conf_pval), function(gene) conf_pval[[gene]] |>
                        as.data.frame() |> mutate(inside_cs = conformal_p_val < thres_in[gene],
                                                  gene = genes[gene]) |> dplyr::select(-conformal_p_val)) |>
    (\(lst) do.call(rbind.data.frame, lst))()
  
  return(label_set)
}


################################################################################
################################################################################
# analysis
################################################################################
################################################################################

set.seed(2)
obs_condition <- "condition"
replicate_id  <- "patient_id"

set.seed(2)
fit <- lemur::lemur(sce, design = ~ condition + patient_id, n_embedding = 60, test_fraction = 0.5)
reducedDim(fit, "fit_umap") <- uwot::umap(t(fit$embedding))
fit <- lemur::align_harmony(fit)
reducedDim(fit, "fit_al_umap") <- uwot::umap(t(fit$embedding))

set.seed(2)

fit <- lemur::test_de(fit, contrast = cond(condition = "panobinostat") - cond(condition = "ctrl"), new_assay_name = "DE_panobinostat")

set.seed(2)
split <- data_processing(fit, "patient_id", "condition")

pred_train <- split$train
pred_cal <- split$cal
pred_test <- split$test

saveRDS(pred_test,paste0(home_dir, "/fit.rds"))
rm(fit,sce)

nei_train <- lemur::find_de_neighborhoods(
  pred_train,
  group_by = vars(patient_id, condition),
  de_mat = assay(pred_train, "DE_panobinostat"),
  test_method = "edgeR"
)

nei_cal <- lemur::find_de_neighborhoods(
  pred_cal,
  group_by = vars(patient_id, condition),
  de_mat = assay(pred_cal, "DE_panobinostat"),
  test_method = "edgeR"
)

nei_test <- lemur::find_de_neighborhoods(
  pred_test,
  group_by = vars(patient_id, condition),
  de_mat = assay(pred_test, "DE_panobinostat"),
  test_method = "edgeR"
)
saveRDS(nei_test, paste0(home_dir, "/nei.rds"))

gene_chunks <- split(all_genes, ceiling(seq_along(all_genes) / gene_batch_size))


set.seed(seed)
param <- MulticoreParam(workers = 7)

chunk_results <- bplapply(
  seq_along(gene_chunks),
  function(i) {
    
    chunk_genes <- gene_chunks[[i]]
    
    pred_train_chunk <- pred_train[chunk_genes, ]
    pred_cal_chunk   <- pred_cal[chunk_genes, ]
    pred_test_chunk  <- pred_test[chunk_genes, ]
    
    nei_train_chunk <- nei_train |> dplyr::filter(name %in% chunk_genes) |>
      mutate(
        neighborhood = map(neighborhood, ~ colnames(pred_train) %in% .x)
      )
    nei_cal_chunk <- nei_cal |> dplyr::filter(name %in% chunk_genes) |>
      mutate(
        neighborhood = map(neighborhood, ~ colnames(pred_cal) %in% .x)
      )
    
    chunk_res <- make_cc(
      pred_train_chunk,
      pred_cal_chunk,
      pred_test_chunk,
      nei_train_chunk,
      nei_cal_chunk
    )
    
    rm(pred_train_chunk, pred_cal_chunk, pred_test_chunk, nei_train_chunk, nei_cal_chunk)
    gc(FALSE)
    
    pred_set_all <- chunk_res |> 
      group_by(gene) |>
      dplyr::summarize(cs_list=list(cell[inside_cs]),
                       .groups="drop") 
    saveRDS(pred_set_all, paste0(home_dir, "/results",i,".rds"))
    
  },
  BPPARAM = param
)
