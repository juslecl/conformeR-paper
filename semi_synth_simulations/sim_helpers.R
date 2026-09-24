ground_truth <- function(fit,output){
  gene_info <- rowData(fit)[rowData(fit)$is_simulated,] |> 
    as.data.frame() |>
    left_join(output, by=c("name"="gene")) |>
    mutate(de_size = map_dbl(is_de_cell, \(x) sum(x != 0))) |>
    mutate(
      is_de_cell = map(is_de_cell, ~ colnames(fit)[.x])
    )
  return(gene_info)
}

label_shifter <- function(fit, gamma = 0.5) {
  
  rd <- rowData(fit)
  idx_is_simulated <- which(rd$is_simulated)
  
  # Initialize columns
  rd$corrup_rate <- NA_real_
  rd$mismatch_rate <- NA_real_
  rd$corrup_labels <- vector("list", nrow(rd))
  
  rd$corrup_rate[idx_is_simulated] <- rep(gamma,length(idx_is_simulated))
  
  # Corrupt labels gene by gene
  for (i in seq_along(idx_is_simulated)) {
    
    id <- idx_is_simulated[i]
    true_labels <- as.vector(rd$is_de_cell[[id]])
    
    is_corrup <- rbinom(
      n = ncol(fit),
      size = 1,
      prob = gamma
    )
    
    random_labels <- rbinom(
      n = ncol(fit),
      size = 1,
      prob = 0.5
    )
    
    corrup_labels <- (1 - is_corrup) * true_labels +
      is_corrup * random_labels
    
    rd$corrup_labels[[id]] <- corrup_labels
    
    # Actual proportion of mismatched labels
    rd$mismatch_rate[id] <- mean(corrup_labels != true_labels)
  }
  
  rowData(fit) <- rd
  
  return(fit)
}


fdr_tpr_rate_cs <- function(output, fit) {
  
  gt <- ground_truth(fit, output)
  res <- lapply(seq_len(nrow(gt)), function(gene_idx) {
    
    truth <- gt$is_de_cell[[gene_idx]]
    pred_adj  <- gt[["cs_list_adj"]][[gene_idx]]
    pred  <- gt[["cs_list"]][[gene_idx]]
    lfc <- gt$lfc[[gene_idx]]
    
    TP <- if (lfc==0) 0 else length(intersect(pred, truth))
    FP <- if (lfc==0) length(pred) else length(setdiff(pred, truth))
    FN <- if (lfc==0) 0 else length(setdiff(truth, pred))
    
    TP_adj <- if (lfc==0) 0 else length(intersect(pred_adj, truth))
    FP_adj <- if (lfc==0) length(pred_adj) else length(setdiff(pred_adj, truth))
    FN_adj <- if (lfc==0) 0 else length(setdiff(truth, pred_adj))
    
    data.frame(
      name = gt$name[gene_idx],
      TPR = if ((TP + FN) == 0) NA_real_ else TP / (TP + FN),
      FDR = if ((TP + FP) == 0) 0 else FP / (TP + FP),
      TPR_adj = if ((TP_adj + FN_adj) == 0) NA_real_ else TP_adj / (TP_adj + FN_adj),
      FDR_adj = if ((TP_adj + FP_adj) == 0) 0 else FP_adj / (TP_adj + FP_adj),
      lfc = gt$lfc[[gene_idx]],
      size= gt$de_size[[gene_idx]]
    )
  })
  
  res <- dplyr::bind_rows(res)
  metric_cols <- c("TPR", "FDR","TPR_adj","FDR_adj")
  names(res)[match(metric_cols, names(res))] <-
    paste0(metric_cols, "_cs")
  
  res |>
    mutate(size=size/ncol(fit),
           FDR_cs = dplyr::coalesce(FDR_cs, 0))
}

fdr_tpr_rate_cc <- function(output, fit) {
  
  gt <- ground_truth(fit, output)
  all_cells <- colnames(fit)
  
  res <- lapply(seq_len(nrow(gt)), function(gene_idx) {
    
    truth <- gt$is_de_cell[[gene_idx]]
    pred_in  <- gt[["cc_list"]][[gene_idx]]
    pred_out <- gt[["cc_list_out"]][[gene_idx]]
    
    non_truth <- setdiff(all_cells, truth)
    
    covered <- length(intersect(truth, pred_in)) + length(intersect(non_truth, pred_out))
    total   <- length(truth) + length(non_truth)
    
    data.frame(
      name    = gt$name[gene_idx],
      covered = covered,
      total   = total,
      coverage = covered / total,
      lfc  = gt$lfc[[gene_idx]],
      size = gt$de_size[[gene_idx]]
    )
  })
  
  res <- dplyr::bind_rows(res)
  
  res |>
    mutate(overall_coverage = sum(covered) / sum(total))
}

conf_layer_shifted <- function(pred_train, pred_cal, pred_test) {
  genes <- rowData(pred_train)$name  
  softies <- lapply(
    genes,
    function(gene_name) {
      idx <- which(rowData(pred_train)$name == gene_name)
      y <- as.factor(as.numeric(rowData(pred_train)$corrup_labels[[idx]]))
      
      expr_de <- assay(pred_train, "DE")[gene_name, ]
      
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
      
      expr_de <- assay(pred_cal, "DE")[gene_name, ]
      
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
      
      expr_de <- assay(pred_test, "DE")[gene_name, ]
      
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
  
  n_cal <- ncol(pred_cal)
  K <- 2
  n_test <- ncol(pred_test)
  alpha <- 0.05
  
  pred_set_list <- lapply(
    seq_along(genes),
    function(row) {
      p_cal <- pred_proba[[row]]$pred
      y_cal <- as.numeric(
        rowData(pred_cal)$corrup_labels[[row]]
      )
      s_cal <- ifelse(
        y_cal == 1,
        1 - p_cal,
        p_cal
      )
      q <- sort(s_cal)[ceiling((n_cal + 1) * (1 - alpha))]
      p_test <- pred_proba_test[[row]]$pred
      cbind.data.frame(
        inside_cc  = (1 - p_test <= q),
        outside_cc = (p_test <= q)
      ) |>
        mutate(
          cell = colnames(pred_test),
          gene = genes[row]
        )
    }
  )
  pred_set <- do.call(rbind.data.frame, pred_set_list)
  pred_proba <- do.call(rbind.data.frame, pred_proba)
  gt_cal <- rowData(pred_cal) |> 
    as.data.frame() |> 
    dplyr::select(name,corrup_labels) |>
    unnest(corrup_labels) |> 
    mutate(cell_id=rep(colnames(pred_cal),nrow(pred_cal))) |>
    dplyr::rename("gene"="name")
  score_cal <- pred_proba |> left_join(gt_cal, by = c("cell_id", "gene"))
  score_cal <- score_cal |> mutate(scores_cal = 1000 *corrup_labels - pred)
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
      corrup_shift <- rowData(pred_cal)$corrup_rate[[gene]]/2
      p1 <- vapply(
        scores_test[[gene]][, 1],
        function(score)
          min(1,((sum(cal_scores < score)+runif(1)*(1+sum(cal_scores == score))) /
                   (ncol(pred_cal) + 1))+corrup_shift),
        numeric(1)
      )
      p2 <- vapply(
        scores_test[[gene]][, 1],
        function(score)
          ((sum(cal_scores < score)+runif(1)*(1+sum(cal_scores == score))) /
             (ncol(pred_cal) + 1)),
        numeric(1)
      )
      cbind.data.frame(
        conformal_p_val_adj = p1,
        conformal_p_val = p2,
        cell = colnames(pred_test)
      )
    }
  )
  scores_test <- do.call(rbind.data.frame, scores_test)
  label_set <- lapply(seq_along(conf_pval), function(gene) {
    
    df <- conf_pval[[gene]] |>
      as.data.frame()
    
    df <- df |>
      mutate(
        inside_cs_adj = p.adjust(
          as.numeric(conformal_p_val_adj),
          method = "BH"
        ) < alpha,
        inside_cs = p.adjust(
          as.numeric(conformal_p_val),
          method = "BH"
        ) < alpha,
        gene = genes[gene]
      ) 
    
    df
  }) |>
    (\(lst) do.call(rbind.data.frame, lst))()
  return(left_join(label_set,pred_set, by=c("gene","cell")))
}


