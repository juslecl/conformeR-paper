sce_simulation <- function(
    sce,
    n_hvgs = 8000,
    n_de_genes = 200,
    cut_at = 5,
    lfc_mean = c(0, .5, 1, 2),
    lfc_sd = 0,
    sample_column = "patient_id",
    seed = 1
) {
  
  set.seed(seed)
  
  sce$fake_condition <- sample(
    c("fake_ctrl", "fake_trt"),
    size = ncol(sce),
    replace = TRUE
  )
  
  sce$sample <- sce$patient_id
  
  # Generate DE
  set.seed(seed)
  
  pca <- lemur:::pca(
    assay(sce, "logcounts"),
    n = min(nrow(sce), ncol(sce), 50)
  )
  
  set.seed(seed)
  
  kmeans_clusterings <- lapply(
    unique(cut_at),
    \(k) stats::kmeans(
      t(pca$embedding),
      centers = k
    )$cluster
  )
  
  names(kmeans_clusterings) <- as.character(unique(cut_at))
  
  set.seed(seed)
  
  de_args <- tibble::tibble(
    name = paste0("simulated_gene-", seq_len(n_de_genes)),
    is_simulated = TRUE,
    cut_at = cut_at,
    base_expr = stats::runif(
      n = n_de_genes,
      min = -7,
      max = 3
    ),
    lfc = rep_len(
      lfc_mean,
      length.out = n_de_genes
    ) *
      sample(
        c(-1, 1),
        size = n_de_genes,
        replace = TRUE
      ),
    is_de_cell = purrr::map(
      n_de_genes,
      \(.) rep(FALSE, ncol(sce))
    )
  )
  
  new_counts <- matrix(
    0,
    nrow = n_de_genes,
    ncol = ncol(sce)
  )
  
  for (idx in seq_len(n_de_genes)) {
    
    gv <- generate_gene_for_cluster(
      n_clusters = de_args$cut_at[idx],
      base_expr = de_args$base_expr[idx],
      lfc_mean = de_args$lfc[idx],
      lfc_sd = 0,
      sample_sd = 0.1,
      overdispersion = 0.2
    )
    
    new_counts[idx, ] <- gv$counts
    de_args$is_de_cell[[idx]] <- gv$is_de_cell
  }
  
  neg_de_args <- tibble::tibble(
    name = rownames(sce),
    is_simulated = FALSE
  )
  
  new_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = rbind(
        unname(SummarizedExperiment::assay(sce, "counts")),
        new_counts
      )
    ),
    colData = SummarizedExperiment::colData(sce),
    rowData = S4Vectors::DataFrame(
      dplyr::bind_rows(
        neg_de_args,
        de_args
      )
    )
  )
  
  colnames(new_sce) <- colnames(sce)
  rownames(new_sce) <- SummarizedExperiment::rowData(new_sce)$name
  
  SummarizedExperiment::rowData(new_sce)$seed <- seed
  
  SummarizedExperiment::assay(
    new_sce,
    "logcounts"
  ) <- transformGamPoi::shifted_log_transform(new_sce)
  
  new_sce
}


generate_gene_for_cluster <- function(
    n_clusters = 10,
    base_expr = -2,
    lfc_mean = lfc_mean,
    lfc_sd = lfc_sd,
    sample_sd = 0.1,
    overdispersion = 0.1,
    ...
) {
  
  cluster_assign <- if (
    as.character(n_clusters) %in% names(kmeans_clusterings)
  ) {
    
    kmeans_clusterings[[as.character(n_clusters)]]
    
  } else {
    
    stats::kmeans(
      t(pca$embedding),
      centers = n_clusters
    )$cluster
  }
  
  set.seed(seed)
  
  sel_cluster <- sample(
    unique(cluster_assign),
    size = 1
  )
  
  is_de_cell <- cluster_assign == sel_cluster
  
  eta_ctrl <- 0
  
  eta_trt <- stats::rnorm(
    1,
    mean = lfc_mean,
    sd = lfc_sd
  )
  
  trt_eff <- colSums(
    lemur:::one_hot_encoding(
      sce$fake_condition
    )[c("fake_ctrl", "fake_trt"), ] *
      c(eta_ctrl, eta_trt)
  )
  
  ctrl_mean <- stats::rnorm(
    length(unique(sce$sample)),
    mean = 0,
    sd = sample_sd
  )
  
  ctrl_eff <- colSums(
    lemur:::one_hot_encoding(sce$sample) *
      ctrl_mean
  )
  
  sf <- colSums(
    SummarizedExperiment::assay(sce, "counts")
  )
  
  sf <- sf / stats::median(sf)
  
  mu <- 2^(
    trt_eff * is_de_cell +
      ctrl_eff +
      base_expr
  ) * sf
  
  counts <- stats::rnbinom(
    n = ncol(sce),
    mu = mu,
    size = 1 / overdispersion
  )
  
  list(
    n_clusters = n_clusters,
    sel_cluster = sel_cluster,
    is_de_cell = is_de_cell,
    lfc = eta_trt - eta_ctrl,
    base_expr = base_expr,
    log_expression_level = unname(
      trt_eff * is_de_cell +
        ctrl_eff +
        base_expr
    ),
    ctrl_mean = ctrl_mean,
    counts = counts
  )
}
