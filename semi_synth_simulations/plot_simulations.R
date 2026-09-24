## --- Results conformal clustering -----------

outputs_cc <- lapply(seq_len(50), function(i) {
  files <- list.files(
    pattern = paste0("results_cc", i, "eps.*\\.rds$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    return(NULL)
  }
  
  lapply(files, function(file) {
    eps <- as.numeric(
      sub(".*eps([0-9.]+)\\.rds$", "\\1", basename(file))
    )
    
    readRDS(file) |>
      dplyr::mutate(
        seed = i,
        eps = eps
      )
  }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(size = size / total)


res_cc <- outputs_cc |>
  dplyr::filter(lfc != 0) |>
  dplyr::group_by(abs(lfc), eps) |>
  dplyr::summarise(
    mean_cov = mean(coverage),
    .groups = "drop"
  )


ggplot2::ggplot(
  res_cc,
  ggplot2::aes(
    x = eps,
    y = mean_cov,
    group = factor(abs(lfc)),
    colour = factor(abs(lfc))
  )
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(
    x = expression(epsilon),
    y = "Marginal coverage",
    colour = "|LFC|"
  ) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(res_cc$eps))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0.95, 1),
    breaks = seq(0.95, 1, 0.01)
  ) +
  ggplot2::geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    colour = "black"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 16),
    axis.title = ggplot2::element_text(size = 18),
    legend.text = ggplot2::element_text(size = 16),
    legend.title = ggplot2::element_text(size = 18)
  )


## --- Results conformal selection (with shifted method) -----------------

outputs_cs <- lapply(seq_len(50), function(i) {
  
  files <- list.files(
    pattern = paste0("results_cs", i, "eps.*\\.rds$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    return(NULL)
  }
  
  lapply(files, function(file) {
    eps <- as.numeric(
      sub(".*eps([0-9.]+)\\.rds$", "\\1", basename(file))
    )
    
    readRDS(file) |>
      dplyr::mutate(
        seed = i,
        eps = eps
      )
  }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows()


res_cs <- outputs_cs |>
  tidyr::pivot_longer(
    cols = c(TPR_cs, FDR_cs, TPR_adj_cs, FDR_adj_cs),
    names_to = c("metric", "method"),
    names_pattern = "^(TPR|FDR)_(.*)$",
    values_to = "value"
  ) |>
  dplyr::mutate(
    method = dplyr::recode(
      method,
      cs = "CS",
      adj_cs = "CS shifted p-val"
    )
  ) |>
  dplyr::filter(lfc != 0) |>
  dplyr::group_by(abs(lfc), eps, method, metric) |>
  dplyr::summarise(
    mean_val = mean(value),
    .groups = "drop"
  )


ggplot2::ggplot(
  res_cs,
  ggplot2::aes(
    x = eps,
    y = mean_val,
    colour = interaction(method, metric),
    group = interaction(method, metric, abs(lfc)),
    linetype = factor(abs(lfc))
  )
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    colour = "black"
  ) +
  ggplot2::labs(
    x = expression(epsilon),
    y = NULL,
    colour = "Method / metric",
    linetype = "|LFC|"
  ) +
  ggplot2::theme_bw()