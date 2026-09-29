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
    group = factor(`abs(lfc)`),
    colour = factor(`abs(lfc)`)
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
    breaks = seq(0, 1, 0.05)
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

res_cc <- res_cc |> dplyr::group_by(eps) |> dplyr::summarize(mean_cov=mean(mean_cov))
ggplot2::ggplot(
  res_cc,
  ggplot2::aes(
    x = eps,
    y = mean_cov
  )
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(
    x = expression(epsilon),
    y = "Marginal coverage"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.05)
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
    cols = c(TPR_cs, FDR_cs),
    names_to = c("metric", "method"),
    names_pattern = "^(TPR|FDR)_(.*)$",
    values_to = "value"
  ) |>
  dplyr::filter(lfc != 0) |>
  dplyr::group_by(abs(lfc), eps, metric) |>
  dplyr::summarise(
    mean_val = mean(value),
    .groups = "drop"
  )

ggplot2::ggplot(
  res_cs,
  ggplot2::aes(
    x = eps,
    y = mean_val,
    colour = metric,
    group = interaction(metric, `abs(lfc)`),
    linetype = factor(`abs(lfc)`)
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
    colour = "Metric",
    linetype = "|LFC|"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 16),
    legend.text = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 15)
  )

res_cs <- res_cs |> dplyr::group_by(eps,metric) |> dplyr::summarize(mean_val=mean(mean_val))

ggplot2::ggplot(
  res_cs,
  ggplot2::aes(
    x = eps,
    y = mean_val,
    colour = metric
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
    colour = "Method / metric"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 16),
    legend.text = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 15)
  )