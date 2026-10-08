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
    sd_cov = sd(coverage),
    .groups = "drop"
  )

ggplot2::ggplot(
  res_cc,
  ggplot2::aes(
    x = eps,
    y = mean_cov,
    group = factor(`abs(lfc)`),
    colour = factor(`abs(lfc)`),
    fill = factor(`abs(lfc)`)
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mean_cov - sd_cov, ymax = mean_cov + sd_cov),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(
    x = expression(epsilon),
    y = "Marginal coverage",
    colour = "|LFC|",
    fill = "|LFC|"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.05)
  ) +
  ggplot2::scale_y_continuous(breaks = seq(0.94, 1, 0.01)) +
  ggplot2::coord_cartesian(ylim = c(0.94, 1)) +
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

res_cc <- outputs_cc |>
  dplyr::filter(lfc != 0) |>
  dplyr::group_by(eps) |>
  dplyr::summarise(
    mean_cov = mean(coverage),
    sd_cov = sd(coverage),
    .groups = "drop"
  )

ggplot2::ggplot(
  res_cc,
  ggplot2::aes(
    x = eps,
    y = mean_cov
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mean_cov - sd_cov, ymax = mean_cov + sd_cov),
    alpha = 0.2,
    fill = "grey40"
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
  ggplot2::scale_y_continuous(breaks = seq(0.94, 1, 0.01)) +
  ggplot2::coord_cartesian(ylim = c(0.94, 1)) +
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
  dplyr::filter(lfc != 0)

res_cs_lfc <- res_cs |>
  dplyr::group_by(abs(lfc), eps, metric) |>
  dplyr::summarise(
    mean_val = mean(value),
    sd_val = sd(value),
    .groups = "drop"
  )

ggplot2::ggplot(
  res_cs_lfc,
  ggplot2::aes(
    x = eps,
    y = mean_val,
    colour = metric,
    fill = metric,
    group = interaction(metric, `abs(lfc)`),
    linetype = factor(`abs(lfc)`)
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    alpha = 0.12,
    colour = NA
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
    fill = "Metric",
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

res_cs <- res_cs |>
  dplyr::group_by(eps, metric) |>
  dplyr::summarise(
    mean_val = mean(value),
    sd_val = sd(value),
    .groups = "drop"
  )

ggplot2::ggplot(
  res_cs,
  ggplot2::aes(
    x = eps,
    y = mean_val,
    colour = metric,
    fill = metric
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    alpha = 0.12,
    colour = NA
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
    fill = "Method / metric"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 16),
    legend.text = ggplot2::element_text(size = 14),
    legend.title = ggplot2::element_text(size = 15)
  )
#################### cc clustering ################################################
## ---- Settings -------------------------------------------------------------- 
cc_pattern <- "results_cc(\\d+)eps([0-9.]+)\\.(rds|RData)$"
data_dir=getwd()
## ---- Helpers ---------------------------------------------------------------
# Load an .rds (returns the object) or an .RData (returns the single object inside)
load_obj <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  e <- new.env()
  nm <- load(path, envir = e)
  get(nm[1], envir = e)
}

# For one gene: how many cells appear 0, 1 or 2 times across the two lists
# (inside_cc = cc_list, outside_cc = cc_list_out). A cell appearing twice is
# in both lists.
count_occurrences <- function(inside, outside, total) {
  inside  <- unique(inside)
  outside <- unique(outside)
  n_twice <- length(intersect(inside, outside))
  n_once  <- length(union(inside, outside)) - n_twice
  n_none  <- total - n_once - n_twice
  c(n_none = n_none, n_once = n_once, n_twice = n_twice)
}

## ---- Process one (seed, epsilon) combination -------------------------------
process_one <- function(cc_file, xl_file, seed, epsilon) {
  cc <- load_obj(cc_file)
  xl <- load_obj(xl_file)
  
  # gene order differs between the two objects -> match by gene name
  idx <- match(xl$gene, cc$name)
  if (anyNA(idx)) warning("Some genes in results_xl not found in results_cc: ", cc_file)
  total <- cc$total[idx]
  
  counts <- t(vapply(seq_len(nrow(xl)), function(i) {
    count_occurrences(xl$cc_list[[i]], xl$cc_list_out[[i]], total[i])
  }, numeric(3)))
  
  data.frame(
    seed          = seed,
    epsilon       = epsilon,
    gene          = xl$gene,
    lfc           = cc$lfc[idx],
    coverage      = cc$coverage[idx],
    conf_quantile = xl$q,
    total_cells   = total,
    n_none        = counts[, "n_none"],
    n_once        = counts[, "n_once"],
    n_twice       = counts[, "n_twice"],
    prop_0     = counts[, "n_none"]  / total,
    prop_1     = counts[, "n_once"]  / total,
    prop_2    = counts[, "n_twice"] / total,
    n_inside_cc   = lengths(xl$cc_list),
    n_outside_cc  = lengths(xl$cc_list_out),
    n_cs          = lengths(xl$cs_list),   # size of cs_list, for reference
    stringsAsFactors = FALSE
  )
}

## ---- Loop over all seeds / epsilons ----------------------------------------
cc_files <- list.files(pattern = cc_pattern, full.names = FALSE)

all_results <- lapply(cc_files, function(f) {
  m       <- regmatches(f, regexec(cc_pattern, f))[[1]]
  seed    <- as.integer(m[2])
  epsilon <- as.numeric(m[3])
  ext     <- m[4]
  xl_file <- file.path(data_dir, sprintf("results_xl%deps%s.%s", seed, m[3], ext))
  if (!file.exists(xl_file)) {
    warning("Missing results_xl file for seed ", seed, ", epsilon ", epsilon)
    return(NULL)
  }
  process_one(file.path(data_dir, f), xl_file, seed, epsilon)
})

summary_df <- bind_rows(all_results) %>% filter(lfc !=0) %>% arrange(epsilon, seed, gene)

## ---- Mean +/- SD across seeds ----------------------------------------------
# Step 1: average over genes within each seed (per epsilon and lfc)
per_seed <- summary_df %>%
  group_by(epsilon, seed, lfc) %>%
  summarise(across(c(prop_0, prop_1, prop_2, conf_quantile), mean),
            .groups = "drop")

# Step 2: mean and SD across seeds
summary_by_lfc <- per_seed %>%
  group_by(epsilon, lfc) %>%
  summarise(across(c(prop_0, prop_1, prop_2, conf_quantile),
                   list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
            n_seeds = n(), .groups = "drop")

## ---- Plot: mean proportion +/- SD -------------------------------------------
plot_df <- per_seed %>%
  select(epsilon, seed, lfc, starts_with("prop_")) %>%
  pivot_longer(starts_with("prop_"), names_to = "Set size",
               names_prefix = "prop_", values_to = "prop") %>%
  mutate(`Set size` = factor(`Set size`)) %>%
  group_by(epsilon, abs(lfc),`Set size`) %>%
  summarise(mean = mean(prop), sd = sd(prop), .groups = "drop")

p_bar <- ggplot(plot_df, aes(x = factor(`abs(lfc)`), y = mean, fill = `Set size`)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.9), width = 0.25) +
  facet_wrap(~ epsilon, labeller = label_both) +
  labs(x = "Absolute log-fold change", y = "Proportion of cells (mean \u00b1 SD across seeds)",
       fill = "Prediction set size") +
  theme_bw()

print(p_bar)

p <- ggplot(plot_df, aes(x = epsilon, y = mean,
                         colour = `Set size`, linetype = as.factor(`abs(lfc)`),
                         group = interaction(`Set size`, as.factor(`abs(lfc)`)))) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
              alpha = 0.12, colour = NA) +
  geom_line() +
  geom_point(size = 1.5) +
  labs(x = expression(epsilon),
       y = "Proportion of cells (mean \u00b1 SD across seeds)",
       colour = "Prediction set size", fill = "Prediction set size",
       linetype = "Absolute log-fold change") +
  theme_bw()

plot(p)
