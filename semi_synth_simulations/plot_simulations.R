library(dplyr)
library(tidyr)
library(ggplot2)

## =============================================================================
## Conformal clustering: marginal coverage
## =============================================================================

outputs_cc <- lapply(seq_len(50), function(i) {
  files <- list.files(
    pattern = paste0("results_cc", i, "eps.*\\.rds$"),
    full.names = TRUE
  )
  if (length(files) == 0) return(NULL)
  
  lapply(files, function(file) {
    eps <- as.numeric(sub(".*eps([0-9.]+)\\.rds$", "\\1", basename(file)))
    readRDS(file) |>
      dplyr::mutate(seed = i, eps = eps)
  }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(size = size / total)

# Step 1: average over genes within each seed (per eps and |lfc|)
cc_per_seed <- outputs_cc |>
  dplyr::filter(lfc != 0) |>
  dplyr::mutate(abs_lfc = abs(lfc)) |>
  dplyr::group_by(seed, eps, abs_lfc) |>
  dplyr::summarise(cov = mean(coverage), .groups = "drop")

# Step 2: mean and SD across seeds
res_cc <- cc_per_seed |>
  dplyr::group_by(eps, abs_lfc) |>
  dplyr::summarise(
    mean_cov = mean(cov),
    sd_cov = sd(cov),
    .groups = "drop"
  )

ggplot2::ggplot(
  res_cc,
  ggplot2::aes(
    x = eps,
    y = mean_cov,
    group = factor(abs_lfc),
    colour = factor(abs_lfc),
    fill = factor(abs_lfc)
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
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.05)) +
  ggplot2::scale_y_continuous(breaks = seq(0.94, 1, 0.01)) +
  ggplot2::coord_cartesian(ylim = c(0.94, 1)) +
  ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", colour = "black") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_line(colour = "grey90"),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()
  ) 

# Overall (all |lfc| pooled): average over genes within seed, then across seeds
cc_per_seed_all <- outputs_cc |>
  dplyr::filter(lfc != 0) |>
  dplyr::group_by(seed, eps) |>
  dplyr::summarise(cov = mean(coverage), .groups = "drop")

res_cc_all <- cc_per_seed_all |>
  dplyr::group_by(eps) |>
  dplyr::summarise(
    mean_cov = mean(cov),
    sd_cov = sd(cov),
    .groups = "drop"
  )

ggplot2::ggplot(
  res_cc_all,
  ggplot2::aes(x = eps, y = mean_cov)
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
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.05)) +
  ggplot2::scale_y_continuous(breaks = seq(0.94, 1, 0.01)) +
  ggplot2::coord_cartesian(ylim = c(0.94, 1)) +
  ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", colour = "black") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_line(colour = "grey90"),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()
  ) 


## =============================================================================
## Conformal selection: TPR / FDR
## =============================================================================

outputs_cs <- lapply(seq_len(50), function(i) {
  files <- list.files(
    pattern = paste0("results_cs", i, "eps.*\\.rds$"),
    full.names = TRUE
  )
  if (length(files) == 0) return(NULL)
  
  lapply(files, function(file) {
    eps <- as.numeric(sub(".*eps([0-9.]+)\\.rds$", "\\1", basename(file)))
    readRDS(file) |>
      dplyr::mutate(seed = i, eps = eps)
  }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows()

res_cs_long <- outputs_cs |>
  tidyr::pivot_longer(
    cols = c(TPR_cs, FDR_cs),
    names_to = c("metric", "method"),
    names_pattern = "^(TPR|FDR)_(.*)$",
    values_to = "value"
  ) |>
  dplyr::filter(lfc != 0) |>
  dplyr::mutate(abs_lfc = abs(lfc))

# Per |lfc|
# Step 1: average within each seed
cs_per_seed_lfc <- res_cs_long |>
  dplyr::group_by(seed, eps, abs_lfc, metric) |>
  dplyr::summarise(val = mean(value), .groups = "drop")

# Step 2: mean and SD across seeds
res_cs_lfc <- cs_per_seed_lfc |>
  dplyr::group_by(eps, abs_lfc, metric) |>
  dplyr::summarise(
    mean_val = mean(val),
    sd_val = sd(val),
    .groups = "drop"
  )

ggplot2::ggplot(
  res_cs_lfc,
  ggplot2::aes(
    x = eps,
    y = mean_val,
    colour = metric,
    fill = metric,
    group = interaction(metric, abs_lfc),
    linetype = factor(abs_lfc)
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", colour = "black") +
  ggplot2::labs(
    x = expression(epsilon),
    y = NULL,
    colour = "Metric",
    fill = "Metric",
    linetype = "|LFC|"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_line(colour = "grey90"),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()
  ) 

# Overall (all |lfc| pooled)
cs_per_seed_all <- res_cs_long |>
  dplyr::group_by(seed, eps, metric) |>
  dplyr::summarise(val = mean(value), .groups = "drop")

res_cs <- cs_per_seed_all |>
  dplyr::group_by(eps, metric) |>
  dplyr::summarise(
    mean_val = mean(val),
    sd_val = sd(val),
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
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", colour = "black") +
  ggplot2::labs(
    x = expression(epsilon),
    y = NULL,
    colour = "Method / metric",
    fill = "Method / metric"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_line(colour = "grey90"),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()
  ) 


## =============================================================================
## cc clustering: prediction set sizes and conformal quantile
## =============================================================================

## ---- Settings --------------------------------------------------------------
cc_pattern <- "results_cc(\\d+)eps([0-9.]+)\\.(rds|RData)$"
data_dir <- getwd()

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
    prop_0        = counts[, "n_none"]  / total,
    prop_1        = counts[, "n_once"]  / total,
    prop_2        = counts[, "n_twice"] / total,
    n_inside_cc   = lengths(xl$cc_list),
    n_outside_cc  = lengths(xl$cc_list_out),
    n_cs          = lengths(xl$cs_list),
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

summary_df <- bind_rows(all_results) |>
  filter(lfc != 0) |>
  mutate(abs_lfc = abs(lfc)) |>
  arrange(epsilon, seed, gene)

## ---- Mean +/- SD across seeds ----------------------------------------------
# Step 1: average over genes within each seed (per epsilon and |lfc|)
per_seed <- summary_df |>
  group_by(epsilon, seed, abs_lfc) |>
  summarise(across(c(prop_0, prop_1, prop_2, conf_quantile), mean),
            .groups = "drop")

# Step 2: mean and SD across seeds
summary_by_lfc <- per_seed |>
  group_by(epsilon, abs_lfc) |>
  summarise(across(c(prop_0, prop_1, prop_2, conf_quantile),
                   list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
            n_seeds = n(), .groups = "drop")

## ---- Plot: mean proportion +/- SD ------------------------------------------
plot_df <- per_seed |>
  select(epsilon, seed, abs_lfc, starts_with("prop_")) |>
  tidyr::pivot_longer(starts_with("prop_"), names_to = "Set size",
               names_prefix = "prop_", values_to = "prop") |>
  mutate(`Set size` = factor(`Set size`)) |>
  group_by(epsilon, abs_lfc, `Set size`) |>
  summarise(mean = mean(prop), sd = sd(prop), .groups = "drop")

# Barplot version (not displayed in the article)
ggplot(plot_df, aes(x = factor(abs_lfc), y = mean, fill = `Set size`)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.9), width = 0.25) +
  facet_wrap(~ epsilon, labeller = label_both) +
  labs(x = "Absolute log-fold change",
       y = "Proportion of cells (mean \u00b1 SD across seeds)",
       fill = "Prediction set size") +
  theme_bw()

# Curves version (proportion vs epsilon, one curve per set size and |LFC|)
ggplot(plot_df, aes(x = epsilon, y = mean,
                    colour = `Set size`, fill = `Set size`,
                    linetype = factor(abs_lfc),
                    group = interaction(`Set size`, factor(abs_lfc)))) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
              alpha = 0.12, colour = NA) +
  geom_line() +
  geom_point(size = 1.5) +
  labs(x = expression(epsilon),
       y = "Proportion of cells (mean \u00b1 SD across seeds)",
       colour = "Set size", fill = "Set size",
       linetype = "|LFC|") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_line(colour = "grey90"),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()
  ) 

## ---- Plot: conformal quantile ----------------------------------------------
ggplot2::ggplot(
  summary_by_lfc,
  ggplot2::aes(
    x = epsilon,
    y = conf_quantile_mean,
    group = factor(abs_lfc),
    colour = factor(abs_lfc),
    fill = factor(abs_lfc)
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = conf_quantile_mean - conf_quantile_sd,
                 ymax = conf_quantile_mean + conf_quantile_sd),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(
    x = expression(epsilon),
    y = "Conformal quantile",
    colour = "|LFC|",
    fill = "|LFC|"
  ) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.05)) +
  ggplot2::theme_bw() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_line(colour = "grey90"),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()
  ) 