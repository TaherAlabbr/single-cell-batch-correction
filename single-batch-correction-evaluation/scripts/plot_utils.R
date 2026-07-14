## Plot functions for technical and biological batch-correction diagnostics.

plot_technical_ratio_heatmap <- function(technical_summary) {
  ggplot(
    technical_summary,
    aes(x = batch_pair, y = cell_type, fill = R_tech)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = round(R_tech, 2)), size = 2.8) +
    facet_wrap(vars(method)) +
    scale_fill_gradient2(
      low = "#2B6CB0",
      mid = "white",
      high = "#C2410C",
      midpoint = 1,
      na.value = "grey90"
    ) +
    labs(
      title = "Technical ratio by cell type and batch pair",
      subtitle = "R_tech = median cross-batch same-cell-type distance / within-batch baseline",
      x = "Batch pair",
      y = "Cell type",
      fill = "R_tech"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_technical_deviation_heatmap <- function(technical_summary) {
  ggplot(
    technical_summary,
    aes(x = batch_pair, y = cell_type, fill = tech_deviation)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = round(tech_deviation, 2)), size = 2.8) +
    facet_wrap(vars(method)) +
    scale_fill_viridis_c(option = "C", direction = -1, na.value = "grey90") +
    ggplot2::coord_fixed(ratio = 1) +
    labs(
      title = "Technical deviation from ideal mixing",
      subtitle = "Lower abs(log(R_tech)) is better",
      x = "Batch pair",
      y = "Cell type",
      fill = "abs(log(R_tech))"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_technical_scatter <- function(technical_scatter) {
  ggplot(
    technical_scatter,
    aes(x = R_tech_before, y = R_tech_after, color = cell_type)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey25") +
    geom_point(size = 2.2, alpha = 0.85) +
    facet_wrap(vars(method)) +
    labs(
      title = "Before-vs-after technical ratio",
      subtitle = "Each point is one cell type and batch pair",
      x = "Uncorrected PCA R_tech",
      y = "Corrected R_tech",
      color = "Cell type"
    ) +
    theme_minimal(base_size = 11)
}

plot_cross_vs_within_overlay <- function(distance_table) {
  ggplot(
    distance_table,
    aes(x = normalized_distance, fill = distance_type, color = distance_type)
  ) +
    geom_density(alpha = 0.25, linewidth = 0.5, na.rm = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey35") +
    facet_grid(
      cell_type + batch_pair ~ method,
      scales = "free_y",
      switch = "y"
    ) +
    labs(
      title = "Cross-batch and within-batch normalized distance distributions",
      subtitle = "Distances are divided by the within-batch same-cell-type baseline inside each representation",
      x = "Normalized distance",
      y = "Density",
      fill = "Distance type",
      color = "Distance type"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "right",
      strip.placement = "outside",
      strip.text.x = element_text(size = 13),
      strip.text.y.left = element_text(angle = 0, size = 9, hjust = 1),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 13),
      panel.spacing.x = grid::unit(1.1, "lines"),
      panel.spacing.y = grid::unit(0.35, "lines")
    )
}

plot_alpha_beta_distribution <- function(biological_distances, beta_cell_type) {
  biological_distances %>%
    filter(target_cell_type == beta_cell_type) %>%
    ggplot(aes(x = relative_distance, color = method, fill = method)) +
    geom_density(alpha = 0.18, linewidth = 0.55, na.rm = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey35") +
    facet_wrap(vars(batch), scales = "free_y") +
    labs(
      title = paste("Relative alpha-to-", beta_cell_type, " distance distributions", sep = ""),
      subtitle = "Values mostly above 1 indicate separation from the alpha baseline",
      x = "Pairwise alpha-to-target distance / alpha-within-batch baseline",
      y = "Density",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal(base_size = 11)
}

plot_mean_relative_heatmap <- function(biological_summary) {
  ggplot(
    biological_summary,
    aes(x = batch, y = target_cell_type, fill = mean_relative_separation)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = round(mean_relative_separation, 2)), size = 2.8) +
    facet_wrap(vars(method)) +
    scale_fill_viridis_c(option = "D", na.value = "grey90") +
    labs(
      title = "Mean alpha-centered relative separation",
      subtitle = "Mean of pairwise relative distances, not a ratio of medians",
      x = "Batch",
      y = "Target cell type",
      fill = "Mean relative\nseparation"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())
}

plot_biological_scatter <- function(biological_scatter) {
  ggplot(
    biological_scatter,
    aes(
      x = mean_relative_before,
      y = mean_relative_after,
      color = target_cell_type
    )
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey25") +
    geom_point(size = 2.2, alpha = 0.85) +
    facet_wrap(vars(method)) +
    labs(
      title = "Biological preservation relative to uncorrected PCA",
      subtitle = "Movement toward 1 can indicate biological collapse",
      x = "Uncorrected PCA mean relative separation",
      y = "Corrected mean relative separation",
      color = "Target cell type"
    ) +
    theme_minimal(base_size = 11)
}

plot_signed_delta_heatmap <- function(biological_delta) {
  ggplot(
    biological_delta,
    aes(x = batch_pair, y = target_cell_type, fill = signed_delta)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = round(signed_delta, 2)), size = 2.8) +
    facet_wrap(vars(method)) +
    scale_fill_gradient2(
      low = "#2B6CB0",
      mid = "white",
      high = "#C2410C",
      midpoint = 0,
      na.value = "grey90"
    ) +
    labs(
      title = "Signed biological geometry difference across batches",
      subtitle = "Delta = mean relative separation in batch 1 minus batch 2",
      x = "Batch pair",
      y = "Target cell type",
      fill = "Signed delta"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
