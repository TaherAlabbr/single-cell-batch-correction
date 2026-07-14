# Utility functions for alpha-centered cellwise distance analysis.

resolve_data_root <- function(data_root) {
  if (dir.exists(data_root)) {
    return(normalizePath(data_root, mustWork = TRUE))
  }

  stop(
    "Could not find data_root: ", data_root, "\n",
    "Expected a folder containing batch directories such as celseq, celseq2, ",
    "fluidigmc1, indrop, and smartseq2.",
    call. = FALSE
  )
}

make_dataset_dirs <- function(data_root, batches) {
  dataset_dirs <- file.path(data_root, batches)
  names(dataset_dirs) <- batches

  missing_dirs <- dataset_dirs[!dir.exists(dataset_dirs)]

  if (length(missing_dirs) > 0) {
    stop(
      "These batch folders were not found: ",
      paste(missing_dirs, collapse = ", "),
      call. = FALSE
    )
  }

  dataset_dirs
}

load_batch_annotations <- function(annotation_file, cell_names) {
  if (!file.exists(annotation_file)) {
    return(rep(NA_character_, length(cell_names)))
  }

  annotations <- read.csv(
    annotation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c("barcode", "cell_type")
  missing_columns <- setdiff(required_columns, colnames(annotations))

  if (length(missing_columns) > 0) {
    warning(
      "Skipping annotation file with missing columns: ",
      annotation_file,
      call. = FALSE
    )
    return(rep(NA_character_, length(cell_names)))
  }

  rownames(annotations) <- annotations$barcode
  annotations[cell_names, "cell_type"]
}

join_layers_if_available <- function(seurat_object) {
  if (!"JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    return(seurat_object)
  }

  tryCatch(
    SeuratObject::JoinLayers(seurat_object),
    error = function(error) seurat_object
  )
}

load_10x_batches <- function(dataset_dirs,
                             min_cells = 3,
                             min_features = 200) {
  objects <- lapply(names(dataset_dirs), function(batch_name) {
    counts <- Seurat::Read10X(data.dir = dataset_dirs[[batch_name]])

    object <- Seurat::CreateSeuratObject(
      counts = counts,
      project = batch_name,
      min.cells = min_cells,
      min.features = min_features
    )

    object$batch <- batch_name
    object$cell_type <- load_batch_annotations(
      annotation_file = file.path(dataset_dirs[[batch_name]], "annotations.csv"),
      cell_names = colnames(object)
    )

    object
  })

  seurat_object <- if (length(objects) == 1) {
    objects[[1]]
  } else {
    merge(objects[[1]], y = objects[-1], add.cell.ids = names(dataset_dirs))
  }

  join_layers_if_available(seurat_object)
}

get_expression_layer <- function(seurat_object,
                                 assay = "RNA",
                                 layer = "data") {
  tryCatch(
    SeuratObject::LayerData(seurat_object, assay = assay, layer = layer),
    error = function(error) {
      Seurat::GetAssayData(seurat_object, assay = assay, slot = layer)
    }
  )
}

select_distance_features <- function(seurat_object,
                                     expression_matrix,
                                     feature_mode = "variable",
                                     n_features = 2000) {
  feature_mode <- match.arg(feature_mode, choices = c("variable", "all"))

  if (identical(feature_mode, "all")) {
    return(rownames(expression_matrix))
  }

  variable_features <- Seurat::VariableFeatures(seurat_object)
  selected_features <- intersect(variable_features, rownames(expression_matrix))

  if (length(selected_features) == 0) {
    warning(
      "No variable features were available in the selected expression layer; ",
      "using all features instead.",
      call. = FALSE
    )
    return(rownames(expression_matrix))
  }

  head(selected_features, n_features)
}

sample_cells <- function(cells, max_cells = NULL) {
  if (is.null(max_cells) || length(cells) <= max_cells) {
    return(cells)
  }

  sample(cells, max_cells)
}

safe_divide <- function(numerator, denominator) {
  result <- numerator / denominator
  denominator <- rep(denominator, length.out = length(result))
  invalid <- !is.finite(denominator) | denominator <= 0
  result[invalid] <- NA_real_
  result
}

calculate_cross_distances <- function(feature_matrix, cells_x, cells_y) {
  x <- Matrix::t(feature_matrix[, cells_x, drop = FALSE])
  y <- Matrix::t(feature_matrix[, cells_y, drop = FALSE])

  squared_distances <-
    outer(Matrix::rowSums(x^2), Matrix::rowSums(y^2), FUN = "+") -
    2 * as.matrix(x %*% Matrix::t(y))

  sqrt(pmax(squared_distances, 0))
}

calculate_within_distances <- function(feature_matrix, cells) {
  if (length(cells) < 2) {
    return(numeric(0))
  }

  x <- Matrix::t(feature_matrix[, cells, drop = FALSE])

  squared_distances <-
    outer(Matrix::rowSums(x^2), Matrix::rowSums(x^2), FUN = "+") -
    2 * as.matrix(x %*% Matrix::t(x))

  distances <- sqrt(pmax(squared_distances, 0))
  distances[upper.tri(distances)]
}

validate_metadata <- function(metadata) {
  required_columns <- c("batch", "cell_type")
  missing_columns <- setdiff(required_columns, colnames(metadata))

  if (length(missing_columns) > 0) {
    stop(
      "The Seurat object is missing required metadata columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

calculate_alpha_to_other_distances <- function(seurat_object,
                                               batches,
                                               alpha_label = "alpha",
                                               assay = "RNA",
                                               layer = "data",
                                               features = NULL,
                                               min_cells_per_type = 20,
                                               max_cells_per_type = 500,
                                               seed = 123) {
  metadata <- seurat_object[[]]
  validate_metadata(metadata)

  expression_matrix <- get_expression_layer(
    seurat_object = seurat_object,
    assay = assay,
    layer = layer
  )

  if (is.null(features)) {
    features <- rownames(expression_matrix)
  }

  features <- intersect(features, rownames(expression_matrix))

  if (length(features) == 0) {
    stop("No requested features are available in the selected layer.", call. = FALSE)
  }

  common_cells <- intersect(colnames(expression_matrix), rownames(metadata))
  expression_matrix <- expression_matrix[features, common_cells, drop = FALSE]
  metadata <- metadata[common_cells, , drop = FALSE]
  metadata$batch <- as.character(metadata$batch)
  metadata$cell_type <- as.character(metadata$cell_type)

  missing_batches <- setdiff(batches, unique(metadata$batch))

  if (length(missing_batches) > 0) {
    stop(
      "These batches were not found in the Seurat object: ",
      paste(missing_batches, collapse = ", "),
      call. = FALSE
    )
  }

  set.seed(seed)

  distance_tables <- list()
  baseline_tables <- list()
  result_index <- 1

  for (current_batch in batches) {
    current_metadata <- metadata[
      metadata$batch == current_batch & !is.na(metadata$cell_type),
      ,
      drop = FALSE
    ]

    alpha_cells <- rownames(current_metadata)[
      current_metadata$cell_type == alpha_label
    ]

    if (length(alpha_cells) < 2) {
      stop(
        "Batch '", current_batch, "' has fewer than two ",
        alpha_label, " cells.",
        call. = FALSE
      )
    }

    alpha_cells <- sample_cells(alpha_cells, max_cells_per_type)
    alpha_distances <- calculate_within_distances(expression_matrix, alpha_cells)
    alpha_baseline <- stats::median(alpha_distances)

    if (!is.finite(alpha_baseline) || alpha_baseline <= 0) {
      stop(
        "The within-batch alpha baseline is zero or invalid for batch '",
        current_batch,
        "'.",
        call. = FALSE
      )
    }

    baseline_tables[[current_batch]] <- data.frame(
      batch = current_batch,
      n_alpha_cells = length(alpha_cells),
      n_alpha_pairs = length(alpha_distances),
      alpha_alpha_mean = mean(alpha_distances),
      alpha_alpha_median = alpha_baseline
    )

    other_cell_types <- setdiff(
      sort(unique(current_metadata$cell_type)),
      alpha_label
    )

    for (other_type in other_cell_types) {
      other_cells <- rownames(current_metadata)[
        current_metadata$cell_type == other_type
      ]

      if (length(other_cells) < min_cells_per_type) {
        next
      }

      other_cells <- sample_cells(other_cells, max_cells_per_type)
      distance_values <- as.vector(
        calculate_cross_distances(expression_matrix, alpha_cells, other_cells)
      )

      distance_tables[[result_index]] <- data.frame(
        batch = current_batch,
        anchor_cell_type = alpha_label,
        other_cell_type = other_type,
        comparison = paste(alpha_label, "vs", other_type),
        n_alpha_cells = length(alpha_cells),
        n_other_cells = length(other_cells),
        n_features = nrow(expression_matrix),
        alpha_baseline_median = alpha_baseline,
        distance = distance_values,
        relative_distance = safe_divide(distance_values, alpha_baseline)
      )

      result_index <- result_index + 1
    }
  }

  if (length(distance_tables) == 0) {
    stop(
      "No alpha-to-other-cell-type distances were calculated. ",
      "Try lowering min_cells_per_type.",
      call. = FALSE
    )
  }

  list(
    distances = dplyr::bind_rows(distance_tables),
    alpha_baselines = dplyr::bind_rows(baseline_tables)
  )
}

summarise_alpha_distances <- function(distance_table) {
  distance_table |>
    dplyr::group_by(batch, comparison, other_cell_type) |>
    dplyr::summarise(
      n_alpha_cells = dplyr::first(n_alpha_cells),
      n_other_cells = dplyr::first(n_other_cells),
      n_pairs = dplyr::n(),
      n_features = dplyr::first(n_features),
      alpha_baseline_median = dplyr::first(alpha_baseline_median),
      mean_distance = mean(distance),
      median_distance = stats::median(distance),
      mean_relative_distance = mean(relative_distance, na.rm = TRUE),
      median_relative_distance = stats::median(relative_distance, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(other_cell_type, batch)
}

filter_shared_cell_types <- function(distance_table, batches) {
  distance_table |>
    dplyr::filter(batch %in% batches) |>
    dplyr::group_by(other_cell_type) |>
    dplyr::filter(dplyr::n_distinct(batch) == length(batches)) |>
    dplyr::ungroup()
}

sample_distance_table <- function(distance_table,
                                  group_cols,
                                  n_per_group = 5000,
                                  seed = 123) {
  if (is.null(n_per_group)) {
    return(distance_table)
  }

  set.seed(seed)

  distance_table |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(
      ~ dplyr::slice_sample(.x, n = min(n_per_group, nrow(.x)))
    ) |>
    dplyr::ungroup()
}

pairwise_batch_differences <- function(summary_table, metric_col) {
  batches <- sort(unique(summary_table$batch))

  if (length(batches) < 2) {
    return(data.frame())
  }

  batch_pairs <- utils::combn(batches, 2, simplify = FALSE)

  dplyr::bind_rows(lapply(batch_pairs, function(current_pair) {
    batch_a <- current_pair[1]
    batch_b <- current_pair[2]

    summary_a <- summary_table |>
      dplyr::filter(batch == batch_a) |>
      dplyr::select(
        other_cell_type,
        metric_a = dplyr::all_of(metric_col)
      )

    summary_b <- summary_table |>
      dplyr::filter(batch == batch_b) |>
      dplyr::select(
        other_cell_type,
        metric_b = dplyr::all_of(metric_col)
      )

    dplyr::inner_join(summary_a, summary_b, by = "other_cell_type") |>
      dplyr::mutate(
        batch_a = batch_a,
        batch_b = batch_b,
        batch_pair = paste(batch_a, "vs", batch_b),
        metric = metric_col,
        signed_difference = metric_a - metric_b,
        absolute_difference = abs(signed_difference)
      )
  }))
}

make_batch_label <- function(batches) {
  paste(batches, collapse = " vs ")
}

plot_relative_distance_distribution <- function(distance_table,
                                                title,
                                                subtitle) {
  ggplot2::ggplot(
    distance_table,
    ggplot2::aes(x = relative_distance, fill = batch)
  ) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.35, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::facet_wrap(~ comparison, ncol = 3, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Alpha-to-cell-type distance / within-batch alpha-to-alpha median",
      y = "Density",
      fill = "Batch"
    ) +
    ggplot2::theme_minimal()
}

plot_raw_distance_distribution <- function(distance_table,
                                           title,
                                           subtitle) {
  ggplot2::ggplot(
    distance_table,
    ggplot2::aes(x = distance, fill = batch)
  ) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.35, na.rm = TRUE) +
    ggplot2::facet_wrap(~ comparison, ncol = 3, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Euclidean distance",
      y = "Density",
      fill = "Batch"
    ) +
    ggplot2::theme_minimal()
}

plot_mean_distance_summary <- function(summary_table,
                                       x_metric,
                                       title,
                                       subtitle,
                                       x_label,
                                       add_baseline = FALSE) {
  plot <- ggplot2::ggplot(
    summary_table,
    ggplot2::aes(
      x = .data[[x_metric]],
      y = stats::reorder(other_cell_type, .data[[x_metric]]),
      color = batch
    )
  ) +
    ggplot2::geom_point(
      size = 3,
      position = ggplot2::position_dodge(width = 0.45)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = "Other cell type",
      color = "Batch"
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(add_baseline)) {
    plot <- plot + ggplot2::geom_vline(xintercept = 1, linetype = "dashed")
  }

  plot
}

plot_raw_vs_relative_summary <- function(distance_table,
                                         title,
                                         subtitle) {
  summary_table <- summarise_alpha_distances(distance_table)

  celltype_order <- summary_table |>
    dplyr::group_by(other_cell_type) |>
    dplyr::summarise(
      average_relative_distance = mean(mean_relative_distance, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(average_relative_distance) |>
    dplyr::pull(other_cell_type)

  summary_table <- summary_table |>
    dplyr::mutate(
      other_cell_type = factor(other_cell_type, levels = celltype_order)
    )

  raw_plot <- ggplot2::ggplot(
    summary_table,
    ggplot2::aes(x = mean_distance, y = other_cell_type)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = other_cell_type),
      color = "grey70",
      linewidth = 0.6
    ) +
    ggplot2::geom_point(ggplot2::aes(color = batch), size = 3) +
    ggplot2::labs(
      title = "Raw distances",
      x = "Mean Euclidean distance",
      y = "Other cell type",
      color = "Batch"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  relative_plot <- ggplot2::ggplot(
    summary_table,
    ggplot2::aes(x = mean_relative_distance, y = other_cell_type)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::geom_line(
      ggplot2::aes(group = other_cell_type),
      color = "grey70",
      linewidth = 0.6
    ) +
    ggplot2::geom_point(ggplot2::aes(color = batch), size = 3) +
    ggplot2::labs(
      title = "Relative distances",
      x = "Mean relative distance",
      y = NULL,
      color = "Batch"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  raw_plot + relative_plot +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = paste(
        "Relative distance = alpha-to-cell-type distance /",
        "within-batch alpha-to-alpha median"
      )
    ) &
    ggplot2::theme(legend.position = "bottom")
}

plot_pairwise_difference_heatmap <- function(difference_table,
                                             title,
                                             subtitle,
                                             fill_label = "Signed difference") {
  ggplot2::ggplot(
    difference_table,
    ggplot2::aes(
      x = batch_pair,
      y = other_cell_type,
      fill = signed_difference
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(signed_difference, 2)),
      size = 3
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      na.value = "grey90"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Batch pair",
      y = "Other cell type",
      fill = fill_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
