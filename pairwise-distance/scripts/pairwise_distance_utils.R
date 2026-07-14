# Utility functions for pairwise same-cell-type batch distance analysis.

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
                                 layer = "scale.data") {
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

get_distance_inputs <- function(seurat_object,
                                assay = "RNA",
                                layer = "scale.data",
                                features = NULL) {
  expression_matrix <- get_expression_layer(
    seurat_object = seurat_object,
    assay = assay,
    layer = layer
  )

  if (nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    stop(
      "The selected expression layer is empty. Run NormalizeData() and ",
      "ScaleData() before calculating distances, or choose another layer.",
      call. = FALSE
    )
  }

  if (is.null(features)) {
    features <- rownames(expression_matrix)
  }

  features <- intersect(features, rownames(expression_matrix))

  if (length(features) == 0) {
    stop("No requested features are available in the selected layer.", call. = FALSE)
  }

  metadata <- seurat_object[[]]
  validate_metadata(metadata)

  common_cells <- intersect(colnames(expression_matrix), rownames(metadata))

  if (length(common_cells) == 0) {
    stop(
      "No cells are shared between the metadata and expression layer.",
      call. = FALSE
    )
  }

  metadata <- metadata[common_cells, , drop = FALSE]
  metadata$batch <- as.character(metadata$batch)
  metadata$cell_type <- as.character(metadata$cell_type)

  list(
    expression_matrix = expression_matrix[features, common_cells, drop = FALSE],
    metadata = metadata,
    features = features
  )
}

get_shared_cell_types <- function(metadata, batch_a, batch_b) {
  cell_types_a <- unique(metadata$cell_type[metadata$batch == batch_a])
  cell_types_b <- unique(metadata$cell_type[metadata$batch == batch_b])
  shared_cell_types <- intersect(cell_types_a, cell_types_b)

  sort(shared_cell_types[!is.na(shared_cell_types) & shared_cell_types != ""])
}

sample_cells <- function(cells, max_cells = NULL) {
  if (is.null(max_cells) || length(cells) <= max_cells) {
    return(cells)
  }

  sample(cells, max_cells)
}

calculate_cross_distances <- function(feature_matrix, cells_x, cells_y) {
  x <- Matrix::t(feature_matrix[, cells_x, drop = FALSE])
  y <- Matrix::t(feature_matrix[, cells_y, drop = FALSE])

  squared_distances <-
    outer(Matrix::rowSums(x^2), Matrix::rowSums(y^2), FUN = "+") -
    2 * as.matrix(x %*% Matrix::t(y))

  sqrt(pmax(squared_distances, 0))
}

calculate_matching_celltype_distances <- function(seurat_object,
                                                  batch_a,
                                                  batch_b,
                                                  assay = "RNA",
                                                  layer = "scale.data",
                                                  features = NULL,
                                                  min_cells_per_type = 20,
                                                  max_cells_per_type_per_batch = 250,
                                                  seed = 123) {
  distance_inputs <- get_distance_inputs(
    seurat_object = seurat_object,
    assay = assay,
    layer = layer,
    features = features
  )

  expression_matrix <- distance_inputs$expression_matrix
  metadata <- distance_inputs$metadata
  available_batches <- sort(unique(metadata$batch))

  missing_batches <- setdiff(c(batch_a, batch_b), available_batches)

  if (length(missing_batches) > 0) {
    stop(
      "These batches were not found: ",
      paste(missing_batches, collapse = ", "),
      call. = FALSE
    )
  }

  shared_cell_types <- get_shared_cell_types(metadata, batch_a, batch_b)

  if (length(shared_cell_types) == 0) {
    stop(
      "No shared cell types were found between ",
      batch_a,
      " and ",
      batch_b,
      ".",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  distance_tables <- lapply(shared_cell_types, function(current_cell_type) {
    cells_a <- rownames(metadata)[
      metadata$batch == batch_a & metadata$cell_type == current_cell_type
    ]

    cells_b <- rownames(metadata)[
      metadata$batch == batch_b & metadata$cell_type == current_cell_type
    ]

    if (
      length(cells_a) < min_cells_per_type ||
        length(cells_b) < min_cells_per_type
    ) {
      return(NULL)
    }

    cells_a <- sample_cells(cells_a, max_cells_per_type_per_batch)
    cells_b <- sample_cells(cells_b, max_cells_per_type_per_batch)

    distances <- calculate_cross_distances(
      feature_matrix = expression_matrix,
      cells_x = cells_a,
      cells_y = cells_b
    )

    data.frame(
      batch_pair = paste(batch_a, "vs", batch_b),
      batch_a = batch_a,
      batch_b = batch_b,
      cell_type = current_cell_type,
      n_cells_a = length(cells_a),
      n_cells_b = length(cells_b),
      n_features = nrow(expression_matrix),
      distance = as.vector(distances)
    )
  })

  result <- dplyr::bind_rows(distance_tables)

  if (nrow(result) == 0) {
    stop(
      "No matching-cell-type distances were calculated for ",
      batch_a,
      " vs ",
      batch_b,
      ". Try lowering min_cells_per_type.",
      call. = FALSE
    )
  }

  result
}

calculate_all_batchpair_celltype_distances <- function(seurat_object,
                                                       batches = NULL,
                                                       assay = "RNA",
                                                       layer = "scale.data",
                                                       features = NULL,
                                                       min_cells_per_type = 20,
                                                       max_cells_per_type_per_batch = 250,
                                                       seed = 123) {
  distance_inputs <- get_distance_inputs(
    seurat_object = seurat_object,
    assay = assay,
    layer = layer,
    features = features
  )

  metadata <- distance_inputs$metadata
  available_batches <- sort(unique(metadata$batch))
  available_batches <- available_batches[!is.na(available_batches) & available_batches != ""]

  if (is.null(batches)) {
    batches <- available_batches
  }

  missing_batches <- setdiff(batches, available_batches)

  if (length(missing_batches) > 0) {
    stop(
      "These batches were not found: ",
      paste(missing_batches, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(batches) < 2) {
    stop("At least two batches are required.", call. = FALSE)
  }

  batch_pairs <- utils::combn(batches, 2, simplify = FALSE)

  distance_tables <- lapply(seq_along(batch_pairs), function(pair_index) {
    current_pair <- batch_pairs[[pair_index]]
    batch_a <- current_pair[1]
    batch_b <- current_pair[2]

    if (length(get_shared_cell_types(metadata, batch_a, batch_b)) == 0) {
      warning(
        "Skipping ",
        batch_a,
        " vs ",
        batch_b,
        ": no shared cell types.",
        call. = FALSE
      )
      return(NULL)
    }

    message("Calculating distances: ", batch_a, " vs ", batch_b)

    calculate_matching_celltype_distances(
      seurat_object = seurat_object,
      batch_a = batch_a,
      batch_b = batch_b,
      assay = assay,
      layer = layer,
      features = distance_inputs$features,
      min_cells_per_type = min_cells_per_type,
      max_cells_per_type_per_batch = max_cells_per_type_per_batch,
      seed = if (is.null(seed)) NULL else seed + pair_index - 1
    )
  })

  dplyr::bind_rows(distance_tables)
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

summarise_distances <- function(distance_table, group_cols) {
  distance_table |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      cells_in_batch_a = dplyr::first(n_cells_a),
      cells_in_batch_b = dplyr::first(n_cells_b),
      number_of_pairs = dplyr::n(),
      number_of_features = dplyr::first(n_features),
      mean_distance = mean(distance),
      median_distance = stats::median(distance),
      sd_distance = stats::sd(distance),
      q25 = stats::quantile(distance, 0.25),
      q75 = stats::quantile(distance, 0.75),
      .groups = "drop"
    )
}

make_batch_label <- function(batches) {
  paste(batches, collapse = " vs ")
}

plot_celltype_histograms <- function(distance_table, title, subtitle = NULL) {
  cell_type_medians <- distance_table |>
    dplyr::group_by(cell_type) |>
    dplyr::summarise(
      median_distance = stats::median(distance),
      .groups = "drop"
    )

  ggplot2::ggplot(distance_table, ggplot2::aes(x = distance)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 50,
      fill = "#77A7A0",
      color = "white",
      linewidth = 0.1
    ) +
    ggplot2::geom_vline(
      data = cell_type_medians,
      ggplot2::aes(xintercept = median_distance),
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~ cell_type, ncol = 3, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Euclidean distance",
      y = "Density"
    ) +
    ggplot2::theme_minimal()
}

plot_celltype_boxplot <- function(distance_table, title, subtitle = NULL) {
  ggplot2::ggplot(
    distance_table,
    ggplot2::aes(
      x = stats::reorder(cell_type, distance, FUN = stats::median),
      y = distance,
      fill = cell_type
    )
  ) +
    ggplot2::geom_boxplot(
      linewidth = 0.35,
      outlier.alpha = 0.25,
      outlier.size = 0.8
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Cell type",
      y = "Euclidean distance"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

plot_celltype_density <- function(distance_table, title, subtitle = NULL) {
  ggplot2::ggplot(
    distance_table,
    ggplot2::aes(x = distance, color = cell_type)
  ) +
    ggplot2::geom_density(linewidth = 0.35, na.rm = TRUE) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Euclidean distance",
      y = "Density",
      color = "Cell type"
    ) +
    ggplot2::theme_minimal()
}

plot_batchpair_density_by_celltype <- function(distance_table,
                                               title,
                                               subtitle) {
  ggplot2::ggplot(
    distance_table,
    ggplot2::aes(x = distance, color = batch_pair)
  ) +
    ggplot2::geom_density(linewidth = 0.35, alpha = 0.85, na.rm = TRUE) +
    ggplot2::facet_wrap(~ cell_type, ncol = 3, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Euclidean distance",
      y = "Density",
      color = "Batch pair"
    ) +
    ggplot2::theme_minimal()
}

plot_batchpair_matrix_density <- function(distance_table,
                                          title,
                                          subtitle) {
  plot_data <- distance_table |>
    dplyr::mutate(
      batch_a = factor(batch_a, levels = unique(batch_a)),
      batch_b = factor(batch_b, levels = unique(batch_b)),
      cell_type = factor(cell_type, levels = sort(unique(as.character(cell_type))))
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = distance, color = cell_type, group = cell_type)
  ) +
    ggplot2::geom_density(linewidth = 0.35, adjust = 1, na.rm = TRUE) +
    ggplot2::facet_grid(batch_a ~ batch_b, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Euclidean distance",
      y = "Density",
      color = "Cell type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 8),
      axis.text = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 9)
    )
}

plot_distance_heatmap <- function(summary_table,
                                  value_column,
                                  title,
                                  subtitle,
                                  fill_label) {
  if (!value_column %in% colnames(summary_table)) {
    stop("Column '", value_column, "' was not found.", call. = FALSE)
  }

  ggplot2::ggplot(
    summary_table,
    ggplot2::aes(
      x = batch_pair,
      y = cell_type,
      fill = .data[[value_column]]
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(.data[[value_column]], 1)),
      size = 3
    ) +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Batch pair",
      y = "Cell type",
      fill = fill_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}
