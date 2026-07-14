## Core evaluation functions extracted from the original R Markdown report.

safe_divide <- function(numerator, denominator) {
  out <- numerator / denominator
  bad <- !is.finite(out) | !is.finite(denominator) | denominator <= 0
  out[bad] <- NA_real_
  out
}

safe_abs_log <- function(x) {
  out <- abs(log(x))
  out[!is.finite(out) | x <= 0] <- NA_real_
  out
}

clean_labels <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

check_required_metadata <- function(seu, batch_col, cell_type_col, assay_name) {
  metadata <- seu[[]]
  missing_cols <- setdiff(c(batch_col, cell_type_col), colnames(metadata))

  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!assay_name %in% names(seu@assays)) {
    stop("Assay `", assay_name, "` was not found in `seu`.", call. = FALSE)
  }

  metadata[[batch_col]] <- clean_labels(metadata[[batch_col]])
  metadata[[cell_type_col]] <- clean_labels(metadata[[cell_type_col]])

  if (sum(!is.na(metadata[[batch_col]])) < 2) {
    stop("The batch column has fewer than two non-missing cells.", call. = FALSE)
  }

  if (sum(!is.na(metadata[[cell_type_col]])) < 2) {
    stop("The cell-type column has fewer than two non-missing cells.", call. = FALSE)
  }

  metadata
}

sample_vector <- function(x, max_n) {
  x <- as.character(x)

  if (is.null(max_n) || is.na(max_n) || length(x) <= max_n) {
    return(x)
  }

  sample(x, max_n)
}

effective_rank <- function(requested_rank, n_cells, n_features) {
  max(1L, min(as.integer(requested_rank), as.integer(n_cells) - 1L, as.integer(n_features)))
}

get_seurat_layer <- function(seu, assay_name, layer_name) {
  tryCatch(
    SeuratObject::LayerData(seu, assay = assay_name, layer = layer_name),
    error = function(e) {
      SeuratObject::GetAssayData(seu, assay = assay_name, slot = layer_name)
    }
  )
}

make_sce_from_seurat <- function(seu, assay_name) {
  logcounts_mat <- get_seurat_layer(seu, assay_name, "data")

  if (nrow(logcounts_mat) == 0 || ncol(logcounts_mat) == 0) {
    stop(
      "No normalized data layer was found. Run NormalizeData() before MNN methods.",
      call. = FALSE
    )
  }

  assays <- list(logcounts = logcounts_mat)

  counts_mat <- tryCatch(
    get_seurat_layer(seu, assay_name, "counts"),
    error = function(e) NULL
  )

  if (!is.null(counts_mat) && identical(dim(counts_mat), dim(logcounts_mat))) {
    assays$counts <- counts_mat[rownames(logcounts_mat), colnames(logcounts_mat), drop = FALSE]
  }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = assays)
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(
    seu[[]][colnames(sce), , drop = FALSE]
  )

  sce
}

align_embedding <- function(embedding, metadata, method_name) {
  embedding <- as.matrix(embedding)

  if (is.null(rownames(embedding))) {
    stop("Embedding for ", method_name, " has no cell rownames.", call. = FALSE)
  }

  if (anyDuplicated(rownames(embedding)) > 0) {
    stop("Embedding for ", method_name, " has duplicated cell rownames.", call. = FALSE)
  }

  missing_cells <- setdiff(rownames(metadata), rownames(embedding))

  if (length(missing_cells) > 0) {
    stop(
      "Embedding for ",
      method_name,
      " is missing ",
      length(missing_cells),
      " metadata cells. First missing cell: ",
      missing_cells[[1]],
      call. = FALSE
    )
  }

  embedding[rownames(metadata), , drop = FALSE]
}

run_uncorrected_pca_embedding <- function(seu,
                                          assay_name,
                                          n_dims,
                                          n_variable_features,
                                          seed) {
  set.seed(seed)
  SeuratObject::DefaultAssay(seu) <- assay_name

  features <- SeuratObject::VariableFeatures(seu)

  if (length(features) == 0) {
    seu <- Seurat::NormalizeData(seu, verbose = FALSE)
    seu <- Seurat::FindVariableFeatures(
      seu,
      selection.method = "vst",
      nfeatures = n_variable_features,
      verbose = FALSE
    )
    features <- SeuratObject::VariableFeatures(seu)
  }

  features <- head(features, n_variable_features)
  rank <- effective_rank(n_dims, ncol(seu), length(features))

  scale_data <- tryCatch(
    get_seurat_layer(seu, assay_name, "scale.data"),
    error = function(e) NULL
  )

  missing_scaled_features <- if (is.null(scale_data) || nrow(scale_data) == 0) {
    features
  } else {
    setdiff(features, rownames(scale_data))
  }

  if (length(missing_scaled_features) > 0) {
    seu <- Seurat::ScaleData(seu, features = features, verbose = FALSE)
  }

  seu <- Seurat::RunPCA(
    seu,
    features = features,
    npcs = rank,
    reduction.name = "pca_uncorrected",
    seed.use = seed,
    verbose = FALSE
  )

  embedding <- SeuratObject::Embeddings(seu, reduction = "pca_uncorrected")
  embedding <- embedding[, seq_len(rank), drop = FALSE]

  list(
    object = seu,
    embedding = embedding,
    features = features
  )
}

run_fastmnn_embedding <- function(seu,
                                  metadata,
                                  batch_col,
                                  assay_name,
                                  features,
                                  n_dims,
                                  seed,
                                  k = 20) {
  set.seed(seed)
  sce <- make_sce_from_seurat(seu, assay_name)
  batch <- factor(metadata[colnames(sce), batch_col])
  features <- intersect(features, rownames(sce))
  rank <- effective_rank(n_dims, ncol(sce), length(features))

  fastmnn_out <- batchelor::fastMNN(
    sce,
    batch = batch,
    subset.row = features,
    d = rank,
    k = k,
    assay.type = "logcounts",
    BPPARAM = BiocParallel::SerialParam()
  )

  embedding <- SingleCellExperiment::reducedDim(fastmnn_out, "corrected")
  rownames(embedding) <- colnames(fastmnn_out)
  colnames(embedding) <- paste0("FastMNN_", seq_len(ncol(embedding)))
  embedding
}

run_mnn_embedding <- function(seu,
                              metadata,
                              batch_col,
                              assay_name,
                              features,
                              n_dims,
                              seed,
                              k = 20,
                              sigma = 0.1) {
  set.seed(seed)
  sce <- make_sce_from_seurat(seu, assay_name)
  batch <- factor(metadata[colnames(sce), batch_col])
  features <- intersect(features, rownames(sce))

  mnn_out <- batchelor::mnnCorrect(
    sce,
    batch = batch,
    subset.row = features,
    k = k,
    sigma = sigma,
    assay.type = "logcounts",
    correct.all = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )

  corrected <- SummarizedExperiment::assay(mnn_out, "corrected")
  rank <- effective_rank(n_dims, ncol(corrected), nrow(corrected))

  pca <- BiocSingular::runPCA(
    t(corrected),
    rank = rank,
    center = TRUE,
    scale = FALSE,
    BSPARAM = BiocSingular::IrlbaParam(),
    BPPARAM = BiocParallel::SerialParam()
  )

  embedding <- pca$x
  rownames(embedding) <- colnames(corrected)
  colnames(embedding) <- paste0("MNN_PC_", seq_len(ncol(embedding)))
  embedding
}

run_seurat_cca_embedding <- function(seu,
                                     batch_col,
                                     assay_name,
                                     n_dims,
                                     n_variable_features,
                                     seed) {
  set.seed(seed)
  future::plan(future::sequential)
  SeuratObject::DefaultAssay(seu) <- assay_name

  object_list <- Seurat::SplitObject(seu, split.by = batch_col)
  object_list <- object_list[vapply(object_list, ncol, numeric(1)) > 1]

  if (length(object_list) < 2) {
    stop("Seurat CCA requires at least two batches with more than one cell.", call. = FALSE)
  }

  object_list <- lapply(object_list, function(obj) {
    SeuratObject::DefaultAssay(obj) <- assay_name
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    Seurat::FindVariableFeatures(
      obj,
      selection.method = "vst",
      nfeatures = n_variable_features,
      verbose = FALSE
    )
  })

  features <- Seurat::SelectIntegrationFeatures(
    object.list = object_list,
    nfeatures = n_variable_features
  )

  min_batch_cells <- min(vapply(object_list, ncol, numeric(1)))
  rank <- effective_rank(n_dims, min_batch_cells, length(features))
  dims <- seq_len(rank)

  object_list <- lapply(object_list, function(obj) {
    Seurat::ScaleData(obj, features = features, verbose = FALSE)
  })

  k_anchor <- min(5L, min_batch_cells - 1L)
  k_filter <- min(200L, min_batch_cells - 1L)
  k_score <- min(30L, min_batch_cells - 1L)
  k_weight <- min(100L, min_batch_cells - 1L)

  anchors <- Seurat::FindIntegrationAnchors(
    object.list = object_list,
    anchor.features = features,
    reduction = "cca",
    dims = dims,
    k.anchor = k_anchor,
    k.filter = k_filter,
    k.score = k_score,
    verbose = FALSE
  )

  integrated <- Seurat::IntegrateData(
    anchorset = anchors,
    new.assay.name = "integrated",
    dims = dims,
    k.weight = k_weight,
    verbose = FALSE
  )

  SeuratObject::DefaultAssay(integrated) <- "integrated"
  integrated <- Seurat::ScaleData(integrated, verbose = FALSE)
  integrated <- Seurat::RunPCA(
    integrated,
    npcs = rank,
    reduction.name = "integrated_pca",
    seed.use = seed,
    verbose = FALSE
  )

  embedding <- SeuratObject::Embeddings(integrated, reduction = "integrated_pca")
  embedding <- embedding[, seq_len(rank), drop = FALSE]
  colnames(embedding) <- paste0("SeuratCCA_PC_", seq_len(ncol(embedding)))

  list(
    object = integrated,
    embedding = embedding,
    features = features
  )
}

make_between_pairs <- function(cells_x, cells_y, max_pairs) {
  n_pairs_total <- length(cells_x) * length(cells_y)
  n_to_sample <- min(max_pairs, n_pairs_total)

  if (n_pairs_total <= max_pairs) {
    grid <- expand.grid(
      cell_x = cells_x,
      cell_y = cells_y,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    return(grid)
  }

  tibble(
    cell_x = sample(cells_x, n_to_sample, replace = TRUE),
    cell_y = sample(cells_y, n_to_sample, replace = TRUE)
  )
}

make_within_pairs <- function(cells, max_pairs) {
  if (length(cells) < 2) {
    return(tibble(cell_x = character(), cell_y = character()))
  }

  n_pairs_total <- length(cells) * (length(cells) - 1) / 2
  n_to_sample <- min(max_pairs, n_pairs_total)

  if (n_pairs_total <= max_pairs) {
    pair_matrix <- utils::combn(cells, 2)
    return(tibble(cell_x = pair_matrix[1, ], cell_y = pair_matrix[2, ]))
  }

  cell_x <- character()
  cell_y <- character()

  while (length(cell_x) < n_to_sample) {
    remaining <- n_to_sample - length(cell_x)
    candidate_x <- sample(cells, remaining * 2, replace = TRUE)
    candidate_y <- sample(cells, remaining * 2, replace = TRUE)
    keep <- candidate_x != candidate_y

    cell_x <- c(cell_x, candidate_x[keep])
    cell_y <- c(cell_y, candidate_y[keep])
  }

  tibble(
    cell_x = cell_x[seq_len(n_to_sample)],
    cell_y = cell_y[seq_len(n_to_sample)]
  )
}

paired_euclidean <- function(embedding, cell_x, cell_y) {
  x <- embedding[cell_x, , drop = FALSE]
  y <- embedding[cell_y, , drop = FALSE]
  sqrt(rowSums((x - y)^2))
}

sample_between_distances <- function(embedding, cells_x, cells_y, max_pairs) {
  pairs <- make_between_pairs(cells_x, cells_y, max_pairs)
  pairs$distance <- paired_euclidean(embedding, pairs$cell_x, pairs$cell_y)
  pairs
}

sample_within_distances <- function(embedding, cells, max_pairs) {
  pairs <- make_within_pairs(cells, max_pairs)

  if (nrow(pairs) == 0) {
    pairs$distance <- numeric()
    return(pairs)
  }

  pairs$distance <- paired_euclidean(embedding, pairs$cell_x, pairs$cell_y)
  pairs
}

cells_for_group <- function(metadata,
                            batch_col,
                            cell_type_col,
                            batch,
                            cell_type,
                            max_cells_per_group) {
  cells <- rownames(metadata)[
    metadata[[batch_col]] == batch &
      metadata[[cell_type_col]] == cell_type &
      !is.na(metadata[[batch_col]]) &
      !is.na(metadata[[cell_type_col]])
  ]

  sample_vector(cells, max_cells_per_group)
}

compute_technical_metrics_for_embedding <- function(embedding,
                                                    metadata,
                                                    batch_col,
                                                    cell_type_col,
                                                    method,
                                                    min_cells_per_group,
                                                    max_cells_per_group,
                                                    max_pairs_per_comparison,
                                                    seed) {
  set.seed(seed)
  metadata <- metadata[rownames(embedding), , drop = FALSE]
  metadata[[batch_col]] <- clean_labels(metadata[[batch_col]])
  metadata[[cell_type_col]] <- clean_labels(metadata[[cell_type_col]])

  batches <- sort(unique(stats::na.omit(metadata[[batch_col]])))
  cell_types <- sort(unique(stats::na.omit(metadata[[cell_type_col]])))
  batch_pairs <- utils::combn(batches, 2, simplify = FALSE)

  summary_rows <- list()
  distance_rows <- list()
  row_index <- 1L
  distance_index <- 1L

  for (batch_pair in batch_pairs) {
    batch_1 <- batch_pair[[1]]
    batch_2 <- batch_pair[[2]]
    batch_pair_label <- paste(batch_1, batch_2, sep = " vs ")

    for (cell_type in cell_types) {
      cells_1 <- cells_for_group(
        metadata,
        batch_col,
        cell_type_col,
        batch_1,
        cell_type,
        max_cells_per_group
      )
      cells_2 <- cells_for_group(
        metadata,
        batch_col,
        cell_type_col,
        batch_2,
        cell_type,
        max_cells_per_group
      )

      if (length(cells_1) < min_cells_per_group || length(cells_2) < min_cells_per_group) {
        message(
          "Skipping technical metric for ",
          method,
          ", ",
          cell_type,
          ", ",
          batch_pair_label,
          ": too few cells."
        )
        next
      }

      cross_distances <- sample_between_distances(
        embedding,
        cells_1,
        cells_2,
        max_pairs_per_comparison
      )
      within_1 <- sample_within_distances(embedding, cells_1, max_pairs_per_comparison)
      within_2 <- sample_within_distances(embedding, cells_2, max_pairs_per_comparison)

      if (nrow(within_1) == 0 || nrow(within_2) == 0) {
        message(
          "Skipping technical metric for ",
          method,
          ", ",
          cell_type,
          ", ",
          batch_pair_label,
          ": no within-batch distances."
        )
        next
      }

      within_baseline <- median(c(within_1$distance, within_2$distance), na.rm = TRUE)
      median_cross <- median(cross_distances$distance, na.rm = TRUE)
      r_tech <- safe_divide(median_cross, within_baseline)

      summary_rows[[row_index]] <- tibble(
        method = method,
        batch_1 = batch_1,
        batch_2 = batch_2,
        batch_pair = batch_pair_label,
        cell_type = cell_type,
        n_cells_batch_1 = length(cells_1),
        n_cells_batch_2 = length(cells_2),
        n_cross_distances = nrow(cross_distances),
        n_within_distances_batch_1 = nrow(within_1),
        n_within_distances_batch_2 = nrow(within_2),
        median_cross_batch_distance = median_cross,
        within_baseline = within_baseline,
        R_tech = r_tech,
        tech_deviation = safe_abs_log(r_tech)
      )
      row_index <- row_index + 1L

      distance_rows[[distance_index]] <- bind_rows(
        cross_distances %>%
          transmute(distance = distance) %>%
          mutate(distance_type = "Cross batch", source_batch = NA_character_),
        within_1 %>%
          transmute(distance = distance) %>%
          mutate(distance_type = "Within batch", source_batch = batch_1),
        within_2 %>%
          transmute(distance = distance) %>%
          mutate(distance_type = "Within batch", source_batch = batch_2)
      ) %>%
        mutate(
          method = method,
          batch_1 = batch_1,
          batch_2 = batch_2,
          batch_pair = batch_pair_label,
          cell_type = cell_type,
          within_baseline = within_baseline,
          normalized_distance = safe_divide(distance, within_baseline),
          .before = 1
        )
      distance_index <- distance_index + 1L
    }
  }

  list(
    summary = bind_rows(summary_rows),
    distances = bind_rows(distance_rows)
  )
}

compute_biological_metrics_for_embedding <- function(embedding,
                                                     metadata,
                                                     batch_col,
                                                     cell_type_col,
                                                     method,
                                                     reference_cell_type,
                                                     min_cells_per_group,
                                                     max_cells_per_group,
                                                     max_pairs_per_comparison,
                                                     seed) {
  set.seed(seed)
  metadata <- metadata[rownames(embedding), , drop = FALSE]
  metadata[[batch_col]] <- clean_labels(metadata[[batch_col]])
  metadata[[cell_type_col]] <- clean_labels(metadata[[cell_type_col]])

  batches <- sort(unique(stats::na.omit(metadata[[batch_col]])))

  summary_rows <- list()
  distance_rows <- list()
  row_index <- 1L
  distance_index <- 1L

  for (batch in batches) {
    batch_metadata <- metadata[
      metadata[[batch_col]] == batch &
        !is.na(metadata[[cell_type_col]]),
      ,
      drop = FALSE
    ]

    alpha_cells <- rownames(batch_metadata)[
      batch_metadata[[cell_type_col]] == reference_cell_type
    ]
    alpha_cells <- sample_vector(alpha_cells, max_cells_per_group)

    if (length(alpha_cells) < max(2L, min_cells_per_group)) {
      message(
        "Skipping biological metric for ",
        method,
        ", ",
        batch,
        ": fewer than ",
        max(2L, min_cells_per_group),
        " ",
        reference_cell_type,
        " cells."
      )
      next
    }

    alpha_within <- sample_within_distances(
      embedding,
      alpha_cells,
      max_pairs_per_comparison
    )
    alpha_baseline <- median(alpha_within$distance, na.rm = TRUE)

    if (!is.finite(alpha_baseline) || alpha_baseline <= 0) {
      message(
        "Skipping biological metric for ",
        method,
        ", ",
        batch,
        ": invalid alpha baseline."
      )
      next
    }

    target_cell_types <- sort(unique(stats::na.omit(batch_metadata[[cell_type_col]])))

    for (target_cell_type in target_cell_types) {
      target_cells <- rownames(batch_metadata)[
        batch_metadata[[cell_type_col]] == target_cell_type
      ]

      if (target_cell_type == reference_cell_type) {
        distances <- alpha_within
        n_target_cells <- length(alpha_cells)
      } else {
        target_cells <- sample_vector(target_cells, max_cells_per_group)

        if (length(target_cells) < min_cells_per_group) {
          message(
            "Skipping biological metric for ",
            method,
            ", ",
            batch,
            ", ",
            reference_cell_type,
            " vs ",
            target_cell_type,
            ": too few target cells."
          )
          next
        }

        distances <- sample_between_distances(
          embedding,
          alpha_cells,
          target_cells,
          max_pairs_per_comparison
        )
        n_target_cells <- length(target_cells)
      }

      distances <- distances %>%
        mutate(relative_distance = safe_divide(distance, alpha_baseline))

      summary_rows[[row_index]] <- tibble(
        method = method,
        batch = batch,
        reference_cell_type = reference_cell_type,
        target_cell_type = target_cell_type,
        n_reference_cells = length(alpha_cells),
        n_target_cells = n_target_cells,
        n_pairwise_distances = nrow(distances),
        alpha_baseline = alpha_baseline,
        mean_relative_separation = mean(distances$relative_distance, na.rm = TRUE),
        median_relative_separation = median(distances$relative_distance, na.rm = TRUE)
      )
      row_index <- row_index + 1L

      distance_rows[[distance_index]] <- distances %>%
        transmute(
          method = method,
          batch = batch,
          reference_cell_type = reference_cell_type,
          target_cell_type = target_cell_type,
          alpha_baseline = alpha_baseline,
          distance = distance,
          relative_distance = relative_distance
        )
      distance_index <- distance_index + 1L
    }
  }

  list(
    summary = bind_rows(summary_rows),
    distances = bind_rows(distance_rows)
  )
}

compute_biological_geometry_delta <- function(biological_summary) {
  split_groups <- biological_summary %>%
    group_by(method, target_cell_type) %>%
    group_split()

  map_dfr(split_groups, function(group_df) {
    batches <- sort(unique(group_df$batch))

    if (length(batches) < 2) {
      return(tibble())
    }

    map_dfr(utils::combn(batches, 2, simplify = FALSE), function(batch_pair) {
      batch_1 <- batch_pair[[1]]
      batch_2 <- batch_pair[[2]]

      value_1 <- group_df$mean_relative_separation[group_df$batch == batch_1][[1]]
      value_2 <- group_df$mean_relative_separation[group_df$batch == batch_2][[1]]

      tibble(
        method = group_df$method[[1]],
        target_cell_type = group_df$target_cell_type[[1]],
        batch_1 = batch_1,
        batch_2 = batch_2,
        batch_pair = paste(batch_1, batch_2, sep = " vs "),
        mean_relative_batch_1 = value_1,
        mean_relative_batch_2 = value_2,
        signed_delta = value_1 - value_2
      )
    })
  })
}

make_technical_scatter_data <- function(technical_summary,
                                        baseline_method = "Uncorrected_PCA") {
  baseline <- technical_summary %>%
    filter(method == baseline_method) %>%
    select(batch_pair, cell_type, R_tech_before = R_tech)

  technical_summary %>%
    filter(method != baseline_method) %>%
    inner_join(baseline, by = c("batch_pair", "cell_type")) %>%
    rename(R_tech_after = R_tech)
}

make_biological_scatter_data <- function(biological_summary,
                                         baseline_method = "Uncorrected_PCA") {
  baseline <- biological_summary %>%
    filter(method == baseline_method) %>%
    select(
      batch,
      target_cell_type,
      mean_relative_before = mean_relative_separation
    )

  biological_summary %>%
    filter(method != baseline_method) %>%
    inner_join(baseline, by = c("batch", "target_cell_type")) %>%
    rename(mean_relative_after = mean_relative_separation)
}

make_final_method_summary <- function(technical_summary,
                                      biological_summary,
                                      biological_delta,
                                      beta_cell_type,
                                      baseline_method = "Uncorrected_PCA") {
  technical_scores <- technical_summary %>%
    group_by(method) %>%
    summarise(
      technical_score = median(tech_deviation, na.rm = TRUE),
      median_R_tech = median(R_tech, na.rm = TRUE),
      n_technical_comparisons = n(),
      .groups = "drop"
    )

  baseline_bio <- biological_summary %>%
    filter(method == baseline_method) %>%
    select(
      batch,
      target_cell_type,
      mean_relative_before = mean_relative_separation
    )

  biological_preservation_scores <- biological_summary %>%
    inner_join(baseline_bio, by = c("batch", "target_cell_type")) %>%
    mutate(
      relative_to_uncorrected = safe_divide(
        mean_relative_separation,
        mean_relative_before
      ),
      biological_preservation_error = safe_abs_log(relative_to_uncorrected)
    ) %>%
    group_by(method) %>%
    summarise(
      biological_preservation_score = median(
        biological_preservation_error,
        na.rm = TRUE
      ),
      n_biological_comparisons = n(),
      .groups = "drop"
    )

  batch_agreement_scores <- biological_delta %>%
    group_by(method) %>%
    summarise(
      biological_batch_agreement_score = median(abs(signed_delta), na.rm = TRUE),
      n_biological_batch_pairs = n(),
      .groups = "drop"
    )

  alpha_beta_scores <- biological_summary %>%
    filter(target_cell_type == beta_cell_type) %>%
    group_by(method) %>%
    summarise(
      alpha_beta_collapse_risk = median(mean_relative_separation, na.rm = TRUE),
      n_beta_batches = n(),
      .groups = "drop"
    )

  technical_scores %>%
    full_join(biological_preservation_scores, by = "method") %>%
    full_join(batch_agreement_scores, by = "method") %>%
    full_join(alpha_beta_scores, by = "method") %>%
    arrange(technical_score)
}

filter_overlay_distances <- function(distance_table,
                                     selected_batch_pairs = NULL,
                                     selected_cell_types = NULL) {
  out <- distance_table

  if (!is.null(selected_batch_pairs)) {
    out <- out %>% filter(batch_pair %in% selected_batch_pairs)
  }

  if (!is.null(selected_cell_types)) {
    out <- out %>% filter(cell_type %in% selected_cell_types)
  }

  out
}

save_plot <- function(plot, filename, width = 10, height = 7, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(figures_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}
