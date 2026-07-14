# Input helpers for loading Seurat objects or 10x batch folders.

find_loaded_seurat_objects <- function() {
  object_names <- ls(envir = .GlobalEnv)
  object_names[
    vapply(
      object_names,
      function(object_name) {
        inherits(get(object_name, envir = .GlobalEnv), "Seurat")
      },
      logical(1)
    )
  ]
}

load_10x_batches_from_root <- function(data_root,
                                       batch_names,
                                       batch_col,
                                       cell_type_col,
                                       min_cells = 3,
                                       min_features = 200) {
  if (is.null(data_root) || !dir.exists(data_root)) {
    stop(
      "The 10x data root does not exist: ",
      data_root,
      ". Set `tenx_data_root` to the folder containing the batch folders.",
      call. = FALSE
    )
  }

  dataset_dirs <- file.path(data_root, batch_names)
  names(dataset_dirs) <- batch_names

  missing_dirs <- dataset_dirs[!dir.exists(dataset_dirs)]

  if (length(missing_dirs) > 0) {
    stop(
      "These 10x batch directories do not exist:\n",
      paste(names(missing_dirs), missing_dirs, sep = ": ", collapse = "\n"),
      call. = FALSE
    )
  }

  objects <- lapply(names(dataset_dirs), function(batch_name) {
    counts <- Seurat::Read10X(data.dir = dataset_dirs[[batch_name]])

    object <- Seurat::CreateSeuratObject(
      counts = counts,
      project = batch_name,
      min.cells = min_cells,
      min.features = min_features
    )

    object[[batch_col]] <- batch_name

    annotation_file <- file.path(dataset_dirs[[batch_name]], "annotations.csv")

    if (file.exists(annotation_file)) {
      annotations <- read.csv(
        annotation_file,
        row.names = "barcode",
        stringsAsFactors = FALSE
      )

      if (!cell_type_col %in% colnames(annotations)) {
        stop(
          "Annotation file for ",
          batch_name,
          " does not contain column `",
          cell_type_col,
          "`.",
          call. = FALSE
        )
      }

      object[[cell_type_col]] <- annotations[colnames(object), cell_type_col]
    } else {
      object[[cell_type_col]] <- NA_character_
      warning("No annotations.csv found for batch: ", batch_name, call. = FALSE)
    }

    object
  })

  merged <- if (length(objects) == 1) {
    objects[[1]]
  } else {
    merge(
      objects[[1]],
      y = objects[-1],
      add.cell.ids = names(dataset_dirs)
    )
  }

  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    merged <- SeuratObject::JoinLayers(merged)
  }

  merged
}

load_input_seurat_object <- function(seurat_rds_path,
                                     seurat_object_name,
                                     tenx_data_root,
                                     tenx_batches,
                                     batch_col,
                                     cell_type_col) {
  if (!is.null(seurat_rds_path) && file.exists(seurat_rds_path)) {
    message("Loaded Seurat object from: ", seurat_rds_path)
    return(readRDS(seurat_rds_path))
  }

  if (
    !is.null(seurat_object_name) &&
      exists(seurat_object_name, envir = .GlobalEnv, inherits = FALSE)
  ) {
    object <- get(seurat_object_name, envir = .GlobalEnv)

    if (inherits(object, "Seurat")) {
      message("Using loaded Seurat object named `", seurat_object_name, "`.")
      return(object)
    }
  }

  if (exists("seu", inherits = TRUE) && inherits(seu, "Seurat")) {
    message("Using loaded Seurat object named `seu`.")
    return(seu)
  }

  if (!is.null(tenx_data_root) && dir.exists(tenx_data_root)) {
    object <- load_10x_batches_from_root(
      data_root = tenx_data_root,
      batch_names = tenx_batches,
      batch_col = batch_col,
      cell_type_col = cell_type_col
    )
    message("Loaded 10x batches from: ", tenx_data_root)
    return(object)
  }

  loaded_seurat_objects <- find_loaded_seurat_objects()
  available_message <- if (length(loaded_seurat_objects) == 0) {
    "No loaded Seurat objects were detected in .GlobalEnv."
  } else {
    paste(
      "Loaded Seurat objects detected:",
      paste(loaded_seurat_objects, collapse = ", ")
    )
  }

  stop(
    "No input Seurat object was found. Load your object as `seu`, ",
    "set `seurat_object_name`, set `seurat_rds`, or set `tenx_data_root` ",
    "to the folder containing batch directories. ",
    available_message,
    call. = FALSE
  )
}
