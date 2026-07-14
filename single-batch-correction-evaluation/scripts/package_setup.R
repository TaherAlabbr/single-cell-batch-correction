# Package checks and local aliases used by the evaluation scripts.

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "Matrix",
  "batchelor",
  "BiocSingular",
  "BiocParallel",
  "dplyr",
  "magrittr",
  "purrr",
  "ggplot2",
  "readr",
  "tibble",
  "future"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them in a fresh R session, then knit again.",
    call. = FALSE
  )
}

`%>%` <- magrittr::`%>%`

arrange <- dplyr::arrange
bind_rows <- dplyr::bind_rows
filter <- dplyr::filter
full_join <- dplyr::full_join
group_by <- dplyr::group_by
group_split <- dplyr::group_split
inner_join <- dplyr::inner_join
mutate <- dplyr::mutate
n <- dplyr::n
rename <- dplyr::rename
select <- dplyr::select
summarise <- dplyr::summarise
transmute <- dplyr::transmute

imap <- purrr::imap
map <- purrr::map
map_dfr <- purrr::map_dfr
tibble <- tibble::tibble

aes <- ggplot2::aes
element_blank <- ggplot2::element_blank
element_text <- ggplot2::element_text
facet_grid <- ggplot2::facet_grid
facet_wrap <- ggplot2::facet_wrap
geom_abline <- ggplot2::geom_abline
geom_density <- ggplot2::geom_density
geom_hline <- ggplot2::geom_hline
geom_point <- ggplot2::geom_point
geom_text <- ggplot2::geom_text
geom_tile <- ggplot2::geom_tile
geom_vline <- ggplot2::geom_vline
ggplot <- ggplot2::ggplot
labs <- ggplot2::labs
scale_fill_gradient2 <- ggplot2::scale_fill_gradient2
scale_fill_viridis_c <- ggplot2::scale_fill_viridis_c
theme <- ggplot2::theme
theme_minimal <- ggplot2::theme_minimal
vars <- ggplot2::vars
