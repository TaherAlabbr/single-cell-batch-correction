# Single-Cell Batch-Correction Evaluation

GitHub-ready R Markdown workflow for comparing single-cell batch-correction
methods with distance-based diagnostics.

## Methods Compared

- Uncorrected PCA
- FastMNN
- MNN
- Seurat CCA

## Project Layout

- `single_cell_batch_correction_evaluation.Rmd`: main report.
- `scripts/package_setup.R`: package checks and local aliases.
- `scripts/data_loading.R`: Seurat object and 10x batch loading helpers.
- `scripts/evaluation_utils.R`: correction, distance, and scoring functions.
- `scripts/plot_utils.R`: plotting functions.
- `data/`: local data folder, ignored by Git except for the README.
- `output/`: generated figures and tables, ignored by Git except for the README.

## Data Layout

The default settings expect 10x-style batch folders under `data/panc8`:

```text
data/panc8/
  celseq/
  celseq2/
  fluidigmc1/
  smartseq2/
  indrop/
```

Each batch folder can include an `annotations.csv` file with a `barcode` column
and a `cell_type` column. You can also bypass the 10x loader by passing a saved
Seurat object with `seurat_rds`.

## Required Packages

Install CRAN packages:

```r
install.packages(c(
  "rmarkdown",
  "knitr",
  "Seurat",
  "Matrix",
  "dplyr",
  "magrittr",
  "purrr",
  "ggplot2",
  "readr",
  "tibble",
  "future",
  "BiocManager"
))
```

Install Bioconductor packages:

```r
BiocManager::install(c(
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "batchelor",
  "BiocSingular",
  "BiocParallel"
))
```

Command-line rendering also requires Pandoc. RStudio usually includes Pandoc,
so knitting from RStudio is often the easiest path.

## Render

From the project folder:

```r
rmarkdown::render("single_cell_batch_correction_evaluation.Rmd")
```

If your data are somewhere else:

```r
rmarkdown::render(
  "single_cell_batch_correction_evaluation.Rmd",
  params = list(tenx_data_root = "/path/to/data/panc8")
)
```

If you already saved a Seurat object:

```r
rmarkdown::render(
  "single_cell_batch_correction_evaluation.Rmd",
  params = list(seurat_rds = "/path/to/seurat_object.rds")
)
```

For a faster preview, reduce the sampled cells and pairs:

```r
rmarkdown::render(
  "single_cell_batch_correction_evaluation.Rmd",
  params = list(
    tenx_data_root = "/path/to/data/panc8",
    max_cells_per_group = 100,
    max_pairs_per_comparison = 2000
  )
)
```

## Outputs

The report writes figures and tables to `output/`. Generated outputs are ignored
by Git so the repository stays small and reproducible.
