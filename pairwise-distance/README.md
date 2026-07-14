# Pairwise Batch Distance Analysis

Clean R Markdown workflow for comparing same-cell-type distances across
single-cell RNA-seq batches.

## Project Layout

- `Pairwise.Rmd`: main analysis report.
- `scripts/pairwise_distance_utils.R`: reusable loading, distance, summary, and plotting helpers.
- `data/`: local data folder, ignored by Git except for the README.
- `output/`: generated summary tables, ignored by Git except for the README.

## Data Layout

The report expects 10x-style batch folders under `data/panc8`:

```text
data/panc8/
  celseq/
  celseq2/
  fluidigmc1/
  smartseq2/
  indrop/
```

Each batch folder can include an optional `annotations.csv` file with:

```text
barcode,cell_type
```

Large data files are intentionally ignored by Git. Keep the analysis code on
GitHub and share the data separately if needed.

## R Packages

Install the required packages before rendering:

```r
install.packages(c(
  "rmarkdown",
  "knitr",
  "Matrix",
  "dplyr",
  "ggplot2"
))

install.packages("Seurat")
```

Rendering from the command line also requires Pandoc. RStudio usually includes
Pandoc, so knitting from RStudio is often the easiest path.

## Render

From the project folder:

```r
rmarkdown::render("Pairwise.Rmd")
```

If your data are somewhere else, pass the path as a parameter:

```r
rmarkdown::render(
  "Pairwise.Rmd",
  params = list(data_root = "/path/to/data/panc8")
)
```

For a faster preview, reduce the cells sampled per cell type in each batch:

```r
rmarkdown::render(
  "Pairwise.Rmd",
  params = list(
    data_root = "/path/to/data/panc8",
    max_cells_per_type_per_batch = 100
  )
)
```

## Main Parameters

Edit the YAML at the top of `Pairwise.Rmd` to change:

- `primary_batches`: the main two-batch comparison shown in detail.
- `batches`: all batches included in the all-pair overview.
- `distance_layer`: usually `"scale.data"` for scaled expression.
- `feature_mode`: `"variable"` for variable genes or `"all"` for all genes.
- `min_cells_per_type`: cell-type filtering threshold.
- `max_cells_per_type_per_batch`: sampling control for pairwise distance size.

## Outputs

When `write_outputs: true`, the report writes CSV summaries to `output/`.
Rendered HTML/PDF files and generated output tables are ignored so the
repository stays small and focused on reproducible code.
