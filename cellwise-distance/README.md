# Cellwise Alpha Distance Analysis

Clean R Markdown workflow for alpha-centered cellwise distance comparisons
across single-cell RNA-seq batches.

## Project Layout

- `Cellwise.Rmd`: main analysis report.
- `scripts/cellwise_distance_utils.R`: reusable loading, distance, summary, and plotting helpers.
- `data/`: local data folder, ignored by Git except for the README.
- `output/`: generated CSV tables, ignored by Git except for the README.

## Data Layout

The report expects 10x-style batch folders under `data/panc8`:

```text
data/panc8/
  celseq/
  celseq2/
  fluidigmc1/
  indrop/
  smartseq2/
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
  "ggplot2",
  "patchwork"
))

install.packages("Seurat")
```

Rendering from the command line also requires Pandoc. RStudio usually includes
Pandoc, so knitting from RStudio is often the easiest path.

## Render

From the project folder:

```r
rmarkdown::render("Cellwise.Rmd")
```

If your data are somewhere else, pass the path as a parameter:

```r
rmarkdown::render(
  "Cellwise.Rmd",
  params = list(data_root = "/path/to/data/panc8")
)
```

## Main Parameters

Edit the YAML at the top of `Cellwise.Rmd` to change:

- `primary_batches`: the main two-batch comparison shown in detail.
- `all_batches`: all batches included in the overview heatmaps.
- `distance_layer`: usually `"data"` for normalized expression or `"scale.data"` for scaled expression.
- `feature_mode`: `"variable"` for variable genes or `"all"` for all genes.
- `min_cells_per_type` and `max_cells_per_type`: filtering and sampling controls.

## Outputs

When `write_outputs: true`, the report writes summary tables to `output/`.
Rendered HTML/PDF files and generated tables are ignored so the repository stays
small and focused on reproducible code.
