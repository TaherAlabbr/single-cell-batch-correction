# single-cell-batch-correction

Evaluation of batch-effect correction methods for pancreatic single-cell RNA-seq data.

This repository contains the analysis code, diagnostic reports, and project poster used to compare batch-correction methods while checking two goals at the same time:

- **Technical alignment:** cells of the same biological type should mix across batches after correction.
- **Biological preservation:** real cell-type structure should not collapse after correction.

## Project Poster

[Open the full poster PDF](Batch_Effect_Correction_in_Single_Cell_Data.pdf)

![Batch Effect Correction in Single-Cell Data poster](assets/poster_preview.png)

## Repository Structure

```text
.
├── cellwise-distance/
├── pairwise-distance/
├── single-batch-correction-evaluation/
├── Batch_Effect_Correction_in_Single_Cell_Data.pdf
├── distance_based_diagnostics_for_batch_effects.pdf
├── results_after_correction.pdf
└── assets/
```

## What Each Part Is For

### `cellwise-distance/`

Contains a focused alpha-centered distance analysis. It compares distances from alpha cells to other cell types within each batch and normalizes those distances by the within-batch alpha-to-alpha baseline.

Use this folder when you want to study whether biological cell-type separation is preserved across batches.

Main files:

- `Cellwise.Rmd`: main report.
- `scripts/cellwise_distance_utils.R`: helper functions for loading data, calculating distances, summarizing results, and plotting.
- `data/`: expected local data location.
- `output/`: generated tables and figures.

### `pairwise-distance/`

Contains pairwise same-cell-type distance diagnostics across batch pairs. It compares cells with the same annotated cell type between two batches, then repeats that for all batch combinations.

Use this folder when you want to see how far matching cell types are from each other across batches.

Main files:

- `Pairwise.Rmd`: main report.
- `scripts/pairwise_distance_utils.R`: helper functions for pairwise distance calculation and plotting.
- `data/`: expected local data location.
- `output/`: generated tables and figures.

### `single-batch-correction-evaluation/`

Contains the full method-comparison workflow. It evaluates:

- Uncorrected PCA
- FastMNN
- MNN
- Seurat CCA

This report combines the technical and biological diagnostics into final method-level summaries.

Main files:

- `single_cell_batch_correction_evaluation.Rmd`: main report.
- `scripts/package_setup.R`: package checks.
- `scripts/data_loading.R`: Seurat object and 10x batch loading helpers.
- `scripts/evaluation_utils.R`: correction methods, distance metrics, and scoring functions.
- `scripts/plot_utils.R`: plotting functions.
- `data/`: expected local data location.
- `output/`: generated figures and tables.

## Additional Project Documents

- [Distance-Based Diagnostics for Batch Effects](distance_based_diagnostics_for_batch_effects.pdf): slide deck explaining the diagnostic framework, including technical mixing and biological preservation.
- [Results After Correction](results_after_correction.pdf): slide deck showing the correction results and method comparisons.
- [Project Poster](Batch_Effect_Correction_in_Single_Cell_Data.pdf): final poster summarizing the project motivation, methods, metrics, results, and conclusions.

## Data

The code expects 10x-style pancreatic single-cell data with batch folders such as:

```text
data/panc8/
  celseq/
  celseq2/
  fluidigmc1/
  smartseq2/
  indrop/
```

Each batch folder can include an `annotations.csv` file with cell barcode and cell-type labels.

Large data files are not included in this repository. Place data locally inside each analysis folder’s `data/` directory, or pass the data path as an R Markdown parameter.

## Running The Analyses

Open any `.Rmd` report in RStudio and click **Knit**.

You can also render from R:

```r
rmarkdown::render("cellwise-distance/Cellwise.Rmd")
rmarkdown::render("pairwise-distance/Pairwise.Rmd")
rmarkdown::render("single-batch-correction-evaluation/single_cell_batch_correction_evaluation.Rmd")
```

If your data are stored elsewhere, pass the path as a parameter. For example:

```r
rmarkdown::render(
  "single-batch-correction-evaluation/single_cell_batch_correction_evaluation.Rmd",
  params = list(tenx_data_root = "/path/to/data/panc8")
)
```

## Output

Each analysis folder writes generated tables and figures to its own `output/` directory. These outputs are intentionally separated from the code so the repository stays organized and reproducible.

## Project Summary

The main idea of the project is that batch correction should be judged by more than visual mixing. A good correction method should:

- reduce technical batch separation for matching cell types;
- preserve meaningful biological separation between different cell types;
- avoid overcorrecting the data so that true biological structure is erased.

The final evaluation found that Seurat CCA gave the strongest technical mixing in this analysis, while FastMNN and MNN provided useful comparison points for understanding the tradeoff between batch alignment and biological preservation.
