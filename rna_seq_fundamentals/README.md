# RNA-seq Fundamentals

RNA-seq mini-projects built around a real GEO count matrix.

## Mini-projects

### `01_count_matrix_summary/01_count_matrix_summary.py`

- summarizes the real GEO count matrix `GSE164073`
- reports library sizes, detected genes, and top expressed genes
- writes a summary plus three plots

### `02_normalization/02_normalization.py`

- computes CPM, FPKM, and TPM
- compares how normalization changes rankings
- includes a two-sample TPM comparison figure

### `03_fold_change/03_fold_change.py`

- compares `Cornea_CoV2` vs `Cornea_mock`
- computes mean counts and log2 fold change
- writes a ranked report and fold-change figure

### `04_volcano_ready_table/04_volcano_ready_table.py`

- builds a volcano-ready table from the same cornea comparison
- adds approximate p-values and adjusted p-values
- writes a TSV plus a preview volcano plot

## Suggested Build Order

1. Start with `01_count_matrix_summary/01_count_matrix_summary.py`
2. Build `02_normalization/02_normalization.py`
3. Continue with `03_fold_change/03_fold_change.py`
4. Finish with `04_volcano_ready_table/04_volcano_ready_table.py`
