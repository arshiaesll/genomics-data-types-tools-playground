# RNA-seq Fundamentals

This section focuses on moving from sequencing reads to expression summaries
and simple differential-style comparisons.

Output convention:

- `outputs/presentable/` contains the reports, figures, and tables meant for
  reading or sharing.

## Mini-projects

### `01_count_matrix_summary/01_count_matrix_summary.py`

Description:

Parse a gene-by-sample count matrix and summarize expression per sample.

Objectives:

- Load a count matrix from a tabular file
- Summarize total counts per sample
- Identify the most highly expressed genes
- Generate simple count-based visualizations

Files we will use:

- `data/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `01_count_matrix_summary/outputs/count_matrix_summary.txt`
- `01_count_matrix_summary/outputs/library_sizes.png`
- `01_count_matrix_summary/outputs/detected_genes_per_sample.png`
- `01_count_matrix_summary/outputs/top_genes_overall.png`

Used public data:

- GEO accession `GSE164073`
- Official NCBI-generated RNA-seq raw count matrix download

### `02_normalization/02_normalization.py`

Description:

Implement and compare common RNA-seq normalization schemes.

Objectives:

- Compute CPM
- Compute FPKM
- Compute TPM
- Compare how normalization changes gene rankings

Files we will use:

- `data/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `data/human_gene_annotations.tsv.gz`
- `02_normalization/outputs/normalization_report.txt`
- `02_normalization/outputs/normalization_comparison.png`
- `02_normalization/outputs/tpm_two_sample_comparison.png`

Used public data:

- GEO accession `GSE164073`
- Official NCBI-generated RNA-seq raw count matrix
- Official GEO human annotation table with gene lengths

### `03_fold_change/03_fold_change.py`

Description:

Compare two conditions and compute simple log2 fold changes on toy data.

Objectives:

- Group samples by condition
- Compute mean expression by condition
- Compute log2 fold change
- Rank genes by change magnitude

Files we will use:

- `data/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `data/human_gene_annotations.tsv.gz`
- `03_fold_change/outputs/fold_change_report.txt`
- `03_fold_change/outputs/top_fold_changes.png`

Used public data:

- GEO accession `GSE164073`
- Real cornea mock vs SARS-CoV-2 sample comparison from the GEO series

### `04_volcano_ready_table/04_volcano_ready_table.py`

Description:

Prepare a volcano-plot-style table for downstream interpretation or plotting.

Objectives:

- Combine fold change with placeholder or toy significance values
- Build a clean output table
- Mark up- and down-regulated genes
- Export a volcano-ready results file

Files we will use:

- `data/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- `data/human_gene_annotations.tsv.gz`
- `04_volcano_ready_table/outputs/volcano_ready_results.tsv`
- `04_volcano_ready_table/outputs/volcano_preview.png`

Used public data:

- GEO accession `GSE164073`
- Real cornea mock vs SARS-CoV-2 comparison from the GEO series

## Suggested Build Order

1. Start with `01_count_matrix_summary/01_count_matrix_summary.py`
2. Build `02_normalization/02_normalization.py`
3. Continue with `03_fold_change/03_fold_change.py`
4. Finish with `04_volcano_ready_table/04_volcano_ready_table.py`
