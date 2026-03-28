# bustools Single-Cell Workflow

This mini-project uses the official bundled kallisto single-cell test reads to
demonstrate the `kallisto bus` plus `bustools` preprocessing workflow.

## What It Does

- builds a transcriptome index from the bundled transcript FASTA
- runs `kallisto bus` with the `10xv2` technology preset
- generates an empirical barcode allowlist with `bustools allowlist`
- corrects barcode errors with `bustools correct`
- sorts BUS records with `bustools sort`
- summarizes barcode and UMI statistics with `bustools inspect`
- collapses transcript evidence into a barcode-by-gene matrix with
  `bustools count --genecounts`

## Inputs

- `mini_projects/rna_seq_tool_workflows/tools/kallisto/test/transcripts.fasta.gz`
- `mini_projects/rna_seq_tool_workflows/tools/kallisto/test/transcripts.gtf.gz`
- `mini_projects/rna_seq_tool_workflows/tools/kallisto/test/sc_reads_1.fastq.gz`
- `mini_projects/rna_seq_tool_workflows/tools/kallisto/test/sc_reads_2.fastq.gz`
- `mini_projects/rna_seq_tool_workflows/tools/kallisto/kallisto`
- `mini_projects/rna_seq_tool_workflows/tools/bustools_pkg/bustools/bustools`

## Outputs

- `outputs/transcripts.idx`
- `outputs/output.bus`
- `outputs/output.sorted.bus`
- `outputs/matrix.ec`
- `outputs/transcripts.txt`
- `outputs/barcode_allowlist.txt`
- `outputs/output.corrected.bus`
- `outputs/output.corrected.sorted.bus`
- `outputs/bustools_inspect.json`
- `outputs/transcript_to_gene.tsv`
- `outputs/cells_x_genes.barcodes.txt`
- `outputs/cells_x_genes.genes.txt`
- `outputs/cells_x_genes.mtx`
- `outputs/gene_by_cell_counts.tsv`
- `outputs/gene_by_top_cells_preview.tsv`
- `outputs/raw_bus_barcode_counts.tsv`
- `outputs/corrected_bus_barcode_counts.tsv`
- `outputs/intermediate_step_summary.txt`
- `outputs/bustools_matrix_summary.txt`
- `outputs/umi_per_barcode.png`
- `outputs/cells_per_gene_distribution.png`
- `outputs/raw_vs_corrected_barcode_rank.png`

## Run

```bash
python3 mini_projects/rna_seq_tool_workflows/06_bustools_single_cell/06_bustools_single_cell.py
```
