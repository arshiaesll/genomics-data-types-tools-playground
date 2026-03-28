# kallisto Quantification

This mini-project uses the official kallisto bundled test data to demonstrate
transcriptome pseudoalignment and abundance estimation.

## What It Does

- uses the official kallisto macOS binary package
- builds a transcriptome index from the bundled transcript FASTA
- quantifies bundled paired-end reads with `kallisto quant`
- parses transcript-level TPM and estimated counts
- summarizes transcript TPM at the gene level using the bundled GTF

## Inputs

- `mini_projects/rna_seq_tool_workflows/tools/kallisto/test/transcripts.fasta.gz`
- `mini_projects/rna_seq_tool_workflows/tools/kallisto/test/transcripts.gtf.gz`
- `mini_projects/rna_seq_tool_workflows/tools/kallisto/test/reads_1.fastq.gz`
- `mini_projects/rna_seq_tool_workflows/tools/kallisto/test/reads_2.fastq.gz`
- `mini_projects/rna_seq_tool_workflows/tools/kallisto/kallisto`

## Outputs

- `outputs/transcripts.idx`
- `outputs/abundance.tsv`
- `outputs/run_info.json`
- `outputs/kallisto_quantification_summary.txt`
- `outputs/kallisto_top_transcripts_tpm.png`
- `outputs/kallisto_top_genes_tpm.png`

## Run

```bash
python3 mini_projects/rna_seq_tool_workflows/05_kallisto_quantification/05_kallisto_quantification.py
```
