# StringTie Quantification

This mini-project runs StringTie on the official HTSeq yeast example
alignment after sorting it by genomic position.

## What It Does

- reuses the official yeast SAM and GTF inputs already downloaded for HTSeq
- creates a coordinate-sorted SAM for StringTie
- runs StringTie in reference-guided mode with `-G`
- writes an assembled transcript GTF and a gene abundance table
- visualizes top genes by TPM and the transcript length distribution

## Inputs

- `mini_projects/rna_seq_tool_workflows/data/reads/yeast_RNASeq_excerpt.sam`
- `mini_projects/rna_seq_tool_workflows/data/reference/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz`
- official StringTie macOS binary unpacked under `mini_projects/rna_seq_tool_workflows/tools/stringtie_pkg/`

## Outputs

- `outputs/yeast_RNASeq_excerpt.sorted.sam`
- `outputs/stringtie_assembled.gtf`
- `outputs/stringtie_gene_abundance.tsv`
- `outputs/stringtie_quantification_summary.txt`
- `outputs/stringtie_top_genes_tpm.png`
- `outputs/stringtie_transcript_lengths.png`

## Run

```bash
python3 mini_projects/rna_seq_tool_workflows/04_stringtie_quantification/04_stringtie_quantification.py
```
