# HTSeq Counting

This mini-project uses the official HTSeq yeast example dataset to produce a
real multi-gene count table.

## What It Does

- downloads the official HTSeq example SAM and yeast GTF if needed
- runs `htseq-count` on the real example alignment
- summarizes the resulting multi-gene count table
- visualizes top-count genes and HTSeq special counters

## Inputs

- `mini_projects/rna_seq_tool_workflows/data/reads/yeast_RNASeq_excerpt.sam`
- `mini_projects/rna_seq_tool_workflows/data/reference/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz`

## Outputs

- `outputs/htseq_counts.tsv`
- `outputs/htseq_counting_summary.txt`
- `outputs/htseq_top_genes.png`
- `outputs/htseq_special_counters.png`

## Run

```bash
python3 mini_projects/rna_seq_tool_workflows/03_htseq_counting/03_htseq_counting.py
```
