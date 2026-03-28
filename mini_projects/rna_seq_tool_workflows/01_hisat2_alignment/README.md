# HISAT2 Alignment

This mini-project uses the official HISAT2 example files to run a real
alignment workflow end to end.

## What It Does

- downloads the example chr22 reference FASTA and paired-end read FASTA files
- builds a small HISAT2 index with `hisat2-build`
- runs `hisat2 --dta`
- saves the SAM output and raw HISAT2 stderr summary
- creates figures for concordant alignment breakdown and MAPQ distribution

## Real Input Files

- `mini_projects/rna_seq_tool_workflows/data/reference/22_20-21M.fa`
- `mini_projects/rna_seq_tool_workflows/data/reads/reads_1.fa`
- `mini_projects/rna_seq_tool_workflows/data/reads/reads_2.fa`

## Outputs

- `outputs/hisat2_alignment_summary.txt`
- `outputs/hisat2_alignment.stderr.txt`
- `outputs/aligned_reads.sam`
- `outputs/hisat2_alignment_breakdown.png`
- `outputs/hisat2_mapq_histogram.png`

## Run

```bash
python3 mini_projects/rna_seq_tool_workflows/01_hisat2_alignment/01_hisat2_alignment.py
```
