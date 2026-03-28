# STAR Alignment

This mini-project aligns the same official chr22 example reads used in the
HISAT2 mini-project, but with STAR.

## What It Does

- reuses the official example reference and paired-end reads
- builds a STAR genome index
- runs STAR on the paired reads
- parses `Log.final.out` for mapping summary metrics
- parses `SJ.out.tab` for splice-junction summaries
- creates figures for mapping percentages and splice-junction motifs

## Real Input Files

- `mini_projects/rna_seq_tool_workflows/data/reference/22_20-21M.fa`
- `mini_projects/rna_seq_tool_workflows/data/reads/reads_1.fa`
- `mini_projects/rna_seq_tool_workflows/data/reads/reads_2.fa`

## Outputs

- `outputs/star_alignment_summary.txt`
- `outputs/star_Log.final.out`
- `outputs/star_SJ.out.tab`
- `outputs/star_Aligned.out.sam`
- `outputs/star_mapping_summary.png`
- `outputs/star_splice_junction_motifs.png`

## Run

```bash
python3 mini_projects/rna_seq_tool_workflows/02_star_alignment/02_star_alignment.py
```
