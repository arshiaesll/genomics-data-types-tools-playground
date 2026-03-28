# HISAT2 Alignment

Uses the official HISAT2 chr22 example to run a small real alignment workflow.

- input: example chr22 FASTA plus paired read FASTA files
- tool: `hisat2` and `hisat2-build` from the Conda environment
- presentable outputs: alignment summary, concordant-alignment figure, MAPQ figure

## Run

```bash
python3 rna_seq_tool_workflows/01_hisat2_alignment/01_hisat2_alignment.py
```
