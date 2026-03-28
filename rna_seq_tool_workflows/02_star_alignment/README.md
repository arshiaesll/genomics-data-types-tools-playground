# STAR Alignment

Runs STAR on the same chr22 example reads used by the HISAT2 project.

- input: chr22 FASTA plus paired read FASTA files
- tool: bundled `STAR` binary under `rna_seq_tool_workflows/tools/STAR/`
- presentable outputs: STAR summary, mapping figure, splice-junction figure

## Run

```bash
python3 rna_seq_tool_workflows/02_star_alignment/02_star_alignment.py
```
