# bustools Single-Cell Workflow

Runs the `kallisto bus` plus `bustools` workflow on bundled single-cell test
reads.

- input: bundled transcript FASTA, GTF, and 10x-style FASTQ files
- tools: bundled `kallisto` and `bustools` binaries under `rna_seq_tool_workflows/tools/`
- presentable outputs: matrix summary, dense gene-by-cell tables, intermediate barcode summaries, and QC-style figures

## Run

```bash
python3 rna_seq_tool_workflows/06_bustools_single_cell/06_bustools_single_cell.py
```
