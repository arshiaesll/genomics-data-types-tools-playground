# kallisto Quantification

Uses the bundled kallisto test dataset to demonstrate pseudoalignment-based
quantification.

- input: bundled transcript FASTA, GTF, and paired FASTQ files
- tool: bundled `kallisto` binary under `rna_seq_tool_workflows/tools/kallisto/`
- presentable outputs: abundance table, summary, transcript TPM figure, gene TPM figure

## Run

```bash
python3 rna_seq_tool_workflows/05_kallisto_quantification/05_kallisto_quantification.py
```
