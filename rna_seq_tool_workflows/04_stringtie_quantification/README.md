# StringTie Quantification

Runs StringTie on the yeast example alignment used in the HTSeq project.

- input: yeast SAM plus matching GTF
- tool: bundled StringTie binary under `rna_seq_tool_workflows/tools/stringtie_pkg/`
- presentable outputs: abundance table, summary, TPM figure, transcript-length figure

## Run

```bash
python3 rna_seq_tool_workflows/04_stringtie_quantification/04_stringtie_quantification.py
```
