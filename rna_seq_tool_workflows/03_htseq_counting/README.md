# HTSeq Counting

Counts reads per gene on the official HTSeq yeast example.

- input: yeast example SAM plus matching GTF
- tool: `htseq-count` from the Conda environment
- presentable outputs: counts table, summary, top-gene figure, special-counter figure

## Run

```bash
python3 rna_seq_tool_workflows/03_htseq_counting/03_htseq_counting.py
```
