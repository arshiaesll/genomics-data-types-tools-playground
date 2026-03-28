# Genomics Data Types & Tools Playground

Small genomics mini-projects for learning core file formats, RNA-seq summaries,
and common RNA-seq tools through runnable examples.

## Setup

```bash
conda env create -f environment.yml
conda activate genomics-tools
```

Most Python-based projects run directly from that environment. A few RNA-seq
tool workflows also rely on small bundled binaries kept under
`rna_seq_tool_workflows/tools/`.

## How To Run

List available projects:

```bash
python3 main.py --list
```

Run one project:

```bash
python3 main.py fasta_parser
python3 main.py normalization star_alignment
```

Run everything:

```bash
python3 main.py
```

## Project Areas

### Sequence Fundamentals

- `01_fasta_parser`: FASTA parsing, GC summary, composition plot
- `02_fastq_parser`: FASTQ parsing and quality summary
- `03_sequence_utilities`: reverse complement, motifs, GC windows, ORFs

### RNA-seq Fundamentals

- `01_count_matrix_summary`: summarize a real GEO count matrix
- `02_normalization`: CPM, FPKM, TPM comparisons
- `03_fold_change`: real two-condition fold-change example
- `04_volcano_ready_table`: volcano-ready results table and preview plot

### RNA-seq Tool Workflows

- `01_hisat2_alignment`
- `02_star_alignment`
- `03_htseq_counting`
- `04_stringtie_quantification`
- `05_kallisto_quantification`
- `06_bustools_single_cell`

## Outputs

- `outputs/presentable/` contains the files meant to be read or shared
- `outputs/technical/` contains raw workflow artifacts written locally

## Reproducibility

- Scripts download public data automatically when practical
- Some RNA-seq tool examples include small bundled binaries or test inputs under
  `rna_seq_tool_workflows/tools/`
- Curated plots and text summaries are committed in `outputs/presentable/`
