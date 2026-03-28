# Genomics Data Types & Tools Playground

This is a portfolio-style genomics learning repo built around small,
understandable mini-projects for the main data types, file formats, and tools
used in bioinformatics.

## Reproducibility

- Public data is downloaded automatically by the scripts whenever practical.
- Presentation-ready outputs that are useful to browse are committed under each
  mini-project's `outputs/presentable/` folder.
- Large technical artifacts such as `SAM`, `BUS`, indexes, logs, and caches are
  written to `outputs/technical/` locally and are not intended to be part of
  the public repo history.
- Small bundled binaries and test files used by the RNA-seq tool workflows are
  kept under `mini_projects/rna_seq_tool_workflows/tools/` so those examples are
  reproducible on a matching environment.

## Environment Setup

Using Conda:

```bash
conda env create -f environment.yml
conda activate genomics-tools
```

Using pip:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Current Python dependency:

- `matplotlib` for plot generation in the FASTA mini-project

## Project Goal

Build hands-on examples for:

- DNA sequence data
- RNA expression data
- Proteomics data
- Common file formats
- Common analysis steps
- Practical genomics tools

Each mini-project should produce something visible and explainable, so the repo becomes both a study tool and a showcase.

## Recommended Roadmap

### Sequence Fundamentals

Start with the simplest objects in genomics: nucleotide strings and their file formats.

Mini-projects:

1. FASTA parser and sequence statistics
2. FASTQ parser and quality score summary
3. Reverse complement, GC content, motif search

Why first:

- Easy to implement in Python
- Builds comfort with real genomics formats
- Creates a foundation for later alignment and expression work

### RNA-seq Fundamentals

Move from raw reads to gene expression tables.

Mini-projects:

1. Parse a count matrix and summarize expression per sample
2. Implement CPM, FPKM, and TPM normalization on toy data
3. Compare two conditions and compute log2 fold change
4. Make a volcano-plot-ready output table

Tools to study alongside:

- HISAT2
- STAR
- HTSeq
- StringTie
- Kallisto
- bustools

### RNA-seq Tool Workflows

Build one mini-project per major RNA-seq tool so the repository includes
workflow-oriented examples in addition to table-based downstream analysis.

Mini-projects:

1. HISAT2 alignment on real example reads
2. STAR alignment on the same real example reads
3. HTSeq gene counting on aligned reads
4. StringTie transcript assembly and quantification
5. Kallisto pseudoalignment-based quantification
6. bustools single-cell preprocessing and count matrix generation

### Differential Expression

Focus on the interpretation side.

Mini-projects:

1. Simulate two groups and compute log2 fold change
2. Parse DE-style result tables
3. Filter genes by adjusted p-value and effect size
4. Generate top-hit reports

Tools to study:

- DESeq2
- edgeR

### Variant and Alignment Concepts

After sequence basics, add the workflow around mapped reads.

Mini-projects:

1. SAM/BAM field explainer on toy examples
2. CIGAR-string parser
3. Coverage calculator from simplified alignments
4. Variant table parser for VCF basics

### Proteomics Basics

Cover the high-level data flow from spectra to proteins.

Mini-projects:

1. Parse a toy peptide identification table
2. Summarize peptide-to-protein mapping
3. Compare protein abundance across conditions
4. Build a simple report for identified proteins

## Suggested Learning Order

1. Sequence strings and FASTA
2. FASTQ and quality scores
3. Read-level utilities
4. Count matrices and normalization
5. Differential expression
6. Alignment concepts
7. Variant formats
8. Proteomics tables

## Notes to Tighten Conceptually

A few places in your notes are worth sharpening as you build:

- The format name is `FASTQ`, not `FASTAQ`.
- In FASTQ, quality scores are encoded as ASCII characters that map to Phred scores. The common modern convention is Phred+33.
- NGS does not always assemble reads into a whole genome from scratch. Often reads are aligned to a reference genome, but de novo assembly is also possible.
- FDR refers to the expected proportion of false positives among the results called significant, not incorrect base pairs.
- RNA-seq expression values may start as raw counts, while TPM/FPKM/CPM are different normalized representations.

## Current Mini-Projects

### Mini-project 1: FASTA and FASTQ Basics

Implemented in Python:

- Parse FASTA records
- Parse FASTQ records
- Compute sequence length and GC content
- Compute reverse complement
- Convert FASTQ quality symbols to Phred scores
- Summarize average read quality

Run it with:

```bash
python3 main.py
```

## Current Project Structure

Sequence fundamentals now lives in:

- `mini_projects/sequence_fundamentals/01_fasta_parser/`
- `mini_projects/sequence_fundamentals/02_fastq_parser/`
- `mini_projects/sequence_fundamentals/03_sequence_utilities/`
- `mini_projects/sequence_fundamentals/data/`

RNA-seq fundamentals is now scaffolded in:

- `mini_projects/rna_seq_fundamentals/01_count_matrix_summary/`
- `mini_projects/rna_seq_fundamentals/02_normalization/`
- `mini_projects/rna_seq_fundamentals/03_fold_change/`
- `mini_projects/rna_seq_fundamentals/04_volcano_ready_table/`
- `mini_projects/rna_seq_fundamentals/data/`

RNA-seq tool workflows are now scaffolded in:

- `mini_projects/rna_seq_tool_workflows/01_hisat2_alignment/`
- `mini_projects/rna_seq_tool_workflows/02_star_alignment/`
- `mini_projects/rna_seq_tool_workflows/03_htseq_counting/`
- `mini_projects/rna_seq_tool_workflows/04_stringtie_quantification/`
- `mini_projects/rna_seq_tool_workflows/05_kallisto_quantification/`
- `mini_projects/rna_seq_tool_workflows/06_bustools_single_cell/`
- `mini_projects/rna_seq_tool_workflows/data/`

## What To Build Next

The best next step after the current sequence work is:

1. Build `01_count_matrix_summary` on a real GEO RNA-seq count matrix
2. Add normalization functions for CPM, FPKM, and TPM
3. Prepare a simple fold-change workflow for two conditions
