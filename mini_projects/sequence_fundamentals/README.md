# Sequence Fundamentals

This section focuses on the core sequence-level objects and formats used throughout genomics. We are using real public data where possible so the mini-projects are portfolio pieces, not just exercises.

## Environment

From the repo root:

```bash
conda env create -f environment.yml
conda activate genomics-tools
```

Generated artifacts now live inside each mini-project folder under its own
`outputs/` directory.

Presentation convention:

- `outputs/presentable/` contains the summary files and figures meant to be
  read or shared.

## Chosen Public Data Sources

### Reference sequence for FASTA-based work

- Organism: SARS-CoV-2 reference genome
- Accession: `NC_045512.2`
- Source: NCBI Nucleotide
- Why this is a good starter file:
  - Small enough to inspect manually
  - Real complete genome
  - Good for motif search, GC content, reverse complement, and ORF-related exercises

Planned local file:

- `data/NC_045512.2.fasta`

### Raw reads for FASTQ-based work

- Dataset: 16S amplicon sequencing of a freshwater bacterial community
- Run accession: `SRR20172987`
- Source: NCBI SRA
- Layout: paired-end
- Size shown by NCBI: `34.8 Mb`
- Why this is a good starter file:
  - Real FASTQ data
  - Much smaller than large whole-genome runs
  - Good for quality parsing and filtering

Planned local files:

- `data/SRR20172987_1.fastq.gz`
- `data/SRR20172987_2.fastq.gz`

## Mini-projects

### `01_fasta_parser/01_fasta_parser.py`

Description:

Parse a real FASTA genome file, extract headers and sequence data, and compute basic sequence statistics.

Objectives:

- Read one or more FASTA records
- Support multi-line sequences
- Compute length and GC content
- Print a clean summary for each record

Files we will use:

- `data/NC_045512.2.fasta`
- `01_fasta_parser/outputs/fasta_summary.txt`
- `01_fasta_parser/outputs/fasta_base_composition.png`

### `02_fastq_parser/02_fastq_parser.py`

Description:

Parse real FASTQ reads and translate ASCII quality characters into Phred scores.

Objectives:

- Read FASTQ records in groups of four lines
- Check that sequence length matches quality length
- Convert quality characters to Phred scores
- Compute average quality per read

Files we will use:

- `data/SRR20172987_1.fastq.gz`
- `data/SRR20172987_2.fastq.gz`
- `02_fastq_parser/outputs/fastq_quality_summary.txt`
- `02_fastq_parser/outputs/fastq_position_quality.png`

### `03_sequence_utilities/03_sequence_utilities.py`

Description:

Build reusable DNA and RNA sequence utilities and run them on a real reference sequence.

Objectives:

- Compute reverse complement
- Compute GC content
- Search for motifs
- Optionally detect open reading frames

Files we will use:

- `data/NC_045512.2.fasta`
- `03_sequence_utilities/outputs/sequence_utilities_report.txt`
- `03_sequence_utilities/outputs/gc_content_windows.png`
- `03_sequence_utilities/outputs/motif_counts.png`
- `03_sequence_utilities/outputs/orf_overview.png`

## Suggested Build Order

1. Start with `01_fasta_parser/01_fasta_parser.py`
2. Build `02_fastq_parser/02_fastq_parser.py`
3. Reuse those helpers in `03_sequence_utilities/03_sequence_utilities.py`

## Notes

- Use a small slice of the FASTQ files at first if full parsing feels slow
- Keep each script runnable on its own
- Promote reusable parsing code into shared helper modules once duplication appears
- `01_fasta_parser.py` now expects `matplotlib` for visualization output
