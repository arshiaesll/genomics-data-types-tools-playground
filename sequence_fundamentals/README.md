# Sequence Fundamentals

Core sequence-format projects using real FASTA and FASTQ data.

Outputs meant for people to read live under each mini-project's
`outputs/presentable/` folder.

## Mini-projects

### `01_fasta_parser/01_fasta_parser.py`

- real FASTA parsing on `NC_045512.2`
- reports GC content, base counts, and repeated k-mers
- writes a summary and composition figure

### `02_fastq_parser/02_fastq_parser.py`

- downloads paired FASTQ files for `SRR20172987`
- converts ASCII quality symbols to Phred scores
- writes a quality summary and per-position quality plot

### `03_sequence_utilities/03_sequence_utilities.py`

- runs reverse complement, motif search, GC windows, and ORF detection
- uses the same SARS-CoV-2 FASTA as project 1
- writes one report plus three figures

## Suggested Build Order

1. Start with `01_fasta_parser/01_fasta_parser.py`
2. Build `02_fastq_parser/02_fastq_parser.py`
3. Reuse those helpers in `03_sequence_utilities/03_sequence_utilities.py`
