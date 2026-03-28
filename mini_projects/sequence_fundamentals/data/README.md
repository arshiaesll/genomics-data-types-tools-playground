# Sequence Fundamentals Data Plan

This folder will hold the real files used by the sequence fundamentals mini-projects.

The FASTA visualization script currently saves a Matplotlib `.png` plot, so
the environment running it should have `matplotlib` installed.

## Planned Files

### FASTA

- `NC_045512.2.fasta`
  - Source: NCBI Nucleotide
  - Accession: `NC_045512.2`
  - Description: SARS-CoV-2 reference genome
  - Use in:
    - `01_fasta_parser.py`
    - `03_sequence_utilities.py`

### FASTQ

- `SRR20172987_1.fastq.gz`
- `SRR20172987_2.fastq.gz`
  - Source: NCBI SRA / ENA
  - Run accession: `SRR20172987`
  - Description: paired-end 16S amplicon reads from a freshwater bacterial community
  - Use in:
    - `02_fastq_parser.py`
    - `04_quality_filtering.py`

The current FASTQ parser script queries ENA for the FASTQ download URLs and
then downloads the gzip-compressed files directly.

## Download Notes

For FASTA:

- NCBI Nucleotide provides the reference sequence view for `NC_045512.2`

For FASTQ:

- NCBI SRA provides metadata for `SRR20172987`
- The SRA Toolkit can export the run to FASTQ
- ENA also exposes FASTQ download locations through its file-report API
