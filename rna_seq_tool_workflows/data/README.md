# RNA-seq Tool Workflow Data Plan

## Bulk RNA-seq tools

Primary lightweight real dataset target:

- HISAT2 official example reads and reference files for a small chromosome 22
  region
- Intended for:
  - HISAT2
  - STAR
  - HTSeq
  - StringTie
  - Kallisto

Why:

- real reads and reference material
- small enough for tool-focused demos
- useful for comparing aligners and downstream tools

Official HISAT2 example files currently used by `01_hisat2_alignment`:

- `https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reference/22_20-21M.fa`
- `https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reads/reads_1.fa`
- `https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reads/reads_2.fa`

## Single-cell tool

Primary real dataset target:

- official 10x PBMC tutorial data used by kallisto|bustools tutorials

Intended for:

- bustools

Why:

- standard public single-cell example
- directly relevant to BUS format and UMI/cell-barcode workflows
