# RNA-seq Tool Workflows

This section extends the RNA-seq phase with one mini-project per major tool.
Each mini-project is tied to a real public dataset or official example dataset
and is designed to produce both a workflow artifact and at least one figure.

Output convention:

- `outputs/presentable/` is for reports, figures, and readable tables.
- `outputs/technical/` is for raw workflow artifacts such as SAM, BUS, index,
  JSON, and other machine-oriented files.

## Shared Strategy

For the bulk RNA-seq tools, we will mostly use the official HISAT2 example
reads and small chromosome 22 reference region because they are real sequencing
reads, lightweight, and suitable for alignment/counting/assembly demos.

For the single-cell workflow, we will use the official 10x PBMC tutorial data
recommended in the kallisto|bustools ecosystem.

## Mini-projects

### `01_hisat2_alignment/01_hisat2_alignment.py`

Goal:

Run HISAT2 on real official paired-end example reads, build the small chr22
index locally, and summarize the alignment outcome.

Real data:

- HISAT2 official example reads
- HISAT2 official chromosome 22 reference region

Current outputs:

- `outputs/hisat2_alignment_summary.txt`
- `outputs/hisat2_alignment_breakdown.png`
- `outputs/hisat2_mapq_histogram.png`
- `outputs/aligned_reads.sam`

What it showcases:

- downloading official example inputs automatically
- building a HISAT2 index with `hisat2-build`
- running `hisat2 --dta` on paired-end reads
- interpreting overall alignment rate and concordant alignment categories
- inspecting mapping quality from the generated SAM file

### `02_star_alignment/02_star_alignment.py`

Goal:

Run STAR on the same real paired-end example reads and summarize mapping and
splice-junction outputs for side-by-side comparison with HISAT2.

Real data:

- same official HISAT2 example reads and small chr22 reference region

Current outputs:

- `outputs/star_alignment_summary.txt`
- `outputs/star_Log.final.out`
- `outputs/star_SJ.out.tab`
- `outputs/star_Aligned.out.sam`
- `outputs/star_mapping_summary.png`
- `outputs/star_splice_junction_motifs.png`

What it showcases:

- building a STAR genome index for the small chr22 reference
- aligning the same read pairs used in the HISAT2 mini-project
- reading STAR's `Log.final.out` summary metrics
- inspecting `SJ.out.tab` splice-junction output
- comparing STAR mapping behavior with HISAT2 on the same example

### `03_htseq_counting/03_htseq_counting.py`

Goal:

Use the official HTSeq yeast example alignment plus a real yeast GTF to count
reads per gene with HTSeq.

Real data:

- official HTSeq example RNA-seq excerpt in SAM format
- official Saccharomyces cerevisiae GTF from the same example set

Current outputs:

- `outputs/htseq_counts.tsv`
- `outputs/htseq_counting_summary.txt`
- `outputs/htseq_top_genes.png`
- `outputs/htseq_special_counters.png`

What it showcases:

- using `htseq-count` on a real multi-gene RNA-seq example
- the role of a feature annotation in assigning reads to genes
- HTSeq special counters such as `__no_feature` and `__ambiguous`
- how aligned read pairs become a count table for downstream RNA-seq analysis

### `04_stringtie_quantification/04_stringtie_quantification.py`

Goal:

Assemble transcripts and summarize abundance with StringTie from a real yeast
RNA-seq example alignment.

Real data:

- official HTSeq yeast RNA-seq excerpt in SAM format
- official Saccharomyces cerevisiae GTF from the same example set

Current outputs:

- `outputs/stringtie_gene_abundance.tsv`
- `outputs/stringtie_assembled.gtf`
- `outputs/stringtie_quantification_summary.txt`
- `outputs/yeast_RNASeq_excerpt.sorted.sam`
- `outputs/stringtie_top_genes_tpm.png`
- `outputs/stringtie_transcript_lengths.png`

What it showcases:

- sorting alignment records into coordinate order for transcript assembly
- running StringTie in reference-guided mode with `-G`
- collecting gene abundance estimates from `-A`
- inspecting assembled transcript structures in the output GTF
- comparing StringTie abundance summaries with HTSeq overlap counts

### `05_kallisto_quantification/05_kallisto_quantification.py`

Goal:

Quantify transcript abundance from real RNA-seq reads using kallisto
pseudoalignment.

Real data:

- official kallisto bundled paired-end test reads
- official kallisto bundled transcriptome FASTA and matching GTF

Current outputs:

- `outputs/kallisto_abundance.tsv`
- `outputs/run_info.json`
- `outputs/transcripts.idx`
- `outputs/kallisto_quantification_summary.txt`
- `outputs/kallisto_top_transcripts_tpm.png`
- `outputs/kallisto_top_genes_tpm.png`

What it showcases:

- building a transcriptome index with `kallisto index`
- quantifying paired-end reads with pseudoalignment
- reading transcript-level TPM and estimated counts from `abundance.tsv`
- rolling transcript TPM up to gene-level summaries for interpretation

### `06_bustools_single_cell/06_bustools_single_cell.py`

Goal:

Process real single-cell RNA-seq data into a cell-by-gene matrix using
kallisto|bustools.

Real data:

- official bundled kallisto single-cell test reads
- official bundled transcriptome FASTA and matching GTF

Current outputs:

- `outputs/bustools_matrix_summary.txt`
- `outputs/bustools_inspect.json`
- `outputs/cells_x_genes.mtx`
- `outputs/cells_x_genes.barcodes.txt`
- `outputs/cells_x_genes.genes.txt`
- `outputs/gene_by_cell_counts.tsv`
- `outputs/intermediate_step_summary.txt`
- `outputs/umi_per_barcode.png`
- `outputs/cells_per_gene_distribution.png`
- `outputs/raw_vs_corrected_barcode_rank.png`

What it showcases:

- running `kallisto bus` on 10x-style single-cell reads
- generating a barcode allowlist directly from BUS records
- correcting barcode errors with `bustools correct`
- turning BUS records into a barcode-by-gene matrix with `bustools count`
- converting the sparse matrix output into a readable dense gene-by-cell table
- inspecting UMI and barcode distributions for candidate cells
