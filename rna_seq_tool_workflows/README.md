# RNA-seq Tool Workflows

One mini-project per common RNA-seq tool, each using real public or official
example data.

Output convention:

- `outputs/presentable/` is for reports, figures, and readable tables.
- `outputs/technical/` is for raw workflow artifacts such as SAM, BUS, index,
  JSON, and other machine-oriented files.

## Tool Availability

- `hisat2` and `htseq-count` are expected from the Conda environment
- `STAR`, `StringTie`, `kallisto`, and `bustools` use the small bundled
  binaries under `rna_seq_tool_workflows/tools/`
- each script checks for its required tool and prints a clear error if it is
  missing

## Mini-projects

### `01_hisat2_alignment/01_hisat2_alignment.py`

- official chr22 example reads
- builds a HISAT2 index and runs `hisat2 --dta`
- summarizes alignment rate and mapping quality

### `02_star_alignment/02_star_alignment.py`

- same chr22 example reads as the HISAT2 project
- builds a STAR index and runs alignment
- summarizes mapping metrics and splice junction output

### `03_htseq_counting/03_htseq_counting.py`

- uses the official HTSeq yeast SAM and GTF
- runs `htseq-count` to produce a multi-gene count table
- highlights both top genes and special counters

### `04_stringtie_quantification/04_stringtie_quantification.py`

- reuses the yeast example alignment from HTSeq
- sorts the alignment and runs StringTie with a guide GTF
- reports gene abundance and assembled transcript structure

### `05_kallisto_quantification/05_kallisto_quantification.py`

- uses bundled kallisto test reads plus transcriptome files
- builds the transcript index and runs `kallisto quant`
- summarizes transcript- and gene-level TPM

### `06_bustools_single_cell/06_bustools_single_cell.py`

- runs `kallisto bus` and `bustools` on bundled single-cell test data
- builds a cell-by-gene matrix plus readable TSV versions
- includes intermediate summaries so the barcode and UMI steps are visible
