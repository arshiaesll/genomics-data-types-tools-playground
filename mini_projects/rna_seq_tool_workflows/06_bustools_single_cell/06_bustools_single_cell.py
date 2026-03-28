"""Mini-project: bustools single-cell workflow.

Description:
Process real bundled single-cell reads with the kallisto | bustools workflow to
produce a barcode-by-gene matrix and summarize single-cell counting outputs.

Official data used:
- bundled kallisto single-cell test FASTQ reads
- bundled kallisto transcriptome FASTA and matching GTF
- official kallisto macOS binary package
- official bustools macOS binary package

Outputs:
- outputs/transcripts.idx
- outputs/output.bus
- outputs/output.sorted.bus
- outputs/matrix.ec
- outputs/transcripts.txt
- outputs/barcode_allowlist.txt
- outputs/output.corrected.bus
- outputs/output.corrected.sorted.bus
- outputs/bustools_inspect.json
- outputs/transcript_to_gene.tsv
- outputs/cells_x_genes.barcodes.txt
- outputs/cells_x_genes.genes.txt
- outputs/cells_x_genes.mtx
- outputs/gene_by_cell_counts.tsv
- outputs/gene_by_top_cells_preview.tsv
- outputs/raw_bus_barcode_counts.tsv
- outputs/corrected_bus_barcode_counts.tsv
- outputs/intermediate_step_summary.txt
- outputs/bustools_matrix_summary.txt
- outputs/umi_per_barcode.png
- outputs/cells_per_gene_distribution.png
- outputs/raw_vs_corrected_barcode_rank.png
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import gzip
import json
import os
import re
import shutil
import subprocess
import textwrap


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = Path(__file__).resolve().parent
TOOLS_DIR = WORKFLOW_DIR / "tools"
OUTPUT_DIR = PROJECT_DIR / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
TECHNICAL_DIR = OUTPUT_DIR / "technical"

KALLISTO_BINARY = TOOLS_DIR / "kallisto" / "kallisto"
BUSTOOLS_BINARY = TOOLS_DIR / "bustools_pkg" / "bustools" / "bustools"
TEST_DIR = TOOLS_DIR / "kallisto" / "test"

TRANSCRIPTS_FASTA_PATH = TEST_DIR / "transcripts.fasta.gz"
TRANSCRIPTS_GTF_PATH = TEST_DIR / "transcripts.gtf.gz"
READ_1_PATH = TEST_DIR / "sc_reads_1.fastq.gz"
READ_2_PATH = TEST_DIR / "sc_reads_2.fastq.gz"

INDEX_PATH = TECHNICAL_DIR / "transcripts.idx"
BUS_PATH = TECHNICAL_DIR / "output.bus"
SORTED_BUS_PATH = TECHNICAL_DIR / "output.sorted.bus"
EC_PATH = TECHNICAL_DIR / "matrix.ec"
TRANSCRIPTS_TXT_PATH = TECHNICAL_DIR / "transcripts.txt"
RUN_INFO_PATH = TECHNICAL_DIR / "run_info.json"
ALLOWLIST_PATH = PRESENTABLE_DIR / "barcode_allowlist.txt"
CORRECTED_BUS_PATH = TECHNICAL_DIR / "output.corrected.bus"
CORRECTED_SORTED_BUS_PATH = TECHNICAL_DIR / "output.corrected.sorted.bus"
INSPECT_JSON_PATH = TECHNICAL_DIR / "bustools_inspect.json"
T2G_PATH = PRESENTABLE_DIR / "transcript_to_gene.tsv"
MATRIX_PREFIX = TECHNICAL_DIR / "cells_x_genes"
MATRIX_PATH = TECHNICAL_DIR / "cells_x_genes.mtx"
TECHNICAL_BARCODES_PATH = TECHNICAL_DIR / "cells_x_genes.barcodes.txt"
TECHNICAL_GENES_PATH = TECHNICAL_DIR / "cells_x_genes.genes.txt"
BARCODES_PATH = PRESENTABLE_DIR / "cells_x_genes.barcodes.txt"
GENES_PATH = PRESENTABLE_DIR / "cells_x_genes.genes.txt"
SUMMARY_PATH = PRESENTABLE_DIR / "bustools_matrix_summary.txt"
INTERMEDIATE_SUMMARY_PATH = PRESENTABLE_DIR / "intermediate_step_summary.txt"
GENE_BY_CELL_TSV_PATH = PRESENTABLE_DIR / "gene_by_cell_counts.tsv"
GENE_BY_TOP_CELLS_PREVIEW_PATH = PRESENTABLE_DIR / "gene_by_top_cells_preview.tsv"
RAW_BARCODE_COUNTS_PATH = PRESENTABLE_DIR / "raw_bus_barcode_counts.tsv"
CORRECTED_BARCODE_COUNTS_PATH = PRESENTABLE_DIR / "corrected_bus_barcode_counts.tsv"
UMI_PER_BARCODE_PLOT_PATH = PRESENTABLE_DIR / "umi_per_barcode.png"
CELLS_PER_GENE_PLOT_PATH = PRESENTABLE_DIR / "cells_per_gene_distribution.png"
RAW_VS_CORRECTED_BARCODE_PLOT_PATH = PRESENTABLE_DIR / "raw_vs_corrected_barcode_rank.png"
THREADS = 2


def ensure_directories() -> None:
    PRESENTABLE_DIR.mkdir(parents=True, exist_ok=True)
    TECHNICAL_DIR.mkdir(parents=True, exist_ok=True)


def configure_matplotlib() -> None:
    mpl_config_dir = TECHNICAL_DIR / ".mplconfig"
    cache_dir = TECHNICAL_DIR / ".cache"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))


def resolve_binary(preferred_path: Path, command_name: str) -> Path | None:
    path_binary = shutil.which(command_name)
    if path_binary:
        return Path(path_binary)
    if preferred_path.exists():
        return preferred_path
    return None


def required_inputs_exist() -> bool:
    return all(
        path.exists()
        for path in [TRANSCRIPTS_FASTA_PATH, TRANSCRIPTS_GTF_PATH, READ_1_PATH, READ_2_PATH]
    )


def run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    return subprocess.run(command, check=False, capture_output=True, text=True, env=env)


def require_success(result: subprocess.CompletedProcess[str], message: str) -> None:
    if result.returncode != 0:
        raise SystemExit(result.stderr or result.stdout or message)


def load_matplotlib():
    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None
    return plt


def build_index_command(binary: Path) -> list[str]:
    return [
        str(binary),
        "index",
        "-i",
        str(INDEX_PATH),
        str(TRANSCRIPTS_FASTA_PATH),
    ]


def build_bus_command(binary: Path) -> list[str]:
    return [
        str(binary),
        "bus",
        "-i",
        str(INDEX_PATH),
        "-x",
        "10xv2",
        "-o",
        str(TECHNICAL_DIR),
        "-t",
        str(THREADS),
        str(READ_1_PATH),
        str(READ_2_PATH),
    ]


def build_sort_command(binary: Path, input_bus: Path, output_bus: Path) -> list[str]:
    return [
        str(binary),
        "sort",
        "-o",
        str(output_bus),
        str(input_bus),
    ]


def build_allowlist_command(binary: Path, sorted_bus: Path) -> list[str]:
    return [
        str(binary),
        "allowlist",
        "-o",
        str(ALLOWLIST_PATH),
        str(sorted_bus),
    ]


def build_correct_command(binary: Path, input_bus: Path) -> list[str]:
    return [
        str(binary),
        "correct",
        "-w",
        str(ALLOWLIST_PATH),
        "-o",
        str(CORRECTED_BUS_PATH),
        str(input_bus),
    ]


def build_inspect_command(binary: Path, sorted_bus: Path) -> list[str]:
    return [
        str(binary),
        "inspect",
        "-o",
        str(INSPECT_JSON_PATH),
        "-e",
        str(EC_PATH),
        "-w",
        str(ALLOWLIST_PATH),
        str(sorted_bus),
    ]


def build_count_command(binary: Path, sorted_bus: Path) -> list[str]:
    return [
        str(binary),
        "count",
        "--genecounts",
        "-o",
        str(MATRIX_PREFIX),
        "-g",
        str(T2G_PATH),
        "-e",
        str(EC_PATH),
        "-t",
        str(TRANSCRIPTS_TXT_PATH),
        str(sorted_bus),
    ]


def build_text_command(binary: Path, bus_path: Path) -> list[str]:
    return [
        str(binary),
        "text",
        "-p",
        str(bus_path),
    ]


def parse_gtf_attributes(attribute_text: str) -> dict[str, str]:
    attribute_pattern = re.compile(r'(\S+) "([^"]+)"')
    return dict(attribute_pattern.findall(attribute_text))


def write_transcript_to_gene_map(gtf_path: Path, txnames_path: Path, destination: Path) -> Path:
    txnames = {
        line.strip()
        for line in txnames_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    seen: set[str] = set()
    rows: list[tuple[str, str]] = []
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            transcript_version = attrs.get("transcript_version", "")
            if transcript_id and transcript_version and not transcript_id.endswith(f".{transcript_version}"):
                transcript_id = f"{transcript_id}.{transcript_version}"
            gene_name = attrs.get("gene_name", attrs.get("gene_id", ""))
            if (
                not transcript_id
                or not gene_name
                or transcript_id in seen
                or transcript_id not in txnames
            ):
                continue
            rows.append((transcript_id, gene_name))
            seen.add(transcript_id)

    destination.write_text(
        "".join(f"{transcript_id}\t{gene_name}\n" for transcript_id, gene_name in rows),
        encoding="utf-8",
    )
    return destination


def parse_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def copy_presentable_matrix_sidecars() -> None:
    shutil.copyfile(TECHNICAL_BARCODES_PATH, BARCODES_PATH)
    shutil.copyfile(TECHNICAL_GENES_PATH, GENES_PATH)


def parse_bus_text_output(stdout_text: str, stderr_text: str = "") -> dict[str, object]:
    log_lines = [line.strip() for line in stderr_text.splitlines() if line.strip()]
    lines = [line.strip() for line in stdout_text.splitlines() if line.strip()]
    records_read = 0
    entries: list[dict[str, object]] = []

    for line in log_lines + lines:
        if line.startswith("Read in "):
            parts = line.split()
            if len(parts) >= 3:
                records_read = int(parts[2])
            continue
        parts = line.split("\t")
        if len(parts) != 4:
            continue
        barcode, umi, ec, count = parts
        entries.append(
            {
                "barcode": barcode,
                "umi": umi,
                "ec": int(ec),
                "count": int(count),
            }
        )

    barcode_records: Counter[str] = Counter()
    barcode_reads: Counter[str] = Counter()
    ec_reads: Counter[int] = Counter()
    umi_reads: Counter[str] = Counter()

    for entry in entries:
        barcode = str(entry["barcode"])
        umi = str(entry["umi"])
        ec = int(entry["ec"])
        count = int(entry["count"])
        barcode_records[barcode] += 1
        barcode_reads[barcode] += count
        ec_reads[ec] += count
        umi_reads[f"{barcode}|{umi}"] += count

    return {
        "records_read": records_read,
        "entries": entries,
        "barcode_records": barcode_records,
        "barcode_reads": barcode_reads,
        "ec_reads": ec_reads,
        "umi_reads": umi_reads,
    }


def parse_matrix_market(
    matrix_path: Path,
    barcodes_path: Path,
    genes_path: Path,
) -> dict[str, object]:
    barcodes = [line.strip() for line in barcodes_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    genes = [line.strip() for line in genes_path.read_text(encoding="utf-8").splitlines() if line.strip()]

    with matrix_path.open(encoding="utf-8") as handle:
        data_lines = [line.strip() for line in handle if line.strip() and not line.startswith("%")]

    dimensions = tuple(int(value) for value in data_lines[0].split())
    entries = [tuple(int(value) for value in line.split()) for line in data_lines[1:]]

    umi_per_barcode: Counter[str] = Counter()
    cells_per_gene: Counter[str] = Counter()
    umi_per_gene: Counter[str] = Counter()
    gene_by_barcode: dict[str, Counter[str]] = {gene: Counter() for gene in genes}

    for row_index, col_index, value in entries:
        barcode = barcodes[row_index - 1]
        gene = genes[col_index - 1]
        umi_per_barcode[barcode] += value
        umi_per_gene[gene] += value
        cells_per_gene[gene] += 1
        gene_by_barcode[gene][barcode] += value

    return {
        "dimensions": dimensions,
        "entries": entries,
        "barcodes": barcodes,
        "genes": genes,
        "umi_per_barcode": umi_per_barcode,
        "umi_per_gene": umi_per_gene,
        "cells_per_gene": cells_per_gene,
        "gene_by_barcode": gene_by_barcode,
    }


def create_umi_per_barcode_plot(
    umi_per_barcode: Counter[str],
    destination: Path = UMI_PER_BARCODE_PLOT_PATH,
) -> Path | None:
    ranked_counts = sorted(umi_per_barcode.values(), reverse=True)
    if not ranked_counts:
        return None

    plt = load_matplotlib()
    if plt is None:
        return None

    figure, axis = plt.subplots(figsize=(9.6, 5.4))
    axis.plot(range(1, len(ranked_counts) + 1), ranked_counts, color="#2b6cb0", linewidth=1.6)
    axis.set_xlabel("Barcode Rank")
    axis.set_ylabel("Gene-level UMI Count")
    axis.set_title("UMIs per Barcode After bustools Correction")
    axis.set_yscale("log")
    axis.grid(alpha=0.25, linestyle=":")
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def create_cells_per_gene_plot(
    cells_per_gene: Counter[str],
    destination: Path = CELLS_PER_GENE_PLOT_PATH,
) -> Path | None:
    items = cells_per_gene.most_common()
    if not items:
        return None

    labels = [gene for gene, _ in items]
    values = [count for _, count in items]

    plt = load_matplotlib()
    if plt is None:
        return None

    figure, axis = plt.subplots(figsize=(10.0, 5.0))
    bars = axis.bar(labels, values, color="#2f855a")
    axis.set_ylabel("Detected Cells")
    axis.set_title("Cells per Gene in the bustools Matrix")
    axis.tick_params(axis="x", rotation=30)
    max_value = max(values)
    for bar, value in zip(bars, values, strict=True):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            value + max_value * 0.01,
            str(value),
            ha="center",
            fontsize=9,
        )
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def create_raw_vs_corrected_barcode_plot(
    raw_barcode_reads: Counter[str],
    corrected_barcode_reads: Counter[str],
    destination: Path = RAW_VS_CORRECTED_BARCODE_PLOT_PATH,
) -> Path | None:
    raw_ranked = sorted(raw_barcode_reads.values(), reverse=True)
    corrected_ranked = sorted(corrected_barcode_reads.values(), reverse=True)
    if not raw_ranked or not corrected_ranked:
        return None

    plt = load_matplotlib()
    if plt is None:
        return None

    figure, axis = plt.subplots(figsize=(9.8, 5.4))
    axis.plot(range(1, len(raw_ranked) + 1), raw_ranked, label="Raw sorted BUS", color="#c05621", linewidth=1.4)
    axis.plot(
        range(1, len(corrected_ranked) + 1),
        corrected_ranked,
        label="Corrected sorted BUS",
        color="#2b6cb0",
        linewidth=1.4,
    )
    axis.set_xlabel("Barcode Rank")
    axis.set_ylabel("Read Count")
    axis.set_yscale("log")
    axis.set_title("Barcode Rank Before and After bustools Correction")
    axis.grid(alpha=0.25, linestyle=":")
    axis.legend()
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def write_barcode_count_table(
    barcode_reads: Counter[str],
    barcode_records: Counter[str],
    destination: Path,
) -> Path:
    lines = ["barcode\tread_count\tbus_record_count"]
    for barcode, read_count in barcode_reads.most_common():
        lines.append(f"{barcode}\t{read_count}\t{barcode_records[barcode]}")
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination


def write_gene_by_cell_tables(
    matrix_info: dict[str, object],
    full_destination: Path = GENE_BY_CELL_TSV_PATH,
    preview_destination: Path = GENE_BY_TOP_CELLS_PREVIEW_PATH,
) -> tuple[Path, Path]:
    genes: list[str] = matrix_info["genes"]  # type: ignore[assignment]
    umi_per_barcode: Counter[str] = matrix_info["umi_per_barcode"]  # type: ignore[assignment]
    gene_by_barcode: dict[str, Counter[str]] = matrix_info["gene_by_barcode"]  # type: ignore[assignment]

    ranked_barcodes = [barcode for barcode, _ in umi_per_barcode.most_common()]
    full_lines = ["gene\t" + "\t".join(ranked_barcodes)]
    for gene in genes:
        counts = [str(gene_by_barcode[gene].get(barcode, 0)) for barcode in ranked_barcodes]
        full_lines.append(f"{gene}\t" + "\t".join(counts))
    full_destination.write_text("\n".join(full_lines) + "\n", encoding="utf-8")

    preview_barcodes = ranked_barcodes[:20]
    preview_lines = ["gene\t" + "\t".join(preview_barcodes)]
    for gene in genes:
        counts = [str(gene_by_barcode[gene].get(barcode, 0)) for barcode in preview_barcodes]
        preview_lines.append(f"{gene}\t" + "\t".join(counts))
    preview_destination.write_text("\n".join(preview_lines) + "\n", encoding="utf-8")

    return full_destination, preview_destination


def explain_scenario() -> str:
    return textwrap.dedent(
        """\
        Scenario
        We use the official bundled kallisto single-cell test reads, which are
        formatted like 10x Genomics v2 data and target a compact human
        transcriptome containing HOXC-region genes plus UGT3A2.

        What bustools adds
        `kallisto bus` turns read pairs into BUS records that capture cell
        barcode, UMI, and transcript equivalence-class information. `bustools`
        then builds an allowlist, corrects barcode errors, sorts BUS records,
        inspects barcode-level statistics, and collapses transcript evidence into
        a cell-by-gene count matrix.

        Decision we want to make
        We want to see whether the dataset yields a sensible barcode-by-gene
        matrix, which genes dominate the signal, and how UMI counts are
        distributed across candidate cells.
        """
    ).strip()


def write_summary(
    run_info: dict[str, object],
    inspect_info: dict[str, object],
    matrix_info: dict[str, object],
    destination: Path = SUMMARY_PATH,
) -> Path:
    umi_per_barcode: Counter[str] = matrix_info["umi_per_barcode"]  # type: ignore[assignment]
    umi_per_gene: Counter[str] = matrix_info["umi_per_gene"]  # type: ignore[assignment]
    cells_per_gene: Counter[str] = matrix_info["cells_per_gene"]  # type: ignore[assignment]
    dimensions = matrix_info["dimensions"]
    total_gene_umis = sum(umi_per_barcode.values())
    nonzero_cells = sum(1 for value in umi_per_barcode.values() if value > 0)

    top_gene_lines = "\n".join(
        f"- {gene}: {count} UMIs across {cells_per_gene[gene]} detected cells"
        for gene, count in umi_per_gene.most_common(8)
    )
    top_barcode_lines = "\n".join(
        f"- {barcode}: {count} UMIs"
        for barcode, count in umi_per_barcode.most_common(8)
    )

    text = (
        f"{explain_scenario()}\n\n"
        "Input files\n"
        f"- {TRANSCRIPTS_FASTA_PATH}\n"
        f"- {TRANSCRIPTS_GTF_PATH}\n"
        f"- {READ_1_PATH}\n"
        f"- {READ_2_PATH}\n\n"
        "Tools\n"
        "- kallisto bus with `-x 10xv2`\n"
        "- bustools allowlist\n"
        "- bustools correct\n"
        "- bustools sort\n"
        "- bustools inspect\n"
        "- bustools count --genecounts\n\n"
        "kallisto bus summary\n"
        f"- Reads processed: {int(run_info.get('n_processed', 0))}\n"
        f"- Reads pseudoaligned: {int(run_info.get('n_pseudoaligned', 0))}\n"
        f"- Pseudoalignment rate: {float(run_info.get('p_pseudoaligned', 0.0)):.1f}%\n"
        f"- Uniquely pseudoaligned reads: {int(run_info.get('n_unique', 0))}\n"
        f"- Transcript targets in index: {int(run_info.get('n_targets', 0))}\n\n"
        "bustools inspect summary\n"
        f"- BUS records after correction/sorting: {int(inspect_info.get('numRecords', 0))}\n"
        f"- Reads represented in corrected sorted BUS: {int(inspect_info.get('numReads', 0))}\n"
        f"- On-list barcodes in corrected BUS: {int(inspect_info.get('numBarcodes', 0))}\n"
        f"- Median reads per barcode: {float(inspect_info.get('medianReadsPerBarcode', 0.0)):.2f}\n"
        f"- Mean reads per barcode: {float(inspect_info.get('meanReadsPerBarcode', 0.0)):.2f}\n"
        f"- Median UMIs per barcode: {float(inspect_info.get('medianUMIsPerBarcode', 0.0)):.2f}\n"
        f"- Mean UMIs per barcode: {float(inspect_info.get('meanUMIsPerBarcode', 0.0)):.2f}\n"
        f"- Allowlist size: {sum(1 for _ in ALLOWLIST_PATH.open(encoding='utf-8'))}\n\n"
        "Cell-by-gene matrix summary\n"
        f"- Matrix dimensions (barcodes x genes x nonzero entries): {dimensions[0]} x {dimensions[1]} x {dimensions[2]}\n"
        f"- Matrix barcode file: {BARCODES_PATH.name}\n"
        f"- Matrix gene file: {GENES_PATH.name}\n"
        f"- Total gene-level UMIs in matrix: {total_gene_umis}\n"
        f"- Nonzero barcodes in matrix: {nonzero_cells}\n"
        f"- Genes represented in matrix: {len(matrix_info['genes'])}\n\n"
        "Top genes by gene-level UMI count\n"
        f"{top_gene_lines}\n\n"
        "Top barcodes by gene-level UMI count\n"
        f"{top_barcode_lines}"
    )

    destination.write_text(text + "\n", encoding="utf-8")
    return destination


def write_intermediate_summary(
    allowlist_result: subprocess.CompletedProcess[str],
    correct_result: subprocess.CompletedProcess[str],
    raw_bus_info: dict[str, object],
    corrected_bus_info: dict[str, object],
    destination: Path = INTERMEDIATE_SUMMARY_PATH,
) -> Path:
    raw_barcode_reads: Counter[str] = raw_bus_info["barcode_reads"]  # type: ignore[assignment]
    corrected_barcode_reads: Counter[str] = corrected_bus_info["barcode_reads"]  # type: ignore[assignment]
    raw_ec_reads: Counter[int] = raw_bus_info["ec_reads"]  # type: ignore[assignment]
    corrected_ec_reads: Counter[int] = corrected_bus_info["ec_reads"]  # type: ignore[assignment]

    raw_top_barcode_lines = "\n".join(
        f"- {barcode}: {count} reads"
        for barcode, count in raw_barcode_reads.most_common(8)
    )
    corrected_top_barcode_lines = "\n".join(
        f"- {barcode}: {count} reads"
        for barcode, count in corrected_barcode_reads.most_common(8)
    )
    ec_lines = "\n".join(
        f"- EC {ec}: raw {raw_ec_reads[ec]} reads, corrected {corrected_ec_reads.get(ec, 0)} reads"
        for ec, _ in raw_ec_reads.most_common(8)
    )

    text = (
        "Intermediate step summary\n"
        "These summaries expose what changes between the raw BUS records and the\n"
        "barcode-corrected BUS records before the final gene-by-cell matrix is built.\n\n"
        "allowlist step\n"
        f"{allowlist_result.stdout.strip() or allowlist_result.stderr.strip()}\n\n"
        "correct step\n"
        f"{correct_result.stdout.strip() or correct_result.stderr.strip()}\n\n"
        "Raw sorted BUS overview\n"
        f"- BUS records read: {int(raw_bus_info['records_read'])}\n"
        f"- Distinct barcodes observed: {len(raw_barcode_reads)}\n"
        f"- Distinct barcode+UMI combinations: {len(raw_bus_info['umi_reads'])}\n\n"
        "Corrected sorted BUS overview\n"
        f"- BUS records read: {int(corrected_bus_info['records_read'])}\n"
        f"- Distinct barcodes observed: {len(corrected_barcode_reads)}\n"
        f"- Distinct barcode+UMI combinations: {len(corrected_bus_info['umi_reads'])}\n\n"
        "Top raw barcodes by read count\n"
        f"{raw_top_barcode_lines}\n\n"
        "Top corrected barcodes by read count\n"
        f"{corrected_top_barcode_lines}\n\n"
        "Top equivalence classes before and after correction\n"
        f"{ec_lines}\n"
    )
    destination.write_text(text, encoding="utf-8")
    return destination


def main() -> None:
    ensure_directories()

    if not required_inputs_exist():
        raise SystemExit(
            "Bundled single-cell kallisto test inputs are missing. "
            "Make sure the kallisto tool bundle is present in "
            "`mini_projects/rna_seq_tool_workflows/tools/kallisto/test`."
        )

    kallisto_binary = resolve_binary(KALLISTO_BINARY, "kallisto")
    if kallisto_binary is None:
        raise SystemExit(
            "Could not find `kallisto`. Add it to PATH or place the binary at "
            "`mini_projects/rna_seq_tool_workflows/tools/kallisto/kallisto`."
        )

    bustools_binary = resolve_binary(BUSTOOLS_BINARY, "bustools")
    if bustools_binary is None:
        raise SystemExit(
            "Could not find `bustools`. Add it to PATH or place the binary at "
            "`mini_projects/rna_seq_tool_workflows/tools/bustools_pkg/bustools/bustools`."
        )

    index_result = run_command(build_index_command(kallisto_binary))
    require_success(index_result, "kallisto index failed.")

    bus_result = run_command(build_bus_command(kallisto_binary))
    require_success(bus_result, "kallisto bus failed.")

    t2g_path = write_transcript_to_gene_map(TRANSCRIPTS_GTF_PATH, TRANSCRIPTS_TXT_PATH, T2G_PATH)

    sort_result = run_command(build_sort_command(bustools_binary, BUS_PATH, SORTED_BUS_PATH))
    require_success(sort_result, "bustools sort failed on raw BUS.")

    allowlist_result = run_command(build_allowlist_command(bustools_binary, SORTED_BUS_PATH))
    require_success(allowlist_result, "bustools allowlist failed.")

    correct_result = run_command(build_correct_command(bustools_binary, BUS_PATH))
    require_success(correct_result, "bustools correct failed.")

    corrected_sort_result = run_command(
        build_sort_command(bustools_binary, CORRECTED_BUS_PATH, CORRECTED_SORTED_BUS_PATH)
    )
    require_success(corrected_sort_result, "bustools sort failed on corrected BUS.")

    inspect_result = run_command(build_inspect_command(bustools_binary, CORRECTED_SORTED_BUS_PATH))
    require_success(inspect_result, "bustools inspect failed.")

    count_result = run_command(build_count_command(bustools_binary, CORRECTED_SORTED_BUS_PATH))
    require_success(count_result, "bustools count failed.")

    if not t2g_path.exists() or not MATRIX_PATH.exists():
        raise SystemExit("The bustools matrix outputs were not created as expected.")

    raw_text_result = run_command(build_text_command(bustools_binary, SORTED_BUS_PATH))
    require_success(raw_text_result, "bustools text failed on raw BUS.")

    corrected_text_result = run_command(build_text_command(bustools_binary, CORRECTED_SORTED_BUS_PATH))
    require_success(corrected_text_result, "bustools text failed on corrected BUS.")

    run_info = parse_json(RUN_INFO_PATH)
    inspect_info = parse_json(INSPECT_JSON_PATH)
    raw_bus_info = parse_bus_text_output(raw_text_result.stdout, raw_text_result.stderr)
    corrected_bus_info = parse_bus_text_output(corrected_text_result.stdout, corrected_text_result.stderr)
    copy_presentable_matrix_sidecars()
    matrix_info = parse_matrix_market(MATRIX_PATH, BARCODES_PATH, GENES_PATH)

    write_barcode_count_table(
        raw_bus_info["barcode_reads"],  # type: ignore[arg-type]
        raw_bus_info["barcode_records"],  # type: ignore[arg-type]
        RAW_BARCODE_COUNTS_PATH,
    )
    write_barcode_count_table(
        corrected_bus_info["barcode_reads"],  # type: ignore[arg-type]
        corrected_bus_info["barcode_records"],  # type: ignore[arg-type]
        CORRECTED_BARCODE_COUNTS_PATH,
    )
    write_intermediate_summary(allowlist_result, correct_result, raw_bus_info, corrected_bus_info)
    write_gene_by_cell_tables(matrix_info)
    write_summary(run_info, inspect_info, matrix_info)
    create_umi_per_barcode_plot(matrix_info["umi_per_barcode"])  # type: ignore[arg-type]
    create_cells_per_gene_plot(matrix_info["cells_per_gene"])  # type: ignore[arg-type]
    create_raw_vs_corrected_barcode_plot(
        raw_bus_info["barcode_reads"],  # type: ignore[arg-type]
        corrected_bus_info["barcode_reads"],  # type: ignore[arg-type]
    )

    print(f"Built transcript-to-gene map: {T2G_PATH}")
    print(f"Wrote dense gene-by-cell table: {GENE_BY_CELL_TSV_PATH}")
    print(f"Wrote preview gene-by-cell table: {GENE_BY_TOP_CELLS_PREVIEW_PATH}")
    print(f"Wrote raw barcode table: {RAW_BARCODE_COUNTS_PATH}")
    print(f"Wrote corrected barcode table: {CORRECTED_BARCODE_COUNTS_PATH}")
    print(f"Wrote intermediate summary: {INTERMEDIATE_SUMMARY_PATH}")
    print(f"Wrote summary: {SUMMARY_PATH}")
    print(f"Wrote matrix: {MATRIX_PATH}")
    print(f"Wrote barcode rank plot: {UMI_PER_BARCODE_PLOT_PATH}")
    print(f"Wrote cells-per-gene plot: {CELLS_PER_GENE_PLOT_PATH}")
    print(f"Wrote raw-vs-corrected barcode plot: {RAW_VS_CORRECTED_BARCODE_PLOT_PATH}")


if __name__ == "__main__":
    main()
