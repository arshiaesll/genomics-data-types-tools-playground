"""Mini-project: count matrix summary.

Description:
Download a real RNA-seq count matrix from GEO, parse it, summarize expression
per sample, and generate a simple visualization of library sizes.

Objectives:
- Download a real gene-by-sample count matrix
- Parse the matrix header and count rows
- Summarize total counts per sample
- Identify the most highly expressed genes overall and per sample
- Generate a library-size plot

Files we will use:
- Input: ../data/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz
- Output: outputs/count_matrix_summary.txt
- Output: outputs/library_sizes.png
- Output: outputs/detected_genes_per_sample.png
- Output: outputs/top_genes_overall.png

Used public data:
- GEO accession: GSE164073
- NCBI-generated RNA-seq raw count matrix
- GEO human annotation table: Human.GRCh38.p13.annot.tsv.gz
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
import sys
from urllib.error import URLError
from urllib.request import urlopen

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


RNA_SEQ_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = RNA_SEQ_DIR / "data"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
COUNT_MATRIX_PATH = DATA_DIR / "GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz"
ANNOTATION_PATH = DATA_DIR / "human_gene_annotations.tsv.gz"
SUMMARY_PATH = PRESENTABLE_DIR / "count_matrix_summary.txt"
LIBRARY_SIZE_PLOT_PATH = PRESENTABLE_DIR / "library_sizes.png"
DETECTED_GENES_PLOT_PATH = PRESENTABLE_DIR / "detected_genes_per_sample.png"
TOP_GENES_PLOT_PATH = PRESENTABLE_DIR / "top_genes_overall.png"
COUNT_MATRIX_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?type=rnaseq_counts&acc=GSE164073&format=file"
    "&file=GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)
ANNOTATION_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?type=rnaseq_counts&format=file"
    "&file=Human.GRCh38.p13.annot.tsv.gz"
)


def ensure_directories() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PRESENTABLE_DIR.mkdir(parents=True, exist_ok=True)


def download_count_matrix(destination: Path = COUNT_MATRIX_PATH) -> tuple[Path, bool]:
    ensure_directories()

    if destination.exists():
        return destination, False

    with urlopen(COUNT_MATRIX_URL, timeout=60) as response:
        destination.write_bytes(response.read())

    return destination, True


def download_gene_annotations(destination: Path = ANNOTATION_PATH) -> tuple[Path, bool]:
    ensure_directories()

    if destination.exists():
        return destination, False

    with urlopen(ANNOTATION_URL, timeout=60) as response:
        destination.write_bytes(response.read())

    return destination, True


def load_gene_annotations(path: Path = ANNOTATION_PATH) -> dict[str, str]:
    if not path.exists():
        return {}

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return {}

        gene_id_key = first_matching_key(reader.fieldnames, ["GeneID", "Gene ID", "gene_id"])
        symbol_key = first_matching_key(
            reader.fieldnames,
            ["Symbol", "GeneSymbol", "gene_symbol", "Official Symbol"],
        )
        name_key = first_matching_key(
            reader.fieldnames,
            ["Name", "GeneName", "description", "Description", "Official Full Name"],
        )

        if gene_id_key is None:
            return {}

        annotations: dict[str, str] = {}
        for row in reader:
            gene_id = row.get(gene_id_key, "").strip()
            if not gene_id:
                continue

            symbol = row.get(symbol_key, "").strip() if symbol_key else ""
            name = row.get(name_key, "").strip() if name_key else ""

            if symbol and name:
                annotations[gene_id] = f"{symbol} [{gene_id}]"
            elif symbol:
                annotations[gene_id] = f"{symbol} [{gene_id}]"
            elif name:
                annotations[gene_id] = f"{gene_id} - {name}"

        return annotations


def first_matching_key(fieldnames: list[str], candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    return None


def display_gene_label(gene_identifier: str, annotations: dict[str, str]) -> str:
    return annotations.get(gene_identifier, gene_identifier)


def parse_count_matrix(path: Path = COUNT_MATRIX_PATH) -> dict[str, object]:
    try:
        download_gene_annotations()
    except URLError:
        pass

    annotations = load_gene_annotations()

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        gene_column_name = header[0]
        sample_names = header[1:]

        library_sizes = [0] * len(sample_names)
        detected_genes_per_sample = [0] * len(sample_names)
        gene_count = 0
        total_counts = 0
        genes_with_nonzero_counts = 0
        top_genes_overall: list[tuple[str, int]] = []
        top_genes_per_sample: list[list[tuple[str, int]]] = [[] for _ in sample_names]

        for row in reader:
            if not row:
                continue

            gene_id = row[0]
            gene_label = display_gene_label(gene_id, annotations)
            counts = [int(value) for value in row[1:]]
            gene_total = sum(counts)
            gene_count += 1
            total_counts += gene_total
            if gene_total > 0:
                genes_with_nonzero_counts += 1

            for index, count in enumerate(counts):
                library_sizes[index] += count
                if count > 0:
                    detected_genes_per_sample[index] += 1
                    update_top_items(top_genes_per_sample[index], (gene_label, count), limit=5)

            if gene_total > 0:
                update_top_items(top_genes_overall, (gene_label, gene_total), limit=10)

    return {
        "series_accession": "GSE164073",
        "gene_column_name": gene_column_name,
        "sample_names": sample_names,
        "sample_count": len(sample_names),
        "gene_count": gene_count,
        "genes_with_nonzero_counts": genes_with_nonzero_counts,
        "library_sizes": library_sizes,
        "detected_genes_per_sample": detected_genes_per_sample,
        "total_counts": total_counts,
        "top_genes_overall": top_genes_overall,
        "top_genes_per_sample": top_genes_per_sample,
    }


def update_top_items(items: list[tuple[str, int]], candidate: tuple[str, int], limit: int) -> None:
    items.append(candidate)
    items.sort(key=lambda item: (-item[1], item[0]))
    del items[limit:]


def build_summary_text(summary: dict[str, object]) -> str:
    sample_names = summary["sample_names"]
    library_sizes = summary["library_sizes"]
    detected_genes_per_sample = summary["detected_genes_per_sample"]
    top_genes_overall = summary["top_genes_overall"]
    top_genes_per_sample = summary["top_genes_per_sample"]

    lines = [
        "Mini-project 1: Count matrix summary",
        f"Dataset: GEO accession {summary['series_accession']}",
        "Source: NCBI-generated RNA-seq raw counts",
        "",
        f"GSE accession: {summary['series_accession']} (the study/series ID)",
        f"GSM accessions: the sample IDs used as matrix columns",
        f"Gene column: {summary['gene_column_name']}",
        f"Samples: {summary['sample_count']}",
        f"Genes: {summary['gene_count']}",
        f"Genes with nonzero counts: {summary['genes_with_nonzero_counts']}",
        f"Total counts across matrix: {summary['total_counts']}",
        "",
        "Library sizes by sample",
    ]

    for sample_name, library_size in zip(sample_names, library_sizes):
        lines.append(f"{sample_name}: {library_size}")

    lines.extend(["", "Detected genes by sample"])
    for sample_name, detected_genes in zip(sample_names, detected_genes_per_sample):
        lines.append(f"{sample_name}: {detected_genes}")

    lines.extend(["", "Top genes overall by total counts"])
    for gene_id, count in top_genes_overall:
        lines.append(f"{gene_id}: {count}")

    lines.extend(["", "Top 5 genes per sample"])
    for sample_name, genes in zip(sample_names, top_genes_per_sample):
        gene_text = ", ".join(f"{gene_id} ({count})" for gene_id, count in genes)
        lines.append(f"{sample_name}: {gene_text}")

    return "\n".join(lines) + "\n"


def save_summary(summary: dict[str, object], destination: Path = SUMMARY_PATH) -> Path:
    ensure_directories()
    destination.write_text(build_summary_text(summary), encoding="utf-8")
    return destination


def get_pyplot():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Matplotlib is required for the RNA-seq count matrix plot. "
            "Install it with: python3 -m pip install matplotlib"
        ) from error
    return plt


def create_library_size_plot(
    summary: dict[str, object],
    destination: Path = LIBRARY_SIZE_PLOT_PATH,
) -> Path:
    plt = get_pyplot()

    sample_names = summary["sample_names"]
    library_sizes = summary["library_sizes"]

    figure, axis = plt.subplots(figsize=(11, 5.5))
    bars = axis.bar(sample_names, library_sizes, color="#2E86AB")
    axis.set_title("GSE164073 Library Sizes by Sample")
    axis.set_xlabel("Sample")
    axis.set_ylabel("Total counts")
    axis.grid(axis="y", alpha=0.3, linestyle="--")
    axis.set_axisbelow(True)
    axis.tick_params(axis="x", rotation=45)

    for bar, value in zip(bars, library_sizes):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            str(value),
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )

    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def create_detected_genes_plot(
    summary: dict[str, object],
    destination: Path = DETECTED_GENES_PLOT_PATH,
) -> Path:
    plt = get_pyplot()

    sample_names = summary["sample_names"]
    detected_genes = summary["detected_genes_per_sample"]

    figure, axis = plt.subplots(figsize=(11, 5.5))
    bars = axis.bar(sample_names, detected_genes, color="#7CB518")
    axis.set_title("GSE164073 Detected Genes by Sample")
    axis.set_xlabel("Sample")
    axis.set_ylabel("Genes with nonzero counts")
    axis.grid(axis="y", alpha=0.3, linestyle="--")
    axis.set_axisbelow(True)
    axis.tick_params(axis="x", rotation=45)

    for bar, value in zip(bars, detected_genes):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            str(value),
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )

    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def create_top_genes_plot(
    summary: dict[str, object],
    destination: Path = TOP_GENES_PLOT_PATH,
) -> Path:
    plt = get_pyplot()

    top_genes = summary["top_genes_overall"]
    gene_ids = [gene_id for gene_id, _ in top_genes]
    counts = [count for _, count in top_genes]

    figure, axis = plt.subplots(figsize=(10, 6))
    bars = axis.barh(gene_ids, counts, color="#9C27B0")
    axis.set_title("Top Genes Overall by Total Counts")
    axis.set_xlabel("Total counts across all samples")
    axis.set_ylabel("Gene")
    axis.grid(axis="x", alpha=0.3, linestyle="--")
    axis.set_axisbelow(True)
    axis.invert_yaxis()

    for bar, value in zip(bars, counts):
        axis.text(
            value,
            bar.get_y() + bar.get_height() / 2,
            f" {value}",
            va="center",
            fontsize=9,
        )

    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def run() -> dict[str, Path]:
    count_matrix_path, _ = download_count_matrix()
    summary = parse_count_matrix(count_matrix_path)
    summary_path = save_summary(summary)
    library_size_plot = create_library_size_plot(summary)
    detected_genes_plot = create_detected_genes_plot(summary)
    top_genes_plot = create_top_genes_plot(summary)
    return {
        "count_matrix_path": count_matrix_path,
        "summary_path": summary_path,
        "library_size_plot": library_size_plot,
        "detected_genes_plot": detected_genes_plot,
        "top_genes_plot": top_genes_plot,
    }


def main() -> None:
    try:
        count_matrix_path, was_downloaded = download_count_matrix()
        summary = parse_count_matrix(count_matrix_path)
        summary_path = save_summary(summary)
        library_size_plot = create_library_size_plot(summary)
        detected_genes_plot = create_detected_genes_plot(summary)
        top_genes_plot = create_top_genes_plot(summary)
    except URLError as error:
        print("Download failed. The script is ready, but it could not reach GEO right now.")
        print(f"Reason: {error}")
        print(f"Expected count matrix URL: {COUNT_MATRIX_URL}")
        return
    except ModuleNotFoundError as error:
        print(error)
        return

    print("RNA-seq mini-project 1 completed.")
    if was_downloaded:
        print(f"Downloaded count matrix: {count_matrix_path}")
    else:
        print(f"Using existing count matrix: {count_matrix_path}")
    print(f"Summary report: {summary_path}")
    print(f"Library size plot: {library_size_plot}")
    print(f"Detected genes plot: {detected_genes_plot}")
    print(f"Top genes plot: {top_genes_plot}")


if __name__ == "__main__":
    main()
