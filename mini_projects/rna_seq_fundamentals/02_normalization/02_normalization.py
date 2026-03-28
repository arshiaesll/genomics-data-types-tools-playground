"""Mini-project: RNA-seq normalization.

Description:
Use a real GEO RNA-seq count matrix and human gene annotation table to compute
CPM, FPKM, and TPM for one sample and compare how these normalizations affect
expression magnitudes and top-ranked genes.

Objectives:
- Download or reuse a real RNA-seq count matrix
- Download or reuse the matching human gene annotation table
- Compute CPM, FPKM, and TPM from raw counts
- Compare top genes under different normalization schemes
- Generate a normalization comparison plot

Files we will use:
- Input: ../data/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz
- Input: ../data/human_gene_annotations.tsv.gz
- Output: outputs/normalization_report.txt
- Output: outputs/normalization_comparison.png
- Output: outputs/tpm_two_sample_comparison.png
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
import sys
from urllib.error import URLError
from urllib.request import urlopen

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


RNA_SEQ_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = RNA_SEQ_DIR / "data"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
COUNT_MATRIX_PATH = DATA_DIR / "GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz"
ANNOTATION_PATH = DATA_DIR / "human_gene_annotations.tsv.gz"
REPORT_PATH = PRESENTABLE_DIR / "normalization_report.txt"
PLOT_PATH = PRESENTABLE_DIR / "normalization_comparison.png"
TPM_COMPARISON_PATH = PRESENTABLE_DIR / "tpm_two_sample_comparison.png"
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


def download_file(url: str, destination: Path) -> tuple[Path, bool]:
    ensure_directories()
    if destination.exists():
        return destination, False

    with urlopen(url, timeout=60) as response:
        destination.write_bytes(response.read())

    return destination, True


def download_inputs() -> dict[str, tuple[Path, bool]]:
    return {
        "count_matrix": download_file(COUNT_MATRIX_URL, COUNT_MATRIX_PATH),
        "annotation": download_file(ANNOTATION_URL, ANNOTATION_PATH),
    }


def first_matching_key(fieldnames: list[str], candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    return None


def load_annotations(path: Path = ANNOTATION_PATH) -> dict[str, dict[str, object]]:
    annotations: dict[str, dict[str, object]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return annotations

        gene_id_key = first_matching_key(reader.fieldnames, ["GeneID", "Gene ID", "gene_id"])
        symbol_key = first_matching_key(reader.fieldnames, ["Symbol", "GeneSymbol", "gene_symbol"])
        description_key = first_matching_key(reader.fieldnames, ["Description", "description", "Name"])
        length_key = first_matching_key(reader.fieldnames, ["Length", "length"])

        if gene_id_key is None or length_key is None:
            return annotations

        for row in reader:
            gene_id = row.get(gene_id_key, "").strip()
            length_text = row.get(length_key, "").strip()
            if not gene_id or not length_text:
                continue

            try:
                length = int(length_text)
            except ValueError:
                continue

            symbol = row.get(symbol_key, "").strip() if symbol_key else ""
            description = row.get(description_key, "").strip() if description_key else ""
            annotations[gene_id] = {
                "symbol": symbol,
                "description": description,
                "length": length,
            }

    return annotations


def gene_label(gene_id: str, annotations: dict[str, dict[str, object]]) -> str:
    annotation = annotations.get(gene_id, {})
    symbol = annotation.get("symbol", "")
    if symbol:
        return f"{symbol} [{gene_id}]"
    return gene_id


def parse_count_matrix(
    count_path: Path = COUNT_MATRIX_PATH,
    annotation_path: Path = ANNOTATION_PATH,
) -> dict[str, object]:
    annotations = load_annotations(annotation_path)

    with gzip.open(count_path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_names = header[1:]
        sample_count = len(sample_names)
        records: list[dict[str, object]] = []
        library_sizes = [0] * sample_count

        for row in reader:
            if not row:
                continue

            gene_id = row[0]
            counts = [int(value) for value in row[1:]]
            annotation = annotations.get(gene_id, {})
            records.append(
                {
                    "gene_id": gene_id,
                    "label": gene_label(gene_id, annotations),
                    "length": int(annotation.get("length", 0)),
                    "counts": counts,
                }
            )

            for index, count in enumerate(counts):
                library_sizes[index] += count

    return {
        "sample_names": sample_names,
        "library_sizes": library_sizes,
        "records": records,
    }


def calculate_cpm(count: int, library_size: int) -> float:
    if library_size == 0:
        return 0.0
    return count / library_size * 1_000_000


def calculate_fpkm(count: int, gene_length: int, library_size: int) -> float:
    if gene_length <= 0 or library_size == 0:
        return 0.0
    gene_length_kb = gene_length / 1_000
    library_size_million = library_size / 1_000_000
    return count / (gene_length_kb * library_size_million)


def calculate_tpm_values(records: list[dict[str, object]], sample_index: int) -> list[float]:
    rpk_values: list[float] = []
    for record in records:
        gene_length = int(record["length"])
        count = int(record["counts"][sample_index])
        if gene_length <= 0:
            rpk_values.append(0.0)
            continue
        rpk_values.append(count / (gene_length / 1_000))

    scaling_factor = sum(rpk_values) / 1_000_000
    if scaling_factor == 0:
        return [0.0] * len(records)

    return [rpk / scaling_factor for rpk in rpk_values]


def normalize_sample(matrix: dict[str, object], sample_index: int = 0) -> dict[str, object]:
    sample_names = matrix["sample_names"]
    library_sizes = matrix["library_sizes"]
    records = matrix["records"]

    sample_name = sample_names[sample_index]
    library_size = library_sizes[sample_index]
    tpm_values = calculate_tpm_values(records, sample_index)

    normalized_rows: list[dict[str, object]] = []
    for record, tpm in zip(records, tpm_values):
        count = int(record["counts"][sample_index])
        length = int(record["length"])
        normalized_rows.append(
            {
                "gene_id": record["gene_id"],
                "label": record["label"],
                "length": length,
                "count": count,
                "cpm": calculate_cpm(count, library_size),
                "fpkm": calculate_fpkm(count, length, library_size),
                "tpm": tpm,
            }
        )

    top_by_count = top_rows(normalized_rows, "count", 10)
    top_by_cpm = top_rows(normalized_rows, "cpm", 10)
    top_by_fpkm = top_rows(normalized_rows, "fpkm", 10)
    top_by_tpm = top_rows(normalized_rows, "tpm", 10)

    return {
        "sample_index": sample_index,
        "sample_name": sample_name,
        "library_size": library_size,
        "normalized_rows": normalized_rows,
        "top_by_count": top_by_count,
        "top_by_cpm": top_by_cpm,
        "top_by_fpkm": top_by_fpkm,
        "top_by_tpm": top_by_tpm,
    }


def top_rows(rows: list[dict[str, object]], key: str, limit: int) -> list[dict[str, object]]:
    nonzero_rows = [row for row in rows if float(row[key]) > 0]
    return sorted(nonzero_rows, key=lambda row: (-float(row[key]), str(row["label"])))[:limit]


def build_report_text(result: dict[str, object]) -> str:
    lines = [
        "Mini-project 2: RNA-seq normalization",
        "Dataset: GEO accession GSE164073",
        f"Selected sample: {result['sample_name']}",
        f"Library size: {result['library_size']}",
        "",
        "Normalization formulas",
        "CPM (Counts Per Million) = count / total_counts * 1,000,000",
        "FPKM (Fragments Per Kilobase Million) = count / (gene_length_kb * library_size_million)",
        "TPM (Transcripts Per Million) = RPK / sum(RPK) * 1,000,000",
        "",
        "Top genes by raw count",
    ]

    for row in result["top_by_count"]:
        lines.append(f"{row['label']}: count={row['count']}")

    lines.extend(["", "Top genes by CPM (Counts Per Million)"])
    for row in result["top_by_cpm"]:
        lines.append(f"{row['label']}: CPM={row['cpm']:.2f}")

    lines.extend(["", "Top genes by FPKM (Fragments Per Kilobase Million)"])
    for row in result["top_by_fpkm"]:
        lines.append(f"{row['label']}: FPKM={row['fpkm']:.2f}")

    lines.extend(["", "Top genes by TPM (Transcripts Per Million)"])
    for row in result["top_by_tpm"]:
        lines.append(f"{row['label']}: TPM={row['tpm']:.2f}")

    return "\n".join(lines) + "\n"


def save_report(result: dict[str, object], destination: Path = REPORT_PATH) -> Path:
    ensure_directories()
    destination.write_text(build_report_text(result), encoding="utf-8")
    return destination


def create_normalization_plot(result: dict[str, object], destination: Path = PLOT_PATH) -> Path:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Matplotlib is required for the normalization comparison plot. "
            "Install it with: python3 -m pip install matplotlib"
        ) from error

    figure, axes = plt.subplots(2, 2, figsize=(14, 10))
    panels = [
        ("Raw count", result["top_by_count"], "count", "#2E86AB"),
        ("CPM\n(Counts Per Million)", result["top_by_cpm"], "cpm", "#7CB518"),
        ("FPKM\n(Fragments Per Kilobase Million)", result["top_by_fpkm"], "fpkm", "#F4A259"),
        ("TPM\n(Transcripts Per Million)", result["top_by_tpm"], "tpm", "#9C27B0"),
    ]

    for axis, (title, rows, key, color) in zip(axes.flat, panels):
        labels = [str(row["label"]) for row in rows]
        values = [float(row[key]) for row in rows]
        axis.barh(labels, values, color=color)
        axis.set_title(title)
        axis.invert_yaxis()
        axis.grid(axis="x", alpha=0.3, linestyle="--")
        axis.set_axisbelow(True)

    figure.suptitle(f"Normalization Comparison for {result['sample_name']}", fontsize=14)
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def create_tpm_two_sample_plot(
    first_result: dict[str, object],
    second_result: dict[str, object],
    destination: Path = TPM_COMPARISON_PATH,
) -> Path:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Matplotlib is required for the TPM comparison plot. "
            "Install it with: python3 -m pip install matplotlib"
        ) from error

    first_rows = {row["gene_id"]: row for row in first_result["normalized_rows"]}
    second_rows = {row["gene_id"]: row for row in second_result["normalized_rows"]}
    selected_gene_ids = []

    for row in first_result["top_by_tpm"][:5]:
        if row["gene_id"] not in selected_gene_ids:
            selected_gene_ids.append(row["gene_id"])
    for row in second_result["top_by_tpm"][:5]:
        if row["gene_id"] not in selected_gene_ids:
            selected_gene_ids.append(row["gene_id"])

    labels = [str(first_rows.get(gene_id, second_rows[gene_id])["label"]) for gene_id in selected_gene_ids]
    first_values = [float(first_rows.get(gene_id, {"tpm": 0.0})["tpm"]) for gene_id in selected_gene_ids]
    second_values = [float(second_rows.get(gene_id, {"tpm": 0.0})["tpm"]) for gene_id in selected_gene_ids]

    figure, axis = plt.subplots(figsize=(12, 6))
    x_positions = list(range(len(labels)))
    width = 0.38
    axis.bar(
        [x - width / 2 for x in x_positions],
        first_values,
        width=width,
        label=first_result["sample_name"],
        color="#2E86AB",
    )
    axis.bar(
        [x + width / 2 for x in x_positions],
        second_values,
        width=width,
        label=second_result["sample_name"],
        color="#F18F01",
    )

    axis.set_title("TPM Comparison Across Two Samples")
    axis.set_ylabel("TPM (Transcripts Per Million)")
    axis.set_xticks(x_positions)
    axis.set_xticklabels(labels, rotation=45, ha="right")
    axis.grid(axis="y", alpha=0.3, linestyle="--")
    axis.set_axisbelow(True)
    axis.legend()
    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def run() -> dict[str, Path]:
    download_inputs()
    matrix = parse_count_matrix()
    first_result = normalize_sample(matrix, sample_index=0)
    second_result = normalize_sample(matrix, sample_index=1)
    report_path = save_report(first_result)
    plot_path = create_normalization_plot(first_result)
    tpm_plot_path = create_tpm_two_sample_plot(first_result, second_result)
    return {
        "report_path": report_path,
        "plot_path": plot_path,
        "tpm_plot_path": tpm_plot_path,
    }


def main() -> None:
    try:
        downloads = download_inputs()
        matrix = parse_count_matrix()
        first_result = normalize_sample(matrix, sample_index=0)
        second_result = normalize_sample(matrix, sample_index=1)
        report_path = save_report(first_result)
        plot_path = create_normalization_plot(first_result)
        tpm_plot_path = create_tpm_two_sample_plot(first_result, second_result)
    except URLError as error:
        print("Download failed. The script is ready, but it could not reach GEO right now.")
        print(f"Reason: {error}")
        print(f"Expected count matrix URL: {COUNT_MATRIX_URL}")
        print(f"Expected annotation URL: {ANNOTATION_URL}")
        return
    except ModuleNotFoundError as error:
        print(error)
        return

    print("RNA-seq mini-project 2 completed.")
    for name, (path, was_downloaded) in downloads.items():
        status = "Downloaded" if was_downloaded else "Using existing"
        print(f"{status} {name}: {path}")
    print(f"Normalization report: {report_path}")
    print(f"Normalization plot: {plot_path}")
    print(f"Two-sample TPM plot: {tpm_plot_path}")


if __name__ == "__main__":
    main()
