"""Mini-project: fold change comparison.

Description:
Use a real GEO RNA-seq count matrix to compare two biological conditions and
compute simple log2 fold changes for each gene.

Objectives:
- Reuse a real GEO RNA-seq count matrix
- Group samples by biological condition
- Compute mean expression for each group
- Compute log2 fold change
- Rank genes by up- and down-regulation
- Generate a fold-change plot

Files we will use:
- Input: ../data/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz
- Input: ../data/human_gene_annotations.tsv.gz
- Output: outputs/fold_change_report.txt
- Output: outputs/top_fold_changes.png
"""

from __future__ import annotations

import csv
import gzip
import math
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
REPORT_PATH = PRESENTABLE_DIR / "fold_change_report.txt"
PLOT_PATH = PRESENTABLE_DIR / "top_fold_changes.png"
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

# Based on the GEO series organization for GSE164073.
GROUPS = {
    "Cornea_mock": ["GSM4996084", "GSM4996085", "GSM4996086"],
    "Cornea_CoV2": ["GSM4996087", "GSM4996088", "GSM4996089"],
    "Limbus_mock": ["GSM4996090", "GSM4996091", "GSM4996092"],
    "Limbus_CoV2": ["GSM4996093", "GSM4996094", "GSM4996095"],
    "Sclera_mock": ["GSM4996096", "GSM4996097", "GSM4996098"],
    "Sclera_CoV2": ["GSM4996099", "GSM4996100", "GSM4996101"],
}
CASE_GROUP = "Cornea_CoV2"
CONTROL_GROUP = "Cornea_mock"
PSEUDOCOUNT = 1.0


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


def load_annotations(path: Path = ANNOTATION_PATH) -> dict[str, str]:
    annotations: dict[str, str] = {}
    if not path.exists():
        return annotations

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return annotations

        gene_id_key = first_matching_key(reader.fieldnames, ["GeneID", "Gene ID", "gene_id"])
        symbol_key = first_matching_key(reader.fieldnames, ["Symbol", "GeneSymbol", "gene_symbol"])
        if gene_id_key is None:
            return annotations

        for row in reader:
            gene_id = row.get(gene_id_key, "").strip()
            if not gene_id:
                continue
            symbol = row.get(symbol_key, "").strip() if symbol_key else ""
            annotations[gene_id] = f"{symbol} [{gene_id}]" if symbol else gene_id

    return annotations


def parse_count_matrix(path: Path = COUNT_MATRIX_PATH, annotation_path: Path = ANNOTATION_PATH) -> dict[str, object]:
    annotations = load_annotations(annotation_path)

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_names = header[1:]

        records: list[dict[str, object]] = []
        for row in reader:
            if not row:
                continue

            gene_id = row[0]
            counts = [int(value) for value in row[1:]]
            records.append(
                {
                    "gene_id": gene_id,
                    "label": annotations.get(gene_id, gene_id),
                    "counts": counts,
                }
            )

    return {
        "sample_names": sample_names,
        "records": records,
    }


def sample_indices(sample_names: list[str], group_name: str) -> list[int]:
    members = GROUPS[group_name]
    return [sample_names.index(sample_name) for sample_name in members]


def mean(values: list[int]) -> float:
    return sum(values) / len(values) if values else 0.0


def compute_fold_changes(matrix: dict[str, object], case_group: str = CASE_GROUP, control_group: str = CONTROL_GROUP) -> dict[str, object]:
    sample_names = matrix["sample_names"]
    records = matrix["records"]
    case_indices = sample_indices(sample_names, case_group)
    control_indices = sample_indices(sample_names, control_group)

    results: list[dict[str, object]] = []
    for record in records:
        counts = record["counts"]
        case_mean = mean([counts[index] for index in case_indices])
        control_mean = mean([counts[index] for index in control_indices])
        log2_fold_change = math.log2((case_mean + PSEUDOCOUNT) / (control_mean + PSEUDOCOUNT))
        results.append(
            {
                "gene_id": record["gene_id"],
                "label": record["label"],
                "case_mean": case_mean,
                "control_mean": control_mean,
                "log2_fold_change": log2_fold_change,
            }
        )

    upregulated = sorted(results, key=lambda row: (-row["log2_fold_change"], str(row["label"])))[:10]
    downregulated = sorted(results, key=lambda row: (row["log2_fold_change"], str(row["label"])))[:10]

    return {
        "case_group": case_group,
        "control_group": control_group,
        "case_samples": [sample_names[index] for index in case_indices],
        "control_samples": [sample_names[index] for index in control_indices],
        "results": results,
        "upregulated": upregulated,
        "downregulated": downregulated,
    }


def build_report_text(result: dict[str, object]) -> str:
    lines = [
        "Mini-project 3: Fold change comparison",
        "Dataset: GEO accession GSE164073",
        f"Comparison: {result['case_group']} vs {result['control_group']}",
        f"Case samples: {', '.join(result['case_samples'])}",
        f"Control samples: {', '.join(result['control_samples'])}",
        f"Pseudocount used in log2 fold change: {PSEUDOCOUNT}",
        "",
        "Top upregulated genes",
    ]

    for row in result["upregulated"]:
        lines.append(
            f"{row['label']}: log2FC={row['log2_fold_change']:.3f}, "
            f"case_mean={row['case_mean']:.2f}, control_mean={row['control_mean']:.2f}"
        )

    lines.extend(["", "Top downregulated genes"])
    for row in result["downregulated"]:
        lines.append(
            f"{row['label']}: log2FC={row['log2_fold_change']:.3f}, "
            f"case_mean={row['case_mean']:.2f}, control_mean={row['control_mean']:.2f}"
        )

    return "\n".join(lines) + "\n"


def save_report(result: dict[str, object], destination: Path = REPORT_PATH) -> Path:
    ensure_directories()
    destination.write_text(build_report_text(result), encoding="utf-8")
    return destination


def create_fold_change_plot(result: dict[str, object], destination: Path = PLOT_PATH) -> Path:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Matplotlib is required for the fold-change plot. "
            "Install it with: python3 -m pip install matplotlib"
        ) from error

    upregulated = list(reversed(result["upregulated"][:8]))
    downregulated = result["downregulated"][:8]
    rows = upregulated + downregulated
    labels = [row["label"] for row in rows]
    values = [row["log2_fold_change"] for row in rows]
    colors = ["#2E7D32" if value > 0 else "#C62828" for value in values]

    figure, axis = plt.subplots(figsize=(12, 7))
    axis.barh(labels, values, color=colors)
    axis.axvline(0, color="black", linewidth=1)
    axis.set_title(f"Top Fold Changes: {result['case_group']} vs {result['control_group']}")
    axis.set_xlabel("log2 fold change")
    axis.set_ylabel("Gene")
    axis.grid(axis="x", alpha=0.3, linestyle="--")
    axis.set_axisbelow(True)
    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def run() -> dict[str, Path]:
    download_inputs()
    matrix = parse_count_matrix()
    result = compute_fold_changes(matrix)
    report_path = save_report(result)
    plot_path = create_fold_change_plot(result)
    return {
        "report_path": report_path,
        "plot_path": plot_path,
    }


def main() -> None:
    try:
        downloads = download_inputs()
        matrix = parse_count_matrix()
        result = compute_fold_changes(matrix)
        report_path = save_report(result)
        plot_path = create_fold_change_plot(result)
    except URLError as error:
        print("Download failed. The script is ready, but it could not reach GEO right now.")
        print(f"Reason: {error}")
        print(f"Expected count matrix URL: {COUNT_MATRIX_URL}")
        print(f"Expected annotation URL: {ANNOTATION_URL}")
        return
    except ModuleNotFoundError as error:
        print(error)
        return

    print("RNA-seq mini-project 3 completed.")
    for name, (path, was_downloaded) in downloads.items():
        status = "Downloaded" if was_downloaded else "Using existing"
        print(f"{status} {name}: {path}")
    print(f"Fold change report: {report_path}")
    print(f"Fold change plot: {plot_path}")


if __name__ == "__main__":
    main()
