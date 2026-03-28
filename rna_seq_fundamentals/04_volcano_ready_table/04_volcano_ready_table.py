"""Mini-project: volcano-ready results table.

Description:
Build a volcano-plot-ready differential expression table from a real GEO
RNA-seq comparison and generate a preview volcano plot.

Objectives:
- Reuse the real GEO count matrix and annotation table
- Compute log2 fold change for a real grouped comparison
- Estimate approximate p-values from replicate variability
- Adjust p-values with Benjamini-Hochberg FDR
- Export a volcano-ready results table
- Generate a preview volcano plot

Files we will use:
- Input: ../data/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz
- Input: ../data/human_gene_annotations.tsv.gz
- Output: outputs/volcano_ready_results.tsv
- Output: outputs/volcano_preview.png
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
RESULTS_PATH = PRESENTABLE_DIR / "volcano_ready_results.tsv"
PLOT_PATH = PRESENTABLE_DIR / "volcano_preview.png"
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
FOLD_CHANGE_THRESHOLD = 1.0
ADJ_P_THRESHOLD = 0.05


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

    return {"sample_names": sample_names, "records": records}


def sample_indices(sample_names: list[str], group_name: str) -> list[int]:
    members = GROUPS[group_name]
    return [sample_names.index(sample_name) for sample_name in members]


def mean(values: list[int]) -> float:
    return sum(values) / len(values) if values else 0.0


def sample_variance(values: list[int]) -> float:
    if len(values) < 2:
        return 0.0
    mean_value = mean(values)
    return sum((value - mean_value) ** 2 for value in values) / (len(values) - 1)


def approximate_two_sided_p_value(case_values: list[int], control_values: list[int]) -> float:
    case_mean = mean(case_values)
    control_mean = mean(control_values)
    case_var = sample_variance(case_values)
    control_var = sample_variance(control_values)
    standard_error = math.sqrt(case_var / len(case_values) + control_var / len(control_values))

    if standard_error == 0:
        return 1.0 if case_mean == control_mean else 1e-12

    z_score = (case_mean - control_mean) / standard_error
    return math.erfc(abs(z_score) / math.sqrt(2))


def benjamini_hochberg(rows: list[dict[str, object]], p_key: str, output_key: str) -> None:
    sorted_rows = sorted(enumerate(rows), key=lambda item: float(item[1][p_key]))
    total = len(rows)
    adjusted = [1.0] * total

    running_min = 1.0
    for reverse_rank, (original_index, row) in enumerate(reversed(sorted_rows), start=1):
        rank = total - reverse_rank + 1
        raw_p = float(row[p_key])
        bh_value = raw_p * total / rank
        running_min = min(running_min, bh_value)
        adjusted[original_index] = min(running_min, 1.0)

    for row, adj_p in zip(rows, adjusted):
        row[output_key] = adj_p


def classify_gene(log2_fold_change: float, adjusted_p_value: float) -> str:
    if adjusted_p_value < ADJ_P_THRESHOLD and log2_fold_change >= FOLD_CHANGE_THRESHOLD:
        return "Upregulated"
    if adjusted_p_value < ADJ_P_THRESHOLD and log2_fold_change <= -FOLD_CHANGE_THRESHOLD:
        return "Downregulated"
    return "Not significant"


def build_volcano_results(matrix: dict[str, object], case_group: str = CASE_GROUP, control_group: str = CONTROL_GROUP) -> dict[str, object]:
    sample_names = matrix["sample_names"]
    records = matrix["records"]
    case_indices = sample_indices(sample_names, case_group)
    control_indices = sample_indices(sample_names, control_group)

    results: list[dict[str, object]] = []
    for record in records:
        counts = record["counts"]
        case_values = [counts[index] for index in case_indices]
        control_values = [counts[index] for index in control_indices]
        case_mean = mean(case_values)
        control_mean = mean(control_values)
        log2_fold_change = math.log2((case_mean + PSEUDOCOUNT) / (control_mean + PSEUDOCOUNT))
        p_value = approximate_two_sided_p_value(case_values, control_values)

        results.append(
            {
                "gene_id": record["gene_id"],
                "gene_label": record["label"],
                "case_mean": case_mean,
                "control_mean": control_mean,
                "log2_fold_change": log2_fold_change,
                "p_value": p_value,
            }
        )

    benjamini_hochberg(results, "p_value", "adjusted_p_value")

    for row in results:
        row["neg_log10_adjusted_p"] = -math.log10(max(float(row["adjusted_p_value"]), 1e-300))
        row["classification"] = classify_gene(float(row["log2_fold_change"]), float(row["adjusted_p_value"]))

    sorted_results = sorted(results, key=lambda row: (float(row["adjusted_p_value"]), -abs(float(row["log2_fold_change"]))))
    return {
        "case_group": case_group,
        "control_group": control_group,
        "case_samples": [sample_names[index] for index in case_indices],
        "control_samples": [sample_names[index] for index in control_indices],
        "results": sorted_results,
    }


def write_results_table(result: dict[str, object], destination: Path = RESULTS_PATH) -> Path:
    ensure_directories()
    fieldnames = [
        "gene_id",
        "gene_label",
        "case_mean",
        "control_mean",
        "log2_fold_change",
        "p_value",
        "adjusted_p_value",
        "neg_log10_adjusted_p",
        "classification",
    ]

    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in result["results"]:
            writer.writerow(
                {
                    "gene_id": row["gene_id"],
                    "gene_label": row["gene_label"],
                    "case_mean": f"{row['case_mean']:.6f}",
                    "control_mean": f"{row['control_mean']:.6f}",
                    "log2_fold_change": f"{row['log2_fold_change']:.6f}",
                    "p_value": f"{row['p_value']:.6e}",
                    "adjusted_p_value": f"{row['adjusted_p_value']:.6e}",
                    "neg_log10_adjusted_p": f"{row['neg_log10_adjusted_p']:.6f}",
                    "classification": row["classification"],
                }
            )

    return destination


def create_volcano_plot(result: dict[str, object], destination: Path = PLOT_PATH) -> Path:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Matplotlib is required for the volcano preview plot. "
            "Install it with: python3 -m pip install matplotlib"
        ) from error

    rows = result["results"]
    colors = {
        "Upregulated": "#C62828",
        "Downregulated": "#1565C0",
        "Not significant": "#9E9E9E",
    }

    figure, axis = plt.subplots(figsize=(10, 7))
    for classification in ["Not significant", "Downregulated", "Upregulated"]:
        subset = [row for row in rows if row["classification"] == classification]
        axis.scatter(
            [row["log2_fold_change"] for row in subset],
            [row["neg_log10_adjusted_p"] for row in subset],
            s=12,
            alpha=0.65,
            color=colors[classification],
            label=classification,
        )

    axis.axvline(FOLD_CHANGE_THRESHOLD, color="black", linestyle="--", linewidth=1)
    axis.axvline(-FOLD_CHANGE_THRESHOLD, color="black", linestyle="--", linewidth=1)
    axis.axhline(-math.log10(ADJ_P_THRESHOLD), color="black", linestyle="--", linewidth=1)
    axis.set_title(f"Volcano Preview: {result['case_group']} vs {result['control_group']}")
    axis.set_xlabel("log2 fold change")
    axis.set_ylabel("-log10 adjusted p-value")
    axis.grid(alpha=0.2, linestyle="--")
    axis.legend()

    top_hits = [
        row for row in rows
        if row["classification"] != "Not significant"
    ][:8]
    for row in top_hits:
        axis.annotate(
            row["gene_label"],
            (row["log2_fold_change"], row["neg_log10_adjusted_p"]),
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )

    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def run() -> dict[str, Path]:
    download_inputs()
    matrix = parse_count_matrix()
    result = build_volcano_results(matrix)
    results_path = write_results_table(result)
    plot_path = create_volcano_plot(result)
    return {"results_path": results_path, "plot_path": plot_path}


def main() -> None:
    try:
        downloads = download_inputs()
        matrix = parse_count_matrix()
        result = build_volcano_results(matrix)
        results_path = write_results_table(result)
        plot_path = create_volcano_plot(result)
    except URLError as error:
        print("Download failed. The script is ready, but it could not reach GEO right now.")
        print(f"Reason: {error}")
        print(f"Expected count matrix URL: {COUNT_MATRIX_URL}")
        print(f"Expected annotation URL: {ANNOTATION_URL}")
        return
    except ModuleNotFoundError as error:
        write_results_table(build_volcano_results(parse_count_matrix()))
        print(error)
        return

    print("RNA-seq mini-project 4 completed.")
    for name, (path, was_downloaded) in downloads.items():
        status = "Downloaded" if was_downloaded else "Using existing"
        print(f"{status} {name}: {path}")
    print(f"Volcano-ready table: {results_path}")
    print(f"Volcano preview plot: {plot_path}")


if __name__ == "__main__":
    main()
