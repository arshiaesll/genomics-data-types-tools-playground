"""Mini-project: HTSeq counting workflow.

Description:
Use the official HTSeq example alignment and annotation files to generate a
real multi-gene count table and visualize the result.

Official example data used:
- https://raw.githubusercontent.com/htseq/htseq/master/example_data/yeast_RNASeq_excerpt.sam
- https://raw.githubusercontent.com/htseq/htseq/master/example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz

Outputs:
- outputs/htseq_counts.tsv
- outputs/htseq_counting_summary.txt
- outputs/htseq_top_genes.png
- outputs/htseq_special_counters.png
"""

from __future__ import annotations

from pathlib import Path
from urllib.error import URLError
from urllib.request import urlretrieve
import gzip
import os
import shutil
import subprocess
import textwrap


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = WORKFLOW_DIR / "data"
READS_DIR = DATA_DIR / "reads"
REFERENCE_DIR = DATA_DIR / "reference"
OUTPUT_DIR = PROJECT_DIR / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
SAM_PATH = READS_DIR / "yeast_RNASeq_excerpt.sam"
GTF_PATH = REFERENCE_DIR / "Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz"
COUNTS_PATH = PRESENTABLE_DIR / "htseq_counts.tsv"
SUMMARY_PATH = PRESENTABLE_DIR / "htseq_counting_summary.txt"
TOP_GENES_PLOT_PATH = PRESENTABLE_DIR / "htseq_top_genes.png"
SPECIAL_PLOT_PATH = PRESENTABLE_DIR / "htseq_special_counters.png"

SAM_URL = "https://raw.githubusercontent.com/htseq/htseq/master/example_data/yeast_RNASeq_excerpt.sam"
GTF_URL = (
    "https://raw.githubusercontent.com/htseq/htseq/master/example_data/"
    "Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz"
)


def ensure_directories() -> None:
    READS_DIR.mkdir(parents=True, exist_ok=True)
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    PRESENTABLE_DIR.mkdir(parents=True, exist_ok=True)


def configure_matplotlib() -> None:
    mpl_config_dir = OUTPUT_DIR / ".mplconfig"
    cache_dir = OUTPUT_DIR / ".cache"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))


def detect_tool(tool_name: str) -> str | None:
    return shutil.which(tool_name)


def download_if_missing(url: str, destination: Path) -> str:
    if destination.exists():
        return "reused"
    try:
        urlretrieve(url, destination)
    except URLError as error:
        raise RuntimeError(f"Could not download {url}: {error}") from error
    return "downloaded"


def ensure_official_example_data() -> dict[str, str]:
    return {
        str(SAM_PATH): download_if_missing(SAM_URL, SAM_PATH),
        str(GTF_PATH): download_if_missing(GTF_URL, GTF_PATH),
    }


def build_htseq_command(sam_path: Path, gtf_path: Path, counts_path: Path) -> list[str]:
    return [
        "htseq-count",
        "-f",
        "sam",
        "-r",
        "pos",
        "-s",
        "no",
        "-t",
        "exon",
        "-i",
        "gene_id",
        "--additional-attr",
        "gene_name",
        "-c",
        str(counts_path),
        str(sam_path),
        str(gtf_path),
    ]


def run_htseq(sam_path: Path, gtf_path: Path, counts_path: Path) -> subprocess.CompletedProcess[str]:
    command = build_htseq_command(sam_path=sam_path, gtf_path=gtf_path, counts_path=counts_path)
    return subprocess.run(command, check=False, capture_output=True, text=True)


def parse_counts_table(path: Path) -> tuple[list[dict[str, object]], dict[str, int]]:
    gene_rows: list[dict[str, object]] = []
    special_counts: dict[str, int] = {}

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue

            gene_id = fields[0]
            count = int(fields[-1])

            if gene_id.startswith("__"):
                special_counts[gene_id] = count
                continue

            gene_name = fields[1] if len(fields) > 2 and fields[1] else gene_id
            gene_rows.append({"gene_id": gene_id, "gene_name": gene_name, "count": count})

    gene_rows.sort(key=lambda row: int(row["count"]), reverse=True)
    return gene_rows, special_counts


def summarize_input_data(sam_path: Path, gtf_path: Path) -> dict[str, object]:
    sam_records = 0
    references: set[str] = set()
    with sam_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@SQ"):
                for field in line.rstrip("\n").split("\t"):
                    if field.startswith("SN:"):
                        references.add(field[3:])
            elif not line.startswith("@"):
                sam_records += 1

    gtf_lines = 0
    exon_features = 0
    genes: set[str] = set()
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            gtf_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] == "exon":
                exon_features += 1
                attr_text = fields[8]
                for chunk in attr_text.split(";"):
                    chunk = chunk.strip()
                    if chunk.startswith("gene_id"):
                        genes.add(chunk.split('"')[1])
                        break

    return {
        "sam_records": sam_records,
        "reference_count": len(references),
        "gtf_lines": gtf_lines,
        "exon_features": exon_features,
        "gene_feature_count": len(genes),
    }


def create_top_genes_plot(gene_rows: list[dict[str, object]], destination: Path = TOP_GENES_PLOT_PATH) -> Path | None:
    top_rows = [row for row in gene_rows if int(row["count"]) > 0][:12]
    if not top_rows:
        return None

    labels = [str(row["gene_name"]) for row in top_rows]
    values = [int(row["count"]) for row in top_rows]

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(9, 5))
    bars = axis.bar(labels, values, color="#2f7d4a")
    axis.set_ylabel("Read counts")
    axis.set_title("Top HTSeq Gene Counts")
    axis.tick_params(axis="x", rotation=30)

    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value + max(values) * 0.01, str(value), ha="center", fontsize=9)

    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def create_special_counter_plot(special_counts: dict[str, int], destination: Path = SPECIAL_PLOT_PATH) -> Path | None:
    if not special_counts:
        return None

    labels = [label.replace("__", "") for label in special_counts]
    values = [special_counts[label] for label in special_counts]

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(8.4, 4.6))
    bars = axis.bar(labels, values, color="#2b6cb0")
    axis.set_ylabel("Alignment records")
    axis.set_title("HTSeq Special Counters")
    axis.tick_params(axis="x", rotation=25)

    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value + max(1, max(values)) * 0.01, str(value), ha="center", fontsize=9)

    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def explain_scenario() -> str:
    return textwrap.dedent(
        """\
        Scenario
        We switch this mini-project to the official HTSeq example dataset so the
        result is a proper multi-gene count table rather than a single-feature
        demo.

        What the files represent
        The alignment file is a real yeast RNA-seq excerpt in SAM format, and
        the annotation is a real Saccharomyces cerevisiae GTF. HTSeq uses the
        annotation to decide which aligned reads overlap which gene exons.

        Decision we want to make
        We want to inspect the resulting gene count table, identify the most
        highly counted genes in the excerpt, and check HTSeq's special counters
        such as `__no_feature`.
        """
    ).strip()


def build_summary_text(
    htseq_path: str | None,
    download_status: dict[str, str] | None,
    input_summary: dict[str, object],
    command: list[str] | None,
    run_result: subprocess.CompletedProcess[str] | None,
    gene_rows: list[dict[str, object]],
    special_counts: dict[str, int],
    top_gene_plot: Path | None,
    special_plot: Path | None,
    error_message: str | None,
) -> str:
    lines = [
        "Mini-project 3: HTSeq counting showcase",
        "",
        explain_scenario(),
        "",
        "Official data sources",
        f"SAM: {SAM_URL}",
        f"GTF: {GTF_URL}",
        "",
        "Workflow inputs",
        f"HTSeq path: {htseq_path or 'not found in PATH'}",
        f"SAM input: {SAM_PATH}",
        f"GTF input: {GTF_PATH}",
        f"Alignment records in SAM: {input_summary.get('sam_records', 'n/a')}",
        f"Reference sequences in SAM header: {input_summary.get('reference_count', 'n/a')}",
        f"GTF feature lines: {input_summary.get('gtf_lines', 'n/a')}",
        f"Exon features in GTF: {input_summary.get('exon_features', 'n/a')}",
        f"Distinct genes in GTF exons: {input_summary.get('gene_feature_count', 'n/a')}",
    ]

    if download_status:
        lines.extend(["", "Download status"])
        for path_text, status in download_status.items():
            lines.append(f"{Path(path_text).name}: {status}")

    if command:
        lines.extend(["", "Command used", " ".join(command)])

    if error_message:
        lines.extend(["", "Scenario outcome", error_message])
        return "\n".join(lines).rstrip() + "\n"

    if run_result is None:
        lines.extend(["", "Scenario outcome", "HTSeq did not run."])
        return "\n".join(lines).rstrip() + "\n"

    nonzero_gene_rows = [row for row in gene_rows if int(row["count"]) > 0]
    total_assigned = sum(int(row["count"]) for row in nonzero_gene_rows)

    lines.extend(
        [
            "",
            "Scenario outcome",
            f"Return code: {run_result.returncode}",
            f"Counts table: {COUNTS_PATH}",
            f"Top genes figure: {top_gene_plot or 'not created'}",
            f"Special counter figure: {special_plot or 'not created'}",
            f"Genes with nonzero counts: {len(nonzero_gene_rows)}",
            f"Total assigned counts across nonzero genes: {total_assigned}",
        ]
    )

    lines.extend(["", "Top counted genes"])
    for row in nonzero_gene_rows[:12]:
        lines.append(f"- {row['gene_name']} ({row['gene_id']}): {row['count']}")

    lines.extend(["", "HTSeq special counters"])
    for key, value in special_counts.items():
        lines.append(f"- {key}: {value}")

    lines.extend(
        [
            "",
            "Interpretation",
            "This is the more typical HTSeq output shape: many genes with counts, plus a small set of special counters that explain where unassigned reads went.",
        ]
    )

    stderr_text = run_result.stderr.strip()
    if stderr_text:
        lines.extend(["", "HTSeq stderr", stderr_text])

    return "\n".join(lines).rstrip() + "\n"


def save_summary(text: str) -> None:
    SUMMARY_PATH.write_text(text, encoding="utf-8")


def main() -> None:
    ensure_directories()
    htseq_path = detect_tool("htseq-count")
    download_status: dict[str, str] | None = None
    input_summary: dict[str, object] = {}
    command: list[str] | None = None
    run_result: subprocess.CompletedProcess[str] | None = None
    gene_rows: list[dict[str, object]] = []
    special_counts: dict[str, int] = {}
    top_gene_plot: Path | None = None
    special_plot: Path | None = None
    error_message: str | None = None

    if not htseq_path:
        error_message = "htseq-count must be installed before this workflow can run."
    else:
        try:
            download_status = ensure_official_example_data()
            input_summary = summarize_input_data(SAM_PATH, GTF_PATH)
            command = build_htseq_command(sam_path=SAM_PATH, gtf_path=GTF_PATH, counts_path=COUNTS_PATH)
            run_result = run_htseq(sam_path=SAM_PATH, gtf_path=GTF_PATH, counts_path=COUNTS_PATH)
            if run_result.returncode != 0:
                error_message = "HTSeq counting failed. Check the terminal output for details."
            else:
                gene_rows, special_counts = parse_counts_table(COUNTS_PATH)
                top_gene_plot = create_top_genes_plot(gene_rows)
                special_plot = create_special_counter_plot(special_counts)
        except RuntimeError as error:
            error_message = str(error)

    summary_text = build_summary_text(
        htseq_path=htseq_path,
        download_status=download_status,
        input_summary=input_summary,
        command=command,
        run_result=run_result,
        gene_rows=gene_rows,
        special_counts=special_counts,
        top_gene_plot=top_gene_plot,
        special_plot=special_plot,
        error_message=error_message,
    )
    save_summary(summary_text)
    print(summary_text)


if __name__ == "__main__":
    main()
