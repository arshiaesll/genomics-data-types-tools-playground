"""Mini-project: HISAT2 alignment workflow.

Description:
Download the official HISAT2 paired-end example reads and chromosome 22
reference region, build a small index locally, run the aligner, and summarize
the results in a way that is easy to showcase.

Official example data used:
- https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reference/22_20-21M.fa
- https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reads/reads_1.fa
- https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reads/reads_2.fa

Outputs:
- outputs/hisat2_alignment_summary.txt
- outputs/hisat2_alignment.stderr.txt
- outputs/aligned_reads.sam
- outputs/hisat2_alignment_breakdown.png
- outputs/hisat2_mapq_histogram.png
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from statistics import mean, median
from urllib.error import URLError
from urllib.request import urlretrieve
import os
import re
import shutil
import subprocess
import textwrap


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = WORKFLOW_DIR / "data"
REFERENCE_DIR = DATA_DIR / "reference"
READS_DIR = DATA_DIR / "reads"
OUTPUT_DIR = PROJECT_DIR / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
TECHNICAL_DIR = OUTPUT_DIR / "technical"
SUMMARY_PATH = PRESENTABLE_DIR / "hisat2_alignment_summary.txt"
ALIGNMENT_STDERR_PATH = TECHNICAL_DIR / "hisat2_alignment.stderr.txt"
SAM_PATH = TECHNICAL_DIR / "aligned_reads.sam"
BREAKDOWN_PLOT_PATH = PRESENTABLE_DIR / "hisat2_alignment_breakdown.png"
MAPQ_PLOT_PATH = PRESENTABLE_DIR / "hisat2_mapq_histogram.png"
REFERENCE_FASTA_PATH = REFERENCE_DIR / "22_20-21M.fa"
READ_1_PATH = READS_DIR / "reads_1.fa"
READ_2_PATH = READS_DIR / "reads_2.fa"
INDEX_PREFIX = REFERENCE_DIR / "22_20-21M"
THREADS = 2

REFERENCE_URL = (
    "https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/"
    "example/reference/22_20-21M.fa"
)
READ_1_URL = (
    "https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/"
    "example/reads/reads_1.fa"
)
READ_2_URL = (
    "https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/"
    "example/reads/reads_2.fa"
)


def ensure_directories() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    READS_DIR.mkdir(parents=True, exist_ok=True)
    PRESENTABLE_DIR.mkdir(parents=True, exist_ok=True)
    TECHNICAL_DIR.mkdir(parents=True, exist_ok=True)


def detect_tool(tool_name: str) -> str | None:
    return shutil.which(tool_name)


def index_files_for_prefix(prefix: Path) -> list[Path]:
    return sorted(prefix.parent.glob(f"{prefix.name}*.ht2"))


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
        str(REFERENCE_FASTA_PATH): download_if_missing(REFERENCE_URL, REFERENCE_FASTA_PATH),
        str(READ_1_PATH): download_if_missing(READ_1_URL, READ_1_PATH),
        str(READ_2_PATH): download_if_missing(READ_2_URL, READ_2_PATH),
    }


def build_hisat2_index(index_prefix: Path, reference_fasta: Path) -> subprocess.CompletedProcess[str]:
    command = ["hisat2-build", str(reference_fasta), str(index_prefix)]
    return subprocess.run(command, check=False, capture_output=True, text=True, env=tool_subprocess_env())


def build_hisat2_command(index_prefix: Path, read_1: Path, read_2: Path, sam_path: Path) -> list[str]:
    return [
        "hisat2",
        "-p",
        str(THREADS),
        "--dta",
        "-f",
        "-x",
        str(index_prefix),
        "-1",
        str(read_1),
        "-2",
        str(read_2),
        "-S",
        str(sam_path),
    ]


def run_alignment(index_prefix: Path, read_1: Path, read_2: Path) -> subprocess.CompletedProcess[str]:
    command = build_hisat2_command(index_prefix=index_prefix, read_1=read_1, read_2=read_2, sam_path=SAM_PATH)
    return subprocess.run(command, check=False, capture_output=True, text=True, env=tool_subprocess_env())


def parse_hisat2_metrics(stderr_text: str) -> dict[str, float | int]:
    metrics: dict[str, float | int] = {}
    patterns: dict[str, str] = {
        "total_reads": r"(\d+)\s+reads; of these:",
        "paired_reads": r"(\d+)\s+\([^)]+\)\s+were paired; of these:",
        "concordant_zero_times": r"(\d+)\s+\([^)]+\)\s+aligned concordantly 0 times",
        "concordant_exactly_once": r"(\d+)\s+\([^)]+\)\s+aligned concordantly exactly 1 time",
        "concordant_multiple_times": r"(\d+)\s+\([^)]+\)\s+aligned concordantly >1 times",
        "overall_alignment_rate": r"([\d.]+)% overall alignment rate",
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, stderr_text)
        if not match:
            continue
        raw_value = match.group(1)
        metrics[key] = float(raw_value) if "." in raw_value else int(raw_value)

    paired_reads = metrics.get("paired_reads")
    if isinstance(paired_reads, int) and paired_reads > 0:
        for numerator_key, percent_key in [
            ("concordant_exactly_once", "concordant_exactly_once_percent"),
            ("concordant_multiple_times", "concordant_multiple_times_percent"),
            ("concordant_zero_times", "concordant_zero_times_percent"),
        ]:
            numerator = metrics.get(numerator_key)
            if isinstance(numerator, int):
                metrics[percent_key] = (numerator / paired_reads) * 100.0

    return metrics


def parse_sam_metrics(sam_path: Path) -> dict[str, object]:
    if not sam_path.exists():
        return {}

    total_records = 0
    mapped_records = 0
    secondary_records = 0
    supplementary_records = 0
    spliced_records = 0
    mapq_values: list[int] = []
    reference_counts: Counter[str] = Counter()

    with sam_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@"):
                continue

            total_records += 1
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            reference_name = fields[2]
            mapq = int(fields[4])
            cigar = fields[5]

            if flag & 0x100:
                secondary_records += 1
            if flag & 0x800:
                supplementary_records += 1
            if flag & 0x4:
                continue

            mapped_records += 1
            reference_counts[reference_name] += 1
            mapq_values.append(mapq)
            if "N" in cigar:
                spliced_records += 1

    result: dict[str, object] = {
        "sam_records": total_records,
        "mapped_records": mapped_records,
        "unmapped_records": total_records - mapped_records,
        "secondary_records": secondary_records,
        "supplementary_records": supplementary_records,
        "spliced_records": spliced_records,
        "top_references": reference_counts.most_common(5),
    }

    if mapq_values:
        result["mean_mapq"] = mean(mapq_values)
        result["median_mapq"] = median(mapq_values)
        result["mapq_values"] = mapq_values

    return result


def alignment_interpretation(metrics: dict[str, float | int], sam_metrics: dict[str, object]) -> str:
    alignment_rate = metrics.get("overall_alignment_rate")
    spliced_records = sam_metrics.get("spliced_records", 0)

    if not isinstance(alignment_rate, float):
        return "Alignment interpretation is unavailable because HISAT2 did not finish cleanly."
    if alignment_rate >= 95.0:
        quality_text = "The alignment rate is excellent for this example run."
    elif alignment_rate >= 80.0:
        quality_text = "The alignment rate is strong and the sample looks usable."
    elif alignment_rate >= 60.0:
        quality_text = "The sample is usable, but the mapping rate should be reviewed before moving downstream."
    else:
        quality_text = "The alignment rate is low, so read quality and reference/sample compatibility should be checked."

    if isinstance(spliced_records, int) and spliced_records > 0:
        return f"{quality_text} The SAM output also includes spliced alignments, which is consistent with transcript-aware mapping."
    return f"{quality_text} This particular demo mostly acts as a lightweight alignment showcase rather than a rich splice-junction example."


def configure_matplotlib() -> None:
    mpl_config_dir = TECHNICAL_DIR / ".mplconfig"
    cache_dir = TECHNICAL_DIR / ".cache"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))


def tool_subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    return env


def create_alignment_breakdown_plot(metrics: dict[str, float | int], destination: Path = BREAKDOWN_PLOT_PATH) -> Path | None:
    paired_reads = metrics.get("paired_reads")
    if not isinstance(paired_reads, int) or paired_reads == 0:
        return None

    values = [
        float(metrics.get("concordant_exactly_once_percent", 0.0)),
        float(metrics.get("concordant_multiple_times_percent", 0.0)),
        float(metrics.get("concordant_zero_times_percent", 0.0)),
    ]
    labels = [
        "Concordant once",
        "Concordant >1",
        "Concordant 0",
    ]

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(7.2, 4.5))
    bars = axis.bar(labels, values, color=["#2f7d4a", "#9ebc3c", "#c85c43"])
    axis.set_ylim(0, 100)
    axis.set_ylabel("Percent of read pairs")
    axis.set_title("HISAT2 Concordant Alignment Breakdown")

    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value + 1.0, f"{value:.1f}%", ha="center")

    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def create_mapq_histogram(sam_metrics: dict[str, object], destination: Path = MAPQ_PLOT_PATH) -> Path | None:
    mapq_values = sam_metrics.get("mapq_values")
    if not isinstance(mapq_values, list) or not mapq_values:
        return None

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(7.2, 4.5))
    max_mapq = max(mapq_values)
    bins = range(0, max_mapq + 5, 5)
    axis.hist(mapq_values, bins=bins, color="#2b6cb0", edgecolor="white")
    axis.set_xlabel("MAPQ")
    axis.set_ylabel("Aligned SAM records")
    axis.set_title("Mapping Quality Distribution")
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def explain_scenario() -> str:
    return textwrap.dedent(
        """\
        Scenario
        We use the official HISAT2 example files: a small chromosome 22
        reference region plus paired-end reads bundled with the HISAT2 project.
        That keeps the run lightweight while still using real public sequencing
        data and the real aligner.

        Why HISAT2 is a good fit here
        HISAT2 is a splice-aware aligner commonly used for RNA-seq workflows.
        We run it with `--dta` so the output is friendly for downstream
        transcript-aware tools such as StringTie.

        Decision we want to make
        After alignment, we want to know whether the example has a strong
        overall mapping rate and what the mapped reads look like in the SAM
        output. Those are the first checks before moving on to counting or
        assembly.
        """
    ).strip()


def build_summary_text(
    hisat2_path: str | None,
    hisat2_build_path: str | None,
    download_status: dict[str, str] | None,
    index_build_result: subprocess.CompletedProcess[str] | None,
    run_result: subprocess.CompletedProcess[str] | None,
    metrics: dict[str, float | int],
    sam_metrics: dict[str, object],
    breakdown_plot: Path | None,
    mapq_plot: Path | None,
    error_message: str | None,
) -> str:
    lines = [
        "Mini-project 1: HISAT2 alignment showcase",
        "",
        explain_scenario(),
        "",
        "Official data sources",
        f"Reference FASTA: {REFERENCE_URL}",
        f"Read 1 FASTA: {READ_1_URL}",
        f"Read 2 FASTA: {READ_2_URL}",
        "",
        "Local workflow inputs",
        f"HISAT2 path: {hisat2_path or 'not found in PATH'}",
        f"HISAT2-build path: {hisat2_build_path or 'not found in PATH'}",
        f"Reference FASTA: {REFERENCE_FASTA_PATH}",
        f"Read 1 path: {READ_1_PATH}",
        f"Read 2 path: {READ_2_PATH}",
        f"Index prefix: {INDEX_PREFIX}",
        f"Index files found: {len(index_files_for_prefix(INDEX_PREFIX))}",
    ]

    if download_status:
        lines.extend(["", "Download status"])
        for path_text, status in download_status.items():
            lines.append(f"{Path(path_text).name}: {status}")

    lines.extend(["", "Command used"])
    lines.append(
        " ".join(
            build_hisat2_command(
                index_prefix=INDEX_PREFIX,
                read_1=READ_1_PATH,
                read_2=READ_2_PATH,
                sam_path=SAM_PATH,
            )
        )
    )

    if index_build_result is not None:
        lines.extend(
            [
                "",
                "Index build",
                f"Return code: {index_build_result.returncode}",
            ]
        )
        if index_build_result.stdout.strip():
            lines.append("Index build stdout was captured.")
        if index_build_result.stderr.strip():
            lines.append("Index build stderr was captured.")

    if error_message:
        lines.extend(["", "Scenario outcome", error_message])
        return "\n".join(lines).rstrip() + "\n"

    if run_result is None:
        lines.extend(
            [
                "",
                "Scenario outcome",
                "The alignment did not run.",
            ]
        )
        return "\n".join(lines).rstrip() + "\n"

    lines.extend(
        [
            "",
            "Scenario outcome",
            f"Return code: {run_result.returncode}",
            f"SAM output: {SAM_PATH}",
            f"HISAT2 stderr log: {ALIGNMENT_STDERR_PATH}",
            f"Breakdown figure: {breakdown_plot or 'not created'}",
            f"MAPQ figure: {mapq_plot or 'not created'}",
        ]
    )

    for metric_name in [
        "total_reads",
        "paired_reads",
        "concordant_exactly_once",
        "concordant_multiple_times",
        "concordant_zero_times",
        "overall_alignment_rate",
    ]:
        value = metrics.get(metric_name)
        if value is None:
            continue
        label = metric_name.replace("_", " ").capitalize()
        if isinstance(value, float):
            lines.append(f"{label}: {value:.2f}")
        else:
            lines.append(f"{label}: {value}")

    lines.extend(
        [
            "",
            "SAM-derived summary",
            f"SAM alignment records: {sam_metrics.get('sam_records', 'n/a')}",
            f"Mapped alignment records: {sam_metrics.get('mapped_records', 'n/a')}",
            f"Unmapped alignment records: {sam_metrics.get('unmapped_records', 'n/a')}",
            f"Secondary alignment records: {sam_metrics.get('secondary_records', 'n/a')}",
            f"Supplementary alignment records: {sam_metrics.get('supplementary_records', 'n/a')}",
            f"Spliced alignment records (CIGAR contains N): {sam_metrics.get('spliced_records', 'n/a')}",
        ]
    )

    mean_mapq = sam_metrics.get("mean_mapq")
    median_mapq = sam_metrics.get("median_mapq")
    if isinstance(mean_mapq, (int, float)):
        lines.append(f"Mean MAPQ: {mean_mapq:.2f}")
    if isinstance(median_mapq, (int, float)):
        lines.append(f"Median MAPQ: {median_mapq:.2f}")

    top_references = sam_metrics.get("top_references", [])
    if isinstance(top_references, list) and top_references:
        lines.extend(["", "Top reference targets in the SAM output"])
        for reference_name, count in top_references:
            lines.append(f"- {reference_name}: {count} records")

    lines.extend(
        [
            "",
            "Interpretation",
            alignment_interpretation(metrics, sam_metrics),
            "",
            "Raw HISAT2 summary",
            run_result.stderr.strip() or "(no stderr output captured)",
        ]
    )

    return "\n".join(lines).rstrip() + "\n"


def save_summary(text: str, destination: Path = SUMMARY_PATH) -> Path:
    destination.write_text(text, encoding="utf-8")
    return destination


def main() -> None:
    ensure_directories()
    hisat2_path = detect_tool("hisat2")
    hisat2_build_path = detect_tool("hisat2-build")

    download_status: dict[str, str] | None = None
    index_build_result: subprocess.CompletedProcess[str] | None = None
    run_result: subprocess.CompletedProcess[str] | None = None
    metrics: dict[str, float | int] = {}
    sam_metrics: dict[str, object] = {}
    breakdown_plot: Path | None = None
    mapq_plot: Path | None = None
    error_message: str | None = None

    if not hisat2_path or not hisat2_build_path:
        error_message = "HISAT2 and hisat2-build must both be installed before this workflow can run."
    else:
        try:
            download_status = ensure_official_example_data()

            if not index_files_for_prefix(INDEX_PREFIX):
                index_build_result = build_hisat2_index(index_prefix=INDEX_PREFIX, reference_fasta=REFERENCE_FASTA_PATH)
                if index_build_result.returncode != 0:
                    error_message = "HISAT2 index building failed. Check the captured build output in the terminal."

            if error_message is None:
                run_result = run_alignment(index_prefix=INDEX_PREFIX, read_1=READ_1_PATH, read_2=READ_2_PATH)
                ALIGNMENT_STDERR_PATH.write_text(run_result.stderr, encoding="utf-8")
                metrics = parse_hisat2_metrics(run_result.stderr)

                if run_result.returncode != 0:
                    error_message = "HISAT2 alignment failed. Check the captured stderr log for details."
                else:
                    sam_metrics = parse_sam_metrics(SAM_PATH)
                    breakdown_plot = create_alignment_breakdown_plot(metrics)
                    mapq_plot = create_mapq_histogram(sam_metrics)

        except RuntimeError as error:
            error_message = str(error)

    summary_text = build_summary_text(
        hisat2_path=hisat2_path,
        hisat2_build_path=hisat2_build_path,
        download_status=download_status,
        index_build_result=index_build_result,
        run_result=run_result,
        metrics=metrics,
        sam_metrics=sam_metrics,
        breakdown_plot=breakdown_plot,
        mapq_plot=mapq_plot,
        error_message=error_message,
    )
    save_summary(summary_text)
    print(summary_text)


if __name__ == "__main__":
    main()
