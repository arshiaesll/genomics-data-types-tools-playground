"""Mini-project: STAR alignment workflow.

Description:
Reuse the official HISAT2 example reads and small chromosome 22 reference
region, align them with STAR, and summarize mapping plus splice-junction
outputs.

Official example data used:
- https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reference/22_20-21M.fa
- https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reads/reads_1.fa
- https://raw.githubusercontent.com/DaehwanKimLab/hisat2/master/example/reads/reads_2.fa

Outputs:
- outputs/star_alignment_summary.txt
- outputs/star_Log.final.out
- outputs/star_SJ.out.tab
- outputs/star_Aligned.out.sam
- outputs/star_mapping_summary.png
- outputs/star_splice_junction_motifs.png
"""

from __future__ import annotations

from collections import Counter
from math import log2
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlretrieve
import os
import platform
import shutil
import subprocess
import textwrap


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = WORKFLOW_DIR / "data"
REFERENCE_DIR = DATA_DIR / "reference"
READS_DIR = DATA_DIR / "reads"
OUTPUT_DIR = PROJECT_DIR / "outputs"
TOOLS_DIR = WORKFLOW_DIR / "tools"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
TECHNICAL_DIR = OUTPUT_DIR / "technical"
SUMMARY_PATH = PRESENTABLE_DIR / "star_alignment_summary.txt"
REFERENCE_FASTA_PATH = REFERENCE_DIR / "22_20-21M.fa"
READ_1_PATH = READS_DIR / "reads_1.fa"
READ_2_PATH = READS_DIR / "reads_2.fa"
STAR_INDEX_DIR = REFERENCE_DIR / "star_chr22_index"
STAR_PREFIX = TECHNICAL_DIR / "star_"
STAR_LOG_PATH = TECHNICAL_DIR / "star_Log.final.out"
STAR_SJ_PATH = TECHNICAL_DIR / "star_SJ.out.tab"
STAR_SAM_PATH = TECHNICAL_DIR / "star_Aligned.out.sam"
MAPPING_PLOT_PATH = PRESENTABLE_DIR / "star_mapping_summary.png"
SPLICE_PLOT_PATH = PRESENTABLE_DIR / "star_splice_junction_motifs.png"
THREADS = 2
LOCAL_STAR_BINARY = TOOLS_DIR / "STAR" / "bin" / "MacOSX_x86_64" / "STAR"

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
    STAR_INDEX_DIR.mkdir(parents=True, exist_ok=True)


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


def detect_tool(tool_name: str) -> str | None:
    return shutil.which(tool_name)


def resolve_star_binary() -> Path | None:
    if platform.system() == "Darwin" and platform.machine() == "arm64" and LOCAL_STAR_BINARY.exists():
        return LOCAL_STAR_BINARY
    star_from_path = detect_tool("STAR")
    if star_from_path:
        return Path(star_from_path)
    if LOCAL_STAR_BINARY.exists():
        return LOCAL_STAR_BINARY
    return None


def star_command_prefix(star_binary: Path) -> list[str]:
    if (
        platform.system() == "Darwin"
        and platform.machine() == "arm64"
        and star_binary.resolve() == LOCAL_STAR_BINARY.resolve()
    ):
        return ["arch", "-x86_64", str(star_binary)]
    return [str(star_binary)]


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


def fasta_total_length(path: Path) -> int:
    total = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def infer_sa_index_nbases(genome_length: int) -> int:
    if genome_length <= 0:
        return 3
    recommended = int(log2(genome_length) / 2 - 1)
    return max(3, min(14, recommended))


def build_star_index_command(star_binary: Path, genome_fasta: Path) -> list[str]:
    genome_length = fasta_total_length(genome_fasta)
    sa_index_nbases = infer_sa_index_nbases(genome_length)
    return star_command_prefix(star_binary) + [
        "--runThreadN",
        str(THREADS),
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        str(STAR_INDEX_DIR),
        "--genomeFastaFiles",
        str(genome_fasta),
        "--genomeSAindexNbases",
        str(sa_index_nbases),
    ]


def build_star_align_command(star_binary: Path) -> list[str]:
    return star_command_prefix(star_binary) + [
        "--runThreadN",
        str(THREADS),
        "--genomeDir",
        str(STAR_INDEX_DIR),
        "--readFilesIn",
        str(READ_1_PATH),
        str(READ_2_PATH),
        "--readFilesType",
        "Fastx",
        "--outSAMtype",
        "SAM",
        "--outSAMattributes",
        "NH",
        "HI",
        "AS",
        "nM",
        "XS",
        "--outFileNamePrefix",
        str(STAR_PREFIX),
    ]


def star_index_exists(genome_dir: Path) -> bool:
    required = [
        genome_dir / "Genome",
        genome_dir / "SA",
        genome_dir / "SAindex",
        genome_dir / "chrLength.txt",
    ]
    return all(path.exists() for path in required)


def build_star_index(star_binary: Path) -> subprocess.CompletedProcess[str]:
    command = build_star_index_command(star_binary, REFERENCE_FASTA_PATH)
    return subprocess.run(command, check=False, capture_output=True, text=True, env=tool_subprocess_env())


def run_star_alignment(star_binary: Path) -> subprocess.CompletedProcess[str]:
    command = build_star_align_command(star_binary)
    return subprocess.run(command, check=False, capture_output=True, text=True, env=tool_subprocess_env())


def parse_star_log(log_path: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    if not log_path.exists():
        return metrics

    with log_path.open(encoding="utf-8") as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = line.split("|", maxsplit=1)
            metrics[key.strip()] = value.strip()
    return metrics


def parse_percent(value: str | None) -> float | None:
    if not value:
        return None
    cleaned = value.replace("%", "").strip()
    try:
        return float(cleaned)
    except ValueError:
        return None


def parse_int(value: str | None) -> int | None:
    if not value:
        return None
    cleaned = value.replace(",", "").strip()
    try:
        return int(cleaned)
    except ValueError:
        return None


def parse_sj_table(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}

    motif_labels = {
        0: "Non-canonical",
        1: "GT/AG",
        2: "CT/AC",
        3: "GC/AG",
        4: "CT/GC",
        5: "AT/AC",
        6: "GT/AT",
    }

    total_rows = 0
    annotated_count = 0
    unique_total = 0
    multimapping_total = 0
    motif_counts: Counter[str] = Counter()

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            total_rows += 1
            motif_code = int(fields[4])
            annotated = int(fields[5])
            unique_reads = int(fields[6])
            multimapping_reads = int(fields[7])

            motif_counts[motif_labels.get(motif_code, f"Motif {motif_code}")] += 1
            if annotated == 1:
                annotated_count += 1
            unique_total += unique_reads
            multimapping_total += multimapping_reads

    return {
        "junction_count": total_rows,
        "annotated_junction_count": annotated_count,
        "unannotated_junction_count": total_rows - annotated_count,
        "unique_junction_reads": unique_total,
        "multimapping_junction_reads": multimapping_total,
        "motif_counts": motif_counts,
    }


def parse_sam_metrics(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}

    total_records = 0
    mapped_records = 0
    spliced_records = 0
    reference_counts: Counter[str] = Counter()

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            total_records += 1
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            reference_name = fields[2]
            cigar = fields[5]
            if flag & 0x4:
                continue
            mapped_records += 1
            reference_counts[reference_name] += 1
            if "N" in cigar:
                spliced_records += 1

    return {
        "sam_records": total_records,
        "mapped_records": mapped_records,
        "unmapped_records": total_records - mapped_records,
        "spliced_records": spliced_records,
        "top_references": reference_counts.most_common(5),
    }


def create_mapping_plot(metrics: dict[str, str], destination: Path = MAPPING_PLOT_PATH) -> Path | None:
    categories = [
        ("Uniquely mapped", parse_percent(metrics.get("Uniquely mapped reads %"))),
        ("Multimapped", parse_percent(metrics.get("% of reads mapped to multiple loci"))),
        ("Too many loci", parse_percent(metrics.get("% of reads mapped to too many loci"))),
        ("Unmapped", parse_percent(metrics.get("% of reads unmapped: too many mismatches")))
    ]
    unmapped_other = [
        parse_percent(metrics.get("% of reads unmapped: too short")) or 0.0,
        parse_percent(metrics.get("% of reads unmapped: other")) or 0.0,
    ]
    if categories[-1][1] is not None:
        categories[-1] = ("Unmapped mismatches", categories[-1][1])
    categories.extend(
        [
            ("Unmapped too short", unmapped_other[0]),
            ("Unmapped other", unmapped_other[1]),
        ]
    )

    labels = [label for label, value in categories if value is not None]
    values = [value for _, value in categories if value is not None]
    if not values:
        return None

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(8, 4.8))
    bars = axis.bar(labels, values, color=["#2f7d4a", "#7aa34b", "#b7c94d", "#d88c4a", "#c85c43", "#9c3d4a"][: len(values)])
    axis.set_ylim(0, max(100, max(values) + 5))
    axis.set_ylabel("Percent of reads")
    axis.set_title("STAR Mapping Summary")
    axis.tick_params(axis="x", rotation=25)
    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value + 0.8, f"{value:.2f}%", ha="center", fontsize=9)
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def create_splice_plot(sj_metrics: dict[str, object], destination: Path = SPLICE_PLOT_PATH) -> Path | None:
    motif_counts = sj_metrics.get("motif_counts")
    if not isinstance(motif_counts, Counter) or not motif_counts:
        return None

    labels = list(motif_counts.keys())
    values = list(motif_counts.values())

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(8, 4.8))
    bars = axis.bar(labels, values, color="#2b6cb0")
    axis.set_ylabel("Observed splice junctions")
    axis.set_title("STAR Splice Junction Motifs")
    axis.tick_params(axis="x", rotation=25)
    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value + max(values) * 0.02, str(value), ha="center", fontsize=9)
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def explain_scenario() -> str:
    return textwrap.dedent(
        """\
        Scenario
        We reuse the same official chr22 example reads from the HISAT2 mini-project
        so the aligner changes but the data do not. That makes the STAR results
        easier to compare directly with HISAT2.

        Why STAR is a good fit here
        STAR is a fast splice-aware RNA-seq aligner. It also reports splice
        junctions explicitly, which makes it useful for inspecting transcript-like
        alignment structure.

        Decision we want to make
        We want to check how well STAR maps the example reads and what kinds of
        splice junctions it reports before moving on to counting or transcript
        assembly.
        """
    ).strip()


def build_summary_text(
    star_path: Path | None,
    download_status: dict[str, str] | None,
    index_result: subprocess.CompletedProcess[str] | None,
    align_result: subprocess.CompletedProcess[str] | None,
    star_metrics: dict[str, str],
    sj_metrics: dict[str, object],
    sam_metrics: dict[str, object],
    mapping_plot: Path | None,
    splice_plot: Path | None,
    error_message: str | None,
) -> str:
    lines = [
        "Mini-project 2: STAR alignment showcase",
        "",
        explain_scenario(),
        "",
        "Official data sources",
        f"Reference FASTA: {REFERENCE_URL}",
        f"Read 1 FASTA: {READ_1_URL}",
        f"Read 2 FASTA: {READ_2_URL}",
        "",
        "Local workflow inputs",
        f"STAR path: {star_path or 'not found'}",
        f"Reference FASTA: {REFERENCE_FASTA_PATH}",
        f"Read 1 path: {READ_1_PATH}",
        f"Read 2 path: {READ_2_PATH}",
        f"STAR genomeDir: {STAR_INDEX_DIR}",
    ]

    if download_status:
        lines.extend(["", "Download status"])
        for path_text, status in download_status.items():
            lines.append(f"{Path(path_text).name}: {status}")

    lines.extend(["", "Commands used"])
    if star_path:
        lines.append(" ".join(build_star_index_command(star_path, REFERENCE_FASTA_PATH)))
        lines.append(" ".join(build_star_align_command(star_path)))

    if index_result is not None:
        lines.extend(["", "Index build", f"Return code: {index_result.returncode}"])

    if error_message:
        lines.extend(["", "Scenario outcome", error_message])
        return "\n".join(lines).rstrip() + "\n"

    if align_result is None:
        lines.extend(["", "Scenario outcome", "The STAR alignment did not run."])
        return "\n".join(lines).rstrip() + "\n"

    lines.extend(
        [
            "",
            "Scenario outcome",
            f"Return code: {align_result.returncode}",
            f"STAR final log: {STAR_LOG_PATH}",
            f"STAR splice junction table: {STAR_SJ_PATH}",
            f"STAR SAM output: {STAR_SAM_PATH}",
            f"Mapping figure: {mapping_plot or 'not created'}",
            f"Splice junction figure: {splice_plot or 'not created'}",
        ]
    )

    for key in [
        "Number of input reads",
        "Average input read length",
        "Uniquely mapped reads number",
        "Uniquely mapped reads %",
        "Number of splices: Total",
        "Number of splices: GT/AG",
        "Mismatch rate per base, %",
        "% of reads mapped to multiple loci",
        "% of reads mapped to too many loci",
        "% of reads unmapped: too short",
        "% of reads unmapped: other",
    ]:
        value = star_metrics.get(key)
        if value is not None:
            lines.append(f"{key}: {value}")

    lines.extend(
        [
            "",
            "SJ.out.tab summary",
            f"Observed junction rows: {sj_metrics.get('junction_count', 'n/a')}",
            f"Annotated junction rows: {sj_metrics.get('annotated_junction_count', 'n/a')}",
            f"Unannotated junction rows: {sj_metrics.get('unannotated_junction_count', 'n/a')}",
            f"Unique reads supporting junctions: {sj_metrics.get('unique_junction_reads', 'n/a')}",
            f"Multimapping reads supporting junctions: {sj_metrics.get('multimapping_junction_reads', 'n/a')}",
        ]
    )

    motif_counts = sj_metrics.get("motif_counts")
    if isinstance(motif_counts, Counter) and motif_counts:
        lines.extend(["Top splice junction motifs"])
        for motif, count in motif_counts.most_common():
            lines.append(f"- {motif}: {count}")

    lines.extend(
        [
            "",
            "SAM-derived summary",
            f"SAM alignment records: {sam_metrics.get('sam_records', 'n/a')}",
            f"Mapped alignment records: {sam_metrics.get('mapped_records', 'n/a')}",
            f"Unmapped alignment records: {sam_metrics.get('unmapped_records', 'n/a')}",
            f"Spliced alignment records (CIGAR contains N): {sam_metrics.get('spliced_records', 'n/a')}",
        ]
    )

    top_references = sam_metrics.get("top_references")
    if isinstance(top_references, list) and top_references:
        lines.extend(["Top reference targets in the SAM output"])
        for ref_name, count in top_references:
            lines.append(f"- {ref_name}: {count} records")

    unique_percent = parse_percent(star_metrics.get("Uniquely mapped reads %"))
    splice_total = parse_int(star_metrics.get("Number of splices: Total"))
    interpretation_parts = []
    if unique_percent is not None:
        if unique_percent >= 90.0:
            interpretation_parts.append("STAR produced a strong unique mapping rate on this example.")
        elif unique_percent >= 70.0:
            interpretation_parts.append("STAR produced a usable unique mapping rate, though it is worth checking the unmapped categories.")
        else:
            interpretation_parts.append("STAR's unique mapping rate is low enough that the run should be reviewed before downstream analysis.")
    if splice_total is not None:
        interpretation_parts.append(f"It reported {splice_total} splice events in the STAR log, which is useful for transcript-aware follow-up.")

    if interpretation_parts:
        lines.extend(["", "Interpretation", " ".join(interpretation_parts)])

    return "\n".join(lines).rstrip() + "\n"


def save_summary(text: str) -> None:
    SUMMARY_PATH.write_text(text, encoding="utf-8")


def main() -> None:
    ensure_directories()
    star_path = resolve_star_binary()
    download_status: dict[str, str] | None = None
    index_result: subprocess.CompletedProcess[str] | None = None
    align_result: subprocess.CompletedProcess[str] | None = None
    star_metrics: dict[str, str] = {}
    sj_metrics: dict[str, object] = {}
    sam_metrics: dict[str, object] = {}
    mapping_plot: Path | None = None
    splice_plot: Path | None = None
    error_message: str | None = None

    if not star_path:
        error_message = "STAR must be installed before this workflow can run."
    else:
        try:
            download_status = ensure_official_example_data()

            if not star_index_exists(STAR_INDEX_DIR):
                index_result = build_star_index(star_path)
                if index_result.returncode != 0:
                    error_message = "STAR genomeGenerate failed. Check the terminal output for details."

            if error_message is None:
                align_result = run_star_alignment(star_path)
                if align_result.returncode != 0:
                    error_message = "STAR alignment failed. Check the terminal output for details."
                else:
                    star_metrics = parse_star_log(STAR_LOG_PATH)
                    sj_metrics = parse_sj_table(STAR_SJ_PATH)
                    sam_metrics = parse_sam_metrics(STAR_SAM_PATH)
                    mapping_plot = create_mapping_plot(star_metrics)
                    splice_plot = create_splice_plot(sj_metrics)

        except RuntimeError as error:
            error_message = str(error)

    summary_text = build_summary_text(
        star_path=star_path,
        download_status=download_status,
        index_result=index_result,
        align_result=align_result,
        star_metrics=star_metrics,
        sj_metrics=sj_metrics,
        sam_metrics=sam_metrics,
        mapping_plot=mapping_plot,
        splice_plot=splice_plot,
        error_message=error_message,
    )
    save_summary(summary_text)
    print(summary_text)


if __name__ == "__main__":
    main()
    if "failed" in SUMMARY_PATH.read_text(encoding="utf-8").lower():
        raise SystemExit(1)
