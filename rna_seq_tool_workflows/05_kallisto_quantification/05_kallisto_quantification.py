"""Mini-project: kallisto quantification workflow.

Description:
Use the official kallisto bundled test data to build a transcriptome index and
run pseudoalignment-based quantification on paired-end reads.

Official data used:
- bundled kallisto test transcriptome: `tools/kallisto/test/transcripts.fasta.gz`
- bundled kallisto test GTF: `tools/kallisto/test/transcripts.gtf.gz`
- bundled kallisto test paired FASTQ reads
- official kallisto macOS binary package

Outputs:
- outputs/transcripts.idx
- outputs/abundance.tsv
- outputs/run_info.json
- outputs/kallisto_quantification_summary.txt
- outputs/kallisto_top_transcripts_tpm.png
- outputs/kallisto_top_genes_tpm.png
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import gzip
import json
import os
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
TEST_DIR = TOOLS_DIR / "kallisto" / "test"
TRANSCRIPTS_FASTA_PATH = TEST_DIR / "transcripts.fasta.gz"
TRANSCRIPTS_GTF_PATH = TEST_DIR / "transcripts.gtf.gz"
READ_1_PATH = TEST_DIR / "reads_1.fastq.gz"
READ_2_PATH = TEST_DIR / "reads_2.fastq.gz"
INDEX_PATH = TECHNICAL_DIR / "transcripts.idx"
TECHNICAL_ABUNDANCE_PATH = TECHNICAL_DIR / "abundance.tsv"
ABUNDANCE_PATH = PRESENTABLE_DIR / "abundance.tsv"
RUN_INFO_PATH = TECHNICAL_DIR / "run_info.json"
SUMMARY_PATH = PRESENTABLE_DIR / "kallisto_quantification_summary.txt"
TOP_TRANSCRIPTS_PLOT_PATH = PRESENTABLE_DIR / "kallisto_top_transcripts_tpm.png"
TOP_GENES_PLOT_PATH = PRESENTABLE_DIR / "kallisto_top_genes_tpm.png"
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


def resolve_kallisto_binary() -> Path | None:
    path_binary = shutil.which("kallisto")
    if path_binary:
        return Path(path_binary)
    if KALLISTO_BINARY.exists():
        return KALLISTO_BINARY
    return None


def required_inputs_exist() -> bool:
    return all(
        path.exists()
        for path in [TRANSCRIPTS_FASTA_PATH, TRANSCRIPTS_GTF_PATH, READ_1_PATH, READ_2_PATH]
    )


def build_index_command(binary: Path) -> list[str]:
    return [
        str(binary),
        "index",
        "-i",
        str(INDEX_PATH),
        str(TRANSCRIPTS_FASTA_PATH),
    ]


def build_quant_command(binary: Path) -> list[str]:
    return [
        str(binary),
        "quant",
        "-i",
        str(INDEX_PATH),
        "-o",
        str(TECHNICAL_DIR),
        "--plaintext",
        "-t",
        str(THREADS),
        str(READ_1_PATH),
        str(READ_2_PATH),
    ]


def run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    return subprocess.run(command, check=False, capture_output=True, text=True, env=env)


def parse_transcript_labels(gtf_path: Path) -> dict[str, dict[str, str]]:
    labels: dict[str, dict[str, str]] = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs: dict[str, str] = {}
            for chunk in fields[8].split(";"):
                chunk = chunk.strip()
                if not chunk or " " not in chunk:
                    continue
                key, rest = chunk.split(" ", 1)
                attrs[key] = rest.strip().strip('"')
            transcript_id = attrs.get("transcript_id", "")
            transcript_version = attrs.get("transcript_version", "")
            if transcript_id:
                key = transcript_id if not transcript_version else f"{transcript_id}.{transcript_version}"
                labels[key] = {
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": attrs.get("gene_name", attrs.get("gene_id", transcript_id)),
                    "transcript_name": attrs.get("transcript_name", transcript_id),
                }
    return labels


def parse_abundance(path: Path, labels: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if len(values) != len(header):
                continue
            row = dict(zip(header, values))
            transcript_id = row.get("target_id", "")
            label = labels.get(transcript_id, {})
            rows.append(
                {
                    "target_id": transcript_id,
                    "transcript_name": label.get("transcript_name", transcript_id),
                    "gene_id": label.get("gene_id", ""),
                    "gene_name": label.get("gene_name", transcript_id),
                    "length": float(row.get("length", "0") or 0),
                    "eff_length": float(row.get("eff_length", "0") or 0),
                    "est_counts": float(row.get("est_counts", "0") or 0),
                    "tpm": float(row.get("tpm", "0") or 0),
                }
            )
    rows.sort(key=lambda row: float(row["tpm"]), reverse=True)
    return rows


def summarize_gene_tpm(transcript_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    gene_totals: dict[tuple[str, str], float] = defaultdict(float)
    for row in transcript_rows:
        gene_key = (str(row["gene_id"]), str(row["gene_name"]))
        gene_totals[gene_key] += float(row["tpm"])
    gene_rows = [
        {"gene_id": gene_id, "gene_name": gene_name, "tpm": tpm}
        for (gene_id, gene_name), tpm in gene_totals.items()
    ]
    gene_rows.sort(key=lambda row: float(row["tpm"]), reverse=True)
    return gene_rows


def parse_run_info(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def copy_presentable_outputs() -> None:
    shutil.copyfile(TECHNICAL_ABUNDANCE_PATH, ABUNDANCE_PATH)


def create_top_transcripts_plot(rows: list[dict[str, object]], destination: Path = TOP_TRANSCRIPTS_PLOT_PATH) -> Path | None:
    top_rows = [row for row in rows if float(row["tpm"]) > 0][:12]
    if not top_rows:
        return None

    labels = [str(row["transcript_name"]) for row in top_rows]
    values = [float(row["tpm"]) for row in top_rows]

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(10, 5.2))
    bars = axis.bar(labels, values, color="#2f7d4a")
    axis.set_ylabel("TPM")
    axis.set_title("Top kallisto Transcripts by TPM")
    axis.tick_params(axis="x", rotation=30)
    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value + max(values) * 0.01, f"{value:.1f}", ha="center", fontsize=9)
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def create_top_genes_plot(rows: list[dict[str, object]], destination: Path = TOP_GENES_PLOT_PATH) -> Path | None:
    top_rows = [row for row in rows if float(row["tpm"]) > 0][:12]
    if not top_rows:
        return None

    labels = [str(row["gene_name"]) for row in top_rows]
    values = [float(row["tpm"]) for row in top_rows]

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(10, 5.2))
    bars = axis.bar(labels, values, color="#2b6cb0")
    axis.set_ylabel("Summed TPM")
    axis.set_title("Top Genes from kallisto Transcript TPM")
    axis.tick_params(axis="x", rotation=30)
    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value + max(values) * 0.01, f"{value:.1f}", ha="center", fontsize=9)
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def explain_scenario() -> str:
    return textwrap.dedent(
        """\
        Scenario
        We use the official kallisto test dataset: a small transcriptome, a
        matching GTF, and paired-end reads bundled with the kallisto source
        distribution.

        Why kallisto is different
        Unlike STAR or HISAT2, kallisto does not perform full genome alignment.
        Instead, it pseudoaligns reads directly to a transcriptome index and
        estimates transcript abundance efficiently.

        Decision we want to make
        We want to inspect which transcripts and genes have the highest inferred
        abundance in the test sample and compare this transcript-centric output
        with the alignment-based workflows from the previous mini-projects.
        """
    ).strip()


def build_summary_text(
    kallisto_binary: Path | None,
    index_result: subprocess.CompletedProcess[str] | None,
    quant_result: subprocess.CompletedProcess[str] | None,
    transcript_rows: list[dict[str, object]],
    gene_rows: list[dict[str, object]],
    run_info: dict[str, object],
    top_transcript_plot: Path | None,
    top_gene_plot: Path | None,
    error_message: str | None,
) -> str:
    lines = [
        "Mini-project 5: kallisto quantification showcase",
        "",
        explain_scenario(),
        "",
        "Workflow inputs",
        f"kallisto binary: {kallisto_binary or 'not found'}",
        f"Transcript FASTA: {TRANSCRIPTS_FASTA_PATH}",
        f"Transcript GTF: {TRANSCRIPTS_GTF_PATH}",
        f"Read 1 FASTQ: {READ_1_PATH}",
        f"Read 2 FASTQ: {READ_2_PATH}",
        f"Index path: {INDEX_PATH}",
    ]

    if index_result is not None:
        lines.extend(["", "Index build", f"Return code: {index_result.returncode}"])
    if quant_result is not None:
        lines.extend(["", "Quantification", f"Return code: {quant_result.returncode}"])

    if error_message:
        lines.extend(["", "Scenario outcome", error_message])
        return "\n".join(lines).rstrip() + "\n"

    if quant_result is None:
        lines.extend(["", "Scenario outcome", "kallisto did not run."])
        return "\n".join(lines).rstrip() + "\n"

    lines.extend(
        [
            "",
            "Scenario outcome",
            f"Abundance table: {ABUNDANCE_PATH}",
            f"Run info JSON: {RUN_INFO_PATH}",
            f"Top transcript figure: {top_transcript_plot or 'not created'}",
            f"Top gene figure: {top_gene_plot or 'not created'}",
            f"Transcript rows quantified: {len(transcript_rows)}",
            f"Gene summaries with nonzero TPM: {sum(1 for row in gene_rows if float(row['tpm']) > 0)}",
        ]
    )

    if run_info:
        lines.extend(
            [
                "",
                "Run info",
                f"Reads processed: {run_info.get('n_processed', 'n/a')}",
                f"Pseudoaligned reads: {run_info.get('n_pseudoaligned', 'n/a')}",
                f"Pseudoalignment rate: {run_info.get('p_pseudoaligned', 'n/a')}",
                f"Kallisto version: {run_info.get('kallisto_version', 'n/a')}",
            ]
        )

    lines.extend(["", "Top transcripts by TPM"])
    for row in transcript_rows[:12]:
        if float(row["tpm"]) <= 0:
            continue
        lines.append(
            f"- {row['transcript_name']} ({row['target_id']}): TPM={float(row['tpm']):.2f}, "
            f"est_counts={float(row['est_counts']):.2f}"
        )

    lines.extend(["", "Top genes by summed transcript TPM"])
    for row in gene_rows[:12]:
        if float(row["tpm"]) <= 0:
            continue
        lines.append(f"- {row['gene_name']} ({row['gene_id']}): TPM={float(row['tpm']):.2f}")

    lines.extend(
        [
            "",
            "Interpretation",
            "kallisto gives transcript-level abundance estimates without full alignment, so the main outputs are TPM and estimated counts per transcript rather than SAM/BAM alignments.",
        ]
    )

    if index_result and index_result.stderr.strip():
        lines.extend(["", "Index stderr", index_result.stderr.strip()])
    if quant_result and quant_result.stderr.strip():
        lines.extend(["", "Quant stderr", quant_result.stderr.strip()])

    return "\n".join(lines).rstrip() + "\n"


def save_summary(text: str) -> None:
    SUMMARY_PATH.write_text(text, encoding="utf-8")


def main() -> None:
    ensure_directories()
    kallisto_binary = resolve_kallisto_binary()
    index_result: subprocess.CompletedProcess[str] | None = None
    quant_result: subprocess.CompletedProcess[str] | None = None
    transcript_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    run_info: dict[str, object] = {}
    top_transcript_plot: Path | None = None
    top_gene_plot: Path | None = None
    error_message: str | None = None

    if not kallisto_binary:
        error_message = (
            "kallisto binary is missing. Place the official macOS package "
            "under rna_seq_tool_workflows/tools."
        )
    elif not required_inputs_exist():
        error_message = "The official kallisto test files are missing from the extracted source tree."
    else:
        index_result = run_command(build_index_command(kallisto_binary))
        if index_result.returncode != 0:
            error_message = "kallisto index failed. Check the terminal output for details."
        else:
            quant_result = run_command(build_quant_command(kallisto_binary))
            if quant_result.returncode != 0:
                error_message = "kallisto quant failed. Check the terminal output for details."
            else:
                labels = parse_transcript_labels(TRANSCRIPTS_GTF_PATH)
                copy_presentable_outputs()
                transcript_rows = parse_abundance(ABUNDANCE_PATH, labels)
                gene_rows = summarize_gene_tpm(transcript_rows)
                run_info = parse_run_info(RUN_INFO_PATH)
                top_transcript_plot = create_top_transcripts_plot(transcript_rows)
                top_gene_plot = create_top_genes_plot(gene_rows)

    summary_text = build_summary_text(
        kallisto_binary=kallisto_binary,
        index_result=index_result,
        quant_result=quant_result,
        transcript_rows=transcript_rows,
        gene_rows=gene_rows,
        run_info=run_info,
        top_transcript_plot=top_transcript_plot,
        top_gene_plot=top_gene_plot,
        error_message=error_message,
    )
    save_summary(summary_text)
    print(summary_text)


if __name__ == "__main__":
    main()
