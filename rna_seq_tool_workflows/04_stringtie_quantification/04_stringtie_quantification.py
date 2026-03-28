"""Mini-project: StringTie quantification workflow.

Description:
Run StringTie on the official HTSeq yeast example alignment after sorting it by
genomic position, then summarize assembled transcripts and gene abundances.

Official data used:
- https://raw.githubusercontent.com/htseq/htseq/master/example_data/yeast_RNASeq_excerpt.sam
- https://raw.githubusercontent.com/htseq/htseq/master/example_data/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz
- official StringTie macOS binary package from the StringTie GitHub releases

Outputs:
- outputs/yeast_RNASeq_excerpt.sorted.sam
- outputs/stringtie_assembled.gtf
- outputs/stringtie_gene_abundance.tsv
- outputs/stringtie_quantification_summary.txt
- outputs/stringtie_top_genes_tpm.png
- outputs/stringtie_transcript_lengths.png
"""

from __future__ import annotations

from pathlib import Path
from urllib.error import URLError
from urllib.request import urlretrieve
import gzip
import os
import platform
import shutil
import subprocess
import textwrap


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = WORKFLOW_DIR / "data"
READS_DIR = DATA_DIR / "reads"
REFERENCE_DIR = DATA_DIR / "reference"
TOOLS_DIR = WORKFLOW_DIR / "tools"
OUTPUT_DIR = PROJECT_DIR / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
TECHNICAL_DIR = OUTPUT_DIR / "technical"
SAM_PATH = READS_DIR / "yeast_RNASeq_excerpt.sam"
GTF_PATH = REFERENCE_DIR / "Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz"
PLAIN_GTF_PATH = REFERENCE_DIR / "Saccharomyces_cerevisiae.SGD1.01.56.gtf"
SORTED_SAM_PATH = TECHNICAL_DIR / "yeast_RNASeq_excerpt.sorted.sam"
ASSEMBLED_GTF_PATH = TECHNICAL_DIR / "stringtie_assembled.gtf"
GENE_ABUNDANCE_PATH = PRESENTABLE_DIR / "stringtie_gene_abundance.tsv"
SUMMARY_PATH = PRESENTABLE_DIR / "stringtie_quantification_summary.txt"
TOP_GENES_PLOT_PATH = PRESENTABLE_DIR / "stringtie_top_genes_tpm.png"
TRANSCRIPT_LENGTHS_PLOT_PATH = PRESENTABLE_DIR / "stringtie_transcript_lengths.png"
THREADS = 2

SAM_URL = "https://raw.githubusercontent.com/htseq/htseq/master/example_data/yeast_RNASeq_excerpt.sam"
GTF_URL = (
    "https://raw.githubusercontent.com/htseq/htseq/master/example_data/"
    "Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz"
)
STRINGTIE_BINARY = (
    TOOLS_DIR
    / "stringtie_pkg"
    / "stringtie-3.0.3.OSX_x86_64"
    / "stringtie"
)
STRINGTIE_RELEASE_URL = (
    "https://github.com/gpertea/stringtie/releases/download/v3.0.3/"
    "stringtie-3.0.3.OSX_x86_64.tar.gz"
)


def ensure_directories() -> None:
    READS_DIR.mkdir(parents=True, exist_ok=True)
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    PRESENTABLE_DIR.mkdir(parents=True, exist_ok=True)
    TECHNICAL_DIR.mkdir(parents=True, exist_ok=True)


def configure_matplotlib() -> None:
    mpl_config_dir = TECHNICAL_DIR / ".mplconfig"
    cache_dir = TECHNICAL_DIR / ".cache"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))


def download_if_missing(url: str, destination: Path) -> str:
    if destination.exists():
        return "reused"
    try:
        urlretrieve(url, destination)
    except URLError as error:
        raise RuntimeError(f"Could not download {url}: {error}") from error
    return "downloaded"


def ensure_input_data() -> dict[str, str]:
    return {
        str(SAM_PATH): download_if_missing(SAM_URL, SAM_PATH),
        str(GTF_PATH): download_if_missing(GTF_URL, GTF_PATH),
    }


def ensure_plain_gtf(source_gzip: Path, destination_plain: Path) -> Path:
    if destination_plain.exists():
        return destination_plain
    with gzip.open(source_gzip, "rb") as source_handle, destination_plain.open("wb") as destination_handle:
        shutil.copyfileobj(source_handle, destination_handle)
    return destination_plain


def resolve_stringtie_binary() -> Path | None:
    path_binary = shutil.which("stringtie")
    if path_binary:
        return Path(path_binary)
    if STRINGTIE_BINARY.exists():
        return STRINGTIE_BINARY
    return None


def stringtie_command_prefix(binary: Path) -> list[str]:
    if (
        platform.system() == "Darwin"
        and platform.machine() == "arm64"
        and binary.resolve() == STRINGTIE_BINARY.resolve()
    ):
        return ["arch", "-x86_64", str(binary)]
    return [str(binary)]


def parse_sam_for_sorting(path: Path) -> tuple[list[str], list[tuple[tuple[int, int, str], str]]]:
    header_lines: list[str] = []
    reference_order: dict[str, int] = {}
    alignment_rows: list[tuple[tuple[int, int, str], str]] = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@"):
                header_lines.append(line)
                if line.startswith("@SQ"):
                    fields = line.rstrip("\n").split("\t")
                    reference_name = ""
                    for field in fields:
                        if field.startswith("SN:"):
                            reference_name = field[3:]
                            break
                    if reference_name:
                        reference_order[reference_name] = len(reference_order)
                continue

            fields = line.rstrip("\n").split("\t")
            reference_name = fields[2]
            position = int(fields[3])
            qname = fields[0]
            ref_rank = reference_order.get(reference_name, len(reference_order) + 1)
            alignment_rows.append(((ref_rank, position, qname), line))

    return header_lines, alignment_rows


def sort_sam_by_coordinate(source: Path, destination: Path) -> Path:
    header_lines, alignment_rows = parse_sam_for_sorting(source)
    alignment_rows.sort(key=lambda item: item[0])

    header_lines = [line for line in header_lines if not line.startswith("@HD")]
    with destination.open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.0\tSO:coordinate\n")
        for line in header_lines:
            handle.write(line)
        for _, line in alignment_rows:
            handle.write(line)

    return destination


def build_stringtie_command(binary: Path, sorted_sam: Path, gtf_path: Path) -> list[str]:
    return stringtie_command_prefix(binary) + [
        "-p",
        str(THREADS),
        "-G",
        str(gtf_path),
        "-A",
        str(GENE_ABUNDANCE_PATH),
        "-o",
        str(ASSEMBLED_GTF_PATH),
        str(sorted_sam),
    ]


def run_stringtie(binary: Path, sorted_sam: Path, gtf_path: Path) -> subprocess.CompletedProcess[str]:
    command = build_stringtie_command(binary, sorted_sam, gtf_path)
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    return subprocess.run(command, check=False, capture_output=True, text=True, env=env)


def parse_gene_abundance(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if len(values) != len(header):
                continue
            row = dict(zip(header, values))
            gene_name = row.get("Gene Name", "") or row.get("Gene ID", "")
            if gene_name == "-":
                continue
            rows.append(
                {
                    "gene_id": row.get("Gene ID", ""),
                    "gene_name": gene_name,
                    "coverage": float(row.get("Coverage", "0") or 0),
                    "fpkm": float(row.get("FPKM", "0") or 0),
                    "tpm": float(row.get("TPM", "0") or 0),
                }
            )
    rows.sort(key=lambda row: float(row["tpm"]), reverse=True)
    return rows


def parse_assembled_gtf(path: Path) -> list[dict[str, object]]:
    transcripts: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            start = int(fields[3])
            end = int(fields[4])
            attr_text = fields[8]
            transcript_id = ""
            reference_gene_id = ""
            for chunk in attr_text.split(";"):
                chunk = chunk.strip()
                if chunk.startswith("transcript_id"):
                    transcript_id = chunk.split('"')[1]
                elif chunk.startswith("ref_gene_id"):
                    reference_gene_id = chunk.split('"')[1]
            transcripts.append(
                {
                    "transcript_id": transcript_id,
                    "reference_gene_id": reference_gene_id,
                    "length": end - start + 1,
                }
            )
    return transcripts


def summarize_inputs() -> dict[str, object]:
    sam_records = 0
    references: set[str] = set()
    with SAM_PATH.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@SQ"):
                for field in line.rstrip("\n").split("\t"):
                    if field.startswith("SN:"):
                        references.add(field[3:])
            elif not line.startswith("@"):
                sam_records += 1

    gtf_lines = 0
    genes: set[str] = set()
    with gzip.open(GTF_PATH, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            gtf_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            for chunk in fields[8].split(";"):
                chunk = chunk.strip()
                if chunk.startswith("gene_id"):
                    genes.add(chunk.split('"')[1])
                    break

    return {
        "sam_records": sam_records,
        "reference_count": len(references),
        "gtf_lines": gtf_lines,
        "gene_count": len(genes),
    }


def create_top_genes_plot(gene_rows: list[dict[str, object]], destination: Path = TOP_GENES_PLOT_PATH) -> Path | None:
    top_rows = [row for row in gene_rows if float(row["tpm"]) > 0][:12]
    if not top_rows:
        return None

    labels = [str(row["gene_name"]) for row in top_rows]
    values = [float(row["tpm"]) for row in top_rows]

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(9, 5))
    bars = axis.bar(labels, values, color="#2f7d4a")
    axis.set_ylabel("TPM")
    axis.set_title("Top StringTie Gene TPM Values")
    axis.tick_params(axis="x", rotation=30)

    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value + max(values) * 0.01, f"{value:.1f}", ha="center", fontsize=9)

    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def create_transcript_length_plot(transcripts: list[dict[str, object]], destination: Path = TRANSCRIPT_LENGTHS_PLOT_PATH) -> Path | None:
    if not transcripts:
        return None

    lengths = [int(item["length"]) for item in transcripts]

    try:
        configure_matplotlib()
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return None

    figure, axis = plt.subplots(figsize=(8.6, 4.8))
    bins = min(20, max(5, len(set(lengths))))
    axis.hist(lengths, bins=bins, color="#2b6cb0", edgecolor="white")
    axis.set_xlabel("Transcript length (bp)")
    axis.set_ylabel("Transcript count")
    axis.set_title("Assembled Transcript Lengths")
    figure.tight_layout()
    figure.savefig(destination, dpi=200)
    plt.close(figure)
    return destination


def explain_scenario() -> str:
    return textwrap.dedent(
        """\
        Scenario
        We reuse the official HTSeq yeast example alignment and annotation, but
        now we ask a different question: instead of just counting overlaps per
        gene, can we assemble transcripts and estimate expression abundance?

        Why StringTie is useful
        StringTie takes coordinate-sorted spliced alignments and produces an
        assembled transcript GTF plus abundance estimates such as FPKM and TPM.
        This makes it a bridge between raw alignments and transcript-level
        interpretation.

        Decision we want to make
        We want to see which genes have the highest StringTie abundance estimates
        in this example and what the assembled transcript set looks like.
        """
    ).strip()


def build_summary_text(
    stringtie_binary: Path | None,
    download_status: dict[str, str] | None,
    input_summary: dict[str, object],
    command: list[str] | None,
    run_result: subprocess.CompletedProcess[str] | None,
    gene_rows: list[dict[str, object]],
    transcripts: list[dict[str, object]],
    top_plot: Path | None,
    length_plot: Path | None,
    error_message: str | None,
) -> str:
    lines = [
        "Mini-project 4: StringTie quantification showcase",
        "",
        explain_scenario(),
        "",
        "Official data sources",
        f"SAM: {SAM_URL}",
        f"GTF: {GTF_URL}",
        f"StringTie release: {STRINGTIE_RELEASE_URL}",
        "",
        "Workflow inputs",
        f"StringTie binary: {stringtie_binary or 'not found'}",
        f"SAM input: {SAM_PATH}",
        f"GTF guide (gz): {GTF_PATH}",
        f"GTF guide used by StringTie: {PLAIN_GTF_PATH}",
        f"Alignment records in SAM: {input_summary.get('sam_records', 'n/a')}",
        f"Reference sequences in SAM header: {input_summary.get('reference_count', 'n/a')}",
        f"GTF feature lines: {input_summary.get('gtf_lines', 'n/a')}",
        f"Distinct genes in guide annotation: {input_summary.get('gene_count', 'n/a')}",
        f"Sorted SAM used for StringTie: {SORTED_SAM_PATH}",
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
        lines.extend(["", "Scenario outcome", "StringTie did not run."])
        return "\n".join(lines).rstrip() + "\n"

    lines.extend(
        [
            "",
            "Scenario outcome",
            f"Return code: {run_result.returncode}",
            f"Assembled GTF: {ASSEMBLED_GTF_PATH}",
            f"Gene abundance table: {GENE_ABUNDANCE_PATH}",
            f"Top genes TPM figure: {top_plot or 'not created'}",
            f"Transcript length figure: {length_plot or 'not created'}",
            f"Assembled transcript records: {len(transcripts)}",
            f"Genes with nonzero TPM: {sum(1 for row in gene_rows if float(row['tpm']) > 0)}",
        ]
    )

    lines.extend(["", "Top StringTie genes by TPM"])
    for row in gene_rows[:12]:
        if float(row["tpm"]) <= 0:
            continue
        lines.append(
            f"- {row['gene_name']} ({row['gene_id']}): TPM={float(row['tpm']):.2f}, "
            f"FPKM={float(row['fpkm']):.2f}, coverage={float(row['coverage']):.2f}"
        )

    if transcripts:
        transcript_lengths = [int(item["length"]) for item in transcripts]
        lines.extend(
            [
                "",
                "Transcript assembly summary",
                f"Minimum transcript length: {min(transcript_lengths)}",
                f"Median-ish length example: {transcript_lengths[len(transcript_lengths) // 2]}",
                f"Maximum transcript length: {max(transcript_lengths)}",
            ]
        )

    lines.extend(
        [
            "",
            "Interpretation",
            "StringTie complements HTSeq by estimating transcript structures and abundances, so the output now includes TPM and FPKM-style expression summaries instead of just overlap counts.",
        ]
    )

    stderr_text = run_result.stderr.strip()
    if stderr_text:
        lines.extend(["", "StringTie stderr", stderr_text])

    return "\n".join(lines).rstrip() + "\n"


def save_summary(text: str) -> None:
    SUMMARY_PATH.write_text(text, encoding="utf-8")


def main() -> None:
    ensure_directories()
    stringtie_binary = resolve_stringtie_binary()
    download_status: dict[str, str] | None = None
    input_summary: dict[str, object] = {}
    command: list[str] | None = None
    run_result: subprocess.CompletedProcess[str] | None = None
    gene_rows: list[dict[str, object]] = []
    transcripts: list[dict[str, object]] = []
    top_plot: Path | None = None
    length_plot: Path | None = None
    error_message: str | None = None

    if not stringtie_binary:
        error_message = (
            "StringTie binary is missing. Install it or place the official "
            "binary package under rna_seq_tool_workflows/tools/stringtie_pkg."
        )
    else:
        try:
            download_status = ensure_input_data()
            input_summary = summarize_inputs()
            ensure_plain_gtf(GTF_PATH, PLAIN_GTF_PATH)
            sort_sam_by_coordinate(SAM_PATH, SORTED_SAM_PATH)
            command = build_stringtie_command(stringtie_binary, SORTED_SAM_PATH, PLAIN_GTF_PATH)
            run_result = run_stringtie(stringtie_binary, SORTED_SAM_PATH, PLAIN_GTF_PATH)
            if run_result.returncode != 0:
                error_message = "StringTie quantification failed. Check the terminal output for details."
            else:
                gene_rows = parse_gene_abundance(GENE_ABUNDANCE_PATH)
                transcripts = parse_assembled_gtf(ASSEMBLED_GTF_PATH)
                top_plot = create_top_genes_plot(gene_rows)
                length_plot = create_transcript_length_plot(transcripts)
        except RuntimeError as error:
            error_message = str(error)

    summary_text = build_summary_text(
        stringtie_binary=stringtie_binary,
        download_status=download_status,
        input_summary=input_summary,
        command=command,
        run_result=run_result,
        gene_rows=gene_rows,
        transcripts=transcripts,
        top_plot=top_plot,
        length_plot=length_plot,
        error_message=error_message,
    )
    save_summary(summary_text)
    print(summary_text)


if __name__ == "__main__":
    main()
    if "failed" in SUMMARY_PATH.read_text(encoding="utf-8").lower():
        raise SystemExit(1)
