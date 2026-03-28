"""Mini-project: FASTQ parser and quality score summary.

Description:
Download a real paired-end FASTQ dataset, parse its reads directly from gzip
files, compute quality-related summaries, and generate a Matplotlib quality
plot.

Objectives:
- Download real FASTQ files for a public run accession
- Parse FASTQ records in blocks of four lines
- Validate sequence and quality lengths
- Convert ASCII quality symbols to Phred scores
- Compute per-read and dataset-level quality summaries
- Visualize average quality by base position

Files we will use:
- Input: data/SRR20172987_1.fastq.gz
- Input: data/SRR20172987_2.fastq.gz
- Output: outputs/fastq_quality_summary.txt
- Output: outputs/fastq_position_quality.png

Used public data:
- ENA / NCBI SRA run accession SRR20172987
- Dataset: paired-end 16S amplicon sequencing reads
"""

from __future__ import annotations

import gzip
from pathlib import Path
import sys
from urllib.error import URLError
from urllib.request import urlopen

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dna_formats import FastqRecord, nucleotide_counts


SEQUENCE_FUNDAMENTALS_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = SEQUENCE_FUNDAMENTALS_DIR / "data"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
RUN_ACCESSION = "SRR20172987"
FASTQ_1_PATH = DATA_DIR / f"{RUN_ACCESSION}_1.fastq.gz"
FASTQ_2_PATH = DATA_DIR / f"{RUN_ACCESSION}_2.fastq.gz"
SUMMARY_PATH = PRESENTABLE_DIR / "fastq_quality_summary.txt"
VISUALIZATION_PATH = PRESENTABLE_DIR / "fastq_position_quality.png"
ENA_REPORT_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    f"?accession={RUN_ACCESSION}"
    "&result=read_run"
    "&fields=run_accession,fastq_ftp"
    "&format=tsv"
)


def ensure_directories() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PRESENTABLE_DIR.mkdir(parents=True, exist_ok=True)


def fetch_fastq_urls() -> list[str]:
    with urlopen(ENA_REPORT_URL, timeout=30) as response:
        tsv = response.read().decode("utf-8").strip().splitlines()

    if len(tsv) < 2:
        raise ValueError(f"No FASTQ download information returned for {RUN_ACCESSION}.")

    headers = tsv[0].split("\t")
    values = tsv[1].split("\t")
    row = dict(zip(headers, values))
    fastq_ftp = row.get("fastq_ftp", "")
    if not fastq_ftp:
        raise ValueError(f"No FASTQ URLs were returned for {RUN_ACCESSION}.")

    urls = [entry.strip() for entry in fastq_ftp.split(";") if entry.strip()]
    if not urls:
        raise ValueError(f"No valid FASTQ URLs were parsed for {RUN_ACCESSION}.")

    normalized_urls: list[str] = []
    for url in urls:
        if url.startswith("http://") or url.startswith("https://"):
            normalized_urls.append(url)
        elif url.startswith("ftp://"):
            normalized_urls.append("https://" + url.removeprefix("ftp://"))
        else:
            normalized_urls.append("https://" + url.lstrip("/"))

    return normalized_urls


def download_file(url: str, destination: Path) -> tuple[Path, bool]:
    ensure_directories()

    if destination.exists():
        return destination, False

    with urlopen(url, timeout=60) as response:
        destination.write_bytes(response.read())

    return destination, True


def download_fastq_files() -> list[tuple[Path, bool]]:
    urls = fetch_fastq_urls()
    destinations = [FASTQ_1_PATH, FASTQ_2_PATH]

    if len(urls) != len(destinations):
        raise ValueError(
            f"Expected {len(destinations)} FASTQ files for paired-end data, got {len(urls)}."
        )

    return [download_file(url, destination) for url, destination in zip(urls, destinations)]


def iterate_fastq_records(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError(f"Incomplete FASTQ record encountered in {path}.")

            header = header.strip()
            sequence = sequence.strip()
            plus = plus.strip()
            quality = quality.strip()

            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header line in {path}: {header}")
            if plus != "+":
                raise ValueError(f"Invalid FASTQ separator line in {path}: {plus}")
            if len(sequence) != len(quality):
                raise ValueError(f"Sequence/quality length mismatch in {path}.")

            yield FastqRecord(
                header=header[1:].strip(),
                sequence=sequence,
                quality=quality,
            )


def summarize_fastq_file(path: Path) -> dict[str, object]:
    read_count = 0
    total_bases = 0
    total_quality = 0
    min_read_length: int | None = None
    max_read_length = 0
    min_avg_quality: float | None = None
    max_avg_quality = 0.0
    base_counts = {base: 0 for base in ["A", "T", "C", "G", "N"]}
    position_quality_sums: list[int] = []
    position_quality_counts: list[int] = []

    for record in iterate_fastq_records(path):
        read_count += 1
        read_length = len(record.sequence)
        scores = record.phred_scores()
        avg_quality = record.average_quality()

        total_bases += read_length
        total_quality += sum(scores)
        min_read_length = read_length if min_read_length is None else min(min_read_length, read_length)
        max_read_length = max(max_read_length, read_length)
        min_avg_quality = avg_quality if min_avg_quality is None else min(min_avg_quality, avg_quality)
        max_avg_quality = max(max_avg_quality, avg_quality)

        for base, count in nucleotide_counts(record.sequence).items():
            base_counts[base] += count

        for index, score in enumerate(scores):
            if index == len(position_quality_sums):
                position_quality_sums.append(0)
                position_quality_counts.append(0)
            position_quality_sums[index] += score
            position_quality_counts[index] += 1

    if read_count == 0:
        raise ValueError(f"No FASTQ reads found in {path}.")

    average_position_quality = [
        quality_sum / count
        for quality_sum, count in zip(position_quality_sums, position_quality_counts)
    ]

    return {
        "path": path,
        "read_count": read_count,
        "total_bases": total_bases,
        "average_read_length": total_bases / read_count,
        "min_read_length": min_read_length or 0,
        "max_read_length": max_read_length,
        "average_quality": total_quality / total_bases if total_bases else 0.0,
        "min_average_read_quality": min_avg_quality or 0.0,
        "max_average_read_quality": max_avg_quality,
        "gc_content": 100 * (base_counts["G"] + base_counts["C"]) / total_bases if total_bases else 0.0,
        "base_counts": base_counts,
        "average_position_quality": average_position_quality,
    }


def build_summary_text(summaries: list[dict[str, object]]) -> str:
    lines = [
        "Mini-project 2: FASTQ parser and quality score summary",
        f"Dataset: ENA / NCBI SRA run accession {RUN_ACCESSION}",
        "",
    ]

    for index, summary in enumerate(summaries, start=1):
        base_counts = summary["base_counts"]
        lines.extend(
            [
                f"FASTQ file {index}",
                f"Path: {summary['path']}",
                f"Read count: {summary['read_count']}",
                f"Total bases: {summary['total_bases']}",
                f"Average read length: {summary['average_read_length']:.2f}",
                f"Read length range: {summary['min_read_length']} - {summary['max_read_length']}",
                f"Average base quality: {summary['average_quality']:.2f}",
                (
                    "Average read quality range: "
                    f"{summary['min_average_read_quality']:.2f} - "
                    f"{summary['max_average_read_quality']:.2f}"
                ),
                f"GC content across reads: {summary['gc_content']:.2f}%",
                (
                    "Base counts: "
                    f"A={base_counts['A']}, T={base_counts['T']}, "
                    f"C={base_counts['C']}, G={base_counts['G']}, N={base_counts['N']}"
                ),
                "",
            ]
        )

    return "\n".join(lines).rstrip() + "\n"


def save_summary(summaries: list[dict[str, object]], destination: Path = SUMMARY_PATH) -> Path:
    destination.write_text(build_summary_text(summaries), encoding="utf-8")
    return destination


def create_quality_plot(summaries: list[dict[str, object]], destination: Path = VISUALIZATION_PATH) -> Path:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Matplotlib is required for the FASTQ visualization step. "
            "Install it with: python3 -m pip install matplotlib"
        ) from error

    figure, axis = plt.subplots(figsize=(10, 5.5))
    colors = ["#1F77B4", "#FF7F0E"]

    for index, summary in enumerate(summaries):
        position_quality = summary["average_position_quality"]
        positions = list(range(1, len(position_quality) + 1))
        axis.plot(
            positions,
            position_quality,
            label=f"Read file {index + 1}",
            color=colors[index % len(colors)],
            linewidth=2,
        )

    axis.set_title(f"{RUN_ACCESSION} Average Quality by Base Position", fontsize=14)
    axis.set_xlabel("Base position")
    axis.set_ylabel("Average Phred quality")
    axis.grid(alpha=0.3, linestyle="--")
    axis.legend()
    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def run() -> dict[str, object]:
    downloads = download_fastq_files()
    summaries = [summarize_fastq_file(FASTQ_1_PATH), summarize_fastq_file(FASTQ_2_PATH)]
    summary_path = save_summary(summaries)
    visualization_path = create_quality_plot(summaries)
    return {
        "downloads": downloads,
        "summary_path": summary_path,
        "visualization_path": visualization_path,
    }


def main() -> None:
    try:
        downloads = download_fastq_files()
        summaries = [summarize_fastq_file(FASTQ_1_PATH), summarize_fastq_file(FASTQ_2_PATH)]
        summary_path = save_summary(summaries)
        visualization_path = create_quality_plot(summaries)
    except URLError as error:
        print("Download failed. The script is ready, but it could not reach ENA right now.")
        print(f"Reason: {error}")
        print(f"Expected ENA API URL: {ENA_REPORT_URL}")
        return
    except ModuleNotFoundError as error:
        print(error)
        return

    print("Mini-project 2 completed.")
    for path, was_downloaded in downloads:
        if was_downloaded:
            print(f"Downloaded FASTQ: {path}")
        else:
            print(f"Using existing FASTQ: {path}")
    print(f"Summary report: {summary_path}")
    print(f"Visualization: {visualization_path}")


if __name__ == "__main__":
    main()
