"""Mini-project: FASTA parser and sequence statistics.

Description:
Download a real FASTA file from NCBI, parse its sequence records, compute
basic statistics, and generate a Matplotlib visualization for the results.

Objectives:
- Download a real FASTA record from NCBI
- Parse FASTA headers and multi-line sequence records
- Compute length, GC count, GC content, and base composition
- Find high-occurring k-mers across increasing k values
- Save a summary report and a Matplotlib plot

Files we will use:
- Input: data/NC_045512.2.fasta
- Output: outputs/fasta_summary.txt
- Output: outputs/fasta_base_composition.png

Used public data:
- NCBI Nucleotide accession NC_045512.2
- Organism: Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)
"""

from __future__ import annotations

from pathlib import Path
import sys
from urllib.error import URLError
from urllib.request import urlopen

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dna_formats import FastaRecord, gc_content, gc_count, nucleotide_counts, parse_fasta


SEQUENCE_FUNDAMENTALS_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = SEQUENCE_FUNDAMENTALS_DIR / "data"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
FASTA_PATH = DATA_DIR / "NC_045512.2.fasta"
SUMMARY_PATH = PRESENTABLE_DIR / "fasta_summary.txt"
VISUALIZATION_PATH = PRESENTABLE_DIR / "fasta_base_composition.png"
NCBI_FASTA_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    "efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text"
)


def ensure_directories() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    PRESENTABLE_DIR.mkdir(parents=True, exist_ok=True)


def download_reference_fasta(destination: Path = FASTA_PATH) -> tuple[Path, bool]:
    ensure_directories()

    if destination.exists():
        return destination, False

    with urlopen(NCBI_FASTA_URL, timeout=30) as response:
        fasta_text = response.read().decode("utf-8")

    destination.write_text(fasta_text, encoding="utf-8")
    return destination, True


def load_fasta_records(path: Path = FASTA_PATH) -> list[FastaRecord]:
    fasta_text = path.read_text(encoding="utf-8")
    records = parse_fasta(fasta_text)
    if not records:
        raise ValueError(f"No FASTA records were found in {path}.")
    return records


def compute_record_summary(record: FastaRecord) -> dict[str, object]:
    sequence = record.sequence.upper()
    counts = nucleotide_counts(sequence)
    sequence_length = len(sequence)
    longest_runs = longest_homopolymer_runs(sequence)
    repeated_kmers = find_repeated_kmers_by_length(sequence)

    return {
        "header": record.header,
        "length": sequence_length,
        "gc_count": gc_count(sequence),
        "gc_content": gc_content(sequence),
        "counts": counts,
        "longest_runs": longest_runs,
        "repeated_kmers": repeated_kmers,
    }


def longest_homopolymer_runs(sequence: str) -> dict[str, int]:
    best = {base: 0 for base in ["A", "T", "C", "G", "N"]}
    current_base = ""
    current_length = 0

    for base in sequence:
        if base == current_base:
            current_length += 1
        else:
            current_base = base
            current_length = 1

        if base not in best:
            best["N"] = max(best["N"], current_length)
        else:
            best[base] = max(best[base], current_length)

    return best


def kmer_counts(sequence: str, k: int) -> dict[str, int]:
    counts: dict[str, int] = {}
    if k <= 0 or len(sequence) < k:
        return counts

    for start in range(len(sequence) - k + 1):
        kmer = sequence[start : start + k]
        if "N" in kmer:
            continue
        counts[kmer] = counts.get(kmer, 0) + 1

    return counts


def find_repeated_kmers_by_length(
    sequence: str,
    start_k: int = 4,
    max_results_per_k: int = 4,
    min_kmer_count: int = 5,
    max_k: int | None = None,
) -> list[dict[str, object]]:
    repeated_groups: list[dict[str, object]] = []

    if max_k is None:
        max_k = max(start_k, min(20, len(sequence)))

    for k in range(start_k, max_k + 1):
        counts = kmer_counts(sequence, k)
        repeated = [
            (kmer, count)
            for kmer, count in counts.items()
            if count >= min_kmer_count
        ]

        if not repeated:
            break

        repeated.sort(key=lambda item: (-item[1], item[0]))
        top_repeated = repeated[:max_results_per_k]
        repeated_groups.append(
            {
                "k": k,
                "results": [
                    {
                        "kmer": kmer,
                        "count": count,
                        "first_position": sequence.find(kmer) + 1,
                    }
                    for kmer, count in top_repeated
                ],
            }
        )

    return repeated_groups


def build_summary_text(records: list[FastaRecord]) -> str:
    lines = [
        "Mini-project 1: FASTA parser and sequence statistics",
        "Dataset: NCBI Nucleotide accession NC_045512.2",
        "Organism: Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)",
        "",
    ]

    for index, record in enumerate(records, start=1):
        summary = compute_record_summary(record)
        counts = summary["counts"]
        longest_runs = summary["longest_runs"]
        repeated_kmers = summary["repeated_kmers"]

        lines.extend(
            [
                f"Record {index}",
                f"Header: {summary['header']}",
                f"Length: {summary['length']} bases",
                f"GC count: {summary['gc_count']}",
                f"GC content: {summary['gc_content']:.2f}%",
                (
                    "Base counts: "
                    f"A={counts['A']}, T={counts['T']}, C={counts['C']}, "
                    f"G={counts['G']}, N={counts['N']}"
                ),
                (
                    "Longest homopolymer runs: "
                    f"A={longest_runs['A']}, T={longest_runs['T']}, "
                    f"C={longest_runs['C']}, G={longest_runs['G']}, "
                    f"N={longest_runs['N']}"
                ),
                "",
            ]
        )

        if repeated_kmers:
            lines.append("High-occurring k-mers by k:")
            for group in repeated_kmers:
                lines.append(f"  k={group['k']}")
                for result in group["results"]:
                    lines.append(
                        (
                            f"    {result['kmer']}: count={result['count']}, "
                            f"first_position={result['first_position']}"
                        )
                    )
            lines.append("")
        else:
            lines.extend(
                [
                    "High-occurring k-mers by k:",
                    "  none found with the current threshold",
                    "",
                ]
            )

    return "\n".join(lines).rstrip() + "\n"


def save_summary(records: list[FastaRecord], destination: Path = SUMMARY_PATH) -> Path:
    summary_text = build_summary_text(records)
    destination.write_text(summary_text, encoding="utf-8")
    return destination


def create_base_composition_plot(record: FastaRecord, destination: Path = VISUALIZATION_PATH) -> Path:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Matplotlib is required for the visualization step. "
            "Install it with: python3 -m pip install matplotlib"
        ) from error

    counts = nucleotide_counts(record.sequence)
    bases = ["A", "T", "C", "G", "N"]
    values = [counts[base] for base in bases]
    colors = {
        "A": "#4CAF50",
        "T": "#FF7043",
        "C": "#42A5F5",
        "G": "#AB47BC",
        "N": "#90A4AE",
    }

    gc_value = gc_content(record.sequence)
    at_value = 100 - gc_value
    figure, axis = plt.subplots(figsize=(9, 5.5))
    bars = axis.bar(bases, values, color=[colors[base] for base in bases], edgecolor="#1F2937")

    axis.set_title("SARS-CoV-2 Reference Genome Base Composition", fontsize=14)
    axis.set_xlabel("Nucleotide")
    axis.set_ylabel("Count")
    axis.grid(axis="y", linestyle="--", alpha=0.35)
    axis.set_axisbelow(True)

    for bar, value in zip(bars, values):
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            f"{value}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    figure.text(0.13, 0.92, f"Length: {len(record.sequence)} bases", fontsize=10)
    figure.text(0.13, 0.88, f"GC content: {gc_value:.2f}% | AT content: {at_value:.2f}%", fontsize=10)
    figure.tight_layout(rect=(0, 0, 1, 0.84))
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def run() -> dict[str, Path]:
    fasta_path, _ = download_reference_fasta()
    records = load_fasta_records(fasta_path)
    summary_path = save_summary(records)
    visualization_path = create_base_composition_plot(records[0])
    return {
        "fasta_path": fasta_path,
        "summary_path": summary_path,
        "visualization_path": visualization_path,
    }


def main() -> None:
    try:
        fasta_path, was_downloaded = download_reference_fasta()
        records = load_fasta_records(fasta_path)
        summary_path = save_summary(records)
        visualization_path = create_base_composition_plot(records[0])
    except URLError as error:
        print("Download failed. The script is ready, but it could not reach NCBI right now.")
        print(f"Reason: {error}")
        print(f"Expected FASTA URL: {NCBI_FASTA_URL}")
        return
    except ModuleNotFoundError as error:
        print(error)
        return

    print("Mini-project 1 completed.")
    if was_downloaded:
        print(f"Downloaded FASTA from NCBI: {fasta_path}")
    else:
        print(f"Using existing FASTA file: {fasta_path}")
    print(f"Summary report: {summary_path}")
    print(f"Visualization: {visualization_path}")


if __name__ == "__main__":
    main()
