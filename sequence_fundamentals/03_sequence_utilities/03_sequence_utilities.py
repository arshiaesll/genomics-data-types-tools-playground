"""Mini-project: sequence utilities for DNA strings.

Description:
Build reusable sequence operations and apply them to a real reference genome.

Objectives:
- Implement reverse complement
- Implement GC content calculation
- Search for sequence motifs
- Add open reading frame detection

Files we will use:
- Input: data/NC_045512.2.fasta
- Output: outputs/sequence_utilities_report.txt
- Output: outputs/gc_content_windows.png
- Output: outputs/motif_counts.png
- Output: outputs/orf_overview.png

Suggested public data:
- NCBI Nucleotide accession NC_045512.2
- Good target motifs: start codon and stop codons
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dna_formats import FastaRecord, gc_content, parse_fasta, reverse_complement


SEQUENCE_FUNDAMENTALS_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = SEQUENCE_FUNDAMENTALS_DIR / "data"
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
PRESENTABLE_DIR = OUTPUT_DIR / "presentable"
FASTA_PATH = DATA_DIR / "NC_045512.2.fasta"
REPORT_PATH = PRESENTABLE_DIR / "sequence_utilities_report.txt"
GC_FIGURE_PATH = PRESENTABLE_DIR / "gc_content_windows.png"
MOTIF_FIGURE_PATH = PRESENTABLE_DIR / "motif_counts.png"
ORF_FIGURE_PATH = PRESENTABLE_DIR / "orf_overview.png"
DEFAULT_MOTIFS = ["ATG", "TAA", "TAG", "TGA", "TTTTAA"]
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class Orf:
    strand: str
    frame: int
    start: int
    end: int
    length_nt: int
    start_codon: str
    stop_codon: str
    sequence_preview: str


def ensure_directories() -> None:
    PRESENTABLE_DIR.mkdir(parents=True, exist_ok=True)


def get_pyplot():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Matplotlib is required for Mini-project 3 figures. "
            "Install it with: python3 -m pip install matplotlib"
        ) from error
    return plt


def load_fasta_records(path: Path = FASTA_PATH) -> list[FastaRecord]:
    fasta_text = path.read_text(encoding="utf-8")
    records = parse_fasta(fasta_text)
    if not records:
        raise ValueError(f"No FASTA records were found in {path}.")
    return records


def find_motif_positions(sequence: str, motif: str) -> list[int]:
    positions: list[int] = []
    start = 0
    while True:
        index = sequence.find(motif, start)
        if index == -1:
            break
        positions.append(index + 1)
        start = index + 1
    return positions


def summarize_motifs(sequence: str, motifs: list[str]) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    for motif in motifs:
        positions = find_motif_positions(sequence, motif)
        summaries.append(
            {
                "motif": motif,
                "count": len(positions),
                "positions": positions[:10],
            }
        )
    return summaries


def find_orfs_in_sequence(sequence: str, strand: str) -> list[Orf]:
    orfs: list[Orf] = []
    original_length = len(sequence)

    for frame_offset in range(3):
        start_index: int | None = None

        for index in range(frame_offset, len(sequence) - 2, 3):
            codon = sequence[index : index + 3]

            if start_index is None:
                if codon in START_CODONS:
                    start_index = index
                continue

            if codon in STOP_CODONS:
                orf_sequence = sequence[start_index : index + 3]
                preview = orf_sequence[:30] + ("..." if len(orf_sequence) > 30 else "")
                if strand == "+":
                    start_position = start_index + 1
                    end_position = index + 3
                else:
                    start_position = original_length - (index + 3) + 1
                    end_position = original_length - start_index
                orfs.append(
                    Orf(
                        strand=strand,
                        frame=frame_offset + 1,
                        start=start_position,
                        end=end_position,
                        length_nt=len(orf_sequence),
                        start_codon=sequence[start_index : start_index + 3],
                        stop_codon=codon,
                        sequence_preview=preview,
                    )
                )
                start_index = None

    return orfs


def find_open_reading_frames(sequence: str) -> list[Orf]:
    forward_orfs = find_orfs_in_sequence(sequence, strand="+")
    reverse_orfs = find_orfs_in_sequence(reverse_complement(sequence), strand="-")
    return sorted(forward_orfs + reverse_orfs, key=lambda orf: (-orf.length_nt, orf.start))


def sliding_window_gc_content(sequence: str, window_size: int = 500) -> tuple[list[int], list[float]]:
    positions: list[int] = []
    values: list[float] = []
    if len(sequence) < window_size:
        return [len(sequence) // 2], [gc_content(sequence)]

    step_size = max(1, window_size // 5)
    for start in range(0, len(sequence) - window_size + 1, step_size):
        window = sequence[start : start + window_size]
        positions.append(start + (window_size // 2) + 1)
        values.append(gc_content(window))

    return positions, values


def create_gc_content_figure(sequence: str, destination: Path = GC_FIGURE_PATH) -> Path:
    plt = get_pyplot()
    positions, values = sliding_window_gc_content(sequence)
    figure, axis = plt.subplots(figsize=(10, 5))
    axis.plot(positions, values, color="#1565C0", linewidth=2)
    axis.axhline(gc_content(sequence), color="#D32F2F", linestyle="--", linewidth=1.5, label="Genome GC%")
    axis.set_title("Sliding-Window GC Content")
    axis.set_xlabel("Genome position")
    axis.set_ylabel("GC content (%)")
    axis.grid(alpha=0.3, linestyle="--")
    axis.legend()
    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def create_motif_figure(
    motif_summaries: list[dict[str, object]],
    destination: Path = MOTIF_FIGURE_PATH,
) -> Path:
    plt = get_pyplot()
    figure, axis = plt.subplots(figsize=(8, 5))
    motifs = [summary["motif"] for summary in motif_summaries]
    counts = [summary["count"] for summary in motif_summaries]
    bars = axis.bar(motifs, counts, color=["#00897B", "#F4511E", "#6D4C41", "#5E35B1", "#546E7A"])
    axis.set_title("Motif Counts in SARS-CoV-2 Reference")
    axis.set_xlabel("Motif")
    axis.set_ylabel("Count")
    axis.grid(axis="y", alpha=0.3, linestyle="--")
    axis.set_axisbelow(True)

    for bar, count in zip(bars, counts):
        axis.text(bar.get_x() + bar.get_width() / 2, count, str(count), ha="center", va="bottom", fontsize=9)

    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def create_orf_figure(orfs: list[Orf], destination: Path = ORF_FIGURE_PATH) -> Path:
    plt = get_pyplot()
    top_orfs = orfs[:15]
    figure, axis = plt.subplots(figsize=(11, 5.5))

    colors = {"+" : "#2E7D32", "-" : "#C62828"}
    for index, orf in enumerate(reversed(top_orfs)):
        axis.barh(
            y=index,
            width=orf.length_nt,
            left=orf.start,
            color=colors[orf.strand],
            alpha=0.85,
        )

    axis.set_yticks(range(len(top_orfs)))
    axis.set_yticklabels([
        f"{orf.strand} frame {orf.frame} ({orf.length_nt} nt)"
        for orf in reversed(top_orfs)
    ])
    axis.set_xlabel("Genome position")
    axis.set_ylabel("Top ORFs")
    axis.set_title("Longest Open Reading Frames")
    axis.grid(axis="x", alpha=0.3, linestyle="--")
    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)
    return destination


def build_report(record: FastaRecord) -> str:
    sequence = record.sequence.upper()
    reverse_comp = reverse_complement(sequence)
    motif_summaries = summarize_motifs(sequence, DEFAULT_MOTIFS)
    orfs = find_open_reading_frames(sequence)

    lines = [
        "Mini-project 3: Sequence utilities",
        f"Header: {record.header}",
        f"Sequence length: {len(sequence)} bases",
        f"GC content: {gc_content(sequence):.2f}%",
        "Figures:",
        f"- {GC_FIGURE_PATH.name}",
        f"- {MOTIF_FIGURE_PATH.name}",
        f"- {ORF_FIGURE_PATH.name}",
        "",
        "Reverse complement",
        f"Length: {len(reverse_comp)} bases",
        f"First 80 bases: {reverse_comp[:80]}",
        f"Last 80 bases: {reverse_comp[-80:]}",
        "",
        "Motif search",
    ]

    for summary in motif_summaries:
        position_text = ", ".join(str(position) for position in summary["positions"]) or "none"
        if summary["count"] > len(summary["positions"]):
            position_text += ", ..."
        lines.append(
            f"Motif {summary['motif']}: count={summary['count']}, first_positions={position_text}"
        )

    lines.extend(
        [
            "",
            "Open reading frames",
            "Legend:",
            "  strand '+' = forward/reference strand",
            "  strand '-' = reverse-complement strand",
            "  frame 1 = codons start at base 1, frame 2 = shifted by 1 base, frame 3 = shifted by 2 bases",
            f"  start codons checked: {', '.join(sorted(START_CODONS))}",
            f"  stop codons checked: {', '.join(sorted(STOP_CODONS))}",
        ]
    )
    if orfs:
        for orf in orfs[:10]:
            lines.append(
                (
                    f"strand={orf.strand}, frame={orf.frame}, start={orf.start}, end={orf.end}, "
                    f"length_nt={orf.length_nt}, start_codon={orf.start_codon}, stop_codon={orf.stop_codon}, "
                    f"preview={orf.sequence_preview}"
                )
            )
        if len(orfs) > 10:
            lines.append(f"... {len(orfs) - 10} additional ORFs not shown")
    else:
        lines.append("No ORFs found.")

    return "\n".join(lines) + "\n"


def save_report(record: FastaRecord, destination: Path = REPORT_PATH) -> Path:
    ensure_directories()
    destination.write_text(build_report(record), encoding="utf-8")
    return destination


def create_figures(record: FastaRecord) -> dict[str, Path]:
    sequence = record.sequence.upper()
    motif_summaries = summarize_motifs(sequence, DEFAULT_MOTIFS)
    orfs = find_open_reading_frames(sequence)
    return {
        "gc_figure": create_gc_content_figure(sequence),
        "motif_figure": create_motif_figure(motif_summaries),
        "orf_figure": create_orf_figure(orfs),
    }


def run() -> dict[str, Path]:
    records = load_fasta_records()
    report_path = save_report(records[0])
    figure_paths = create_figures(records[0])
    return {"report_path": report_path, **figure_paths}


def main() -> None:
    try:
        outputs = run()
    except ModuleNotFoundError as error:
        print(error)
        return

    print("Mini-project 3 completed.")
    print(f"Report: {outputs['report_path']}")
    print(f"GC content figure: {outputs['gc_figure']}")
    print(f"Motif figure: {outputs['motif_figure']}")
    print(f"ORF figure: {outputs['orf_figure']}")


if __name__ == "__main__":
    main()
