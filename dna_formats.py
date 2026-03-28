from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


DNA_COMPLEMENT = str.maketrans("ATCGNatcgn", "TAGCNtagcn")


@dataclass
class FastaRecord:
    header: str
    sequence: str


@dataclass
class FastqRecord:
    header: str
    sequence: str
    quality: str

    def phred_scores(self) -> list[int]:
        return [ord(char) - 33 for char in self.quality]

    def average_quality(self) -> float:
        scores = self.phred_scores()
        return sum(scores) / len(scores) if scores else 0.0


def parse_fasta(text: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(FastaRecord(header=header, sequence="".join(chunks)))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)

    if header is not None:
        records.append(FastaRecord(header=header, sequence="".join(chunks)))

    return records


def parse_fastq(text: str) -> list[FastqRecord]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) % 4 != 0:
        raise ValueError("FASTQ input must contain groups of 4 non-empty lines.")

    records: list[FastqRecord] = []
    for index in range(0, len(lines), 4):
        header, sequence, plus, quality = lines[index : index + 4]

        if not header.startswith("@"):
            raise ValueError(f"Invalid FASTQ header line: {header}")
        if plus != "+":
            raise ValueError(f"Invalid FASTQ separator line: {plus}")
        if len(sequence) != len(quality):
            raise ValueError("FASTQ sequence and quality lengths must match.")

        records.append(
            FastqRecord(
                header=header[1:].strip(),
                sequence=sequence,
                quality=quality,
            )
        )

    return records


def gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0

    return 100 * gc_count(sequence) / len(sequence)


def gc_count(sequence: str) -> int:
    return sum(1 for base in sequence.upper() if base in {"G", "C"})


def nucleotide_counts(sequence: str) -> dict[str, int]:
    counts = {base: 0 for base in ["A", "T", "C", "G", "N"]}

    for base in sequence.upper():
        if base in counts:
            counts[base] += 1
        else:
            counts["N"] += 1

    return counts


def total_sequence_length(records: Iterable[FastaRecord]) -> int:
    return sum(len(record.sequence) for record in records)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def summarize_fasta(records: list[FastaRecord]) -> list[str]:
    output: list[str] = []
    for record in records:
        output.append(
            (
                f"{record.header}: length={len(record.sequence)}, "
                f"GC={gc_content(record.sequence):.2f}%, "
                f"reverse_complement={reverse_complement(record.sequence)}"
            )
        )
    return output


def summarize_fastq(records: list[FastqRecord]) -> list[str]:
    output: list[str] = []
    for record in records:
        output.append(
            (
                f"{record.header}: length={len(record.sequence)}, "
                f"avg_quality={record.average_quality():.2f}, "
                f"phred_scores={record.phred_scores()}"
            )
        )
    return output


SAMPLE_FASTA = """>seq1 human_example
ATGCGTACGTAGCTAG
>seq2 gc_rich_example
GGGCCCGGTTAA
"""


SAMPLE_FASTQ = """@read1
ATGCATGC
+
IIIIIIII
@read2
GGTTCCAA
+
!''*(())
"""


def demo_report() -> str:
    fasta_records = parse_fasta(SAMPLE_FASTA)
    fastq_records = parse_fastq(SAMPLE_FASTQ)

    lines = [
        "Mini-project 1: FASTA and FASTQ basics",
        "",
        "FASTA summary",
        *summarize_fasta(fasta_records),
        "",
        "FASTQ summary",
        *summarize_fastq(fastq_records),
        "",
        "Next good extension: add motif search and open reading frame detection.",
    ]
    return "\n".join(lines)
