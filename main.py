from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
MPLCONFIG_DIR = ROOT / ".mplconfig"
XDG_CACHE_DIR = ROOT / ".cache"


@dataclass(frozen=True)
class Project:
    slug: str
    group: str
    script: Path
    summary: str
    aliases: tuple[str, ...] = ()


PROJECTS: tuple[Project, ...] = (
    Project(
        slug="fasta_parser",
        group="sequence_fundamentals",
        script=ROOT / "sequence_fundamentals/01_fasta_parser/01_fasta_parser.py",
        summary="FASTA parsing, GC summary, and base-composition plot",
        aliases=("01_fasta_parser", "sequence_fundamentals/01_fasta_parser"),
    ),
    Project(
        slug="fastq_parser",
        group="sequence_fundamentals",
        script=ROOT / "sequence_fundamentals/02_fastq_parser/02_fastq_parser.py",
        summary="FASTQ parsing and quality-by-position summary",
        aliases=("02_fastq_parser", "sequence_fundamentals/02_fastq_parser"),
    ),
    Project(
        slug="sequence_utilities",
        group="sequence_fundamentals",
        script=ROOT / "sequence_fundamentals/03_sequence_utilities/03_sequence_utilities.py",
        summary="Reverse complement, motif search, GC windows, and ORFs",
        aliases=("03_sequence_utilities", "sequence_fundamentals/03_sequence_utilities"),
    ),
    Project(
        slug="count_matrix_summary",
        group="rna_seq_fundamentals",
        script=ROOT / "rna_seq_fundamentals/01_count_matrix_summary/01_count_matrix_summary.py",
        summary="Real GEO count matrix summary and sample-level plots",
        aliases=("01_count_matrix_summary", "rna_seq_fundamentals/01_count_matrix_summary"),
    ),
    Project(
        slug="normalization",
        group="rna_seq_fundamentals",
        script=ROOT / "rna_seq_fundamentals/02_normalization/02_normalization.py",
        summary="CPM, FPKM, and TPM comparison on GEO counts",
        aliases=("02_normalization", "rna_seq_fundamentals/02_normalization"),
    ),
    Project(
        slug="fold_change",
        group="rna_seq_fundamentals",
        script=ROOT / "rna_seq_fundamentals/03_fold_change/03_fold_change.py",
        summary="Grouped RNA-seq fold-change example on real GEO samples",
        aliases=("03_fold_change", "rna_seq_fundamentals/03_fold_change"),
    ),
    Project(
        slug="volcano_ready_table",
        group="rna_seq_fundamentals",
        script=ROOT / "rna_seq_fundamentals/04_volcano_ready_table/04_volcano_ready_table.py",
        summary="Volcano-ready table and preview plot",
        aliases=("04_volcano_ready_table", "rna_seq_fundamentals/04_volcano_ready_table"),
    ),
    Project(
        slug="hisat2_alignment",
        group="rna_seq_tool_workflows",
        script=ROOT / "rna_seq_tool_workflows/01_hisat2_alignment/01_hisat2_alignment.py",
        summary="HISAT2 alignment workflow on official example reads",
        aliases=("01_hisat2_alignment", "rna_seq_tool_workflows/01_hisat2_alignment"),
    ),
    Project(
        slug="star_alignment",
        group="rna_seq_tool_workflows",
        script=ROOT / "rna_seq_tool_workflows/02_star_alignment/02_star_alignment.py",
        summary="STAR alignment workflow and splice-junction summary",
        aliases=("02_star_alignment", "rna_seq_tool_workflows/02_star_alignment"),
    ),
    Project(
        slug="htseq_counting",
        group="rna_seq_tool_workflows",
        script=ROOT / "rna_seq_tool_workflows/03_htseq_counting/03_htseq_counting.py",
        summary="HTSeq gene counting on a yeast example alignment",
        aliases=("03_htseq_counting", "rna_seq_tool_workflows/03_htseq_counting"),
    ),
    Project(
        slug="stringtie_quantification",
        group="rna_seq_tool_workflows",
        script=ROOT / "rna_seq_tool_workflows/04_stringtie_quantification/04_stringtie_quantification.py",
        summary="StringTie transcript assembly and gene abundance summary",
        aliases=("04_stringtie_quantification", "rna_seq_tool_workflows/04_stringtie_quantification"),
    ),
    Project(
        slug="kallisto_quantification",
        group="rna_seq_tool_workflows",
        script=ROOT / "rna_seq_tool_workflows/05_kallisto_quantification/05_kallisto_quantification.py",
        summary="kallisto pseudoalignment and transcript quantification",
        aliases=("05_kallisto_quantification", "rna_seq_tool_workflows/05_kallisto_quantification"),
    ),
    Project(
        slug="bustools_single_cell",
        group="rna_seq_tool_workflows",
        script=ROOT / "rna_seq_tool_workflows/06_bustools_single_cell/06_bustools_single_cell.py",
        summary="kallisto|bustools single-cell workflow and cell-by-gene matrix",
        aliases=("06_bustools_single_cell", "rna_seq_tool_workflows/06_bustools_single_cell"),
    ),
)


def build_lookup() -> dict[str, Project]:
    lookup: dict[str, Project] = {}
    for project in PROJECTS:
        lookup[project.slug] = project
        for alias in project.aliases:
            lookup[alias] = project
    return lookup


def resolve_projects(requested: list[str]) -> list[Project]:
    if not requested:
        return list(PROJECTS)

    lookup = build_lookup()
    resolved: list[Project] = []
    seen: set[str] = set()

    for item in requested:
        if item in {"all", "*"}:
            for project in PROJECTS:
                if project.slug not in seen:
                    resolved.append(project)
                    seen.add(project.slug)
            continue

        project = lookup.get(item)
        if project is None:
            valid = ", ".join(project.slug for project in PROJECTS)
            raise SystemExit(f"Unknown project '{item}'. Use one of: {valid}")

        if project.slug not in seen:
            resolved.append(project)
            seen.add(project.slug)

    return resolved


def print_project_list() -> None:
    current_group = ""
    for project in PROJECTS:
        if project.group != current_group:
            current_group = project.group
            print(f"\n{current_group}:")
        print(f"  - {project.slug}: {project.summary}")


def run_project(project: Project) -> int:
    print(f"\n=== Running {project.slug} ===", flush=True)
    MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
    XDG_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["MPLCONFIGDIR"] = str(MPLCONFIG_DIR)
    env["XDG_CACHE_HOME"] = str(XDG_CACHE_DIR)
    result = subprocess.run(
        [sys.executable, str(project.script)],
        cwd=ROOT,
        check=False,
        env=env,
    )
    if result.returncode == 0:
        print(f"Completed {project.slug}", flush=True)
    else:
        print(f"Failed {project.slug} with exit code {result.returncode}", flush=True)
    return result.returncode


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one genomics mini-project or the full collection."
    )
    parser.add_argument(
        "projects",
        nargs="*",
        help="Optional project slug(s) to run. If omitted, all projects run in order.",
    )
    parser.add_argument(
        "-p",
        "--project",
        dest="project_flags",
        action="append",
        default=[],
        help="Project slug to run. Repeat the flag to run more than one project.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Show available project slugs and exit.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.list:
        print_project_list()
        return

    requested_projects = list(args.projects) + list(args.project_flags)
    selected_projects = resolve_projects(requested_projects)
    failures = [project.slug for project in selected_projects if run_project(project) != 0]

    if failures:
        raise SystemExit(f"Some projects failed: {', '.join(failures)}")


if __name__ == "__main__":
    main()
