#!/usr/bin/env python3
"""Create concise README files for tracked directories that do not have one.

Run from the repository root after adding a new calculation directory.
Existing README files are never modified.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def tracked_directories() -> list[Path]:
    files = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True)
    directories: set[Path] = set()
    for name in files.splitlines():
        directory = (ROOT / name).parent
        while directory != ROOT:
            directories.add(directory)
            directory = directory.parent
    return sorted(directories)


def description(relative: Path) -> str:
    parts = set(relative.parts)
    name = relative.name.replace("_", " ")
    if relative.as_posix() == "materials-modeling/report_figures_part1":
        return "the two Part I figures embedded in the final English and Russian reports"
    if relative.as_posix() == "molecular-modeling/software_performance/quantum_espresso":
        return "Quantum ESPRESSO MPI-scaling benchmark inputs and measured timings"
    if relative.as_posix().endswith("software_checks/basic_functionality_checks"):
        return "end-to-end software-verification scripts and their saved outputs"
    if relative.as_posix().endswith("5_Electronic_Properties_Al/Challenge_Al"):
        return "the extended aluminium convergence, relaxation, and electronic-properties challenge"
    if relative.as_posix().endswith("8_Adsorption/electronic"):
        return "electronic-structure outputs for the hydrogen-on-graphene adsorption calculation"
    if "reference_outputs" in parts:
        return "reference outputs retained for comparison with the reproduced calculation"
    if "outputs" in parts:
        return "saved program outputs produced by the calculation in the parent directory"
    if "scf_files" in parts:
        return "self-consistent-field files used as input to the subsequent analysis"
    if "pdos_results" in parts:
        return "projected-density-of-states data generated from the parent calculation"
    if "potential_results" in parts:
        return "electrostatic-potential and planar-averaging outputs"
    if "vib" in name.lower():
        return "finite-displacement vibrational-analysis cache files"
    if "stored_files" in parts:
        return "archived copies of inputs and outputs retained for provenance"
    if "benchmark" in name.lower() or "performance" in parts:
        return "benchmark inputs, scripts, and measured performance outputs"
    if "mopac" in parts or "nwchem" in parts or "pyscf" in parts or "xtb" in parts:
        return "method-specific inputs, scripts, and output logs for this calculation"
    return "inputs, scripts, and output artifacts for this calculation step"


def write_readme(directory: Path) -> None:
    relative = directory.relative_to(ROOT)
    title = f"{relative.name.replace('_', ' ')}"
    underline = "-" * len(title)
    text = f"""{title}\n{underline}\n\nThis directory contains {description(relative)}.\n\nContents\n--------\n\nFiles are retained as part of the reproducible course-work record. Inputs and\nscripts define the calculation; log, trajectory, and data files record its\nresults. Consult the nearest parent ``readme.rst`` for the exercise-level\nworkflow and interpretation.\n"""
    (directory / "readme.rst").write_text(text, encoding="utf-8")


def main() -> None:
    for directory in tracked_directories():
        if directory == ROOT or (directory / "readme.rst").exists():
            continue
        write_readme(directory)


if __name__ == "__main__":
    main()
