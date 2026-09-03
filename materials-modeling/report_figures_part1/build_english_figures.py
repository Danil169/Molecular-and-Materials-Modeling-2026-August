#!/usr/bin/env python3
"""Render the English-language Part I figures from the saved calculation data."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent


def build_qe_scaling() -> None:
    data_path = ROOT / "molecular-modeling/software_performance/quantum_espresso/benchmark_results.json"
    records = json.loads(data_path.read_text())
    nproc = np.array([1, 2, 4, 8, 16, 24])
    runs = []
    for offset in (0, 6):
        runs.append(np.array([records[offset + i]["wall_time"] for i in range(6)]))

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(nproc, runs[0], "o-", color="#3b6ea8", lw=2.2, ms=7, label="Run 1")
    ax.plot(nproc, runs[1], "s--", color="#c9534c", lw=2.0, ms=6.5, label="Run 2")
    ax.plot(nproc, runs[0][0] / nproc, ":", color="black", lw=1.8, label="Ideal scaling")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(nproc)
    ax.set_xticklabels([str(x) for x in nproc])
    ax.set_xlabel("Number of MPI Processes")
    ax.set_ylabel("Runtime (s)")
    ax.set_title("Quantum ESPRESSO: Tl$_2$ SCF Scaling\n(50 k-points, spin-orbit coupling, OMP=1)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "qe_mpi_scaling.png", dpi=300)


def build_h2_curve() -> None:
    data_path = ROOT / "molecular-modeling/ase_mace/ase_tests/ase_exercises/H2_molecule/mace_calculator/H2_potential_curve_MACE-MP-0_small.txt_SAVED"
    data = np.loadtxt(data_path, skiprows=7)
    r, e = data[:, 0], data[:, 1]
    re, de = 0.75399, 4.12746

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(r, e, "o-", color="#3b6ea8", lw=2.0, ms=5.5)
    ax.axvline(re, color="#c9534c", ls="--", lw=1.6, label=r"Experimental $r_e = 0.7414$ Å")
    idx = np.argmin(np.abs(r - re))
    ax.plot(r[idx], e[idx], "o", color="#3b6ea8")
    ax.annotate(
        rf"$r_e = {re:.4f}$ Å" + "\n" + rf"$D_e = {de:.3f}$ eV" + "\n" + r"(experimental $D_e = 4.52$ eV)",
        xy=(r[idx], e[idx]), xytext=(1.18, -6.2),
        arrowprops={"arrowstyle": "->", "lw": 1.2}, fontsize=11,
    )
    ax.set_xlabel("Internuclear Distance r(H–H) (Å)")
    ax.set_ylabel("Energy, E (eV)")
    ax.set_title("H$_2$ Potential Energy Curve (MACE-MP-0)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "h2_potential_curve_mace.png", dpi=300)


if __name__ == "__main__":
    build_qe_scaling()
    build_h2_curve()
