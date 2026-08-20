Recalculation of ionization potential of Hg (v0.1)
==================================================

This computational exercise builds on our previous paper: https://doi.org/10.1016/S0009-2614(99)00665-X.
We focused on recalculating the ionization potential (IP) of the mercury atom using several software packages installed via the conda‑forge channel, namely MOPAC, xTB, NWChem, and PySCF. Python scripts were developed within the Atomic Simulation Environment (ASE) with assistance from online AI tools.

The experimental ionization energy of mercury is 10.437 eV. The semi‑empirical PM7 method in MOPAC, using a UHF reference, yields an excellent value of 10.518 eV, indicating well‑fitted parameters even for heavy elements like Hg. In contrast, the GFN1‑xTB extended tight‑binding semi‑empirical DFT method significantly overestimates the IP, giving 16.768 eV.

We then performed ab‑initio calculations based on ECP‑SCF single‑determinant wave functions, where the ECP corresponds to the DK scalar relativistic treatment used in the reference paper. NWChem CCSD(T) results show systematic improvement with basis set size: def2‑SVP gives 10.018 eV, def2‑TZVP yields 10.139 eV, and def2‑QZVP provides the best value at 10.325 eV. Similarly, PySCF produces an ECP‑UHF‑CCSD(T) value of 10.125 eV with the def2‑TZVP basis set.

All working files are available at:
https://github.com/miroi/Molecular-and-Materials-Modeling-2026-August/tree/main/molecular-modeling/abstract_ionization_potentials
