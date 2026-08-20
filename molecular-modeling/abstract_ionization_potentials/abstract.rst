Recalculation of ionization potential of Hg (v0.1)
==================================================

This computational molecular modeling exercise follows our previous paper https://doi.org/10.1016/S0009-2614(99)00665-X.

We focused on recalculating the ionization potential of the mercury atom.

For that, we use available software, installed via conda-forge channel, namely MOPAC, xTB, NWChem and PySCF. 
Python scripts were developed withing the Atomic Simulation Environment (ASE), with the help of online AI.

The experimental value of mercury ionization energy is 10.437 eV.

Semiemprical PM7 method of MOPAC gives on UHF reference very good value of 10.518 eV, what indicated well fitted parameters also for heavy element like Hg.
The GFN1-xTB (extended tight-binding) semiemprical DFT method strongly exceeds IP, 16.768 eV.

Next group of calculations are based on ECP-SCF single determinant wave-function.
NWChem CCSD(T) values show proper improvement with the size of basis set.
The def2-svp gives 10.018 eV, larger basis def2-tzvp gives 10.139 eV and the largest basis def2-qzvp brings 10.325 eV, what is the best value.

Similarly PySCF code  gives ECP-UHF-CCSD(T) value  10.125 eV  in def2-tzvp basis set. 


All working files are in https://github.com/miroi/Molecular-and-Materials-Modeling-2026-August/tree/main/molecular-modeling/abstract_ionization_potentials .
