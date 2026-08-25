===============================
Ethane C-C bond energy with ASE
===============================

The C-C bond dissociation energy in ethane is cca 348-377 kJ/mol.

It is calculated by comparing the total energy of optimized ethane with two methyl radicals,
often using methods like DFT (e.g., B3LYP) or CCSD(T) with basis sets like cc-pVTZ. 

The C-C equilibrium bond length is 1.53-1.54 A.

run
---
python ethane_CC_bond_en.py
Ethane Energy: 1.4622 eV
Optimized C-C bond length in ethane: 1.0096 Å
Methyl Radical Energy: 1.8106 eV
------------------------------
Calculated C-C Bond Energy: 2.1591 eV
Calculated C-C Bond Energy: 208.32 kJ/mol
------------------------------


Avogadro try
~~~~~~~~~~~~
kcal/mol ??
          Ethane    CH3.
MMFF94    -19.8086  none
UFF        0.811     0   

challenge
---------
modify the program to compute C-H bond energy in methane (form the reaction of CH4 -> CH3. + .H)


