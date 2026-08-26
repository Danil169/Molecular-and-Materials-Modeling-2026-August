===============================
Ethane C-C bond energy with ASE
===============================

see codes construction in https://share.google/aimode/nJUzk4Gqx6YMnO99N

The C-C bond dissociation energy in ethane is cca 348-377 kJ/mol.
Experimental BDE: Approximately 376 kJ/mol

The C-C equilibrium bond length is 1.53-1.54 A.
Experimental C-C Bond Length: Typically cited at 1.522–1.54 Å.

run
---
python code3.py

Method       | BDE (kJ/mol)    | C-C Length (Å)
------------------------------------------------
EMT          | 208.32          | 1.0089
CHGNet       | 314.79          | 1.5351
MACE         | 360.37          | 1.5203


Avogadro try
~~~~~~~~~~~~
kcal/mol ??
          Ethane    CH3.
MMFF94    -19.8086  none
UFF        0.811     0   

challenge
---------
modify the program to compute C-H bond energy in methane (form the reaction of CH4 -> CH3. + .H)

