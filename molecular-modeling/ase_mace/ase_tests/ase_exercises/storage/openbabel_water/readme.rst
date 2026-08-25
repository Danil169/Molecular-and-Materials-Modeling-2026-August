====================
ASE driven OpenBabel 
====================

ASE driven OpenBabel optimization of water molecule.

install & run
-------------
pip install openbabel-wheel (Version: 3.1.1.23)

python ase_openbabel_water_optim.py

 Selected FF : mmff94
Energy before optimization: 10.33855133442349 kJ/mol
Energy after optimization: 8.729929234125912e-07 kJ/mol
 Optimized geometry:
O-H1: 0.9690
O-H2: 0.9690
H-O-H Angle (deg): 103.98

Experimental water geometry (Google AI) :
 bond angle: cca 104.45 deg - 104.0 deg
 bond length: 0.9578 Ang


Challenge
---------
Modify the Python program: make a loop of all available force fields and print geometry coordinates of the molecule.

Which ff does give the best geometry for the chosen molecule ?

Try this approach for another simple molecule (ethane, benzene, SF6, BF3, ..).

