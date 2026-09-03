============================================
ASE driven Q.E. relaxation of water molecule
============================================

python qe_h2o_optimize.py

Optimized Bond Lengths (A):
O-H1: 0.9734
O-H2: 0.9734
H-O-H Angle (deg): 104.51
Final Energy (eV): -598.7717

Experimental water geometry (Google AI) :
 bond angle: cca 104.45 deg - 104.0 deg
 bond length: 0.9578 Ang


Challenge
---------
Adapt the code for the ASE-driven Q.E. carried out geometry relaxation of CH4 molecule.
Compare computed data with experimental values.

**Status: ✓ Completed.** See ``challenge_ch4_qe/`` for script, output, and results.

Key CH₄ results (PBE, 46 Ry, Gamma point, 12 Å box):

- Total energy: −315.702 eV
- Mean C–H bond: 1.0957 Å (exp. 1.087 Å, error 0.8 %)
- H–C–H angle: 109.47° (exp. 109.47°, perfect tetrahedral)
- Max force: 0.006 eV/Å
