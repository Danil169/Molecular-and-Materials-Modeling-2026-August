Challenge I.7 — QE relaxation of methane
=========================================

``relax_ch4_qe.py`` adapts the water workflow in the parent directory to
methane.  It uses ASE to drive Quantum ESPRESSO (PBE PAW, 46 Ry, Gamma point,
12 Å cubic box) and BFGS to ``fmax = 0.01 eV/Å``.

Run::

   python3 relax_ch4_qe.py > relax_ch4_qe.out

The calculation log, final XYZ structure and stdout record are retained in
this directory.  The analysis compares the mean C-H length and H-C-H angle
with gas-phase experimental values (1.087 Å and 109.47 degrees).
