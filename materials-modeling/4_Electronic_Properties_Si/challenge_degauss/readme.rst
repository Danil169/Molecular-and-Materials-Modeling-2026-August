Challenge II.4 — Gaussian broadening in Si DOS
===============================================

This directory reruns the Si SCF/NSCF/DOS workflow at Gaussian ``degauss``
values of 0.005, 0.010 and 0.020 Ry.  Each calculation uses PBE, the
repository ONCV ``Si.upf``, 65 Ry, a 15x15x15 SCF mesh and a 24x24x24 NSCF
mesh.  ``degauss_summary.csv`` records energy, Fermi level and DOS at Fermi;
the three ``dos_*.dat`` files preserve the broadened spectra.

Run::

   /home/danil/chem_venv/bin/python compare_degauss.py > compare_degauss.out
