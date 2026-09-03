Challenges II.1 and II.2.4 — Si pseudopotential comparison
============================================================

``compare_si_pseudopotentials.py`` repeats short, explicit PBE convergence
series for the repository ONCV ``Si.upf`` and the QE-distributed RRKJ
``Si.pbe-rrkj.UPF``.  The identical two-atom Si cell is used in both cases.
The stdout log and CSV retain total energies, cutoff changes and k-grid
changes; absolute energies are not compared between different pseudopotentials.

Run with the ASE/QE environment::

   /home/danil/chem_venv/bin/python compare_si_pseudopotentials.py > compare_si_pseudopotentials.out
