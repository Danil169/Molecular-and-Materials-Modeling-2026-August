Challenge II.7 — effect of ``assume_isolated = '2D'``
======================================================

``compare_isolation.py`` repeats the graphene SCF → ``pp.x`` → ``average.x``
workflow twice, once with QE's ``assume_isolated = '2D'`` and once without it.
The two planar/macroscopic potential profiles are saved as ``avg_2d.dat`` and
``avg_none.dat``.  ``isolation_summary.csv`` records the total energy, Fermi
level, a central-vacuum plateau estimate and the resulting work function.

Run::

   /home/danil/chem_venv/bin/python compare_isolation.py > compare_isolation.out
