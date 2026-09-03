==============================
Measuring software performance
==============================

Exercise I.2: benchmark NWChem, Quantum ESPRESSO, and MOPAC on the local machine
to understand MPI/thread scaling.  Results below are for the machine used during the school.
See subdirectories for full input/output files.

MOPAC — "mini-DNA" PM7 thread scaling
--------------------------------------

System: DNA fragment (~300 atoms), PM7, various thread counts.

Notebook (Intel Core i5-12450H)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Nthr", "Wall time (s)"

   1,  27.50
   6,  18.49
   12, 18.34

Desktop PC (Intel Core i7-12700K)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Nthr", "Wall time (s)"

   3,  30.58
   6,  10.82
   10, 15.09
   12, 15.61

Optimal performance at ~6 threads; beyond that, returns diminish.
See also: http://openmopac.net/Manual/Reducing_computation_time.html

NWChem — CH₃(•) radical ZORA-B3LYP MPI scaling
------------------------------------------------

System: CH3 radical, B3LYP/6-311G with ZORA scalar-relativistic corrections.
NWChem 7.2.2 (Ubuntu 24.04.2), Intel Core i5-12450H.

.. csv-table::
   :header: "Nproc (MPI)", "Wall time (s)"

   2, 11.5
   4,  7.5
   6,  7.7

Details in ``nwchem/primitive/readme.rst``.

Quantum ESPRESSO — Si SCF MPI scaling (Exercise I.2 / Fig. 1)
--------------------------------------------------------------

System: Si unit cell, PBE, ecutwfc = 60 Ry, 50 k-points, 34 bands.
Measured 2026-08-26, 12th Gen Intel Core i5-12450H.

.. csv-table::
   :header: "Nproc (MPI)", "Wall time (s)", "CPU time (s)"

   1,  45.19, 44.64
   2,  31.60, 30.55
   4,  20.67, 19.62
   8,  12.09, 11.75
   16, 11.22, 10.86
   24, 10.08,  9.78

Speed-up from 1 → 8 MPI: ×3.7.  Saturation visible at 16–24 MPI.
Full data in ``quantum_espresso/benchmark_results.json``; plot: ``../report_figures_part1/qe_mpi_scaling.png``.

Challenge
=========
On your personal computer, measure software performance for all three codes and save
results in this or a child ``readme.rst`` file.  Compare conda-installed vs
system-package binaries (MOPAC and NWChem show measurable differences due to
linked OpenMP runtime; see ``mopac/readme.rst`` for details).
