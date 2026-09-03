===========================
NWChem software performance
===========================

Exercise I.2: NWChem MPI-scaling benchmark using the CH₃(•) radical
(ZORA-B3LYP/6-311G, scalar-relativistic DFT with property calculations).

Two subdirectories are provided:

``primitive/``
    Direct bash-driven NWChem runs; raw timing data.

``python-driven/``
    Python-driven benchmark script with AI-assisted scripting example.

Performance data (Intel Core i5-12450H)
-----------------------------------------

NWChem 7.0.2 (Ubuntu 22.04):

.. csv-table::
   :header: "Nproc (MPI)", "Wall time (s)"

   2,  9.9
   4,  6.8
   6,  6.6

NWChem 7.2.2 (Ubuntu 24.04.2):

.. csv-table::
   :header: "Nproc (MPI)", "Wall time (s)"

   2, 11.5
   4,  7.5
   6,  7.7

The newer 7.2.2 package is slightly slower here; both saturate at 4–6 MPI
for this small system.  Full inputs and outputs are in ``primitive/``.

