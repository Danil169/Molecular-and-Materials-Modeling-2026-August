Challenge II.2.2 — final Si force/stress calculation
=====================================================

``Si_force_stress.in`` is a direct QE SCF calculation using the local
converged settings (65 Ry and 15x15x15 k-points).  ``tstress`` and ``tprnfor``
are enabled, so the saved ``Si_force_stress.out`` contains the final total
energy, stress tensor and forces.

Run::

   mpirun -np 4 pw.x -in Si_force_stress.in > Si_force_stress.out
