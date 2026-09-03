Challenges II.2.3 and II.3 — shifted Si cell with tight thresholds
===================================================================

Both workflows start from the same deliberately displaced Si atom and the
experimental primitive cell.  The direct QE ``vc-relax`` uses
``etot_conv_thr = 1e-8 Ry``, ``forc_conv_thr = 1e-5 Ry/bohr`` and
``press_conv_thr = 0.1 kbar``.  The ASE-driven QE BFGS calculation uses
``conv_thr = 1e-12`` and ``fmax = 0.005 eV/Å``.

Run the direct and ASE variants separately::

   mpirun -np 4 pw.x -in si_vc_relax_tight.in > si_vc_relax_tight.out
   /home/danil/chem_venv/bin/python relax_si_tight_ase.py > relax_si_tight_ase.out

The two output logs and final structures provide an independent convergence
check after the manual displacement.
