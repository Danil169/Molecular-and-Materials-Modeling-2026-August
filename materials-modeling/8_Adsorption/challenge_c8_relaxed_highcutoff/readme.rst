Challenge II.8 — relaxed C8 and 100 Ry adsorption energy
==========================================================

``relax_c8_highcutoff.py`` relaxes the clean C8 flake and H@C8 at 100 Ry,
then recomputes the isolated H energy at the same cutoff.  This makes the
reported binding energy internally consistent:
``E_bind = E(C8) + E(H) - E(H@C8)``.  The workflow otherwise retains the
parent calculation choices (PBE, DFT-D3, spin-polarization, Gamma point and
``fmax = 0.01 eV/Å``).

Run::

   /home/danil/chem_venv/bin/python relax_c8_highcutoff.py > relax_c8_highcutoff.out
