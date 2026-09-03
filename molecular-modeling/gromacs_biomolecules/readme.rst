================
Gromacs exercise
================

Exercise: MD simulation of crambin protein (PDB: 1CRN) with GROMACS.

Installation
------------

::

    conda install -c conda-forge gromacs

Verify::

    gmx --version > gmx_version.log

GROMACS 2026.3 is used (see ``gmx_version.log``).

Workflow (``protein_1crn/``)
-----------------------------

Three progressive Python scripts (ASE-driven GROMACS interface):

``gromacs_ase_example.py``
    Basic setup: read 1CRN.pdb, generate topology (CHARMM27/TIP3P), write .tpr.

``gromacs_ase_example_02.py``
    Full pipeline: solvation → energy minimisation → NVT equilibration.

``gromacs_ase_example_03.py``
    Extended analysis: RMSD, radius of gyration, temperature monitoring; 
    generates ``rmsd.png``, ``gyrate.png``, ``energy_profile.png``.

Key results (1CRN, 46 residues, 642 atoms, CHARMM27 + TIP3P)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Energy after minimisation: **−2582 kJ/mol**
- NVT equilibration (100 ps): T = 298.1 ± 11.2 K
- RMSD (Cα, vs crystal): **0.139 nm**
- Radius of gyration: **0.953 nm**

Challenge
---------
Use a different protein from the RCSB PDB and repeat the workflow.
Try the AMBER99sb-ildn force field instead of CHARMM27.

