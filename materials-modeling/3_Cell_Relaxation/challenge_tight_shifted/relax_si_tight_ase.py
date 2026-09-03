#!/usr/bin/env python3
"""Challenge II.2.3/II.3: shifted Si cell, tight ASE-driven QE relaxation."""
import os

import numpy as np
from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.filters import UnitCellFilter
from ase.io import write
from ase.optimize import BFGS

os.environ['OMP_NUM_THREADS'] = '1'
atoms = Atoms(
    'Si2',
    positions=[[1.40745, 4.04235, 4.09235], [0.0, 0.0, 0.0]],
    cell=[[0.0, 2.7149, 2.7149], [2.7149, 0.0, 2.7149], [2.7149, 2.7149, 0.0]],
    pbc=True,
)
profile = EspressoProfile('mpirun -np 4 pw.x', pseudo_dir='../..')
atoms.calc = Espresso(
    profile=profile,
    pseudopotentials={'Si': 'Si.upf'},
    input_data={
        'control': {'calculation': 'scf', 'prefix': 'si_tight_ase', 'outdir': './tmp_ase',
                    'tstress': True, 'tprnfor': True, 'verbosity': 'low'},
        'system': {'ecutwfc': 65, 'ibrav': 0, 'occupations': 'smearing', 'smearing': 'gauss', 'degauss': 0.01},
        'electrons': {'conv_thr': 1.0e-12},
    },
    kpts=(15, 15, 15),
)
filter_atoms = UnitCellFilter(atoms, scalar_pressure=0.0)
opt = BFGS(filter_atoms, trajectory='si_tight_ase.traj', logfile='si_tight_ase.log')
opt.run(fmax=0.005, steps=100)

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress(voigt=False)
a_cubic = np.linalg.norm(atoms.cell[0]) * np.sqrt(2.0)
write('si_tight_ase_final.vasp', atoms, direct=False)
print('Challenge II.2.3 / II.3 — shifted-cell tight ASE→QE relaxation')
print(f'Total energy: {energy:.6f} eV')
print(f'Equivalent cubic lattice parameter: {a_cubic:.6f} Å')
print(f'Max |force component|: {np.abs(forces).max():.6f} eV/Å')
print(f'Mean normal stress: {np.trace(stress) / 3.0 * 160.21766208:.6f} GPa')
print('Convergence target: fmax = 0.005 eV/Å; QE SCF conv_thr = 1e-12')
print('Final structure: si_tight_ase_final.vasp')
