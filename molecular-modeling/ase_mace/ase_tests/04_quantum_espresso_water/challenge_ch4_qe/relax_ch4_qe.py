#!/usr/bin/env python3
"""Challenge I.7: ASE-driven PBE/QE relaxation of methane.

The workflow intentionally follows ``qe_h2o_optimize.py``: Gamma-only
calculation in a 12 Å box, PBE PAW pseudopotentials, and BFGS to 0.01 eV/Å.
"""
import os

from ase.build import molecule
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.io import write
from ase.optimize import BFGS

os.environ['OMP_NUM_THREADS'] = '1'

atoms = molecule('CH4')
atoms.set_cell([12.0, 12.0, 12.0])
atoms.center()

profile = EspressoProfile('mpirun -np 4 pw.x', pseudo_dir='/usr/share/espresso/pseudo')
calc = Espresso(
    profile=profile,
    pseudopotentials={
        'C': 'C.pbe-n-kjpaw_psl.0.1.UPF',
        'H': 'H.pbe-kjpaw.UPF',
    },
    input_data={
        'control': {
            'calculation': 'scf', 'prefix': 'ch4_pbe', 'outdir': './tmp',
            'verbosity': 'low', 'tstress': True, 'tprnfor': True,
        },
        'system': {
            'ecutwfc': 46.0, 'ecutrho': 184.0, 'ibrav': 0, 'nosym': True, 'noinv': True,
            'occupations': 'smearing', 'smearing': 'gaussian', 'degauss': 0.02,
        },
        'electrons': {'conv_thr': 1.0e-8, 'mixing_beta': 0.7},
    },
    kpts=(1, 1, 1),
)
atoms.calc = calc

opt = BFGS(atoms, trajectory='ch4_relax.traj', logfile='ch4_relax.log')
opt.run(fmax=0.01, steps=100)

bond_lengths = [atoms.get_distance(0, i) for i in range(1, 5)]
angles = [atoms.get_angle(i, 0, j) for i in range(1, 5) for j in range(i + 1, 5)]
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
fmax = abs(forces).max()
write('ch4_final.xyz', atoms)

print('Challenge I.7 — QE/ASE methane relaxation')
print(f'Total energy: {energy:.6f} eV')
print(f'Mean C-H bond: {sum(bond_lengths) / len(bond_lengths):.4f} Å')
print(f'C-H range: {min(bond_lengths):.4f}–{max(bond_lengths):.4f} Å')
print(f'Mean H-C-H angle: {sum(angles) / len(angles):.2f}°')
print(f'Angle range: {min(angles):.2f}–{max(angles):.2f}°')
print(f'Max force component: {fmax:.6f} eV/Å')
print('Experimental gas-phase reference: C-H = 1.087 Å; H-C-H = 109.47°')
print('Final structure: ch4_final.xyz')
