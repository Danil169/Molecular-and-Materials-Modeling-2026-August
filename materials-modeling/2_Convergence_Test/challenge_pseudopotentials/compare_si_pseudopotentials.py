#!/usr/bin/env python3
"""Challenge II.1/II.2.4: repeat Si convergence with two PBE potentials."""
import csv
import os

from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile

os.environ['OMP_NUM_THREADS'] = '1'

atoms = Atoms(
    'Si2',
    positions=[[1.35745, 4.07235, 4.07235], [0.0, 0.0, 0.0]],
    cell=[[0.0, 2.7149, 2.7149], [2.7149, 0.0, 2.7149], [2.7149, 2.7149, 0.0]],
    pbc=True,
)
potentials = {
    'ONCV (repository Si.upf)': ('../../', 'Si.upf'),
    'RRKJ (QE distribution)': ('/usr/share/espresso/pseudo', 'Si.pbe-rrkj.UPF'),
}
ecut_values = (45, 55, 65, 75)
k_values = (10, 12, 14, 15)

rows = []
for label, (pseudo_dir, pseudo) in potentials.items():
    print(f'\n=== {label} ===', flush=True)
    profile = EspressoProfile('mpirun -np 4 pw.x', pseudo_dir=pseudo_dir)

    def run(ecut, kpts):
        calc = Espresso(
            profile=profile,
            pseudopotentials={'Si': pseudo},
            input_data={
                'control': {'calculation': 'scf', 'prefix': 'si_pp', 'outdir': './tmp', 'verbosity': 'low'},
                'system': {'ecutwfc': ecut, 'ibrav': 0, 'occupations': 'smearing', 'smearing': 'gauss', 'degauss': 0.01},
                'electrons': {'conv_thr': 1.0e-10},
            },
            kpts=(kpts, kpts, kpts),
        )
        atoms.calc = calc
        return atoms.get_potential_energy()

    previous = None
    chosen_ecut = None
    for ecut in ecut_values:
        energy = run(ecut, 15)
        delta = None if previous is None else (energy - previous) * 500.0
        rows.append([label, 'ecut', ecut, '15x15x15', energy, delta])
        print(f'ecut={ecut:2d} Ry  E={energy:.6f} eV  dE={delta if delta is not None else float("nan"):.4f} meV/atom')
        if delta is not None and abs(delta) < 0.05 and chosen_ecut is None:
            chosen_ecut = ecut
        previous = energy
    chosen_ecut = chosen_ecut or ecut_values[-1]

    previous = None
    for k in k_values:
        energy = run(chosen_ecut, k)
        delta = None if previous is None else (energy - previous) * 500.0
        rows.append([label, 'k-grid', chosen_ecut, f'{k}x{k}x{k}', energy, delta])
        print(f'k={k:2d}^3 at {chosen_ecut} Ry  E={energy:.6f} eV  dE={delta if delta is not None else float("nan"):.4f} meV/atom')
        previous = energy

with open('si_pseudopotential_convergence.csv', 'w', newline='') as handle:
    writer = csv.writer(handle)
    writer.writerow(['potential', 'series', 'ecutwfc_Ry', 'k_grid', 'energy_eV', 'delta_meV_per_atom'])
    writer.writerows(rows)
print('\nSaved si_pseudopotential_convergence.csv')
