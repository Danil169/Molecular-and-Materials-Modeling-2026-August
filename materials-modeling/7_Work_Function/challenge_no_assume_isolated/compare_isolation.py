#!/usr/bin/env python3
"""Challenge II.7: compare graphene potentials with and without 2D isolation."""
import os
import subprocess

import numpy as np
from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile

os.environ['OMP_NUM_THREADS'] = '1'
atoms = Atoms(
    'C2', positions=[[0., 0., 7.5], [0., 1.4239328281, 7.5]],
    cell=[[2.466324005, 0., 0.], [-1.2331620025, 2.1358992423, 0.], [0., 0., 15.]],
    pbc=True,
)
profile = EspressoProfile('mpirun -np 4 pw.x', pseudo_dir='../..')
summary = []
for tag, isolated in (('2d', True), ('none', False)):
    system = {'ecutwfc': 100, 'ibrav': 0, 'occupations': 'smearing', 'smearing': 'gauss', 'degauss': 0.01}
    if isolated:
        system['assume_isolated'] = '2D'
    prefix = f'graphene_{tag}'
    atoms.calc = Espresso(
        profile=profile, pseudopotentials={'C': 'C.upf'},
        input_data={'control': {'calculation': 'scf', 'prefix': prefix, 'outdir': './tmp', 'verbosity': 'low',
                                'tstress': True, 'tprnfor': True},
                    'system': system, 'electrons': {'conv_thr': 1.0e-10}},
        kpts=(27, 27, 1),
    )
    energy = atoms.get_potential_energy()
    fermi = atoms.calc.get_fermi_level()
    pp_in = f"""&INPUTPP\n prefix = '{prefix}',\n outdir = './tmp',\n filplot = 'potential_{tag}',\n plot_num = 11,\n/\n&PLOT\n iflag = 3,\n output_format = 6,\n/\n"""
    with open(f'pp_{tag}.in', 'w') as handle:
        handle.write(pp_in)
    with open(f'pp_{tag}.out', 'w') as handle:
        subprocess.run(['pp.x', '-in', f'pp_{tag}.in'], stdout=handle, stderr=subprocess.STDOUT, check=True)
    with open(f'average_{tag}.in', 'w') as handle:
        handle.write(f'1\npotential_{tag}\n1.0\n200\n3\n3.8149\n')
    with open(f'average_{tag}.out', 'w') as handle:
        subprocess.run(['average.x'], stdin=open(f'average_{tag}.in'), stdout=handle, stderr=subprocess.STDOUT, check=True)
    os.replace('avg.dat', f'avg_{tag}.dat')
    data = np.loadtxt(f'avg_{tag}.dat')
    macroscopic = data[:, 2] * 13.605693
    vacuum_slice = macroscopic[(data[:, 0] > 5.0) & (data[:, 0] < 9.0)]
    vacuum_ev = float(np.mean(vacuum_slice))
    summary.append((tag, energy, fermi, vacuum_ev, vacuum_ev - fermi))
    print(f'{tag}: E={energy:.6f} eV; E_F={fermi:.6f} eV; vacuum plateau={vacuum_ev:.6f} eV; W={vacuum_ev-fermi:.6f} eV')

with open('isolation_summary.csv', 'w') as handle:
    handle.write('mode,total_energy_eV,fermi_eV,vacuum_plateau_eV,work_function_eV\n')
    for row in summary:
        handle.write(','.join([row[0], *(f'{value:.8f}' for value in row[1:])]) + '\n')
print('Saved avg_2d.dat, avg_none.dat and isolation_summary.csv')
