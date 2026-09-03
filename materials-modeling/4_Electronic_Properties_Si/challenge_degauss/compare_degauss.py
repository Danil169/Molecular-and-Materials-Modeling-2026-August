#!/usr/bin/env python3
"""Challenge II.4: quantify Gaussian-smearing effects for Si DOS."""
import csv
import os
import subprocess

import numpy as np
from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.calculators.calculator import PropertyNotImplementedError

os.environ['OMP_NUM_THREADS'] = '1'
atoms = Atoms(
    'Si2', positions=[[1.3669840666, 4.1009521997, 4.1009521997], [0., 0., 0.]],
    cell=[[0., 2.7339681331, 2.7339681331], [2.7339681331, 0., 2.7339681331],
          [2.7339681331, 2.7339681331, 0.]], pbc=True,
)
profile = EspressoProfile('mpirun -np 4 pw.x', pseudo_dir='../..')
summary = []
for degauss in (0.005, 0.010, 0.020):
    tag = f'dg{int(degauss * 1000):03d}'
    common = {
        'control': {'prefix': tag, 'outdir': './tmp', 'wf_collect': True, 'verbosity': 'low'},
        'system': {'ecutwfc': 65, 'ibrav': 0, 'occupations': 'smearing', 'smearing': 'gauss', 'degauss': degauss},
        'electrons': {'conv_thr': 1.0e-10},
    }
    scf_data = {key: dict(value) for key, value in common.items()}
    scf_data['control']['calculation'] = 'scf'
    atoms.calc = Espresso(profile=profile, pseudopotentials={'Si': 'Si.upf'}, input_data=scf_data, kpts=(15, 15, 15))
    energy = atoms.get_potential_energy()
    fermi = atoms.calc.get_fermi_level()

    nscf_data = {key: dict(value) for key, value in common.items()}
    nscf_data['control']['calculation'] = 'nscf'
    nscf_data['system']['nosym'] = True
    atoms.calc = Espresso(profile=profile, pseudopotentials={'Si': 'Si.upf'}, input_data=nscf_data, kpts=(24, 24, 24))
    # ASE executes QE when an energy is requested; QE successfully completes
    # the NSCF wavefunction calculation but does not expose a variational
    # energy for this calculation type.  Catch that expected post-execution
    # exception and continue to ``dos.x``.
    try:
        atoms.get_potential_energy()
    except PropertyNotImplementedError:
        pass

    input_name = f'dos_{tag}.in'
    with open(input_name, 'w') as handle:
        handle.write(f"&DOS\n prefix = '{tag}',\n outdir = './tmp',\n fildos = 'dos_{tag}.dat',\n Emin = -15.0, Emax = 15.0, DeltaE = 0.02,\n ngauss = 0, degauss = {degauss},\n/\n")
    with open(f'dos_{tag}.out', 'w') as handle:
        subprocess.run(['dos.x', '-in', input_name], stdout=handle, stderr=subprocess.STDOUT, check=True)
    data = np.loadtxt(f'dos_{tag}.dat')
    dos_at_fermi = data[np.argmin(abs(data[:, 0] - fermi)), 1]
    summary.append([degauss, energy, fermi, dos_at_fermi])
    print(f'degauss={degauss:.3f} Ry  E={energy:.6f} eV  E_F={fermi:.6f} eV  DOS(E_F)={dos_at_fermi:.6f}')

with open('degauss_summary.csv', 'w', newline='') as handle:
    writer = csv.writer(handle)
    writer.writerow(['degauss_Ry', 'scf_energy_eV', 'fermi_eV', 'DOS_at_fermi_states_per_eV'])
    writer.writerows(summary)
print('Saved degauss_summary.csv and three DOS files.')
