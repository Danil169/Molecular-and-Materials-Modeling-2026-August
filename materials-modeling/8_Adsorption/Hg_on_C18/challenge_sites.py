#!/usr/bin/env python
"""Challenge II.8: adsorption of Hg on bridge and hollow sites of C18.

Top-site result comes from the main exercise (relaxation_Hg_C18.py).
Energies of the isolated slab and atom are site-independent, so the
values from energy_c18.py and energy_Hg.py are reused.
"""
import os
import sys

from ase import Atom
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS

os.environ['OMP_NUM_THREADS'] = '1'

E_C18 = -2947.320598  # eV, energy_c18.py
E_HG = -4533.787985   # eV, energy_Hg.py
E_HGC18_TOP = -7481.346431  # eV, relaxation_Hg_C18.py

pw_profile = EspressoProfile(command='mpirun -np 4 pw.x', pseudo_dir='./')

input_data = {
    'control': {
        'outdir': './tmp',
        'verbosity': 'low',
        'tstress': True,
        'tprnfor': True
    },
    'system': {
        'ecutwfc': 80,
        'occupations': 'smearing',
        'smearing': 'gauss',
        'degauss': 0.01,
        'ibrav': 0,
        'nat': 19,
        'ntyp': 2,
        'vdw_corr': 'DFT-D3',
        'dftd3_version': 4
    },
    'electrons': {
        'conv_thr': 1.0e-10
    }
}

def make_calculator(prefix):
    data = dict(input_data)
    data['control'] = dict(input_data['control'], prefix=prefix)
    return Espresso(
        profile=pw_profile,
        pseudopotentials={'C': 'C.upf', 'Hg': 'Hg.upf'},
        input_data=data,
        kpts=(1, 1, 1)
    )

slab = read('c18.vasp', format='vasp')

# High-symmetry adsorption sites above the sheet (Cartesian x, y in A)
sites = {
    'bridge': (3.084655006, 1.780926402),  # midpoint of a C-C bond
    'hollow': (2.467724005, 2.849482241),  # center of a carbon hexagon
}
Z_START = 10.9  # A above the sheet, as in c18_hg.vasp

results = {'top': E_HGC18_TOP}

for name, (x, y) in sites.items():
    print(f"\n=== {name} site ===", flush=True)
    atoms = slab.copy()
    atoms.append(Atom('Hg', (x, y, Z_START)))

    d_min = min(((atoms.positions[:-1] - [x, y, Z_START])**2)
                .sum(axis=1)**0.5)
    print(f"Initial Hg-nearest-C distance: {d_min:.3f} A", flush=True)

    vasp_name = f'c18_hg_{name}.vasp'
    write(vasp_name, atoms, format='vasp', direct=False)
    print(f"Initial structure saved to {vasp_name}", flush=True)

    atoms.calc = make_calculator(f'chg_{name}')
    traj = Trajectory(f'relaxation_{name}.traj', 'w', atoms)
    opt = BFGS(atoms, trajectory=traj, logfile=f'relaxation_{name}.log')
    try:
        opt.run(fmax=0.01)
        e = atoms.get_potential_energy()
        results[name] = e
        write(f'final_relaxed_structure_{name}.vasp', atoms,
              format='vasp', direct=False)
        print(f"E(Hg@C18, {name}) = {e:.6f} eV", flush=True)
    except Exception as exc:
        print(f"Relaxation failed for {name}: {exc}", flush=True)
    finally:
        traj.close()

print("\n=== Adsorption energies ===", flush=True)
print(f"E_ads = E(Hg@C18) - E(C18) - E(Hg)")
for name, e in results.items():
    e_ads = e - E_C18 - E_HG
    print(f"  {name:>7}: E(Hg@C18) = {e:>13.6f} eV   E_ads = {e_ads:+.4f} eV",
          flush=True)
