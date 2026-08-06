import torch
import torchani
from ase import Atoms
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.units import mol, kJ, J, kB
import numpy as np
import os
import shutil

# --- 1. Setup & Optimize ---
atoms = Atoms('SH2', positions=[[0.0, 0.0, 0.0], [0.0, 0.9, 1.0], [0.0, -0.9, 1.0]])
atoms.calc = torchani.models.ANI2x().ase()
BFGS(atoms, logfile=None).run(fmax=0.0001)

# --- 2. Vibrations ---
vib = Vibrations(atoms, name='h2s_vib')
vib.run()
# Get energies and filter for real vibrational modes (highest 3)
vib_energies = sorted([e.real for e in vib.get_energies()])[-3:]

# --- 3. Thermochemistry Calculation ---
temp = 298.15
pres = 101325 
thermo = IdealGasThermo(vib_energies=vib_energies, geometry='nonlinear',
                        atoms=atoms, symmetrynumber=2, spin=0)

# A. Standard Properties (Supported by almost all ASE versions)
zpe_ev = sum(vib_energies) / 2.0
s_ev = thermo.get_entropy(temperature=temp, pressure=pres, verbose=False)

# B. Manual Heat Capacity Calculation (to bypass missing ASE methods)
# Cv = Cv_trans + Cv_rot + Cv_vib
# For a non-linear molecule: Cv_trans = 3/2 k, Cv_rot = 3/2 k
cv_trans_rot = 3.0 * kB 

# Cv_vib = k * sum( (hv/kT)^2 * exp(hv/kT) / (exp(hv/kT) - 1)^2 )
cv_vib = 0.0
for hv in vib_energies:
    x = hv / (kB * temp)
    cv_vib += kB * (x**2 * np.exp(x)) / (np.exp(x) - 1)**2

cv_ev = cv_trans_rot + cv_vib
cp_ev = cv_ev + kB  # Cp = Cv + R

# --- 4. Units & Comparison ---
s_std = s_ev * mol / J         
cp_std = cp_ev * mol / J       
zpe_std = zpe_ev * mol / kJ    

exp_s, exp_cp, exp_zpe = 205.81, 34.23, 41.71

# --- 5. Print Comparison ---
print("\n" + "="*70)
print(f" THERMOCHEMICAL COMPARISON (T={temp}K, P=1.0 atm)")
print("="*70)
print(f"{'Property':<25} | {'ANI-2x':<12} | {'Exp.':<12} | {'% Error'}")
print("-" * 70)

def print_row(name, calc, ref, unit):
    err = abs(calc - ref) / ref * 100
    print(f"{name:<25} | {calc:>8.2f} {unit:<3} | {ref:>8.2f} {unit:<3} | {err:>6.2f}%")

print_row("Entropy (S)", s_std, exp_s, "J/K")
print_row("Heat Capacity (Cp)", cp_std, exp_cp, "J/K")
print_row("Zero Point Energy (ZPE)", zpe_std, exp_zpe, "kJ ")
print("="*70)

# Cleanup
vib.clean()
if os.path.exists('h2s_vib'): shutil.rmtree('h2s_vib')

