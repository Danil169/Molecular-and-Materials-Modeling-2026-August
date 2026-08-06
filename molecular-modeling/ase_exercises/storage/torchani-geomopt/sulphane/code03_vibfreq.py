import torch
import torchani
from ase import Atoms
from ase.optimize import BFGS
from ase.vibrations import Vibrations
import os
import shutil

# --- Setup Molecule ---
# H2S initial coordinates
atoms = Atoms('SH2', positions=[[0.0, 0.0, 0.0], 
                                [0.0, 0.9, 1.0], 
                                [0.0, -0.9, 1.0]])

# Use ANI-2x which supports Sulfur
calculator = torchani.models.ANI2x().ase()
atoms.calc = calculator

# --- 1. High-precision Optimization ---
# Required for accurate frequencies (Hessian calculation)
opt = BFGS(atoms, logfile=None)
opt.run(fmax=0.0001)

# --- 2. Frequency Calculation ---
# Reference for Experimental H2S (Gas Phase):
# [1] Shimanouchi, T., Tables of Molecular Vibrational Frequencies, 
#     Consolidated Volume I, NSRDS-NBS 39, 1972.
# [2] NIST Chemistry WebBook (https://nist.gov)
exp_freqs = [1183.0, 2615.0, 2626.0]  # [v2, v1, v3]
mode_names = ["Bending (v2)", "Symm. Stretch (v1)", "Asymm. Stretch (v3)"]

vib = Vibrations(atoms, name='h2s_vib')
vib.run()
all_freqs = vib.get_frequencies() 

# Get the 3 highest frequencies (vibrational modes)
ani_vibs = sorted([f.real for f in all_freqs])[-3:]

# --- 3. Print Results & Relative Differences ---
print("-" * 75)
print(f"{'Mode':<20} | {'ANI-2x':<12} | {'Exp.':<10} | {'Relat. Diff. (%)'}")
print("-" * 75)

for i in range(3):
    calc = ani_vibs[i]
    ref = exp_freqs[i]
    rel_diff = abs(calc - ref) / ref * 100
    
    print(f"{mode_names[i]:<20} | {calc:<12.1f} | {ref:<10.1f} | {rel_diff:.2f}%")

print("-" * 75)

# --- Cleanup ---
vib.clean()
if os.path.exists('h2s_vib'):
    shutil.rmtree('h2s_vib')

