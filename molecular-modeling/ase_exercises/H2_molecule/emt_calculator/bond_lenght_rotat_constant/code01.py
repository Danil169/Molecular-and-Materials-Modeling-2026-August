import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
import ase.units as units

# 1. Accepted Experimental Reference Values (for H2 Ground State)
R_E_EXP = 0.7416  # Angstroms (Å)
B_E_EXP = 60.853  # Wavenumbers (cm^-1)

# 2. Construct and Optimize H2 Molecule with ASE
h2 = Atoms('H2', positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]])
h2.calc = EMT()

# Run geometry optimization
opt = BFGS(h2, logfile=None)
opt.run(fmax=0.01)

# Extract optimized distance
r_e_calc = h2.get_distance(0, 1)

# 3. Call get_moments_of_inertia() directly on the Atoms object
moments = h2.get_moments_of_inertia()
I_amu_ang2 = np.max(moments)  # Pick the non-zero principal moment

# 4. Standard Unit Conversion to cm^-1 (using explicit values)
# Convert Moment of Inertia from (amu * Å^2) to (kg * m^2)
I_kg_m2 = I_amu_ang2 * units._amu * (1e-10)**2

# Define Planck constant (h) and Speed of light (c) in SI units
h_si = 6.62607015e-34  # J * s
c_cm_s = 29979245800.0  # cm / s

# Calculate B_e using the spectroscopic formula: h / (8 * pi^2 * c * I)
b_e_calc = h_si / (8 * np.pi**2 * c_cm_s * I_kg_m2)

# 5. Compute Error Percentages
r_error = ((r_e_calc - R_E_EXP) / R_E_EXP) * 100
b_error = ((b_e_calc - B_E_EXP) / B_E_EXP) * 100

# 6. Format and Print Comparison Summary
print("==================================================")
print("     H2 ROTATIONAL METRICS COMPARISON SUMMARY     ")
print("==================================================")
print(f"Equilibrium Bond Length (r_e):")
print(f"  - Calculated (EMT):   {r_e_calc:.4f} Å")
print(f"  - Experimental:       {R_E_EXP:.4f} Å")
print(f"  - Relative Error:     {r_error:+.2f}%")
print("--------------------------------------------------")
print(f"Rotational Constant (B_e):")
print(f"  - Calculated (EMT):   {b_e_calc:.3f} cm^-1")
print(f"  - Experimental:       {B_E_EXP:.3f} cm^-1")
print(f"  - Relative Error:     {b_error:+.2f}%")
print("==================================================")

