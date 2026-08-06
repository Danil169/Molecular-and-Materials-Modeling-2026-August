from ase import Atoms
from ase.calculators.emt import EMT
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

# Define bond lengths (wider range, more points)
bond_lengths = np.linspace(0.3, 4.0, 80)
energies = []

# Calculate H2 energies
for d in bond_lengths:
    h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, d]])
    h2.calc = EMT()
    energies.append(h2.get_potential_energy())

# Energy of isolated H atom
H = Atoms('H', positions=[[0, 0, 0]])
H.calc = EMT()
E_H = H.get_potential_energy()

# Normalize energies relative to separated atoms (2E_H = 0 reference)
energies_norm = np.array([E - 2*E_H for E in energies])

# Morse potential function
def morse(r, De, a, re):
    return De * (1 - np.exp(-a*(r-re)))**2 - De

# Fit Morse potential to computed data
popt, _ = curve_fit(morse, bond_lengths, energies_norm, p0=[4.0, 1.5, 0.75])
De_fit, a_fit, re_fit = popt

# Binding energy at equilibrium (positive value)
binding_energy = De_fit
eq_bond_length = re_fit

# Experimental values
exp_bond_length = 0.741  # Å
exp_binding_energy = 4.52  # eV

# Print results in a table
data = {
    "Quantity": ["Equilibrium bond length (Å)", "Binding energy (eV)", "Morse parameter a (Å⁻¹)"],
    "Computed (EMT/Morse)": [f"{eq_bond_length:.3f}", f"{binding_energy:.3f}", f"{a_fit:.3f}"],
    "Experimental": [f"{exp_bond_length:.3f}", f"{exp_binding_energy:.2f}", "—"]
}
df = pd.DataFrame(data)
print(df.to_string(index=False))

# Plot normalized curve with Morse fit
plt.plot(bond_lengths, energies_norm, 'o', label='ASE (EMT, normalized)')
plt.plot(bond_lengths, morse(bond_lengths, *popt), '-', label='Morse fit')
plt.axvline(exp_bond_length, color='r', linestyle='--', label='Exp bond length')
plt.axhline(-exp_binding_energy, color='g', linestyle='--', label='Exp binding energy')

plt.xlabel('Bond length (Å)')
plt.ylabel('Energy relative to 2H (eV)')
plt.title('H2 Potential Energy Curve (Morse Fit)')
plt.legend()
plt.show()

