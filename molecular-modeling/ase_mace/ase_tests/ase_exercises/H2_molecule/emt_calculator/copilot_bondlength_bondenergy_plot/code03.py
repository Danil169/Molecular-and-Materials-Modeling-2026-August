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

# Physical constants
eV_to_J = 1.602e-19
amu_to_kg = 1.6605e-27
c = 2.998e10  # cm/s

# Reduced mass of H2 (amu -> kg)
m_H = 1.00784 * amu_to_kg
mu = m_H / 2

# Convert De to Joules
De_J = De_fit * eV_to_J

# Convert a from Å⁻¹ to m⁻¹
a_m = a_fit * 1e10

# Vibrational frequency (Hz)
omega = a_m / (2*np.pi) * np.sqrt(2*De_J/mu)

# Convert to wavenumber (cm⁻¹)
omega_cm = omega / c

# Experimental values
exp_bond_length = 0.741  # Å
exp_binding_energy = 4.52  # eV
exp_vibrational_freq = 4401  # cm^-1

# Print results in a table
data = {
    "Quantity": ["Equilibrium bond length (Å)", "Binding energy (eV)", "Morse parameter a (Å⁻¹)", "Vibrational frequency (cm⁻¹)"],
    "Computed (EMT/Morse)": [f"{re_fit:.3f}", f"{De_fit:.3f}", f"{a_fit:.3f}", f"{omega_cm:.0f}"],
    "Experimental": [f"{exp_bond_length:.3f}", f"{exp_binding_energy:.2f}", "—", f"{exp_vibrational_freq}"]
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

