import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.calculators.emt import EMT

# 1. Calculate the energy of a single isolated Hydrogen atom
isolated_atom = Atoms('H', positions=[(0, 0, 0)])
isolated_atom.calc = EMT()
e_single_atom = isolated_atom.get_potential_energy()
e_two_isolated_atoms = 2 * e_single_atom

# 2. Define the range of interatomic distances (in Angstroms)
distances = np.linspace(0.4, 2.5, 50)
energies = []

# 3. Loop over distances to calculate molecular potential energy
for r in distances:
    molecule = Atoms('H2', positions=[(0, 0, 0), (0, 0, r)])
    molecule.calc = EMT()
    energies.append(molecule.get_potential_energy())

energies = np.array(energies)

# 4. Find the equilibrium bond length and minimum molecular energy
min_idx = np.argmin(energies)
r_eq = distances[min_idx]
e_min = energies[min_idx]

# 5. Compute Simulated Dissociation Energy (De)
# De = E(separated atoms) - E(equilibrium molecule)
simulated_de = e_two_isolated_atoms - e_min
experimental_de = 4.52  # Precise experimental De for H2 in eV

# 6. Save data points and metrics to a text file
output_filename = "potential_curve_with_dissociation.txt"
data_to_save = np.column_stack((distances, energies))

header_text = (
    f"H2 Potential Energy Curve and Dissociation Data (ASE/EMT)\n"
    f"Equilibrium Bond Length (r_e): {r_eq:.3f} A\n"
    f"Simulated Dissociation Energy (De): {simulated_de:.3f} eV\n"
    f"Experimental Dissociation Energy (De): {experimental_de:.2f} eV\n"
    f"{'Distance(A)':<15} {'Energy(eV)':<15}"
)
np.savetxt(output_filename, data_to_save, fmt="%-15.6f", header=header_text)
print(f"Data saved to {output_filename}")

# 7. Print the comparison metrics
print("\n--- Summary of Results ---")
print(f"Equilibrium Bond Length: {r_eq:.3f} Å")
print(f"Simulated Dissociation Energy (De): {simulated_de:.3f} eV")
print(f"Experimental Dissociation Energy (De): {experimental_de:.3f} eV")
print(f"Absolute Error: {abs(simulated_de - experimental_de):.3f} eV")

# 8. Plot the potential energy curve with discrete simulation steps
plt.figure(figsize=(7, 5))
plt.plot(distances, energies, 'b-o', markersize=4, lw=1.5, label='ASE Data Points')
plt.axvline(r_eq, color='r', linestyle='--', label=f'r_e = {r_eq:.2f} Å')
plt.axhline(e_min, color='g', linestyle='--', label=f'E_min = {e_min:.2f} eV')
# Reference line indicating separated atom energy limit
plt.axhline(e_two_isolated_atoms, color='m', linestyle=':', label='2x Isolated H Atoms')

plt.title(r'Potential Energy Curve & Dissociation Energy of $\mathrm{H_2}$')
plt.xlabel(r'Interatomic Distance $r$ (Å)')
plt.ylabel('Potential Energy $E$ (eV)')
plt.legend()
plt.grid(True, linestyle=':')
plt.show()

