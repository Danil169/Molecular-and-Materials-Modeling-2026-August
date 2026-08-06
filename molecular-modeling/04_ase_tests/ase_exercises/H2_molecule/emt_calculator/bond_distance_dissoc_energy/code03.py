import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.calculators.emt import EMT

# 1. Define the range of interatomic distances (in Angstroms)
distances = np.linspace(0.4, 2.5, 50)
energies = []

# 2. Loop over distances to calculate potential energy
for r in distances:
    # Create H2 molecule aligned along the Z-axis
    molecule = Atoms('H2', positions=[(0, 0, 0), (0, 0, r)])
    
    # Attach the EMT calculator
    molecule.calc = EMT()
    
    # Calculate and store potential energy (eV)
    energies.append(molecule.get_potential_energy())

# Convert lists to a numpy array for file operations
energies = np.array(energies)

# 3. Find the equilibrium bond length and minimum energy
min_idx = np.argmin(energies)
r_eq = distances[min_idx]
e_min = energies[min_idx]

# 4. Write data points to a text file
output_filename = "potential_curve_data.txt"

# Stack data into two columns: Distance (A) and Energy (eV)
data_to_save = np.column_stack((distances, energies))

# Save with a descriptive file header
header_text = (
    f"H2 Potential Energy Curve Data (ASE/EMT)\n"
    f"Equilibrium Bond Length: {r_eq:.3f} Angstroms\n"
    f"Minimum Energy: {e_min:.3f} eV\n"
    f"{'Distance(A)':<15} {'Energy(eV)':<15}"
)

np.savetxt(output_filename, data_to_save, fmt="%-15.6f", header=header_text)
print(f"Data successfully saved to {output_filename}")

# 5. Plot the potential energy curve with explicit points
plt.figure(figsize=(7, 5))
plt.plot(distances, energies, 'b-o', markersize=4, lw=1.5, label='Calculated Points')
plt.axvline(r_eq, color='r', linestyle='--', label=f'r_e = {r_eq:.2f} Å')
plt.axhline(e_min, color='g', linestyle='--', label=f'E_e = {e_min:.2f} eV')
plt.title(r'Potential Energy Curve of $\mathrm{H_2}$ Molecule')
plt.xlabel(r'Interatomic Distance $r$ (Å)')
plt.ylabel('Potential Energy $E$ (eV)')
plt.legend()
plt.grid(True, linestyle=':')
plt.show()

