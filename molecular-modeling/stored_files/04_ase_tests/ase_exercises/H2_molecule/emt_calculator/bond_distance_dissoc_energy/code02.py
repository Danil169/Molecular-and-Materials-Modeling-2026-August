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

# 3. Find the equilibrium bond length and minimum energy
min_idx = np.argmin(energies)
r_eq = distances[min_idx]
e_min = energies[min_idx]

print(f"Equilibrium Bond Length: {r_eq:.3f} Å")
print(f"Minimum Energy: {e_min:.3f} eV")

# 4. Plot the potential energy curve with explicit points
plt.figure(figsize=(7, 5))
# Added 'o-' format to draw both lines and scatter marker points
plt.plot(distances, energies, 'b-o', markersize=4, lw=1.5, label='Calculated Points')
plt.axvline(r_eq, color='r', linestyle='--', label=f'r_e = {r_eq:.2f} Å')
plt.axhline(e_min, color='g', linestyle='--', label=f'E_e = {e_min:.2f} eV')
plt.title(r'Potential Energy Curve of $\mathrm{H_2}$ Molecule')
plt.xlabel(r'Interatomic Distance $r$ (Å)')
plt.ylabel('Potential Energy $E$ (eV)')
plt.legend()
plt.grid(True, linestyle=':')
plt.show()

