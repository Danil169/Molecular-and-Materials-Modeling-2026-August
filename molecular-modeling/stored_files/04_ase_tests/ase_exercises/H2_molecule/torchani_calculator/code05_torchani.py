import sys
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms

# Ensure TorchANI is installed before proceeding
try:
    import torchani
except ImportError:
    print("\n" + "="*60)
    print("ERROR: 'torchani' library is not installed.")
    print("To run this script, please install TorchANI and its dependencies:")
    print("    pip install torch torchani")
    print("="*60 + "\n")
    sys.exit(1)

# 1. Calculate the energy of a single isolated Hydrogen atom
isolated_atom_ani = Atoms('H', positions=[(0, 0, 0)])
isolated_atom_ani.calc = torchani.models.ANI2x().ase()
e_single_atom_ani = isolated_atom_ani.get_potential_energy()
e_two_isolated_atoms_ani = 2 * e_single_atom_ani

# 2. Define the range of interatomic distances (in Angstroms)
distances = np.linspace(0.4, 2.5, 50)
energies_ani = []

# 3. Loop over distances to calculate molecular potential energy
for r in distances:
    molecule = Atoms('H2', positions=[(0, 0, 0), (0, 0, r)])
    
    # TorchANI Calculation
    molecule.calc = torchani.models.ANI2x().ase()
    energies_ani.append(molecule.get_potential_energy())

energies_ani = np.array(energies_ani)

# 4. Find the equilibrium bond length and minimum molecular energy
min_idx_ani = np.argmin(energies_ani)
r_eq_ani = distances[min_idx_ani]
e_min_ani = energies_ani[min_idx_ani]

# 5. Compute Simulated Dissociation Energy (De)
simulated_de_ani = e_two_isolated_atoms_ani - e_min_ani
experimental_de = 4.52  # Precise experimental De for H2 in eV

# 6. Save data points and metrics to a text file
output_filename = "potential_curve_with_dissociation.txt"
data_to_save = np.column_stack((distances, energies_ani))

header_text = (
    f"H2 Potential Energy Curve and Dissociation Data (TorchANI)\n"
    f"TorchANI Equilibrium Bond Length (r_e): {r_eq_ani:.3f} A\n"
    f"TorchANI Simulated Dissociation Energy (De): {simulated_de_ani:.3f} eV\n"
    f"Experimental Dissociation Energy (De): {experimental_de:.2f} eV\n"
    f"{'Distance(A)':<15} {'Energy_ANI(eV)':<15}"
)

np.savetxt(output_filename, data_to_save, fmt="%-15.6f", header=header_text)
print(f"Data saved to {output_filename}")

# 7. Print the comparison metrics
print("\n--- Summary of Results (TorchANI ANI-2x) ---")
print(f"Equilibrium Bond Length: {r_eq_ani:.3f} Å")
print(f"Simulated Dissociation Energy (De): {simulated_de_ani:.3f} eV")
print(f"Experimental Dissociation Energy (De): {experimental_de:.3f} eV")
print(f"Absolute Error: {abs(simulated_de_ani - experimental_de):.3f} eV")

# 8. Plot the potential energy curve
plt.figure(figsize=(8, 6))

# Plot TorchANI Data
plt.plot(distances, energies_ani, 'g-s', markersize=4, lw=1.5, label='TorchANI Data Points')
plt.axvline(r_eq_ani, color='g', linestyle='--', alpha=0.5, label=f'r_e (ANI) = {r_eq_ani:.2f} Å')
plt.axhline(e_two_isolated_atoms_ani, color='g', linestyle=':', alpha=0.5, label='2x Isolated H (ANI)')

plt.title(r'Potential Energy Curve & Dissociation Energy of $\mathrm{H_2}$ (TorchANI)')
plt.xlabel(r'Interatomic Distance $r$ (Å)')
plt.ylabel('Potential Energy $E$ (eV)')
plt.legend()
plt.grid(True, linestyle=':')
plt.show()

