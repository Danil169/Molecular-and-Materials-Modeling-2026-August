from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
import numpy as np
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Experimental values for N2 at 298.15 K and 1 atm
EXP_FREQ = 2358.57  # cm^-1 (harmonic vibrational frequency)
EXP_H = -0.524  # eV (formation enthalpy relative to atoms, approximate)
EXP_S = 0.001916  # eV/K (191.6 J/mol·K)
EXP_BOND_LENGTH = 1.0975  # Angstroms

print("=" * 70)
print("N2 MOLECULE THERMOCHEMISTRY CALCULATION (ASE with EMT)")
print("=" * 70)

# 1. Setup and Optimize
print("\n1. Structure Optimization")
print("-" * 40)
atoms = molecule('N2')
print(f"Initial bond length: {atoms.get_distance(0, 1):.4f} Å")
print(f"Experimental bond length: {EXP_BOND_LENGTH:.4f} Å")

atoms.calc = EMT()
dyn = QuasiNewton(atoms)
dyn.run(fmax=0.01)

final_bond_length = atoms.get_distance(0, 1)
print(f"Optimized bond length: {final_bond_length:.4f} Å")
print(f"Error in bond length: {(final_bond_length - EXP_BOND_LENGTH)*1000:.2f} mÅ")
print(f"Relative error: {abs(final_bond_length - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:.2f}%")

# 2. Get Potential Energy and Vibrations
print("\n2. Vibrational Analysis")
print("-" * 40)
potentialenergy = atoms.get_potential_energy()
vib = Vibrations(atoms)
vib.run()
vib_energies = vib.get_energies()  # in eV
vib_freqs_cm1 = vib.get_frequencies()  # in cm^-1

print("Calculated vibrational frequencies (cm⁻¹):")
for i, freq in enumerate(vib_freqs_cm1):
    print(f"  Mode {i+1}: {freq:.2f} cm⁻¹")

# Only keep positive frequencies (remove imaginary/zero modes)
positive_freqs = vib_freqs_cm1[vib_freqs_cm1 > 0]
if len(positive_freqs) > 0:
    calc_freq = positive_freqs[0]  # First positive frequency for linear molecule
    print(f"\nFundamental vibrational frequency: {calc_freq:.2f} cm⁻¹")
    print(f"Experimental frequency: {EXP_FREQ:.2f} cm⁻¹")
    print(f"Error: {calc_freq - EXP_FREQ:.2f} cm⁻¹ ({abs(calc_freq - EXP_FREQ)/EXP_FREQ*100:.2f}%)")
    print(f"Note: EMT is a simple effective medium theory, so deviations are expected")

# 3. Thermochemistry Calculation
print("\n3. Thermochemistry at 298.15 K, 1 atm")
print("-" * 40)

# Use vib_energies for thermochemistry (in eV)
thermo = IdealGasThermo(
    vib_energies=vib_energies,
    potentialenergy=potentialenergy,
    atoms=atoms,
    geometry='linear',
    symmetrynumber=2,
    spin=0
)

temperature = 298.15
pressure = 101325.0  # 1 atm in Pa

# Calculate thermodynamic properties
H = thermo.get_enthalpy(temperature=temperature)
S = thermo.get_entropy(temperature=temperature, pressure=pressure)
G = thermo.get_gibbs_energy(temperature=temperature, pressure=pressure)
Cv = thermo.get_heat_capacity(temperature=temperature)
Cp = Cv + 8.314462618 / 96485.33212  # Add R to convert to Cp (eV/K)

# Get individual contributions (for analysis)
H_trans = thermo.get_enthalpy(temperature=temperature, verbose=False) - potentialenergy
# We need to compute contributions separately
# For rotational enthalpy: RT for linear molecule
R_ev = 8.314462618 / 96485.33212  # R in eV/K
H_rot = R_ev * temperature  # For linear molecule
H_vib = thermo.get_enthalpy(temperature=temperature) - thermo.potentialenergy - H_trans - H_rot

print(f"Total Enthalpy (H):")
print(f"  Calculated: {H:.6f} eV")
print(f"  Experimental: {EXP_H:.6f} eV (approximate)")
print(f"  Error: {H - EXP_H:.6f} eV")

print(f"\nTotal Entropy (S):")
print(f"  Calculated: {S:.6f} eV/K")
print(f"  Experimental: {EXP_S:.6f} eV/K")
print(f"  Error: {S - EXP_S:.6f} eV/K ({abs(S - EXP_S)/EXP_S*100:.2f}%)")

print(f"\nGibbs Free Energy (G): {G:.6f} eV")
print(f"Heat Capacity (Cv): {Cv:.6f} eV/K")
print(f"Heat Capacity (Cp): {Cp:.6f} eV/K")

# 4. Comparison Summary
print("\n" + "=" * 70)
print("SUMMARY OF COMPARISON WITH EXPERIMENT")
print("=" * 70)

properties = [
    ("Bond Length", final_bond_length, EXP_BOND_LENGTH, "Å"),
    ("Vibrational Frequency", calc_freq, EXP_FREQ, "cm⁻¹"),
    ("Entropy", S, EXP_S, "eV/K")
]

print("\n{:<20} {:<15} {:<15} {:<15} {:<10}".format(
    "Property", "Calculated", "Experimental", "Error", "Error %"
))
print("-" * 75)

for name, calc, exp, unit in properties:
    error = calc - exp
    error_pct = abs(error)/exp * 100
    print("{:<20} {:<15.4f} {:<15.4f} {:<+15.4f} {:<9.2f}%".format(
        name, calc, exp, error, error_pct
    ))

print("\nDISCUSSION:")
print("-" * 40)
print("• EMT is a fast but approximate method, not designed for high accuracy")
print("• For quantitative results, use DFT (e.g., GPAW, VASP) or quantum chemistry")
print("• The N-N triple bond is stiff, making vibrational frequency sensitive to the method")
print("• Zero-point energy contribution is not explicitly included in this simple calculation")
print("• Experimental entropy includes all degrees of freedom (translational, rotational, vibrational)")

# 5. Additional Analysis - Zero-point energy
print("\n5. Zero-Point Energy (ZPE)")
print("-" * 40)
zpe = 0.5 * sum(vib_energies[vib_energies > 0])  # Sum over positive frequencies only
print(f"Zero-point energy: {zpe:.6f} eV ({zpe*96485.33212:.2f} kJ/mol)")
print(f"ZPE as fraction of bond energy: {zpe/abs(potentialenergy)*100:.2f}%")
print(f"Note: The true ZPE for N2 is approximately 0.145 eV (14.0 kJ/mol)")

print("\n" + "=" * 70)
