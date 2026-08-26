from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
import numpy as np
import warnings

warnings.filterwarnings('ignore')

# Experimental values
EXP_FREQ = 2358.57  # cm^-1
EXP_S = 0.001916  # eV/K
EXP_BOND_LENGTH = 1.0975  # Å
EXP_ZPE = 0.145  # eV

print("=" * 70)
print("N2 MOLECULE THERMOCHEMISTRY CALCULATION (ASE with EMT)")
print("=" * 70)

# 1. Structure Optimization
print("\n1. Structure Optimization")
print("-" * 40)
atoms = molecule('N2')
print(f"Initial bond length: {atoms.get_distance(0, 1):.4f} Å")

atoms.calc = EMT()
dyn = QuasiNewton(atoms)
dyn.run(fmax=0.01)

final_bond_length = atoms.get_distance(0, 1)
print(f"Optimized bond length: {final_bond_length:.4f} Å")
print(f"Error: {(final_bond_length - EXP_BOND_LENGTH)*1000:.2f} mÅ")

# 2. Vibrational Analysis
print("\n2. Vibrational Analysis")
print("-" * 40)
potentialenergy = atoms.get_potential_energy()
print(f"Potential energy: {potentialenergy:.6f} eV")

vib = Vibrations(atoms)
vib.run()
vib_energies = vib.get_energies()
vib_freqs_cm1 = vib.get_frequencies()

print("Calculated vibrational frequencies (cm⁻¹):")
positive_modes = []
for i, freq in enumerate(vib_freqs_cm1):
    freq_real = np.real(freq)
    if abs(freq_real) > 1.0:
        positive_modes.append(freq_real)
        print(f"  Mode {i+1}: {freq_real:.2f} cm⁻¹")
    else:
        print(f"  Mode {i+1}: {freq_real:.2f} cm⁻¹ (translational/rotational)")

if positive_modes:
    calc_freq = max(positive_modes)
    print(f"\nN-N stretching frequency: {calc_freq:.2f} cm⁻¹")
    print(f"Experimental: {EXP_FREQ:.2f} cm⁻¹")
    print(f"Error: {calc_freq - EXP_FREQ:.2f} cm⁻¹ ({abs(calc_freq - EXP_FREQ)/EXP_FREQ*100:.2f}%)")

# 3. Thermochemistry
print("\n3. Thermochemistry at 298.15 K, 1 atm")
print("-" * 40)

# Filter vib_energies for thermochemistry
vib_energies_filtered = np.array([np.real(e) for e in vib_energies if np.real(e) > 0.001])

# Manual ZPE calculation (since it's not a separate method)
zpe = 0.5 * np.sum(vib_energies_filtered)
print(f"Zero-Point Energy (ZPE): {zpe:.6f} eV")
print(f"Experimental ZPE: {EXP_ZPE:.4f} eV")
print(f"ZPE Error: {zpe - EXP_ZPE:.4f} eV ({abs(zpe - EXP_ZPE)/EXP_ZPE*100:.2f}%)")

# Create thermo object
thermo = IdealGasThermo(
    vib_energies=vib_energies_filtered,
    potentialenergy=potentialenergy,
    atoms=atoms,
    geometry='linear',
    symmetrynumber=2,
    spin=0
)

temperature = 298.15
pressure = 101325.0

# Calculate properties silently
H = thermo.get_enthalpy(temperature=temperature, verbose=False)
S = thermo.get_entropy(temperature=temperature, pressure=pressure, verbose=False)
G = thermo.get_gibbs_energy(temperature=temperature, pressure=pressure, verbose=False)
U = thermo.get_internal_energy(temperature=temperature, verbose=False)

print(f"Internal Energy (U): {U:.6f} eV")
print(f"Enthalpy (H): {H:.6f} eV")
print(f"Entropy (S): {S:.6f} eV/K")
print(f"Gibbs Free Energy (G): {G:.6f} eV")

# Manual heat capacity calculation (since it's not a method)
def calculate_Cv(thermo, T, deltaT=0.1):
    U1 = thermo.get_internal_energy(temperature=T - deltaT/2, verbose=False)
    U2 = thermo.get_internal_energy(temperature=T + deltaT/2, verbose=False)
    return (U2 - U1) / deltaT

Cv = calculate_Cv(thermo, temperature)
kB = 8.617333262145e-5  # Boltzmann constant in eV/K
Cp = Cv + kB

print(f"Heat Capacity (Cv): {Cv:.6f} eV/K")
print(f"Heat Capacity (Cp): {Cp:.6f} eV/K")

# Print final summary (this will show the formatted output once)
print("\n" + "=" * 70)
print("FINAL THERMODYNAMIC SUMMARY")
print("=" * 70)
thermo.get_enthalpy(temperature=temperature, verbose=True)

# 4. Comparison Summary
print("\n" + "=" * 70)
print("COMPARISON WITH EXPERIMENT")
print("=" * 70)

properties = [
    ("Bond Length", final_bond_length, EXP_BOND_LENGTH, "Å"),
    ("Vibrational Freq", calc_freq, EXP_FREQ, "cm⁻¹"),
    ("Entropy", S, EXP_S, "eV/K"),
    ("ZPE", zpe, EXP_ZPE, "eV")
]

print("\n{:<20} {:<15} {:<15} {:<15} {:<10}".format(
    "Property", "Calculated", "Experimental", "Error", "Error %"
))
print("-" * 75)

for name, calc, exp, unit in properties:
    error = calc - exp
    error_pct = abs(error)/exp * 100 if exp != 0 else 0
    print("{:<20} {:<15.4f} {:<15.4f} {:<+15.4f} {:<9.2f}%".format(
        name, calc, exp, error, error_pct
    ))

print("\n" + "=" * 70)
