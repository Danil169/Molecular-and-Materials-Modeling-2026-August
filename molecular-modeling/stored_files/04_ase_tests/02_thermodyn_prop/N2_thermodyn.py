from ase.build import molecule
from ase.calculators.emt import EMT # Example calculator
from ase.optimize import QuasiNewton
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations

# 1. Setup and Optimize
atoms = molecule('N2')
atoms.calc = EMT()
dyn = QuasiNewton(atoms)
dyn.run(fmax=0.01)

# 2. Get Potential Energy and Vibrations
potentialenergy = atoms.get_potential_energy()
vib = Vibrations(atoms)
vib.run()
vib_energies = vib.get_energies()

# print the results - vibrational frequencies
vib.summary()

# 3. Thermochemistry Calculation
thermo = IdealGasThermo(
    vib_energies=vib_energies,
    potentialenergy=potentialenergy,
    atoms=atoms,
    geometry='linear', # or 'nonlinear'
    symmetrynumber=2,
    spin=0
)

# 4. Get H and S at 298.15K
H = thermo.get_enthalpy(temperature=298.15)
S = thermo.get_entropy(temperature=298.15, pressure=101325.0)
print(f"N2 molecule Enthalpy: {H} eV")
print(f"N2 molecule Entropy: {S} eV/K")

