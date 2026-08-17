from ase.build import molecule
from ase.calculators.psi4 import Psi4

print("==================================================")
# 1. Initialize a target molecular configuration
print("Initializing test water molecule...")
atoms = molecule('H2O')

# 2. Attach the native Psi4 Calculator
# Note: Lower memory limits prevent allocating large default chunks on shared systems
print("Configuring Psi4 calculator (HF/3-21g)...")
calc = Psi4(
    atoms=atoms, 
    method='hf', 
    basis='3-21g', 
    memory='500MB'
)
atoms.calc = calc

# 3. Request evaluation streams 
try:
    print("\n--- Running Psi4 Potential Energy Sweep ---")
    energy = atoms.get_potential_energy()
    print(f"✔ Success! Potential Energy: {energy:.4f} eV")
    
    print("\n--- Running Psi4 Gradient/Force Sweep ---")
    forces = atoms.get_forces()
    print("✔ Success! Forces Matrix Recycled:")
    print(forces)
    
except Exception as e:
    print(f"\n❌ Execution Check Failed: {e}")

print("==================================================")

