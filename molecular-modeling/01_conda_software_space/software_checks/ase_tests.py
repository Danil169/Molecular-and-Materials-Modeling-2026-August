from ase import Atoms
from ase.build import molecule
from ase.calculators.mopac import MOPAC
from ase.calculators.pyscf_calc import PySCF

# 1. Define a simple Water molecule
atoms = molecule('H2O')

# 2. Quick test for PySCF Calculator
print("Testing ASE + PySCF integration...")
atoms.calc = PySCF(method='RHF', basis='6-31g')
energy_pyscf = atoms.get_potential_energy()
print(f"PySCF Potential Energy: {energy_pyscf:.4f} eV\n")

# 3. Quick test for MOPAC Calculator (Semi-empirical PM7)
print("Testing ASE + MOPAC integration...")
atoms.calc = MOPAC(method='PM7')
energy_mopac = atoms.get_potential_energy()
print(f"MOPAC Potential Energy: {energy_mopac:.4f} eV")

print("\nAll quick tests passed successfully!")

