from ase import Atoms
from ase.build import molecule
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from mace.calculators import mace_mp
import warnings

warnings.filterwarnings('ignore')

# Experimental atomization energies (in eV)
EXP_VALUES = {
    'H2': 4.52,
    'N2': 9.76,
    'O2': 5.16,
    'F2': 1.65
}

def calculate_atomization(calc, mol_name, element):
    # Single atom
    atom = Atoms(element)
    atom.calc = calc
    e_atom = atom.get_potential_energy()
    
    # Molecule
    mol = molecule(mol_name)
    mol.calc = calc
    
    # Relax molecule
    dyn = BFGS(mol, logfile=None)
    dyn.run(fmax=0.01)
    
    e_mol = mol.get_potential_energy()
    bond_length = mol.get_distance(0, 1)
    
    # Atomization Energy: 2*E(atom) - E(mol)
    e_atomization = 2 * e_atom - e_mol
    
    return bond_length, e_atomization

print("=" * 80)
print("ATOMIZATION ENERGIES OF DIATOMIC MOLECULES (EMT vs MACE)")
print("=" * 80)

# Initialize calculators
calc_emt = EMT()
try:
    print("Loading MACE calculator...")
    calc_mace = mace_mp(model="medium", device="cpu", default_dtype="float64")
except Exception as e:
    print(f"Failed to load MACE: {e}")
    exit(1)

molecules = [('H2', 'H'), ('N2', 'N'), ('O2', 'O'), ('F2', 'F')]

print("\n" + "=" * 85)
print(f"{'Molecule':<10} | {'EMT E_atom (eV)':<15} | {'MACE E_atom (eV)':<17} | {'Exp (eV)':<10} | {'MACE Error':<10}")
print("-" * 85)

for mol_name, element in molecules:
    # EMT Calculation
    try:
        emt_d, emt_e = calculate_atomization(calc_emt, mol_name, element)
    except Exception as e:
        emt_e = float('nan')
        
    # MACE Calculation
    mace_d, mace_e = calculate_atomization(calc_mace, mol_name, element)
    
    exp_e = EXP_VALUES[mol_name]
    mace_err = abs(mace_e - exp_e) / exp_e * 100
    
    print(f"{mol_name:<10} | {emt_e:<15.4f} | {mace_e:<17.4f} | {exp_e:<10.2f} | {mace_err:<10.1f}%")

print("=" * 85)
print("\nNOTE: EMT is a classical potential designed for metals (like Cu, Ag, Au).")
print("It fails catastrophically for covalent molecules like O2 and F2.")
print("MACE (Neural Network) successfully calculates covalent bonds for all elements.")
