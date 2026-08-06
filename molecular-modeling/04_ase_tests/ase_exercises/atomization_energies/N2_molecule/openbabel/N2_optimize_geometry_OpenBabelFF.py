#
#
#
from ase import Atoms
from openbabel import openbabel as ob
from openbabel import pybel

from ase import Atoms
from ase.io import write
from ase.optimize import BFGS

def pybel_to_ase(mol):
    # 1. Get atomic numbers
    numbers = [a.atomicnum for a in mol.atoms]
    # 2. Get 3D coordinates
    coords = [a.coords for a in mol.atoms]
    # 3. Create ASE Atoms object
    return Atoms(numbers=numbers, positions=coords)

def ase_to_pybel(ase_atoms):
    """Converts an ASE Atoms object to a Pybel Molecule object."""
    # 1. Create an XYZ string from ASE
    # Format: NumAtoms \n Comment \n Symbol X Y Z
    n_atoms = len(ase_atoms)
    symbols = ase_atoms.get_chemical_symbols()
    positions = ase_atoms.get_positions()
    
    xyz_lines = [f"{n_atoms}", "Converted from ASE"]
    for sym, pos in zip(symbols, positions):
        xyz_lines.append(f"{sym} {pos[0]} {pos[1]} {pos[2]}")
    
    xyz_string = "\n".join(xyz_lines)
    
    # 2. Read the string into Pybel
    return pybel.readstring("xyz", xyz_string)


d = 1.1 # this is the N2 experimental bond length (in Angs)
N2_molecule_ase = Atoms('2N', [(0., 0., 0.), (0., 0., d)])
N2_molecule_pybel = ase_to_pybel(N2_molecule_ase)

# List available force fields
print("\nList of openbabel FFs :",pybel.forcefields,"\n")

#thisff="uff"
thisff="gaff"
print("\n Selected FF :", thisff)
ff = pybel._forcefields[thisff]

success_molecule = ff.Setup(N2_molecule_pybel.OBMol)
if success_molecule:
    print(f"Energy of N2 molecule before optimization: {ff.Energy()} kJ/mol")

ob_log = pybel.ob.obErrorLog
# Set level to capture everything (4 = Debug, 2 = Info)
ob_log.SetOutputLevel(4)
N2_molecule_pybel.localopt(forcefield=thisff, steps=500)

all_messages = ob_log.GetMessagesOfLevel(2)
for msg in all_messages:
    print(f"Log: {msg}")

# 5. Clear the log so it doesn't grow indefinitely
ob_log.ClearLog()

success = ff.Setup(N2_molecule_pybel.OBMol)
if success:
    print(f"Energy of N2 molecule after optimization: {ff.Energy()} kJ/mol")


# 1. Sync the optimized coordinates back to ASE
N2_molecule_optimized_ase = pybel_to_ase(N2_molecule_pybel)

# 2. Check the new distance
opt_dist = N2_molecule_optimized_ase.get_distance(0, 1)

print(f"\nOptimized N-N distance (ASE): {opt_dist:.4f} Ang  (experiment is 1.098 Ang)")

# 3. Final potential energy in eV (ASE standard)
# 1 kJ/mol = 0.010364 eV
#e_ev = ff.Energy() * 0.010364
#print(f"Final Energy: {e_ev:.4f} eV")

