#
# with Google AI help
#

from openbabel import openbabel as ob
from openbabel import pybel

from ase import Atoms
from ase.io import write
from ase.visualize import view

def pybel_to_ase(mol):
    # 1. Get atomic numbers
    numbers = [a.atomicnum for a in mol.atoms]
    # 2. Get 3D coordinates
    coords = [a.coords for a in mol.atoms]
    # 3. Create ASE Atoms object
    return Atoms(numbers=numbers, positions=coords)

def print_water_geometry(mol,title):
    ase_mol = pybel_to_ase(mol)
    print(title)
    print(f"O-H1: {ase_mol.get_distance(0, 1):.4f}")
    print(f"O-H2: {ase_mol.get_distance(0, 2):.4f}")
    print(f"H-O-H Angle (deg): {ase_mol.get_angle(1, 0, 2):.3f}")

# List available force fields
print("\nList of openbabel FFs :",pybel.forcefields,"\n")

#mol = pybel.readstring("smi", "O")
#mol.addh()
#mol.make3D()
#mol.write("xyz", "water.xyz", overwrite=True)

mol = next(pybel.readfile("xyz", "water_uglygeom.xyz"))
print(f"Molecular Formula: {mol.formula}")

print_water_geometry(mol," Crude geometry:")

thisff="mmff94"
print("\n Selected FF :", thisff)

ff = pybel._forcefields[thisff]
success = ff.Setup(mol.OBMol)
if success:
    print(f"Energy before optimization: {ff.Energy()} kJ/mol")

# Redirect Open Babel internal logs to a file
#ob.obErrorLog.SetOutputStream("optimization_log.txt")
#ob.obErrorLog.SetOutputLevel(5) # Adjust level as needed

ob_log = pybel.ob.obErrorLog 

# Set level to capture everything (4 = Debug, 2 = Info)
ob_log.SetOutputLevel(2) 

mol.localopt(forcefield=thisff, steps=500)

all_messages = ob_log.GetMessagesOfLevel(2)

for msg in all_messages:
    print(f"Log: {msg}")

# 5. Clear the log so it doesn't grow indefinitely
ob_log.ClearLog()


success = ff.Setup(mol.OBMol)
if success:
    print(f"Energy after optimization: {ff.Energy()} kJ/mol")

#mol = pybel.readstring("smi", "C1=NC2=C(N1)C(=NC=N2)N") # Adenine
#mol.draw(show=True, filename=None)

# Convert, and write as ase geometry
ase_mol = pybel_to_ase(mol)
write("water_optgem_ase.xyz",ase_mol)

# This opens a separate GUI window (ase-gui)
#view(ase_mol)

mol.title = "FF optimized H2O"
mol.write("xyz", "water_optgem.xyz", overwrite=True)

# 4. Print optimized geometry
print_water_geometry(mol," Optimized geometry:")


