from openbabel import pybel
from ase import Atoms

def ase_to_pybel(ase_atoms):
    """Converts ASE Atoms to Pybel Molecule via XYZ string."""
    n_atoms = len(ase_atoms)
    symbols = ase_atoms.get_chemical_symbols()
    positions = ase_atoms.get_positions()
    
    xyz_lines = [f"{n_atoms}", "N2_Optimization_Test"]
    for sym, pos in zip(symbols, positions):
        xyz_lines.append(f"{sym} {pos[0]} {pos[1]} {pos[2]}")
    
    xyz_string = "\n".join(xyz_lines)
    return pybel.readstring("xyz", xyz_string)

# 1. Create a starting N2 molecule (slightly stretched at 1.2 A)
d_start = 1.2
n2_ase_initial = Atoms('2N', [(0., 0., 0.), (0., 0., d_start)])

print(f"{'Force Field':<12} | {'Bond Length (Å)':<15} | {'Energy (kJ/mol)':<15}")
print("-" * 50)

# 2. Iterate through all available force fields
for ff_name in pybel.forcefields:
    # Always start from the same geometry
    mol = ase_to_pybel(n2_ase_initial)
    ff = pybel._forcefields[ff_name]
    
    # 3. Setup the force field for this molecule
    # Setup can fail if the FF doesn't have parameters for the atom types
    if not ff.Setup(mol.OBMol):
        print(f"{ff_name:<12} | {'Setup Failed':<15} | {'N/A':<15}")
        continue
    
    # 4. Perform Optimization (500 steps of Steepest Descent)
    ff.SteepestDescent(500)
    ff.GetCoordinates(mol.OBMol) # Sync optimized coords back to the OBMol
    
    # 5. Extract Results
    # Get distance between the two Nitrogen atoms (indices 0 and 1)
    a1_coords = mol.atoms[0].coords
    a2_coords = mol.atoms[1].coords
    dist = ((a1_coords[0]-a2_coords[0])**2 + 
            (a1_coords[1]-a2_coords[1])**2 + 
            (a1_coords[2]-a2_coords[2])**2)**0.5
    
    energy = ff.Energy()
    
    print(f"{ff_name:<12} | {dist:<15.4f} | {energy:<15.4f}")

print("-" * 50)
print("Note: Experimental N-N distance is ~1.098 Å.")

