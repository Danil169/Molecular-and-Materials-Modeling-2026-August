from ase.build import molecule
from ase.optimize import LBFGS
from ase.calculators.emt import EMT
from ase.io import write
from ase import Atoms
try:
    from openbabel import pybel
except ImportError:
    pybel = None

# ML Calculators
try:
    from chgnet.model.dynamics import CHGNetCalculator
except ImportError:
    CHGNetCalculator = None
try:
    from mace.calculators import mace_mp
except ImportError:
    mace_mp = None

def pybel_to_ase(mol):
    """Convert Pybel molecule to ASE Atoms object."""
    numbers = [a.atomicnum for a in mol.atoms]
    coords = [a.coords for a in mol.atoms]
    return Atoms(numbers=numbers, positions=coords)

def run_openbabel_ff(ff_name):
    """Calculates BDE and Bond Length using OpenBabel Force Fields."""
    results = []
    for species in ['C2H6', 'CH3']:
        # Create molecule from SMILES for best bond perception
        smi = "CC" if species == 'C2H6' else "[CH3]"
        mol = pybel.readstring("smi", smi)
        mol.addh()
        mol.make3D()
        
        # Optimize
        mol.localopt(forcefield=ff_name, steps=500)
        
        # Get Energy (converted from kcal/mol to kJ/mol)
        ff = pybel._forcefields[ff_name]
        ff.Setup(mol.OBMol)
        energy_kj = ff.Energy() * 4.184 
        
        # Convert to ASE for distance and saving
        ase_mol = pybel_to_ase(mol)
        write(f"{species}_{ff_name}_final.xyz", ase_mol)
        
        results.append((energy_kj, ase_mol))
    
    e_ethane, mol_e = results[0]
    e_methyl, mol_m = results[1]
    
    bde = (2 * e_methyl) - e_ethane
    dist = mol_e.get_distance(0, 1)
    return bde, dist

def run_ase_calc(calculator, name):
    """Calculates BDE and Bond Length using ASE Calculators."""
    # Ethane
    c2h6 = molecule('C2H6')
    c2h6.set_cell((15, 15, 15))
    c2h6.center()
    c2h6.pbc = True
    c2h6.rattle(stdev=0.01)
    c2h6.calc = calculator
    
    dyn_e = LBFGS(c2h6, trajectory=f"ethane_{name}.traj", logfile=None)
    dyn_e.run(fmax=0.05)
    write(f"ethane_{name}_final.xyz", c2h6)
    
    # Methyl
    ch3 = molecule('CH3')
    ch3.set_cell((15, 15, 15))
    ch3.center()
    ch3.pbc = True
    ch3.rattle(stdev=0.01)
    ch3.calc = calculator
    
    dyn_m = LBFGS(ch3, trajectory=f"methyl_{name}.traj", logfile=None)
    dyn_m.run(fmax=0.05)
    write(f"methyl_{name}_final.xyz", ch3)
    
    bde = (2 * ch3.get_potential_energy() - c2h6.get_potential_energy()) * 96.485
    dist = c2h6.get_distance(0, 1)
    return bde, dist

def main():
    table_data = []

    # 1. Run OpenBabel Force Fields
    if pybel:
        for ff in ["mmff94", "uff", "ghemical"]:
            print(f"Running OpenBabel {ff}...")
            try:
                bde, dist = run_openbabel_ff(ff)
                table_data.append((ff.upper(), bde, dist))
            except Exception as e:
                table_data.append((ff.upper(), f"Error: {e}", "N/A"))

    # 2. Run ASE Calculators
    ase_calcs = [("EMT", EMT())]
    if CHGNetCalculator:
        ase_calcs.append(("CHGNet", CHGNetCalculator()))
    if mace_mp:
        ase_calcs.append(("MACE", mace_mp(model="medium", device='cpu')))

    for name, calc in ase_calcs:
        print(f"Running ASE {name}...")
        try:
            bde, dist = run_ase_calc(calc, name)
            table_data.append((name, bde, dist))
        except Exception as e:
            table_data.append((name, "Error", "N/A"))

    # 3. Print Final Results Table
    header = f"{'Method':<12} | {'BDE (kJ/mol)':<15} | {'C-C Length (Å)':<15}"
    print(f"\n{header}")
    print("-" * len(header))
    for name, bde, dist in table_data:
        if isinstance(bde, float):
            print(f"{name:<12} | {bde:<15.2f} | {dist:<15.4f}")
        else:
            print(f"{name:<12} | {bde:<15} | {dist:<15}")

if __name__ == "__main__":
    main()

