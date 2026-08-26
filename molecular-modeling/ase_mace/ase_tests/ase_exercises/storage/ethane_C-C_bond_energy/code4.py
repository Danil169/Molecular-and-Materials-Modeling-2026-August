from ase.build import molecule
from ase.optimize import LBFGS
from ase.calculators.emt import EMT
from ase.io import write
import numpy as np

# Import ML Calculators
try:
    from chgnet.model.dynamics import CHGNetCalculator
except ImportError:
    CHGNetCalculator = None

try:
    from mace.calculators import mace_mp
except ImportError:
    mace_mp = None

def get_molecular_properties(calculator, name):
    """Calculates C-C BDE and equilibrium bond length, saving trajectories."""
    
    # 1. Setup Ethane
    c2h6 = molecule('C2H6')
    c2h6.set_cell([15, 15, 15])
    c2h6.center()
    c2h6.pbc = True 
    c2h6.rattle(stdev=0.01) # Avoid symmetry-based singular matrix errors
    
    c2h6.calc = calculator
    # Save optimization path to .traj file
    traj_file_e = f"ethane_{name}.traj"
    dyn_e = LBFGS(c2h6, trajectory=traj_file_e, logfile=None)
    dyn_e.run(fmax=0.05)
    
    # Save final converged structure
    write(f"ethane_{name}_final.xyz", c2h6)
    
    e_ethane = c2h6.get_potential_energy()
    bond_length = c2h6.get_distance(0, 1)

    # 2. Setup Methyl Radical
    ch3 = molecule('CH3')
    ch3.set_cell([15, 15, 15])
    ch3.center()
    ch3.pbc = True
    ch3.rattle(stdev=0.01)
    
    ch3.calc = calculator
    traj_file_m = f"methyl_{name}.traj"
    dyn_m = LBFGS(ch3, trajectory=traj_file_m, logfile=None)
    dyn_m.run(fmax=0.05)
    
    # Save final converged structure
    write(f"methyl_{name}_final.xyz", ch3)
    
    e_methyl = ch3.get_potential_energy()

    # 3. BDE Calculation (kJ/mol)
    bde_kj = (2 * e_methyl - e_ethane) * 96.485
    
    print(f"  [Files saved: {traj_file_e}, {traj_file_m}, and .xyz finals]")
    return bde_kj, bond_length

def main():
    calcs = [("EMT", EMT())]
    
    if CHGNetCalculator:
        calcs.append(("CHGNet", CHGNetCalculator()))
    
    if mace_mp:
        calcs.append(("MACE", mace_mp(model="medium", device='cpu')))

    results = []
    for name, calc in calcs:
        try:
            print(f"Running {name}...")
            bde, dist = get_molecular_properties(calc, name)
            results.append((name, bde, dist))
        except Exception as e:
            results.append((name, f"Error: {str(e)[:15]}", "N/A"))

    # Output Table
    header = f"{'Method':<12} | {'BDE (kJ/mol)':<15} | {'C-C Length (Å)':<15}"
    print(f"\n{header}")
    print("-" * len(header))
    for name, bde, dist in results:
        if isinstance(bde, float):
            print(f"{name:<12} | {bde:<15.2f} | {dist:<15.4f}")
        else:
            print(f"{name:<12} | {bde:<15} | {dist:<15}")

if __name__ == "__main__":
    main()

