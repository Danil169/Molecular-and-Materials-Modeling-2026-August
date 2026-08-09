from ase.build import molecule
from ase.optimize import LBFGS
from ase.calculators.emt import EMT
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
    """Calculates C-C BDE and equilibrium bond length."""
    
    # 1. Setup Ethane
    c2h6 = molecule('C2H6')
    # Put in a 15A periodic box (required for CHGNet/MACE)
    c2h6.set_cell([15, 15, 15])
    c2h6.center()
    c2h6.pbc = True 
    
    # Tiny random displacement to prevent 'Singular Matrix' due to perfect symmetry
    c2h6.rattle(stdev=0.01)
    
    c2h6.calc = calculator
    dyn_e = LBFGS(c2h6, logfile=None)
    dyn_e.run(fmax=0.05)
    
    e_ethane = c2h6.get_potential_energy()
    bond_length = c2h6.get_distance(0, 1)

    # 2. Setup Methyl Radical
    ch3 = molecule('CH3')
    ch3.set_cell([15, 15, 15])
    ch3.center()
    ch3.pbc = True
    ch3.rattle(stdev=0.01)
    
    ch3.calc = calculator
    dyn_m = LBFGS(ch3, logfile=None)
    dyn_m.run(fmax=0.05)
    e_methyl = ch3.get_potential_energy()

    # 3. BDE Calculation (kJ/mol)
    bde_kj = (2 * e_methyl - e_ethane) * 96.485
    
    return bde_kj, bond_length

def main():
    # Define Calculators
    calcs = [("EMT", EMT())]
    
    if CHGNetCalculator:
        calcs.append(("CHGNet", CHGNetCalculator()))
    
    if mace_mp:
        # 'medium' model is a good balance of speed and accuracy
        calcs.append(("MACE", mace_mp(model="medium", device='cpu')))

    results = []
    for name, calc in calcs:
        try:
            print(f"Running {name}...")
            bde, dist = get_molecular_properties(calc, name)
            results.append((name, bde, dist))
        except Exception as e:
            results.append((name, f"Error: {str(e)[:20]}...", "N/A"))

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

