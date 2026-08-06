from ase.build import molecule
from ase.optimize import LBFGS
from ase.calculators.emt import EMT
import numpy as np

# Attempt to import CHGNet
try:
    from chgnet.model.dynamics import CHGNetCalculator
except ImportError:
    CHGNetCalculator = None

def get_bde(calculator, name):
    """Calculates C-C BDE with periodic boundary conditions for ML potentials."""
    
    # 1. Setup Ethane
    c2h6 = molecule('C2H6')
    # CHGNet needs a box. 15A is enough to avoid self-interaction.
    c2h6.set_cell([15, 15, 15])
    c2h6.center()
    c2h6.pbc = True  # Enable periodic boundary conditions
    
    c2h6.calc = calculator
    dyn_e = LBFGS(c2h6, logfile=None)
    dyn_e.run(fmax=0.05)
    e_ethane = c2h6.get_potential_energy()

    # 2. Setup Methyl Radical
    ch3 = molecule('CH3')
    ch3.set_cell([15, 15, 15])
    ch3.center()
    ch3.pbc = True
    
    ch3.calc = calculator
    dyn_m = LBFGS(ch3, logfile=None)
    dyn_m.run(fmax=0.05)
    e_methyl = ch3.get_potential_energy()

    # 3. Calculate BDE (2*CH3 - C2H6)
    bde_ev = (2 * e_methyl) - e_ethane
    return bde_ev * 96.485

def main():
    calcs_to_run = [("EMT", EMT())]
    
    if CHGNetCalculator:
        # stress=None is often safer for isolated molecules in fixed boxes
        calcs_to_run.append(("CHGNet", CHGNetCalculator()))

    results = []
    for name, calc in calcs_to_run:
        try:
            bde_kj = get_bde(calc, name)
            results.append((name, bde_kj))
        except Exception as e:
            results.append((name, f"Error: {e}"))

    # Print Table Output
    print(f"\n{'Method':<15} | {'BDE (kJ/mol)':<15}")
    print("-" * 33)
    for name, val in results:
        if isinstance(val, float):
            print(f"{name:<15} | {val:<15.2f}")
        else:
            print(f"{name:<15} | {val:<15}")

if __name__ == "__main__":
    main()

