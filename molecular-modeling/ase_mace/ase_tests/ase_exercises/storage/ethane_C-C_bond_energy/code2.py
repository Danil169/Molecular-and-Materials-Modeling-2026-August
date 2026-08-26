from ase.build import molecule
from ase.optimize import LBFGS
from ase.calculators.emt import EMT

# Attempt to import CHGNet
try:
    from chgnet.model.dynamics import CHGNetCalculator
except ImportError:
    CHGNetCalculator = None

def get_molecular_properties(calculator):
    """Calculates C-C BDE and equilibrium bond length for a given calculator."""
    
    # 1. Setup and Relax Ethane
    # CHGNet requires a periodic box to define the neighbor graph correctly
    c2h6 = molecule('C2H6')
    c2h6.set_cell([15, 15, 15]) # 15A box to avoid self-interaction
    c2h6.center()
    c2h6.pbc = True 
    
    c2h6.calc = calculator
    dyn_e = LBFGS(c2h6, logfile=None)
    dyn_e.run(fmax=0.05)
    
    e_ethane = c2h6.get_potential_energy()
    # In ASE's C2H6 molecule, atoms 0 and 1 are the Carbons
    bond_length = c2h6.get_distance(0, 1)

    # 2. Setup and Relax Methyl Radical
    ch3 = molecule('CH3')
    ch3.set_cell([15, 15, 15])
    ch3.center()
    ch3.pbc = True
    
    ch3.calc = calculator
    dyn_m = LBFGS(ch3, logfile=None)
    dyn_m.run(fmax=0.05)
    e_methyl = ch3.get_potential_energy()

    # 3. Calculate BDE (2*CH3 - C2H6)
    # Conversion: 1 eV = 96.485 kJ/mol
    bde_kj = (2 * e_methyl - e_ethane) * 96.485
    
    return bde_kj, bond_length

def main():
    # Initialize calculators
    calcs_to_run = [("EMT", EMT())]
    
    if CHGNetCalculator:
        calcs_to_run.append(("CHGNet", CHGNetCalculator()))

    results = []
    for name, calc in calcs_to_run:
        try:
            bde, dist = get_molecular_properties(calc)
            results.append((name, bde, dist))
        except Exception as e:
            results.append((name, f"Error: {e}", "N/A"))

    # Print Table Output
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

