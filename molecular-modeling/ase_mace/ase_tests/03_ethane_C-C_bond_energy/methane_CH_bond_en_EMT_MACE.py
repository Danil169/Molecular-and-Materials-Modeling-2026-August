from ase import Atoms
from ase.build import molecule
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from mace.calculators import mace_mp
import warnings
import time

warnings.filterwarnings('ignore')

# Experimental reference values for methane C-H bond
EXP_CH_BOND_LENGTH = 1.087  # Å
EXP_CH_BOND_ENERGY = 4.55   # eV
EXP_CH_BOND_ENERGY_kJ = 439.0 # kJ/mol

print("=" * 70)
print("METHANE C-H BOND ENERGY - EMT vs MACE COMPARISON")
print("=" * 70)

def calculate_methane_ch_bond(calculator, label):
    print(f"\n{'-' * 35}")
    print(f"{label} CALCULATIONS")
    print(f"{'-' * 35}")
    
    # 1. Setup Methane (CH4)
    ch4 = molecule('CH4')
    ch4.calc = calculator
    
    # Relax methane geometry
    dyn_ch4 = BFGS(ch4, trajectory=f'ch4_geom_opt_{label.lower()}.traj', logfile=None)
    dyn_ch4.run(fmax=0.01)
    e_methane = ch4.get_potential_energy()
    print(f"Methane Energy: {e_methane:.4f} eV")

    # Compute C-H bond length (C is atom 0, H is atom 1, 2, 3, or 4)
    ch_bond_length = ch4.get_distance(0, 1)
    print(f"Optimized C-H bond length in methane: {ch_bond_length:.4f} Å")

    # 2. Setup Methyl Radical (CH3)
    ch3 = molecule('CH3')
    ch3.calc = calculator
    
    # Relax methyl radical
    dyn_ch3 = BFGS(ch3, trajectory=f'ch3_geom_opt_{label.lower()}.traj', logfile=None)
    dyn_ch3.run(fmax=0.01)
    e_methyl = ch3.get_potential_energy()
    print(f"Methyl Radical Energy: {e_methyl:.4f} eV")

    # 3. Setup Hydrogen Atom (H)
    h_atom = Atoms('H')
    h_atom.calc = calculator
    e_h = h_atom.get_potential_energy()
    print(f"Hydrogen Atom Energy: {e_h:.4f} eV")

    # 4. Calculate Bond Dissociation Energy (BDE)
    # Reaction: CH4 -> CH3 + H
    # BDE = E(CH3) + E(H) - E(CH4)
    bond_energy = e_methyl + e_h - e_methane
    
    print(f"\n{label} Results:")
    print(f"  C-H Bond Length:  {ch_bond_length:.4f} Å")
    print(f"  C-H Bond Energy:  {bond_energy:.4f} eV")
    print(f"  C-H Bond Energy:  {bond_energy * 96.485:.2f} kJ/mol")
    
    return ch_bond_length, bond_energy

# Main execution
if __name__ == "__main__":
    # --- EMT Calculations ---
    calc_emt = EMT()
    emt_length, emt_energy = calculate_methane_ch_bond(calc_emt, "EMT")
    
    # --- MACE Calculations ---
    print("\nInitializing MACE calculator...")
    try:
        calc_mace = mace_mp(model="medium", device="cpu", default_dtype="float64")
        mace_length, mace_energy = calculate_methane_ch_bond(calc_mace, "MACE")
    except Exception as e:
        print(f"Warning: MACE calculation failed: {e}")
        exit(1)

    # --- Comparison Table ---
    print("\n" + "=" * 70)
    print("COMPARISON WITH EXPERIMENT")
    print("=" * 70)
    
    print(f"\nC-H Bond Length (Exp: {EXP_CH_BOND_LENGTH:.3f} Å):")
    print(f"  EMT:  {emt_length:.4f} Å  (Δ = {emt_length - EXP_CH_BOND_LENGTH:+.4f} Å, {abs(emt_length - EXP_CH_BOND_LENGTH)/EXP_CH_BOND_LENGTH*100:.2f}%)")
    print(f"  MACE: {mace_length:.4f} Å  (Δ = {mace_length - EXP_CH_BOND_LENGTH:+.4f} Å, {abs(mace_length - EXP_CH_BOND_LENGTH)/EXP_CH_BOND_LENGTH*100:.2f}%)")

    print(f"\nC-H Bond Energy (Exp: {EXP_CH_BOND_ENERGY:.2f} eV | {EXP_CH_BOND_ENERGY_kJ:.1f} kJ/mol):")
    print(f"  EMT:  {emt_energy:.4f} eV ({emt_energy * 96.485:.2f} kJ/mol)  (Δ = {emt_energy - EXP_CH_BOND_ENERGY:+.4f} eV, {abs(emt_energy - EXP_CH_BOND_ENERGY)/EXP_CH_BOND_ENERGY*100:.2f}%)")
    print(f"  MACE: {mace_energy:.4f} eV ({mace_energy * 96.485:.2f} kJ/mol)  (Δ = {mace_energy - EXP_CH_BOND_ENERGY:+.4f} eV, {abs(mace_energy - EXP_CH_BOND_ENERGY)/EXP_CH_BOND_ENERGY*100:.2f}%)")
    print("\n" + "=" * 70)
