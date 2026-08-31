from ase.build import molecule
from ase.optimize import QuasiNewton, BFGS
from ase.calculators.emt import EMT
from mace.calculators import mace_mp
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
import numpy as np
import warnings
import time

warnings.filterwarnings('ignore')

# Experimental reference values for O2 at 298.15 K and 1 atm
EXP_FREQ = 1580.0          # cm^-1
EXP_S = 0.002126           # eV/K (205.2 J/mol*K)
EXP_BOND_LENGTH = 1.207    # Å
EXP_ZPE = 0.098            # eV

print("=" * 80)
print("O2 MOLECULE THERMOCHEMISTRY: EMT vs MACE")
print("=" * 80)

def calculate_o2_properties(calculator, label):
    print(f"\n{'-'*40}")
    print(f"{label} CALCULATOR")
    print(f"{'-'*40}")
    
    atoms = molecule('O2')
    atoms.calc = calculator
    
    # 1. Optimization
    print("1. Optimizing geometry...")
    dyn = BFGS(atoms, logfile=None)
    dyn.run(fmax=0.01)
    
    bond_length = atoms.get_distance(0, 1)
    potentialenergy = atoms.get_potential_energy()
    print(f"Optimized bond length: {bond_length:.4f} Å")
    print(f"Potential energy: {potentialenergy:.6f} eV")
    
    # 2. Vibrational Analysis
    print("2. Calculating vibrations...")
    vib = Vibrations(atoms, name=f"vib_o2_{label.lower()}")
    vib.clean() # Prevent cache collision
    vib.run()
    
    vib_energies = vib.get_energies()
    vib_freqs = vib.get_frequencies()
    
    # Extract stretching mode
    positive_modes = [np.real(f) for f in vib_freqs if np.real(f) > 50] # filter translation/rotation
    calc_freq = max(positive_modes) if positive_modes else 0.0
    print(f"O-O stretching frequency: {calc_freq:.2f} cm⁻¹")
    
    # Filter energies for thermochemistry (>0.001 eV)
    vib_energies_filtered = np.array([np.real(e) for e in vib_energies if np.real(e) > 0.001])
    
    # 3. Thermochemistry
    print("3. Calculating thermochemistry (T=298.15 K)...")
    zpe = 0.5 * np.sum(vib_energies_filtered)
    
    # O2 is a triplet ground state! So spin=1 (which means 2S+1 = 3)
    thermo = IdealGasThermo(
        vib_energies=vib_energies_filtered,
        potentialenergy=potentialenergy,
        atoms=atoms,
        geometry='linear',
        symmetrynumber=2,
        spin=1  # Important for Oxygen!
    )
    
    H = thermo.get_enthalpy(temperature=298.15, verbose=False)
    S = thermo.get_entropy(temperature=298.15, pressure=101325.0, verbose=False)
    
    print(f"Zero-Point Energy (ZPE): {zpe:.4f} eV")
    print(f"Enthalpy (H): {H:.4f} eV")
    print(f"Entropy (S): {S:.6f} eV/K")
    
    return bond_length, calc_freq, zpe, S

if __name__ == "__main__":
    # EMT
    calc_emt = EMT()
    emt_d, emt_freq, emt_zpe, emt_S = calculate_o2_properties(calc_emt, "EMT")
    
    # MACE
    print("\nInitializing MACE...")
    try:
        calc_mace = mace_mp(model="medium", device="cpu", default_dtype="float64")
        mace_d, mace_freq, mace_zpe, mace_S = calculate_o2_properties(calc_mace, "MACE")
    except Exception as e:
        print(f"MACE Failed: {e}")
        exit(1)
        
    # Summary Table
    print("\n" + "=" * 80)
    print("SUMMARY COMPARISON VS EXPERIMENT")
    print("=" * 80)
    print(f"{'Property':<20} | {'EMT':<12} | {'MACE':<12} | {'Experiment':<12}")
    print("-" * 65)
    print(f"{'Bond Length (Å)':<20} | {emt_d:<12.4f} | {mace_d:<12.4f} | {EXP_BOND_LENGTH:<12.4f}")
    print(f"{'Frequency (cm⁻¹)':<20} | {emt_freq:<12.1f} | {mace_freq:<12.1f} | {EXP_FREQ:<12.1f}")
    print(f"{'ZPE (eV)':<20} | {emt_zpe:<12.4f} | {mace_zpe:<12.4f} | {EXP_ZPE:<12.4f}")
    print(f"{'Entropy (eV/K)':<20} | {emt_S:<12.6f} | {mace_S:<12.6f} | {EXP_S:<12.6f}")
    print("=" * 80)
