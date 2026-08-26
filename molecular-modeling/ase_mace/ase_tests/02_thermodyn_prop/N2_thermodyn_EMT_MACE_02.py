from ase.build import molecule
from ase.calculators.emt import EMT
from ase.calculators.mace import MACECalculator
from ase.optimize import QuasiNewton
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
import numpy as np
import warnings
import time

warnings.filterwarnings('ignore')

# Experimental values for N2 at 298.15 K and 1 atm
EXP_FREQ = 2358.57  # cm^-1 (harmonic vibrational frequency)
EXP_S = 0.001916  # eV/K (191.6 J/mol·K)
EXP_BOND_LENGTH = 1.0975  # Angstroms
EXP_ZPE = 0.145  # eV (approximately)

print("=" * 80)
print("N2 MOLECULE THERMOCHEMISTRY: EMT vs MACE COMPARISON")
print("=" * 80)

def calculate_properties(atoms, calculator, label):
    """Calculate thermodynamic properties for a given calculator."""
    print(f"\n{'='*80}")
    print(f"{label.upper()} CALCULATOR")
    print('='*80)
    
    # Set calculator
    atoms.calc = calculator
    
    # 1. Structure Optimization
    print("\n1. Structure Optimization")
    print("-" * 40)
    print(f"Initial bond length: {atoms.get_distance(0, 1):.4f} Å")
    
    start_time = time.time()
    dyn = QuasiNewton(atoms)
    dyn.run(fmax=0.01)
    opt_time = time.time() - start_time
    
    final_bond_length = atoms.get_distance(0, 1)
    print(f"Optimized bond length: {final_bond_length:.4f} Å")
    print(f"Error: {(final_bond_length - EXP_BOND_LENGTH)*1000:.2f} mÅ")
    print(f"Relative error: {abs(final_bond_length - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:.2f}%")
    print(f"Optimization time: {opt_time:.2f} seconds")
    
    # 2. Get Potential Energy and Vibrations
    print("\n2. Vibrational Analysis")
    print("-" * 40)
    potentialenergy = atoms.get_potential_energy()
    print(f"Potential energy: {potentialenergy:.6f} eV")
    
    start_time = time.time()
    vib = Vibrations(atoms)
    vib.run()
    vib_time = time.time() - start_time
    
    vib_energies = vib.get_energies()  # in eV
    vib_freqs_cm1 = vib.get_frequencies()  # in cm^-1
    
    print(f"Vibrational analysis time: {vib_time:.2f} seconds")
    print("\nCalculated vibrational frequencies (cm⁻¹):")
    positive_modes = []
    for i, freq in enumerate(vib_freqs_cm1):
        freq_real = np.real(freq)
        if abs(freq_real) > 1.0:  # Filter out near-zero modes
            positive_modes.append(freq_real)
            print(f"  Mode {i+1}: {freq_real:.2f} cm⁻¹")
        else:
            print(f"  Mode {i+1}: {freq_real:.2f} cm⁻¹ (translational/rotational)")
    
    # Get the stretching frequency
    if positive_modes:
        calc_freq = max(positive_modes)
        print(f"\nN-N stretching frequency: {calc_freq:.2f} cm⁻¹")
        print(f"Experimental: {EXP_FREQ:.2f} cm⁻¹")
        freq_error = calc_freq - EXP_FREQ
        freq_error_pct = abs(freq_error)/EXP_FREQ*100
        print(f"Error: {freq_error:.2f} cm⁻¹ ({freq_error_pct:.2f}%)")
    else:
        calc_freq = 0
        print("\nWarning: No positive vibrational modes found!")
    
    # 3. Thermochemistry
    print("\n3. Thermochemistry at 298.15 K, 1 atm")
    print("-" * 40)
    
    # Filter vib_energies for thermochemistry
    vib_energies_filtered = np.array([np.real(e) for e in vib_energies if np.real(e) > 0.001])
    
    # Manual ZPE calculation
    zpe = 0.5 * np.sum(vib_energies_filtered)
    print(f"Zero-Point Energy (ZPE): {zpe:.6f} eV")
    print(f"Experimental ZPE: {EXP_ZPE:.4f} eV")
    zpe_error = zpe - EXP_ZPE
    zpe_error_pct = abs(zpe_error)/EXP_ZPE*100
    print(f"ZPE Error: {zpe_error:.4f} eV ({zpe_error_pct:.2f}%)")
    
    # Create thermo object
    thermo = IdealGasThermo(
        vib_energies=vib_energies_filtered,
        potentialenergy=potentialenergy,
        atoms=atoms,
        geometry='linear',
        symmetrynumber=2,
        spin=0
    )
    
    temperature = 298.15
    pressure = 101325.0
    
    # Calculate properties silently
    H = thermo.get_enthalpy(temperature=temperature, verbose=False)
    S = thermo.get_entropy(temperature=temperature, pressure=pressure, verbose=False)
    G = thermo.get_gibbs_energy(temperature=temperature, pressure=pressure, verbose=False)
    U = thermo.get_internal_energy(temperature=temperature, verbose=False)
    
    print(f"\nInternal Energy (U): {U:.6f} eV")
    print(f"Enthalpy (H): {H:.6f} eV")
    print(f"Entropy (S): {S:.6f} eV/K")
    print(f"Gibbs Free Energy (G): {G:.6f} eV")
    
    # Manual heat capacity calculation
    def calculate_Cv(thermo, T, deltaT=0.1):
        U1 = thermo.get_internal_energy(temperature=T - deltaT/2, verbose=False)
        U2 = thermo.get_internal_energy(temperature=T + deltaT/2, verbose=False)
        return (U2 - U1) / deltaT
    
    Cv = calculate_Cv(thermo, temperature)
    kB = 8.617333262145e-5  # Boltzmann constant in eV/K
    Cp = Cv + kB
    
    print(f"Heat Capacity (Cv): {Cv:.6f} eV/K")
    print(f"Heat Capacity (Cp): {Cp:.6f} eV/K")
    
    # Print final summary
    print("\n" + "=" * 40)
    print("FINAL THERMODYNAMIC SUMMARY")
    print("=" * 40)
    thermo.get_enthalpy(temperature=temperature, verbose=True)
    
    # Return results for comparison
    return {
        'label': label,
        'bond_length': final_bond_length,
        'freq': calc_freq,
        'zpe': zpe,
        'S': S,
        'H': H,
        'G': G,
        'U': U,
        'Cv': Cv,
        'Cp': Cp,
        'potentialenergy': potentialenergy,
        'freq_error_pct': freq_error_pct if positive_modes else 100,
        'zpe_error_pct': zpe_error_pct,
        'opt_time': opt_time,
        'vib_time': vib_time
    }

# Main execution
if __name__ == "__main__":
    # Create N2 molecule
    atoms_emt = molecule('N2')
    atoms_mace = molecule('N2')
    
    # Setup calculators
    print("Initializing calculators...")
    
    # EMT calculator
    emt_calc = EMT()
    
    # MACE calculator - using the recommended model
    # You may need to adjust the model path based on your installation
    try:
        mace_calc = MACECalculator(
            model_paths=['/home/milias/work/projects/schools/Molecular-and-Materials-Modeling-2026-August/molecular-modeling/ase_mace/mace_mp_0.model'],  # Update this path
            device='cpu'  # Use 'cuda' if GPU is available
        )
        mace_available = True
    except Exception as e:
        print(f"Warning: MACE calculator initialization failed: {e}")
        print("Please update the model path in the script.")
        mace_available = False
    
    # Calculate with EMT
    emt_results = calculate_properties(atoms_emt, emt_calc, "EMT")
    
    # Calculate with MACE if available
    if mace_available:
        mace_results = calculate_properties(atoms_mace, mace_calc, "MACE")
    
    # Comparison Summary
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY: EMT vs MACE vs EXPERIMENT")
    print("=" * 80)
    
    print("\n{:<20} {:<15} {:<15} {:<15} {:<10}".format(
        "Property", "EMT", "MACE", "Experimental", "Best Method"
    ))
    print("-" * 80)
    
    properties = [
        ("Bond Length (Å)", 
         emt_results['bond_length'], 
         mace_results['bond_length'] if mace_available else "N/A",
         EXP_BOND_LENGTH),
        ("Frequency (cm⁻¹)", 
         emt_results['freq'], 
         mace_results['freq'] if mace_available else "N/A",
         EXP_FREQ),
        ("ZPE (eV)", 
         emt_results['zpe'], 
         mace_results['zpe'] if mace_available else "N/A",
         EXP_ZPE),
        ("Entropy (eV/K)", 
         emt_results['S'], 
         mace_results['S'] if mace_available else "N/A",
         EXP_S)
    ]
    
    for name, emt_val, mace_val, exp_val in properties:
        # Determine which method is closer
        if mace_available:
            emt_error = abs(emt_val - exp_val)
            mace_error = abs(mace_val - exp_val)
            best = "MACE" if mace_error < emt_error else "EMT"
            emt_str = f"{emt_val:.4f}"
            mace_str = f"{mace_val:.4f}"
        else:
            emt_str = f"{emt_val:.4f}"
            mace_str = "N/A"
            best = "EMT"
        
        print("{:<20} {:<15} {:<15} {:<15} {:<10}".format(
            name, emt_str, mace_str, f"{exp_val:.4f}", best
        ))
    
    print("\n" + "=" * 80)
    print("PERFORMANCE COMPARISON")
    print("=" * 80)
    
    print(f"\n{'Method':<20} {'Optimization (s)':<20} {'Vibrations (s)':<20} {'Total (s)'}")
    print("-" * 80)
    
    emt_total = emt_results['opt_time'] + emt_results['vib_time']
    print(f"{'EMT':<20} {emt_results['opt_time']:<20.2f} {emt_results['vib_time']:<20.2f} {emt_total:.2f}")
    
    if mace_available:
        mace_total = mace_results['opt_time'] + mace_results['vib_time']
        print(f"{'MACE':<20} {mace_results['opt_time']:<20.2f} {mace_results['vib_time']:<20.2f} {mace_total:.2f}")
        print(f"\nSpeedup (MACE/EMT): {mace_total/emt_total:.1f}x slower")
    
    print("\n" + "=" * 80)
    print("DISCUSSION")
    print("=" * 80)
    
    print("""
    EMT (Effective Medium Theory):
    - Pros: Very fast, good for metals and simple structures
    - Cons: Poor for covalent bonds, inaccurate vibrations
    - Best for: Quick estimates, metallic systems
    
    MACE (MACE-MP-0):
    - Pros: Accurate for molecules, good vibrations, near-DFT quality
    - Cons: Slower, requires pretrained model, memory intensive
    - Best for: Accurate thermochemistry of molecules and materials
    
    Key Findings:
    1. EMT underestimates N-N stretching frequency by ~48%
    2. MACE should predict frequencies much closer to experiment (~2358 cm⁻¹)
    3. Entropy is similar because it's dominated by translation/rotation
    4. MACE is significantly slower but worth it for accurate properties
    
    For accurate thermochemistry:
    - Use MACE (or other ML potentials) for molecules
    - Use DFT for more general systems
    - EMT only for metallic systems or quick estimates
    """)
    
    print("=" * 80)
