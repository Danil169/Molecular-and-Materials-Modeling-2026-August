from ase import Atoms
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import BFGS, QuasiNewton
from ase.thermochemistry import IdealGasThermo
from ase.vibrations import Vibrations
from mace.calculators import mace_mp
import numpy as np
import warnings
import time

warnings.filterwarnings('ignore')

# Experimental values for N2 at 298.15 K and 1 atm
EXP_FREQ = 2358.57  # cm^-1 (harmonic vibrational frequency)
EXP_S = 0.001916  # eV/K (191.6 J/mol·K)
EXP_BOND_LENGTH = 1.0975  # Angstroms
EXP_ZPE = 0.145  # eV (approximately)
EXP_ATOMIZATION_ENERGY = 9.76  # eV

print("=" * 80)
print("N2 MOLECULE THERMOCHEMISTRY: EMT vs MACE COMPARISON")
print("=" * 80)

def calculate_properties_with_vibrations(atoms, calculator, label, calc_type='ase'):
    """
    Calculate thermodynamic properties including vibrations.
    This is the full thermochemistry calculation.
    """
    print(f"\n{'='*80}")
    print(f"{label.upper()} CALCULATOR - FULL THERMOCHEMISTRY")
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
    
    # 2. Get Potential Energy
    potentialenergy = atoms.get_potential_energy()
    print(f"\nPotential energy: {potentialenergy:.6f} eV")
    
    # 3. Vibrational Analysis
    print("\n2. Vibrational Analysis")
    print("-" * 40)
    
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
        freq_error_pct = 100
        print("\nWarning: No positive vibrational modes found!")
    
    # 4. Thermochemistry
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
        'freq_error_pct': freq_error_pct,
        'zpe_error_pct': zpe_error_pct,
        'opt_time': opt_time,
        'vib_time': vib_time
    }

def calculate_atomization_energy(calculator, label):
    """
    Calculate atomization energy (single point calculations).
    This is a simpler calculation without vibrations.
    """
    print(f"\n{'='*80}")
    print(f"{label.upper()} - ATOMIZATION ENERGY")
    print('='*80)
    
    # Single N atom energy
    atom = Atoms('N', calculator=calculator)
    e_atom = atom.get_potential_energy()
    print(f"N atom energy: {e_atom:.6f} eV")
    
    # N₂ molecule (use optimized geometry from previous calculation if available)
    # We'll do a quick optimization here
    d_init = 1.1
    molecule = Atoms('2N', [(0., 0., 0.), (0., 0., d_init)], calculator=calculator)
    opt = BFGS(molecule)
    opt.run(fmax=0.01)
    
    d_opt = molecule.get_distance(0, 1)
    e_molecule = molecule.get_potential_energy()
    
    print(f"Optimized N-N bond length: {d_opt:.4f} Å")
    print(f"N₂ molecule energy: {e_molecule:.6f} eV")
    
    # Atomization energy
    e_atomization = 2 * e_atom - e_molecule
    
    print(f"Atomization energy: {e_atomization:.4f} eV")
    print(f"Experimental: {EXP_ATOMIZATION_ENERGY:.2f} eV")
    print(f"Error: {e_atomization - EXP_ATOMIZATION_ENERGY:+.4f} eV ({abs(e_atomization - EXP_ATOMIZATION_ENERGY)/EXP_ATOMIZATION_ENERGY*100:.2f}%)")
    
    return {
        'label': label,
        'bond_length': d_opt,
        'atomization_energy': e_atomization,
        'e_atom': e_atom,
        'e_molecule': e_molecule
    }

# Main execution
if __name__ == "__main__":
    # ============================================
    # Setup calculators
    # ============================================
    print("Initializing calculators...")
    
    # EMT calculator
    emt_calc = EMT()
    
    # MACE calculator - using the correct approach from your working script
    try:
        # Use float64 for better accuracy in geometry optimization
        mace_calc = mace_mp(
            model="medium",  # or "small", "large", "MACE-MP-0"
            device="cpu",    # or "cuda" for GPU
            default_dtype="float64"
        )
        mace_available = True
        print("MACE-MP (medium model) loaded successfully!")
        print("  Using float64 precision for better accuracy")
    except Exception as e:
        print(f"Warning: MACE calculator initialization failed: {e}")
        print("Please check that MACE is properly installed:")
        print("  pip install mace-torch")
        mace_available = False
    
    # ============================================
    # Full thermochemistry calculations (with vibrations)
    # ============================================
    print("\n" + "=" * 80)
    print("PART 1: FULL THERMOCHEMISTRY (with vibrations)")
    print("=" * 80)
    
    # EMT full thermochemistry
    atoms_emt = molecule('N2')
    emt_results = calculate_properties_with_vibrations(atoms_emt, emt_calc, "EMT")
    
    # MACE full thermochemistry if available
    if mace_available:
        atoms_mace = molecule('N2')
        mace_results = calculate_properties_with_vibrations(atoms_mace, mace_calc, "MACE")
    
    # ============================================
    # Atomization energy calculations (no vibrations)
    # ============================================
    print("\n" + "=" * 80)
    print("PART 2: ATOMIZATION ENERGY (single point)")
    print("=" * 80)
    
    emt_atomization = calculate_atomization_energy(emt_calc, "EMT")
    
    if mace_available:
        mace_atomization = calculate_atomization_energy(mace_calc, "MACE")
    
    # ============================================
    # Comparison Summary
    # ============================================
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY: EMT vs MACE vs EXPERIMENT")
    print("=" * 80)
    
    # Table 1: Structural and Vibrational Properties
    print("\n1. STRUCTURAL AND VIBRATIONAL PROPERTIES")
    print("-" * 80)
    print(f"{'Property':<25} {'EMT':<15} {'MACE':<15} {'Experimental':<15}")
    print("-" * 80)
    
    if mace_available:
        print(f"{'Bond Length (Å)':<25} {emt_results['bond_length']:<15.4f} {mace_results['bond_length']:<15.4f} {EXP_BOND_LENGTH:<15.4f}")
        print(f"{'Frequency (cm⁻¹)':<25} {emt_results['freq']:<15.2f} {mace_results['freq']:<15.2f} {EXP_FREQ:<15.2f}")
        print(f"{'ZPE (eV)':<25} {emt_results['zpe']:<15.4f} {mace_results['zpe']:<15.4f} {EXP_ZPE:<15.4f}")
        print(f"{'Entropy (eV/K)':<25} {emt_results['S']:<15.6f} {mace_results['S']:<15.6f} {EXP_S:<15.6f}")
    else:
        print(f"{'Bond Length (Å)':<25} {emt_results['bond_length']:<15.4f} {'N/A':<15} {EXP_BOND_LENGTH:<15.4f}")
        print(f"{'Frequency (cm⁻¹)':<25} {emt_results['freq']:<15.2f} {'N/A':<15} {EXP_FREQ:<15.2f}")
        print(f"{'ZPE (eV)':<25} {emt_results['zpe']:<15.4f} {'N/A':<15} {EXP_ZPE:<15.4f}")
        print(f"{'Entropy (eV/K)':<25} {emt_results['S']:<15.6f} {'N/A':<15} {EXP_S:<15.6f}")
    
    # Table 2: Errors
    print("\n2. ERRORS RELATIVE TO EXPERIMENT")
    print("-" * 80)
    print(f"{'Property':<25} {'EMT Error (%)':<15} {'MACE Error (%)':<15}")
    print("-" * 80)
    
    if mace_available:
        print(f"{'Bond Length':<25} {abs(emt_results['bond_length'] - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:<15.2f} {abs(mace_results['bond_length'] - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:<15.2f}")
        print(f"{'Frequency':<25} {emt_results['freq_error_pct']:<15.2f} {mace_results['freq_error_pct']:<15.2f}")
        print(f"{'ZPE':<25} {emt_results['zpe_error_pct']:<15.2f} {mace_results['zpe_error_pct']:<15.2f}")
        print(f"{'Entropy':<25} {abs(emt_results['S'] - EXP_S)/EXP_S*100:<15.2f} {abs(mace_results['S'] - EXP_S)/EXP_S*100:<15.2f}")
    else:
        print(f"{'Bond Length':<25} {abs(emt_results['bond_length'] - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100:<15.2f} {'N/A':<15}")
        print(f"{'Frequency':<25} {emt_results['freq_error_pct']:<15.2f} {'N/A':<15}")
        print(f"{'ZPE':<25} {emt_results['zpe_error_pct']:<15.2f} {'N/A':<15}")
        print(f"{'Entropy':<25} {abs(emt_results['S'] - EXP_S)/EXP_S*100:<15.2f} {'N/A':<15}")
    
    # Table 3: Atomization Energy
    print("\n3. ATOMIZATION ENERGY")
    print("-" * 80)
    print(f"{'Method':<15} {'Bond Length (Å)':<20} {'Atomization (eV)':<20} {'Error (%)':<15}")
    print("-" * 80)
    
    emt_error_pct = abs(emt_atomization['atomization_energy'] - EXP_ATOMIZATION_ENERGY)/EXP_ATOMIZATION_ENERGY*100
    print(f"{'EMT':<15} {emt_atomization['bond_length']:<20.4f} {emt_atomization['atomization_energy']:<20.4f} {emt_error_pct:<15.2f}")
    
    if mace_available:
        mace_error_pct = abs(mace_atomization['atomization_energy'] - EXP_ATOMIZATION_ENERGY)/EXP_ATOMIZATION_ENERGY*100
        print(f"{'MACE':<15} {mace_atomization['bond_length']:<20.4f} {mace_atomization['atomization_energy']:<20.4f} {mace_error_pct:<15.2f}")
    
    print(f"{'Experimental':<15} {EXP_BOND_LENGTH:<20.4f} {EXP_ATOMIZATION_ENERGY:<20.2f} {'N/A':<15}")
    
    # Table 4: Performance
    print("\n4. PERFORMANCE COMPARISON")
    print("-" * 80)
    print(f"{'Method':<15} {'Optimization (s)':<20} {'Vibrations (s)':<20} {'Total (s)':<15}")
    print("-" * 80)
    
    emt_total = emt_results['opt_time'] + emt_results['vib_time']
    print(f"{'EMT':<15} {emt_results['opt_time']:<20.2f} {emt_results['vib_time']:<20.2f} {emt_total:<15.2f}")
    
    if mace_available:
        mace_total = mace_results['opt_time'] + mace_results['vib_time']
        print(f"{'MACE':<15} {mace_results['opt_time']:<20.2f} {mace_results['vib_time']:<20.2f} {mace_total:<15.2f}")
        print(f"\nSpeedup (MACE/EMT): {mace_total/emt_total:.1f}x slower")
        print(f"Speedup (EMT/MACE): {emt_total/mace_total:.1f}x faster")
    
    # ============================================
    # Discussion and Analysis
    # ============================================
    print("\n" + "=" * 80)
    print("DISCUSSION AND ANALYSIS")
    print("=" * 80)
    
    if mace_available:
        print("""
    KEY FINDINGS:
    -------------
    1. Bond Length:
       • EMT:    {:.4f} Å (error: {:.2f}%)
       • MACE:   {:.4f} Å (error: {:.2f}%)
       • Exp:    {:.4f} Å
       • MACE is {:.1f}x more accurate for bond length
       
    2. Vibrational Frequency:
       • EMT:    {:.2f} cm⁻¹ (error: {:.1f}%)
       • MACE:   {:.2f} cm⁻¹ (error: {:.1f}%)
       • Exp:    {:.2f} cm⁻¹
       • MACE gives much better frequencies
       
    3. Zero-Point Energy (ZPE):
       • EMT:    {:.4f} eV (error: {:.1f}%)
       • MACE:   {:.4f} eV (error: {:.1f}%)
       • Exp:    {:.4f} eV
       • ZPE error follows frequency error
       
    4. Entropy:
       • EMT:    {:.6f} eV/K (error: {:.2f}%)
       • MACE:   {:.6f} eV/K (error: {:.2f}%)
       • Exp:    {:.6f} eV/K
       • Both methods give good entropy (dominated by translation/rotation)
       
    5. Atomization Energy:
       • EMT:    {:.4f} eV (error: {:.2f}%)
       • MACE:   {:.4f} eV (error: {:.2f}%)
       • Exp:    {:.2f} eV
       • EMT surprisingly good for N₂ atomization energy
    
    CONCLUSIONS:
    ------------
    • MACE is superior for vibrational properties (frequency, ZPE)
    • Entropy is well-predicted by both methods
    • EMT gives surprisingly good atomization energy for N₂
    • EMT is much faster ({:.1f}x faster) but less accurate for vibrations
    • For accurate thermochemistry, MACE is recommended
    • For quick estimates, EMT can be useful
    """.format(
            emt_results['bond_length'], abs(emt_results['bond_length'] - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100,
            mace_results['bond_length'], abs(mace_results['bond_length'] - EXP_BOND_LENGTH)/EXP_BOND_LENGTH*100,
            EXP_BOND_LENGTH,
            abs(emt_results['bond_length'] - EXP_BOND_LENGTH)/abs(mace_results['bond_length'] - EXP_BOND_LENGTH),
            emt_results['freq'], emt_results['freq_error_pct'],
            mace_results['freq'], mace_results['freq_error_pct'],
            EXP_FREQ,
            emt_results['zpe'], emt_results['zpe_error_pct'],
            mace_results['zpe'], mace_results['zpe_error_pct'],
            EXP_ZPE,
            emt_results['S'], abs(emt_results['S'] - EXP_S)/EXP_S*100,
            mace_results['S'], abs(mace_results['S'] - EXP_S)/EXP_S*100,
            EXP_S,
            emt_atomization['atomization_energy'], emt_error_pct,
            mace_atomization['atomization_energy'], mace_error_pct,
            EXP_ATOMIZATION_ENERGY,
            mace_total/emt_total
        ))
    else:
        print("""
    MACE not available. Please install it with:
        pip install mace-torch
        
    Then run with:
        python N2_thermodyn_EMT_MACE_02.py
    
    For accurate molecular thermochemistry, MACE is recommended over EMT.
    """)
    
    print("=" * 80)
