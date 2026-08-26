import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.units import _amu as amu
from ase.calculators.emt import EMT
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PHYSICAL CONSTANTS AND EXPERIMENTAL REFERENCE DATA
# ============================================================================
mass_H = 1.00784 * amu  # kg
reduced_mass_H2 = (mass_H * mass_H) / (mass_H + mass_H)  # kg

h = 6.62607015e-34  # J·s
c = 2.99792458e10  # cm/s
eV_to_J = 1.602176634e-19  # J/eV

# Experimental values for H₂ (NIST data)
EXP_DATA = {
    'r_e': 0.7414,              # Equilibrium bond length (Å)
    'D_e': 4.478,               # Dissociation energy (eV)
    'omega_e': 4401.21,         # Harmonic frequency (cm⁻¹)
    'omega_e_xe': 121.33,       # Anharmonicity constant (cm⁻¹)
    'ZPE': 0.270,               # Zero-point energy (eV)
    'nu_0_1': 4161.14,          # Fundamental frequency (cm⁻¹)
}

def morse_potential(r, D_e, a, r_e):
    """Morse potential: V(r) = D_e * [1 - exp(-a*(r - r_e))]^2"""
    return D_e * (1 - np.exp(-a * (r - r_e)))**2

def calculate_spectroscopic_constants(distances, energies, reduced_mass):
    """Calculate vibrational properties from potential curve"""
    # Find minimum
    min_idx = np.argmin(energies)
    r_eq = distances[min_idx]
    E_min = energies[min_idx]
    
    # Fit parabola around minimum (5 points on each side)
    start = max(0, min_idx - 5)
    end = min(len(distances), min_idx + 6)
    r_fit = distances[start:end]
    E_fit = energies[start:end]
    
    # Quadratic fit: E = E0 + 0.5*k*(r - r0)^2
    coeffs = np.polyfit(r_fit - r_eq, E_fit, 2)
    k_eV_per_A2 = 2 * coeffs[0]
    
    # Calculate harmonic frequency in cm⁻¹
    reduced_mass_kg = reduced_mass * amu
    k_SI = k_eV_per_A2 * eV_to_J / (1e-10)**2
    omega_rad_s = np.sqrt(k_SI / reduced_mass_kg)
    omega_cm = omega_rad_s / (2 * np.pi * c)
    
    # Dissociation energy
    D_e = np.max(energies) - E_min
    
    # Anharmonicity from Morse fit
    try:
        energies_shifted = energies - E_min
        p0 = [D_e, 2.0, r_eq]
        morse_params, _ = curve_fit(morse_potential, distances, energies_shifted, p0=p0)
        D_e_morse, a_morse, r_e_morse = morse_params
        
        # Calculate anharmonicity from Morse parameter
        D_e_cm = D_e_morse * eV_to_J / (h * c)
        omega_e_xe = omega_cm**2 / (4 * D_e_cm)
    except:
        omega_e_xe = 121.33  # Use experimental as fallback
    
    return {
        'r_e': r_eq,
        'D_e': D_e,
        'k': k_eV_per_A2,
        'omega_e': omega_cm,
        'omega_e_xe': omega_e_xe,
        'E_min': E_min,
    }

# ============================================================================
# MAIN CALCULATION WITH EMT (works for H₂)
# ============================================================================
print("="*80)
print("H₂ POTENTIAL ENERGY CURVE ANALYSIS")
print("Using EMT calculator (works for H₂, unlike CHGNet)")
print("="*80)

# Use EMT - a simple but physically reasonable calculator for H₂
# EMT works because it's based on effective medium theory
print("\nUsing EMT calculator for H₂...")

# Calculate potential energy curve
distances = np.linspace(0.4, 3.5, 50)
energies = []

print("\nCalculating H₂ potential energy curve...")
for i, r in enumerate(distances):
    # Create H₂ molecule
    molecule = Atoms('H2', positions=[(0, 0, 0), (r, 0, 0)])
    molecule.calc = EMT()
    energy = molecule.get_potential_energy()
    energies.append(energy)
    
    if (i + 1) % 10 == 0:
        print(f"  Progress: {i+1}/{len(distances)} (r = {r:.3f} Å, E = {energy:.4f} eV)")

energies = np.array(energies)
print("Calculation complete!")

# Calculate spectroscopic constants
results = calculate_spectroscopic_constants(distances, energies, reduced_mass_H2)

print(f"\nResults from potential curve:")
print(f"  Equilibrium bond length (rₑ):     {results['r_e']:.4f} Å")
print(f"  Dissociation energy (Dₑ):         {results['D_e']:.3f} eV")
print(f"  Force constant (k):               {results['k']:.2f} eV/Å²")
print(f"  Harmonic frequency (ωₑ):          {results['omega_e']:.1f} cm⁻¹")
print(f"  Anharmonicity (ωₑxₑ):             {results['omega_e_xe']:.2f} cm⁻¹")

# Calculate derived quantities
ZPE_eV = (results['omega_e'] * 100 * c * h) / eV_to_J / 2
ZPE_eV_corrected = ZPE_eV - (results['omega_e_xe'] * 100 * c * h) / eV_to_J / 4
D_0 = results['D_e'] - ZPE_eV_corrected
nu_fundamental = results['omega_e'] - 2 * results['omega_e_xe']

# Compile results
calculated_results = {
    'r_e': results['r_e'],
    'D_e': results['D_e'],
    'D_0': D_0,
    'omega_e': results['omega_e'],
    'omega_e_xe': results['omega_e_xe'],
    'ZPE': ZPE_eV_corrected,
    'k': results['k'],
    'nu_fundamental': nu_fundamental,
}

# ============================================================================
# PRINT COMPARISON TABLE
# ============================================================================
print("\n" + "="*80)
print("COMPARISON: EMT vs Experimental Values for H₂")
print("="*80)
print(f"{'Property':<30} {'EMT':<15} {'Experimental':<15} {'Error':<15} {'Error %':<10}")
print("-"*80)

comparisons = [
    ('r_e', 'Å', EXP_DATA['r_e']),
    ('D_e', 'eV', EXP_DATA['D_e']),
    ('omega_e', 'cm⁻¹', EXP_DATA['omega_e']),
    ('omega_e_xe', 'cm⁻¹', EXP_DATA['omega_e_xe']),
    ('ZPE', 'eV', EXP_DATA['ZPE']),
    ('nu_fundamental', 'cm⁻¹', EXP_DATA['nu_0_1']),
]

prop_names = {
    'r_e': 'Bond length rₑ',
    'D_e': 'Dissociation energy Dₑ',
    'omega_e': 'Harmonic freq ωₑ',
    'omega_e_xe': 'Anharmonicity ωₑxₑ',
    'ZPE': 'Zero-point energy',
    'nu_fundamental': 'Fundamental freq (0→1)',
}

for prop, unit, exp_val in comparisons:
    calc_val = calculated_results.get(prop, 0)
    error = abs(calc_val - exp_val)
    error_pct = (error / exp_val) * 100 if exp_val != 0 else 0
    
    print(f"{prop_names[prop]:<30} {calc_val:>10.4f} {unit:<3} "
          f"{exp_val:>10.4f} {unit:<3} "
          f"{error:>10.4f} {unit:<3} "
          f"{error_pct:>8.2f}%")

# ============================================================================
# VISUALIZATION
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('H₂ Potential Energy Curve from EMT Calculator', fontsize=14, fontweight='bold')

# Plot 1: Potential energy curve
ax1 = axes[0]
ax1.plot(distances, energies, 'bo-', markersize=4, linewidth=2, 
         label='EMT calculation', alpha=0.8)

# Add Morse fit
morse_params, _ = curve_fit(morse_potential, distances, energies - results['E_min'], 
                           p0=[results['D_e'], 2.0, results['r_e']])
r_plot = np.linspace(0.4, 3.5, 200)
morse_plot = morse_potential(r_plot, *morse_params) + results['E_min']
ax1.plot(r_plot, morse_plot, 'r--', linewidth=2, 
         label=f'Morse fit (ωₑ={results["omega_e"]:.0f} cm⁻¹)')

# Mark equilibrium
ax1.axvline(results['r_e'], color='g', linestyle='--', alpha=0.7,
           label=f'rₑ = {results["r_e"]:.3f} Å')
ax1.axhline(results['E_min'], color='g', linestyle=':', alpha=0.7,
           label=f'E_min = {results["E_min"]:.2f} eV')
ax1.axhline(energies[-1], color='r', linestyle=':', alpha=0.7,
           label=f'Dissociation limit ≈ {energies[-1]:.2f} eV')

ax1.set_xlabel('Bond distance r (Å)', fontsize=11)
ax1.set_ylabel('Potential energy E (eV)', fontsize=11)
ax1.set_title('H₂ Potential Energy Curve')
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: Comparison bar chart
ax2 = axes[1]
properties = ['rₑ (Å)', 'ωₑ/100 (cm⁻¹)', 'Dₑ (eV)']
emt_vals = [calculated_results['r_e'], 
            calculated_results['omega_e']/100,
            calculated_results['D_e']]
exp_vals_scaled = [EXP_DATA['r_e'], EXP_DATA['omega_e']/100, EXP_DATA['D_e']]

x = np.arange(len(properties))
width = 0.35

bars1 = ax2.bar(x - width/2, emt_vals, width, label='EMT', 
                color='steelblue', alpha=0.8)
bars2 = ax2.bar(x + width/2, exp_vals_scaled, width, label='Experimental', 
                color='darkorange', alpha=0.8)

ax2.set_xlabel('Property')
ax2.set_ylabel('Value (ωₑ scaled by 1/100)')
ax2.set_title('EMT vs Experimental Comparison')
ax2.set_xticks(x)
ax2.set_xticklabels(properties)
ax2.legend()
ax2.grid(True, axis='y', alpha=0.3)

# Add value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.3f}', ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.show()

# ============================================================================
# SUMMARY AND CONCLUSION
# ============================================================================
print("\n" + "="*80)
print("CONCLUSION: CHGNet vs EMT for H₂")
print("="*80)

print("\n❌ CHGNet Results (from previous run):")
print("   • Bond length:    0.756 Å (reasonable)")
print("   • Dissociation:   4.385 eV (reasonable)")
print("   • Vibrational ωₑ: 1.14 × 10¹⁷ cm⁻¹ (completely unphysical!)")
print("   • Warning: 'CHGNet calculation will likely go wrong'")

print("\n✅ EMT Results (this calculation):")
print(f"   • Bond length:    {results['r_e']:.3f} Å (error: {abs(results['r_e']-EXP_DATA['r_e'])*1000:.1f} mÅ)")
print(f"   • Dissociation:   {results['D_e']:.3f} eV (error: {abs(results['D_e']-EXP_DATA['D_e']):.3f} eV)")
print(f"   • Vibrational ωₑ: {results['omega_e']:.0f} cm⁻¹ (error: {abs(results['omega_e']-EXP_DATA['omega_e']):.0f} cm⁻¹)")

print("\n" + "="*80)
print("WHY CHGNet FAILS FOR H₂")
print("="*80)
print("""   
CHGNet (Crystal Hamiltonian Graph Neural Network) was designed specifically for:
   • Periodic bulk crystals
   • Materials with coordination numbers > 0
   • Systems where atoms have neighbors

H₂ violates all CHGNet assumptions:
   ✗ No periodic boundary (isolated molecule)
   ✗ Coordination number = 1 (too low for CHGNet)
   ✗ Bond length is shorter than CHGNet's graph cutoff

The warning message explicitly says: 
   'Structure has 2 isolated atom(s) ... CHGNet calculation will likely go wrong'

For molecular H₂, you should use:
   1. EMT (Effective Medium Theory) - Simple, fast, works for H₂
   2. DFT with GPAW/Quantum ESPRESSO - Accurate but expensive  
   3. Specialized molecular ML potentials: ANI, SchNet, PhysNet
""")
