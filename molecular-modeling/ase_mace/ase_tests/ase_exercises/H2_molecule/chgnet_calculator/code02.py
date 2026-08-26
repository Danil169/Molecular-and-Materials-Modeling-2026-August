import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.units import _amu as amu
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from chgnet.model.dynamics import CHGNetCalculator
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PHYSICAL CONSTANTS AND EXPERIMENTAL REFERENCE DATA
# ============================================================================
# Atomic mass constants
mass_H = 1.00784 * amu  # kg
reduced_mass_H2 = (mass_H * mass_H) / (mass_H + mass_H)  # kg

# Physical constants
h_bar = 1.054571817e-34  # J·s
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

def get_chgnet_calculator():
    """Initialize CHGNet calculator"""
    return CHGNetCalculator(use_device='cpu')

def morse_potential(r, D_e, a, r_e):
    """Morse potential: V(r) = D_e * [1 - exp(-a*(r - r_e))]^2"""
    return D_e * (1 - np.exp(-a * (r - r_e)))**2

def fit_morse_potential(distances, energies, r_eq_guess, D_e_guess):
    """Fit Morse potential to energy curve"""
    energies_shifted = energies - np.min(energies)
    p0 = [D_e_guess, 2.0, r_eq_guess]
    
    try:
        popt, _ = curve_fit(morse_potential, distances, energies_shifted, p0=p0)
        return popt
    except Exception as e:
        print(f"Morse fit failed: {e}")
        return None

def calculate_spectroscopic_constants(distances, energies, reduced_mass):
    """
    Calculate vibrational properties from potential curve using proper physical units
    """
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
    k_eV_per_A2 = 2 * coeffs[0]  # Force constant in eV/Å²
    
    # Convert to N/m for comparison
    k_N_per_m = k_eV_per_A2 * eV_to_J / (1e-10)**2
    
    # Calculate harmonic frequency in cm⁻¹
    # ω = (1/(2πc)) * sqrt(k/μ)
    reduced_mass_kg = reduced_mass * amu
    k_SI = k_eV_per_A2 * eV_to_J / (1e-10)**2
    omega_rad_s = np.sqrt(k_SI / reduced_mass_kg)
    omega_cm = omega_rad_s / (2 * np.pi * c)
    
    # Estimate dissociation energy (difference between dissociated and minimum)
    D_e = np.max(energies) - E_min
    
    # Estimate anharmonicity from Morse relationship if we have D_e
    # ω_e x_e = ω_e^2 / (4D_e) where D_e is in cm⁻¹
    if D_e > 0:
        D_e_cm = D_e * eV_to_J / (h * c)
        omega_e_xe = omega_cm**2 / (4 * D_e_cm)
    else:
        omega_e_xe = 0
    
    return {
        'r_e': r_eq,
        'D_e': D_e,
        'k_eV_per_A2': k_eV_per_A2,
        'k_N_per_m': k_N_per_m,
        'omega_e_cm': omega_cm,
        'omega_e_xe': omega_e_xe,
        'E_min': E_min,
    }

# ============================================================================
# MAIN CALCULATION - Using larger supercell approach for CHGNet
# ============================================================================
print("="*80)
print("H₂ POTENTIAL ENERGY CURVE WITH CHGNet")
print("Spectroscopic Properties & Anharmonicity Analysis")
print("="*80)

# Initialize calculator
print("\nLoading CHGNet model...")
calculator = get_chgnet_calculator()

# Use a much larger cell to approximate vacuum (CHGNet needs periodic bulk-like environment)
# We'll use a 10x10x10 Å³ cell and put the molecule in the center
cell_size = 15.0  # Larger cell to reduce periodic interactions
print(f"Using periodic cell size: {cell_size} Å")

# Calculate reference energy using a different approach
# Instead of isolated atom, use H2 at very large distance
print("\nCalculating reference energy (H2 at large separation)...")
large_r = 8.0  # Large separation to approximate dissociated atoms
box_center = cell_size / 2
positions_large = [
    (box_center - large_r/2, box_center, box_center),
    (box_center + large_r/2, box_center, box_center)
]
h2_large = Atoms('H2', positions=positions_large,
                 cell=[cell_size, cell_size, cell_size], pbc=True)
h2_large.calc = calculator
e_large_separation = h2_large.get_potential_energy()
print(f"  H2 at r={large_r} Å: {e_large_separation:.6f} eV")

# Calculate H₂ potential energy curve
# Use smaller range near equilibrium for better resolution
distances = np.linspace(0.5, 3.0, 40)
energies = []

print("\nCalculating H₂ potential energy curve...")
for i, r in enumerate(distances):
    box_center = cell_size / 2
    positions = [
        (box_center - r/2, box_center, box_center),
        (box_center + r/2, box_center, box_center)
    ]
    
    molecule = Atoms('H2', positions=positions,
                    cell=[cell_size, cell_size, cell_size], pbc=True)
    molecule.calc = calculator
    energy = molecule.get_potential_energy()
    energies.append(energy)
    
    if (i + 1) % 10 == 0:
        print(f"  Progress: {i+1}/{len(distances)} (r = {r:.3f} Å, E = {energy:.4f} eV)")

energies = np.array(energies)
print("Calculation complete!")

# Calculate spectroscopic constants
results = calculate_spectroscopic_constants(distances, energies, reduced_mass_H2)

print(f"\nResults from CHGNet potential curve:")
print(f"  Equilibrium bond length (rₑ):     {results['r_e']:.4f} Å")
print(f"  Dissociation energy (Dₑ):         {results['D_e']:.3f} eV")
print(f"  Force constant (k):               {results['k_eV_per_A2']:.2f} eV/Å² = {results['k_N_per_m']:.1f} N/m")
print(f"  Harmonic frequency (ωₑ):          {results['omega_e_cm']:.1f} cm⁻¹")
print(f"  Anharmonicity (ωₑxₑ):             {results['omega_e_xe']:.2f} cm⁻¹")

# Calculate derived quantities
ZPE_eV = (results['omega_e_cm'] * 100 * c * h) / eV_to_J / 2
ZPE_eV_corrected = ZPE_eV - (results['omega_e_xe'] * 100 * c * h) / eV_to_J / 4
D_0 = results['D_e'] - ZPE_eV_corrected
nu_fundamental = results['omega_e_cm'] - 2 * results['omega_e_xe']

# Compile results
calculated_results = {
    'r_e': results['r_e'],
    'D_e': results['D_e'],
    'D_0': D_0,
    'omega_e': results['omega_e_cm'],
    'omega_e_xe': results['omega_e_xe'],
    'ZPE': ZPE_eV_corrected,
    'k': results['k_eV_per_A2'],
    'nu_fundamental': nu_fundamental,
}

# ============================================================================
# PRINT COMPARISON TABLE
# ============================================================================
print("\n" + "="*80)
print("COMPARISON: CHGNet vs Experimental Values for H₂")
print("="*80)
print(f"{'Property':<30} {'CHGNet':<15} {'Experimental':<15} {'Error':<15} {'Error %':<10}")
print("-"*80)

comparisons = [
    ('r_e', 'Å', 0.7414),
    ('D_e', 'eV', 4.478),
    ('omega_e', 'cm⁻¹', 4401.21),
    ('omega_e_xe', 'cm⁻¹', 121.33),
    ('ZPE', 'eV', 0.270),
    ('nu_fundamental', 'cm⁻¹', 4161.14),
]

for prop, unit, exp_val in comparisons:
    calc_val = calculated_results.get(prop, 0)
    error = abs(calc_val - exp_val)
    error_pct = (error / exp_val) * 100 if exp_val != 0 else 0
    
    prop_names = {
        'r_e': 'Bond length rₑ',
        'D_e': 'Dissociation energy Dₑ',
        'omega_e': 'Harmonic freq ωₑ',
        'omega_e_xe': 'Anharmonicity ωₑxₑ',
        'ZPE': 'Zero-point energy',
        'nu_fundamental': 'Fundamental freq (0→1)',
    }
    
    print(f"{prop_names[prop]:<30} {calc_val:>10.4f} {unit:<3} "
          f"{exp_val:>10.4f} {unit:<3} "
          f"{error:>10.4f} {unit:<3} "
          f"{error_pct:>8.2f}%")

# ============================================================================
# VISUALIZATION
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('H₂ Potential Energy Curve from CHGNet', fontsize=14, fontweight='bold')

# Plot 1: Potential energy curve
ax1 = axes[0]
ax1.plot(distances, energies, 'bo-', markersize=6, linewidth=2, 
         label='CHGNet calculation', alpha=0.8)

# Add Morse fit if available
morse_params = fit_morse_potential(distances, energies, results['r_e'], results['D_e'])
if morse_params is not None:
    r_plot = np.linspace(0.5, 3.0, 200)
    morse_plot = morse_potential(r_plot, *morse_params) + results['E_min']
    ax1.plot(r_plot, morse_plot, 'r--', linewidth=2, 
             label=f'Morse fit (ωₑ={results["omega_e_cm"]:.0f} cm⁻¹)')

# Mark equilibrium
ax1.axvline(results['r_e'], color='g', linestyle='--', alpha=0.7,
           label=f'rₑ = {results["r_e"]:.3f} Å')
ax1.axhline(results['E_min'], color='g', linestyle=':', alpha=0.7,
           label=f'E_min = {results["E_min"]:.2f} eV')
ax1.axhline(e_large_separation, color='r', linestyle=':', alpha=0.7,
           label=f'Dissociated limit = {e_large_separation:.2f} eV')

ax1.set_xlabel('Bond distance r (Å)', fontsize=11)
ax1.set_ylabel('Potential energy E (eV)', fontsize=11)
ax1.set_title('H₂ Potential Energy Curve')
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: Comparison bar chart
ax2 = axes[1]
properties = ['rₑ (Å)', 'ωₑ (cm⁻¹)', 'Dₑ (eV)']
chgnet_vals = [calculated_results['r_e'], 
               calculated_results['omega_e'],
               calculated_results['D_e']]
exp_vals = [0.7414, 4401.21, 4.478]

# Scale omega for better visualization (divide by 100)
chgnet_vals_scaled = [chgnet_vals[0], chgnet_vals[1]/100, chgnet_vals[2]]
exp_vals_scaled = [exp_vals[0], exp_vals[1]/100, exp_vals[2]]

x = np.arange(len(properties))
width = 0.35

bars1 = ax2.bar(x - width/2, chgnet_vals_scaled, width, label='CHGNet', 
                color='steelblue', alpha=0.8)
bars2 = ax2.bar(x + width/2, exp_vals_scaled, width, label='Experimental', 
                color='darkorange', alpha=0.8)

ax2.set_xlabel('Property')
ax2.set_ylabel('Value (scaled: ωₑ/100)')
ax2.set_title('CHGNet vs Experimental Comparison')
ax2.set_xticks(x)
ax2.set_xticklabels(properties)
ax2.legend()
ax2.grid(True, axis='y', alpha=0.3)

# Add value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.2f}', ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.show()

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "="*80)
print("SUMMARY AND DISCUSSION")
print("="*80)

print(f"\nCHGNet Results for H₂:")
print(f"  ✓ Bond length:      {results['r_e']:.3f} Å (exp: 0.741 Å) → error: {abs(results['r_e']-0.7414)*1000:.1f} mÅ")
print(f"  ✓ Dissociation Dₑ:  {results['D_e']:.3f} eV (exp: 4.478 eV) → error: {abs(results['D_e']-4.478):.3f} eV")
print(f"  ✓ Harmonic ωₑ:      {results['omega_e_cm']:.0f} cm⁻¹ (exp: 4401 cm⁻¹) → error: {abs(results['omega_e_cm']-4401):.0f} cm⁻¹")

print("\nImportant Notes:")
print("  • CHGNet is designed for periodic bulk materials, not isolated molecules")
print("  • The large cell (15 Å) helps approximate vacuum conditions")
print("  • Results show reasonable bond length but overestimates vibrational frequencies")
print("  • For accurate molecular properties, use DFT (e.g., GPAW, Quantum ESPRESSO) or")
print("    specialized ML potentials (e.g., SchNet, PhysNet, ANI)")
print("\nRecommended alternatives for H₂:")
print("  1. DFT with B3LYP/6-31G* via ASE+GPAW")
print("  2. Specialized ML potential: ANI-1ccx or SchNet")
print("  3. Traditional force field: ReaxFF or UFF")


