import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PHYSICAL CONSTANTS - USING CORRECT UNITS
# ============================================================================
# Atomic masses
mass_H_kg = 1.6735575e-27  # kg (1.00784 amu)
reduced_mass_H2_kg = mass_H_kg / 2  # kg for H2

# Physical constants
h = 6.62607015e-34      # J·s
hbar = 1.054571817e-34  # J·s
c_cm_per_s = 2.99792458e10  # cm/s
c_m_per_s = 2.99792458e8    # m/s
eV_to_J = 1.602176634e-19   # J/eV
J_to_eV = 1 / eV_to_J
eV_to_cm = 8065.54429   # 1 eV = 8065.54429 cm⁻¹
cm_to_eV = 1 / eV_to_cm

# Experimental values for H₂ (NIST data)
EXP_DATA = {
    'r_e': 0.7414,              # Å
    'D_e': 4.478,               # eV
    'omega_e': 4401.21,         # cm⁻¹
    'omega_e_xe': 121.33,       # cm⁻¹
    'ZPE': 0.270,               # eV
    'nu_0_1': 4161.14,          # cm⁻¹
}

def morse_potential(r, D_e, a, r_e):
    """Morse potential: V(r) = D_e * [1 - exp(-a*(r - r_e))]^2"""
    return D_e * (1 - np.exp(-a * (r - r_e)))**2

def calculate_spectroscopic_constants(D_e, a, r_e):
    """
    Calculate spectroscopic constants from Morse parameters
    Using CORRECT formulas with proper unit conversions
    """
    # Force constant: k = 2 * D_e * a^2 (eV/Å²)
    k_eV_per_A2 = 2 * D_e * a**2
    
    # Convert to N/m for comparison
    # 1 eV/Å² = 1.602e-19 J / (1e-10 m)² = 1.602e-19 / 1e-20 = 16.02 J/m²
    k_N_per_m = k_eV_per_A2 * 16.02
    
    # ========== HARMONIC FREQUENCY CALCULATION (CRITICAL FIX) ==========
    # ω = (1/(2πc)) * sqrt(k/μ)
    # First convert k to SI units (J/m²)
    k_SI = k_eV_per_A2 * eV_to_J / (1e-10)**2  # J/m²
    
    # Angular frequency in rad/s
    omega_rad_s = np.sqrt(k_SI / reduced_mass_H2_kg)
    
    # Convert to wavenumbers: ν̃ = ω/(2πc) where c in cm/s
    omega_cm = omega_rad_s / (2 * np.pi * c_cm_per_s)
    
    # ========== ANHARMONICITY CALCULATION ==========
    # For Morse oscillator: ωₑxₑ = ωₑ²/(4Dₑ) where Dₑ in cm⁻¹
    D_e_cm = D_e * eV_to_cm
    omega_e_xe_cm = omega_cm**2 / (4 * D_e_cm)
    
    # Alternative formula (should match):
    # ωₑxₑ = h_bar * a² / (8π² μ c)
    a_m = a * 1e10  # Convert Å⁻¹ to m⁻¹
    omega_e_xe_alt = (hbar * a_m**2) / (8 * np.pi**2 * reduced_mass_H2_kg * c_m_per_s)
    omega_e_xe_alt_cm = omega_e_xe_alt / (100 * c_cm_per_s)  # Convert to cm⁻¹
    
    # ========== ZERO-POINT ENERGY ==========
    ZPE_cm = omega_cm/2 - omega_e_xe_cm/4
    ZPE_eV = ZPE_cm * cm_to_eV
    
    # ========== FUNDAMENTAL FREQUENCY ==========
    nu_fundamental_cm = omega_cm - 2 * omega_e_xe_cm
    
    # ========== VIBRATIONAL ENERGY LEVELS ==========
    vibrational_levels = []
    v = 0
    while True:
        E_v_cm = omega_cm * (v + 0.5) - omega_e_xe_cm * (v + 0.5)**2
        if E_v_cm < D_e_cm:
            vibrational_levels.append((v, E_v_cm))
            v += 1
        else:
            break
    
    return {
        'r_e': r_e,
        'D_e': D_e,
        'D_e_cm': D_e_cm,
        'k_eV_per_A2': k_eV_per_A2,
        'k_N_per_m': k_N_per_m,
        'omega_cm': omega_cm,
        'omega_e_xe_cm': omega_e_xe_cm,
        'omega_e_xe_alt_cm': omega_e_xe_alt_cm,
        'ZPE_eV': ZPE_eV,
        'ZPE_cm': ZPE_cm,
        'nu_fundamental_cm': nu_fundamental_cm,
        'a': a,
        'vibrational_levels': vibrational_levels,
    }

# ============================================================================
# MAIN CALCULATION
# ============================================================================
print("="*80)
print("H₂ MORSE POTENTIAL - CORRECTED UNIT CONVERSIONS")
print("="*80)

# Use experimental values
D_e_exp = EXP_DATA['D_e']
r_e_exp = EXP_DATA['r_e']
omega_e_exp = EXP_DATA['omega_e']

# Calculate a from ωₑ using the correct formula
# From ω_e = (a/(2πc)) * sqrt(2D_e/μ)
# Solve for a: a = ω_e * 2πc / sqrt(2D_e/μ)
mu_kg = reduced_mass_H2_kg
a_correct = omega_e_exp * (2 * np.pi * c_cm_per_s) / np.sqrt(2 * D_e_exp * eV_to_J / mu_kg)
a_correct = a_correct / 1e10  # Convert to Å⁻¹

print(f"\n✓ Morse parameter a from experimental ωₑ = {omega_e_exp} cm⁻¹")
print(f"  a = {a_correct:.4f} Å⁻¹")

# Calculate all spectroscopic constants
results = calculate_spectroscopic_constants(D_e_exp, a_correct, r_e_exp)

print(f"\n{'='*80}")
print("MORSE POTENTIAL SPECTROSCOPIC CONSTANTS")
print(f"{'='*80}")
print(f"  Morse parameter a:           {results['a']:.4f} Å⁻¹")
print(f"  Bond length rₑ:              {results['r_e']:.4f} Å")
print(f"  Dissociation Dₑ:             {results['D_e']:.3f} eV ({results['D_e_cm']:.0f} cm⁻¹)")
print(f"  Force constant k:            {results['k_eV_per_A2']:.2f} eV/Å² = {results['k_N_per_m']:.0f} N/m")
print(f"  Harmonic freq ωₑ:            {results['omega_cm']:.2f} cm⁻¹")
print(f"  Anharmonicity ωₑxₑ:          {results['omega_e_xe_cm']:.2f} cm⁻¹")
print(f"  Anharmonicity (alt formula): {results['omega_e_xe_alt_cm']:.2f} cm⁻¹")
print(f"  Zero-point energy:           {results['ZPE_cm']:.2f} cm⁻¹ = {results['ZPE_eV']:.3f} eV")
print(f"  Fundamental freq (0→1):      {results['nu_fundamental_cm']:.2f} cm⁻¹")
print(f"  Number of bound states:      {len(results['vibrational_levels'])}")

# ============================================================================
# GENERATE POTENTIAL ENERGY CURVE
# ============================================================================
distances = np.linspace(0.4, 3.5, 500)
energies = morse_potential(distances, D_e_exp, a_correct, r_e_exp)

# ============================================================================
# COMPARISON TABLE
# ============================================================================
print("\n" + "="*80)
print("COMPARISON: Morse Potential vs Experimental Values")
print("="*80)
print(f"{'Property':<30} {'Morse':<15} {'Experimental':<15} {'Error':<15} {'Error %':<10}")
print("-"*80)

comparisons = [
    ('r_e', 'Å', EXP_DATA['r_e'], results['r_e']),
    ('D_e', 'eV', EXP_DATA['D_e'], results['D_e']),
    ('omega_e', 'cm⁻¹', EXP_DATA['omega_e'], results['omega_cm']),
    ('omega_e_xe', 'cm⁻¹', EXP_DATA['omega_e_xe'], results['omega_e_xe_cm']),
    ('ZPE', 'eV', EXP_DATA['ZPE'], results['ZPE_eV']),
    ('nu_fundamental', 'cm⁻¹', EXP_DATA['nu_0_1'], results['nu_fundamental_cm']),
]

prop_names = {
    'r_e': 'Bond length rₑ',
    'D_e': 'Dissociation energy Dₑ',
    'omega_e': 'Harmonic freq ωₑ',
    'omega_e_xe': 'Anharmonicity ωₑxₑ',
    'ZPE': 'Zero-point energy',
    'nu_fundamental': 'Fundamental freq (0→1)',
}

for prop, unit, exp_val, morse_val in comparisons:
    error = abs(morse_val - exp_val)
    error_pct = (error / exp_val) * 100 if exp_val != 0 else 0
    print(f"{prop_names[prop]:<30} {morse_val:>10.4f} {unit:<3} "
          f"{exp_val:>10.4f} {unit:<3} "
          f"{error:>10.4f} {unit:<3} "
          f"{error_pct:>8.2f}%")

# ============================================================================
# VISUALIZATION
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle('H₂ Morse Potential - Accurate Quantum Mechanical Model', 
             fontsize=14, fontweight='bold')

# Plot 1: Potential energy curve with vibrational levels
ax1 = axes[0, 0]
ax1.plot(distances, energies, 'b-', linewidth=2, label='Morse potential')

# Mark equilibrium and dissociation
ax1.axvline(r_e_exp, color='g', linestyle='--', alpha=0.7, 
           label=f'rₑ = {r_e_exp} Å')
ax1.axhline(0, color='g', linestyle=':', alpha=0.7, label='Minimum')
ax1.axhline(D_e_exp, color='r', linestyle=':', alpha=0.7, 
           label=f'Dissociation Dₑ = {D_e_exp} eV')

# Plot vibrational levels
for v, E_cm in results['vibrational_levels'][:10]:
    E_eV = E_cm * cm_to_eV
    ax1.axhline(E_eV, xmin=0.2, xmax=0.8, color='purple', alpha=0.3, linestyle=':')
    if v <= 5:
        ax1.text(3.2, E_eV, f'v={v}', fontsize=8, verticalalignment='bottom')

ax1.set_xlabel('Bond distance r (Å)', fontsize=11)
ax1.set_ylabel('Potential energy E (eV)', fontsize=11)
ax1.set_title('H₂ Potential Energy Curve with Vibrational Bound States')
ax1.set_xlim(0.4, 3.5)
ax1.set_ylim(-0.2, 5.0)
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: Vibrational energy level diagram
ax2 = axes[0, 1]
v_states = [v for v, _ in results['vibrational_levels']]
v_energies = [E_cm for _, E_cm in results['vibrational_levels']]

colors = plt.cm.plasma(np.linspace(0, 1, len(v_states)))
ax2.barh(v_states, v_energies, color=colors, edgecolor='navy', alpha=0.8)
ax2.set_xlabel('Energy (cm⁻¹)', fontsize=11)
ax2.set_ylabel('Vibrational quantum number v', fontsize=11)
ax2.set_title(f'H₂ Vibrational Energy Levels\nωₑ = {results["omega_cm"]:.1f} cm⁻¹, ωₑxₑ = {results["omega_e_xe_cm"]:.2f} cm⁻¹')

# Add energy labels for lowest levels
for v, E in zip(v_states[:8], v_energies[:8]):
    ax2.text(E + 30, v, f'{E:.0f} cm⁻¹', va='center', fontsize=7)

# Mark dissociation limit
ax2.axvline(results['D_e_cm'], color='red', linestyle='--', alpha=0.7, 
           label=f'Dissociation: {results["D_e_cm"]:.0f} cm⁻¹')
ax2.legend(loc='upper right', fontsize=8)
ax2.grid(True, axis='x', alpha=0.3)

# Plot 3: Comparison bar chart
ax3 = axes[1, 0]
properties = ['rₑ', 'Dₑ', 'ωₑ', 'ωₑxₑ', 'ZPE', 'ν₀→₁']
morse_vals = [results['r_e'], results['D_e'], results['omega_cm'], 
              results['omega_e_xe_cm'], results['ZPE_eV'], results['nu_fundamental_cm']]
exp_vals = [EXP_DATA['r_e'], EXP_DATA['D_e'], EXP_DATA['omega_e'], 
            EXP_DATA['omega_e_xe'], EXP_DATA['ZPE'], EXP_DATA['nu_0_1']]

# Scale values for better visualization
scaled_morse = [morse_vals[0], morse_vals[1], morse_vals[2]/100, 
                morse_vals[3], morse_vals[4], morse_vals[5]/1000]
scaled_exp = [exp_vals[0], exp_vals[1], exp_vals[2]/100, 
              exp_vals[3], exp_vals[4], exp_vals[5]/1000]

x = np.arange(len(properties))
width = 0.35

bars1 = ax3.bar(x - width/2, scaled_morse, width, label='Morse potential', 
                color='steelblue', alpha=0.8)
bars2 = ax3.bar(x + width/2, scaled_exp, width, label='Experimental', 
                color='darkorange', alpha=0.8)

ax3.set_xlabel('Property', fontsize=11)
ax3.set_ylabel('Value (scaled: ωₑ/100, ν₀→₁/1000)', fontsize=11)
ax3.set_title('Morse Potential vs Experimental Comparison')
ax3.set_xticks(x)
ax3.set_xticklabels(properties, fontsize=10)
ax3.legend()
ax3.grid(True, axis='y', alpha=0.3)

# Plot 4: Error analysis
ax4 = axes[1, 1]
errors = []
prop_labels = ['rₑ', 'Dₑ', 'ωₑ', 'ωₑxₑ', 'ZPE']

for i, prop in enumerate(prop_labels):
    exp_val = exp_vals[i]
    morse_val = morse_vals[i]
    if exp_val != 0:
        error_pct = abs(morse_val - exp_val) / exp_val * 100
        errors.append(error_pct)

colors = ['green' if e < 2 else 'orange' if e < 10 else 'red' for e in errors]
bars = ax4.bar(prop_labels, errors, color=colors, alpha=0.7, edgecolor='black')
ax4.set_ylabel('Relative Error (%)', fontsize=11)
ax4.set_title('Morse Potential Accuracy')
ax4.axhline(2, color='green', linestyle='--', alpha=0.5, label='2%')
ax4.axhline(10, color='orange', linestyle='--', alpha=0.5, label='10%')
ax4.legend()
ax4.grid(True, axis='y', alpha=0.3)

for bar, err in zip(bars, errors):
    height = bar.get_height()
    ax4.text(bar.get_x() + bar.get_width()/2., height,
            f'{err:.1f}%', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.show()

# ============================================================================
# PRINT VIBRATIONAL ENERGY LEVELS
# ============================================================================
print("\n" + "="*80)
print("VIBRATIONAL ENERGY LEVELS (cm⁻¹)")
print("="*80)
print(f"{'v':<5} {'E(v) (cm⁻¹)':<20} {'E(v) (eV)':<15} {'ΔE (cm⁻¹)':<15}")
print("-"*60)

prev_E = 0
for v, E_cm in results['vibrational_levels'][:15]:
    E_eV = E_cm * cm_to_eV
    delta = E_cm - prev_E if v > 0 else 0
    print(f"{v:<5} {E_cm:>10.1f} {'':>8} {E_eV:>8.4f} {'':>4} {delta:>10.1f}")
    prev_E = E_cm

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("\n" + "="*80)
print("CONCLUSION: Correct H₂ Morse Potential")
print("="*80)

print(f"""
✅ SUCCESS: Morse potential now correctly reproduces ωₑ = {results['omega_cm']:.2f} cm⁻¹

OPTIMAL MORSE PARAMETERS:
   Dₑ = {D_e_exp} eV
   a = {a_correct:.4f} Å⁻¹
   rₑ = {r_e_exp} Å

MORSE POTENTIAL FUNCTION:
   V(r) = {D_e_exp:.3f} * [1 - exp(-{a_correct:.4f} * (r - {r_e_exp:.4f}))]²

ACCURACY:
   • ωₑ error:      {abs(results['omega_cm']-EXP_DATA['omega_e']):.2f} cm⁻¹ ({abs(results['omega_cm']-EXP_DATA['omega_e'])/EXP_DATA['omega_e']*100:.2f}%)
   • ωₑxₑ error:    {abs(results['omega_e_xe_cm']-EXP_DATA['omega_e_xe']):.2f} cm⁻¹
   • ZPE error:     {abs(results['ZPE_eV']-EXP_DATA['ZPE']):.3f} eV

FINAL ANSWER TO YOUR ORIGINAL QUESTION:
   
   ❌ CHGNet CANNOT calculate H₂ properly - it's designed for periodic crystals
   ❌ The warnings from CHGNet explicitly state calculation will "go wrong"
   
   ✅ For H₂, use the Morse potential above with a = {a_correct:.4f} Å⁻¹
   ✅ This gives the correct potential energy surface and spectroscopic constants
   
   The Morse potential is the standard analytical model for diatomic molecules
   and will serve as your reference for H₂ calculations.
""")
