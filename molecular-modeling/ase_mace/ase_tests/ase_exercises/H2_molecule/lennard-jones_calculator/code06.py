import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.calculators.emt import EMT

# Lennard-Jones potential calculator for H2
class LennardJonesCalculator:
    """Lennard-Jones potential calculator for H2 molecule"""
    def __init__(self, epsilon=0.24, sigma=2.93):
        """
        Parameters:
        epsilon: well depth in eV (for H2: 0.24 eV)
        sigma: distance at zero potential in Angstroms (for H2: 2.93 Å)
        """
        self.epsilon = epsilon
        self.sigma = sigma
        self.name = "Lennard-Jones"
    
    def get_potential_energy(self, atoms):
        """Calculate Lennard-Jones potential energy for H2 molecule"""
        positions = atoms.get_positions()
        if len(positions) != 2:
            raise ValueError("Lennard-Jones calculator only works for H2 (2 atoms)")
        
        # Calculate distance between the two H atoms
        r = np.linalg.norm(positions[1] - positions[0])
        
        # Lennard-Jones potential: V(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
        if r > 0:
            term = (self.sigma / r) ** 6
            energy = 4 * self.epsilon * (term**2 - term)
        else:
            energy = np.inf
        
        return energy

# Morse potential calculator for H2
class MorseCalculator:
    """Morse potential calculator for H2 molecule"""
    def __init__(self, de=4.52, a=1.94, r0=0.74):
        """
        Parameters:
        de: dissociation energy in eV (for H2: 4.52 eV)
        a: width parameter in Angstroms^-1 (for H2: 1.94 Å^-1)
        r0: equilibrium bond length in Angstroms (for H2: 0.74 Å)
        """
        self.de = de
        self.a = a
        self.r0 = r0
        self.name = "Morse"
    
    def get_potential_energy(self, atoms):
        """Calculate Morse potential energy for H2 molecule"""
        positions = atoms.get_positions()
        if len(positions) != 2:
            raise ValueError("Morse calculator only works for H2 (2 atoms)")
        
        # Calculate distance between the two H atoms
        r = np.linalg.norm(positions[1] - positions[0])
        
        # Morse potential: V(r) = De * [1 - exp(-a*(r - r0))]^2
        energy = self.de * (1 - np.exp(-self.a * (r - self.r0))) ** 2
        
        return energy

# Function to calculate energy for a single atom (for reference)
def calculate_single_atom_energy(calculator_type="EMT"):
    """Calculate energy of a single isolated H atom"""
    atom = Atoms('H', positions=[(0, 0, 0)])
    
    if calculator_type == "EMT":
        atom.calc = EMT()
    elif calculator_type == "Lennard-Jones":
        atom.calc = LennardJonesCalculator()
    elif calculator_type == "Morse":
        atom.calc = MorseCalculator()
    else:
        raise ValueError(f"Unknown calculator type: {calculator_type}")
    
    return atom.get_potential_energy()

# Main calculation function
def calculate_h2_potential_curve(distances, calculator_type="EMT"):
    """Calculate potential energies for H2 molecule at given distances"""
    energies = []
    
    for r in distances:
        # Create H2 molecule with atoms along z-axis
        molecule = Atoms('H2', positions=[(0, 0, 0), (0, 0, r)])
        
        # Set the appropriate calculator
        if calculator_type == "EMT":
            molecule.calc = EMT()
        elif calculator_type == "Lennard-Jones":
            molecule.calc = LennardJonesCalculator()
        elif calculator_type == "Morse":
            molecule.calc = MorseCalculator()
        else:
            raise ValueError(f"Unknown calculator type: {calculator_type}")
        
        energy = molecule.get_potential_energy()
        energies.append(energy)
    
    return np.array(energies)

# Function to find equilibrium properties
def find_equilibrium(distances, energies):
    """Find equilibrium bond length and minimum energy"""
    min_idx = np.argmin(energies)
    r_eq = distances[min_idx]
    e_min = energies[min_idx]
    
    # Optional: Fit parabola around minimum for more accurate r_eq
    if min_idx > 0 and min_idx < len(distances) - 1:
        # Simple parabolic fit using three points around minimum
        x = distances[min_idx-1:min_idx+2]
        y = energies[min_idx-1:min_idx+2]
        coeffs = np.polyfit(x, y, 2)
        r_eq_fit = -coeffs[1] / (2 * coeffs[0])
        e_min_fit = np.polyval(coeffs, r_eq_fit)
        return r_eq_fit, e_min_fit
    
    return r_eq, e_min

# Main execution
def main():
    print("="*70)
    print("H2 Potential Energy Curve Analysis")
    print("Comparing EMT, Lennard-Jones, and Morse potentials")
    print("="*70)
    
    # Define distance range (extended to 3.5 Å)
    distances = np.linspace(0.4, 3.5, 150)  # 150 points for smooth curves
    
    # Calculate potential energies for each model
    print("\nCalculating potential energies...")
    
    print("  - EMT calculator...")
    energies_emt = calculate_h2_potential_curve(distances, "EMT")
    
    print("  - Lennard-Jones calculator...")
    energies_lj = calculate_h2_potential_curve(distances, "Lennard-Jones")
    
    print("  - Morse calculator...")
    energies_morse = calculate_h2_potential_curve(distances, "Morse")
    
    # Calculate single atom energies for each model
    print("\nCalculating single atom reference energies...")
    e_isolated_emt = calculate_single_atom_energy("EMT")
    e_isolated_lj = calculate_single_atom_energy("Lennard-Jones")
    e_isolated_morse = calculate_single_atom_energy("Morse")
    
    e_two_isolated_emt = 2 * e_isolated_emt
    e_two_isolated_lj = 2 * e_isolated_lj
    e_two_isolated_morse = 2 * e_isolated_morse
    
    # Shift LJ and Morse potentials to match separated atom limit of EMT for fair comparison
    # This aligns all curves at the dissociation limit
    shift_lj = e_two_isolated_emt - energies_lj[-1]
    shift_morse = e_two_isolated_emt - energies_morse[-1]
    energies_lj_shifted = energies_lj + shift_lj
    energies_morse_shifted = energies_morse + shift_morse
    
    # Find equilibrium properties
    print("\nFinding equilibrium properties...")
    r_eq_emt, e_min_emt = find_equilibrium(distances, energies_emt)
    r_eq_lj, e_min_lj = find_equilibrium(distances, energies_lj_shifted)
    r_eq_morse, e_min_morse = find_equilibrium(distances, energies_morse_shifted)
    
    # Calculate dissociation energies
    de_emt = e_two_isolated_emt - e_min_emt
    de_lj = e_two_isolated_emt - e_min_lj
    de_morse = e_two_isolated_emt - e_min_morse
    
    # Experimental values
    exp_de = 4.52  # eV
    exp_r_eq = 0.74  # Angstroms
    
    # Print results
    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)
    print(f"{'Model':<15} {'r_eq (Å)':<12} {'De (eV)':<12} {'De Error (eV)':<15} {'r_eq Error (Å)':<12}")
    print("-"*70)
    print(f"{'EMT':<15} {r_eq_emt:<12.3f} {de_emt:<12.3f} {abs(de_emt - exp_de):<15.3f} {abs(r_eq_emt - exp_r_eq):<12.3f}")
    print(f"{'Lennard-Jones':<15} {r_eq_lj:<12.3f} {de_lj:<12.3f} {abs(de_lj - exp_de):<15.3f} {abs(r_eq_lj - exp_r_eq):<12.3f}")
    print(f"{'Morse':<15} {r_eq_morse:<12.3f} {de_morse:<12.3f} {abs(de_morse - exp_de):<15.3f} {abs(r_eq_morse - exp_r_eq):<12.3f}")
    print("="*70)
    print(f"Experimental reference: r_eq = {exp_r_eq} Å, De = {exp_de} eV")
    
    # Save data to file
    output_file = "H2_potential_curves_data.txt"
    header = f"""H2 Potential Energy Curves Comparison
{'='*70}
Experimental reference: r_eq = {exp_r_eq} Å, De = {exp_de} eV

Results:
Model           r_eq (Å)    De (eV)     De Error (eV)
EMT             {r_eq_emt:.3f}       {de_emt:.3f}        {abs(de_emt - exp_de):.3f}
Lennard-Jones   {r_eq_lj:.3f}       {de_lj:.3f}        {abs(de_lj - exp_de):.3f}
Morse           {r_eq_morse:.3f}       {de_morse:.3f}        {abs(de_morse - exp_de):.3f}

Data columns: Distance(Å), EMT(eV), LJ_shifted(eV), Morse_shifted(eV)
{'='*70}
"""
    
    data = np.column_stack((distances, energies_emt, energies_lj_shifted, energies_morse_shifted))
    np.savetxt(output_file, data, fmt='%12.6f', header=header, comments='')
    print(f"\nData saved to {output_file}")
    
    # Create plots
    print("\nGenerating plots...")
    
    # Figure 1: Full potential energy curves
    fig1, ax1 = plt.subplots(figsize=(10, 7))
    
    ax1.plot(distances, energies_emt, 'b-', linewidth=2, label='EMT (ASE)', alpha=0.8)
    ax1.plot(distances, energies_lj_shifted, 'g-', linewidth=2, label='Lennard-Jones (shifted)', alpha=0.8)
    ax1.plot(distances, energies_morse_shifted, 'r-', linewidth=2, label='Morse (shifted)', alpha=0.8)
    
    # Mark equilibrium points
    ax1.plot(r_eq_emt, e_min_emt, 'bo', markersize=10, label=f'EMT min: r = {r_eq_emt:.3f} Å, E = {e_min_emt:.2f} eV')
    ax1.plot(r_eq_lj, e_min_lj, 'gs', markersize=10, label=f'LJ min: r = {r_eq_lj:.3f} Å, E = {e_min_lj:.2f} eV')
    ax1.plot(r_eq_morse, e_min_morse, 'r^', markersize=10, label=f'Morse min: r = {r_eq_morse:.3f} Å, E = {e_min_morse:.2f} eV')
    
    # Reference lines
    ax1.axhline(e_two_isolated_emt, color='m', linestyle='--', linewidth=1.5, 
                label=f'Dissociation limit: {e_two_isolated_emt:.2f} eV')
    ax1.axvline(exp_r_eq, color='k', linestyle=':', linewidth=1.5, 
                label=f'Experimental r_eq = {exp_r_eq} Å')
    
    ax1.set_xlabel('Interatomic Distance r (Å)', fontsize=12)
    ax1.set_ylabel('Potential Energy E (eV)', fontsize=12)
    ax1.set_title('H₂ Potential Energy Curves: EMT vs Lennard-Jones vs Morse', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3, linestyle=':')
    ax1.set_xlim(0.4, 3.5)
    
    plt.tight_layout()
    
    # Figure 2: Zoomed view near equilibrium
    fig2, ax2 = plt.subplots(figsize=(10, 7))
    
    ax2.plot(distances, energies_emt, 'b-', linewidth=2.5, label='EMT', alpha=0.8)
    ax2.plot(distances, energies_lj_shifted, 'g-', linewidth=2.5, label='Lennard-Jones', alpha=0.8)
    ax2.plot(distances, energies_morse_shifted, 'r-', linewidth=2.5, label='Morse', alpha=0.8)
    
    # Mark equilibrium points in zoomed view
    ax2.plot(r_eq_emt, e_min_emt, 'bo', markersize=12, markeredgecolor='black', markeredgewidth=1)
    ax2.plot(r_eq_lj, e_min_lj, 'gs', markersize=12, markeredgecolor='black', markeredgewidth=1)
    ax2.plot(r_eq_morse, e_min_morse, 'r^', markersize=12, markeredgecolor='black', markeredgewidth=1)
    
    ax2.set_xlabel('Interatomic Distance r (Å)', fontsize=12)
    ax2.set_ylabel('Potential Energy E (eV)', fontsize=12)
    ax2.set_title('H₂ Potential Energy Curves - Zoomed View Near Equilibrium', fontsize=14, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3, linestyle=':')
    ax2.set_xlim(0.5, 1.5)
    
    # Set y-limits to focus on the well region
    y_min = min(e_min_emt, e_min_lj, e_min_morse)
    y_max = max(e_min_emt + 3, e_min_lj + 3, e_min_morse + 3)
    ax2.set_ylim(y_min - 0.5, y_max)
    
    plt.tight_layout()
    
    # Figure 3: Error comparison bar chart
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    
    models = ['EMT', 'Lennard-Jones', 'Morse']
    de_errors = [abs(de_emt - exp_de), abs(de_lj - exp_de), abs(de_morse - exp_de)]
    r_eq_errors = [abs(r_eq_emt - exp_r_eq), abs(r_eq_lj - exp_r_eq), abs(r_eq_morse - exp_r_eq)]
    
    x = np.arange(len(models))
    width = 0.35
    
    bars1 = ax3.bar(x - width/2, de_errors, width, label='Dissociation Energy Error (eV)', color='skyblue', alpha=0.7)
    bars2 = ax3.bar(x + width/2, r_eq_errors, width, label='Bond Length Error (Å)', color='lightcoral', alpha=0.7)
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height + 0.01, f'{height:.3f}', ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height + 0.01, f'{height:.3f}', ha='center', va='bottom', fontsize=9)
    
    ax3.set_xlabel('Model', fontsize=12)
    ax3.set_ylabel('Error', fontsize=12)
    ax3.set_title('Comparison of Model Errors Relative to Experimental Values', fontsize=14, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(models)
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, alpha=0.3, axis='y', linestyle=':')
    
    plt.tight_layout()
    
    # Show all plots
    plt.show()
    
    # Print additional analysis
    print("\n" + "="*70)
    print("ADDITIONAL ANALYSIS")
    print("="*70)
    
    # Calculate force constants (spring constant at equilibrium)
    def calculate_force_constant(distances, energies, r_eq):
        """Calculate force constant k = d²E/dr² at equilibrium"""
        # Find indices near equilibrium
        idx = np.argmin(np.abs(distances - r_eq))
        if idx > 0 and idx < len(distances) - 1:
            # Second derivative using central difference
            h = distances[idx+1] - distances[idx]
            k = (energies[idx+1] - 2*energies[idx] + energies[idx-1]) / (h**2)
            # Convert to eV/Å²
            return k
        return 0.0
    
    k_emt = calculate_force_constant(distances, energies_emt, r_eq_emt)
    k_lj = calculate_force_constant(distances, energies_lj_shifted, r_eq_lj)
    k_morse = calculate_force_constant(distances, energies_morse_shifted, r_eq_morse)
    
    # Convert to N/m for comparison (1 eV/Å² = 160.2176 N/m)
    conversion = 160.2176
    print(f"\nForce Constants (spring constants) at equilibrium:")
    print(f"  EMT:           {k_emt:.2f} eV/Å² ({k_emt*conversion:.1f} N/m)")
    print(f"  Lennard-Jones: {k_lj:.2f} eV/Å² ({k_lj*conversion:.1f} N/m)")
    print(f"  Morse:         {k_morse:.2f} eV/Å² ({k_morse*conversion:.1f} N/m)")
    
    # Harmonic frequency approximation (ω = sqrt(k/μ))
    # Reduced mass of H2: μ = m_H/2 = 1.00784/2 = 0.50392 u = 8.367×10⁻²⁸ kg
    mu_kg = 8.367e-28  # kg
    eV_to_J = 1.602e-19
    hbar = 1.0546e-34  # J·s
    
    def calculate_frequency(k_eV_per_A2):
        """Calculate vibrational frequency in cm⁻¹"""
        k_N_per_m = k_eV_per_A2 * conversion
        omega_rad_per_s = np.sqrt(k_N_per_m / mu_kg)
        nu_cm1 = omega_rad_per_s / (2 * np.pi * 3e10)  # Convert to wavenumbers
        return nu_cm1
    
    if k_emt > 0:
        freq_emt = calculate_frequency(k_emt)
        print(f"\nVibrational frequencies (harmonic approximation):")
        print(f"  EMT:           {freq_emt:.0f} cm⁻¹")
        print(f"  Lennard-Jones: {calculate_frequency(k_lj):.0f} cm⁻¹")
        print(f"  Morse:         {calculate_frequency(k_morse):.0f} cm⁻¹")
        print(f"  Experimental:  4401 cm⁻¹ (fundamental)")  # Experimental H2 vibrational frequency
    
    # Determine best model
    de_errors_array = np.array([abs(de_emt - exp_de), abs(de_lj - exp_de), abs(de_morse - exp_de)])
    r_eq_errors_array = np.array([abs(r_eq_emt - exp_r_eq), abs(r_eq_lj - exp_r_eq), abs(r_eq_morse - exp_r_eq)])
    
    best_de_model = models[np.argmin(de_errors_array)]
    best_r_eq_model = models[np.argmin(r_eq_errors_array)]
    
    print(f"\nBest model for dissociation energy: {best_de_model}")
    print(f"Best model for bond length: {best_r_eq_model}")
    
    print("\n" + "="*70)
    print("Analysis complete!")
    print("="*70)

# Run the main function
if __name__ == "__main__":
    main()
