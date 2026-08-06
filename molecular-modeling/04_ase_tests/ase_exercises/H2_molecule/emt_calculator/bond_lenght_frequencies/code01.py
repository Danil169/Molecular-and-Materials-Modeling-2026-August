import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, u, eV
from ase import Atoms
from ase.calculators.emt import EMT

# =========================================================================
# 1. POTENTIAL ENERGY SURFACE (PES) MAPPING
# =========================================================================
# Wide range to plot the asymmetric well shape
distances_wide = np.linspace(0.4, 2.5, 100)
energies_wide = []

for r in distances_wide:
    atoms = Atoms('H2', positions=[(0, 0, 0), (0, 0, r)])
    atoms.calc = EMT()
    energies_wide.append(atoms.get_potential_energy())

# Find precise equilibrium r_e using a polynomial minimum locator
poly_coarse = np.polyfit(distances_wide, energies_wide, 6)
deriv1_coarse = np.polyder(poly_coarse, 1)
roots = np.roots(deriv1_coarse)
r_e_list = [r.real for r in roots if np.isreal(r) and 0.6 < r < 1.0]
r_e = float(r_e_list[0])

# =========================================================================
# 2. SEPARATED SCALAR SAMPLING FOR HIGH-RES DERIVATIVES
# =========================================================================
# Standard displacement step size in Angstroms
dx = 0.005  

def get_energy(bond_length):
    atoms = Atoms('H2', positions=[(0, 0, 0), (0, 0, bond_length)])
    atoms.calc = EMT()
    return float(atoms.get_potential_energy())

# Sample individual flat scalar variables to block array/matrix parsing bugs
E_minus2 = get_energy(r_e - 2 * dx)
E_minus1 = get_energy(r_e - dx)
E_zero   = get_energy(r_e)
E_plus1  = get_energy(r_e + dx)
E_plus2  = get_energy(r_e + 2 * dx)

# 5-Point Central Finite Difference formulas (units: eV/Å^n)
f2 = (E_minus1 - 2 * E_zero + E_plus1) / (dx ** 2)
f3 = (-0.5 * E_minus2 + E_minus1 - E_plus1 + 0.5 * E_plus2) / (dx ** 3)
f4 = (E_minus2 - 4 * E_minus1 + 6 * E_zero - 4 * E_plus1 + E_plus2) / (dx ** 4)

# =========================================================================
# 3. UNIT CONVERSION & SPECTROSCOPIC PHYSICAL MATH (VPT2)
# =========================================================================
# Convert derivatives to SI units (Joules and Meters)
k2 = f2 * eV / (1e-10) ** 2   # Force constant (N/m)
k3 = f3 * eV / (1e-10) ** 3   # Cubic constant (N/m^2)
k4 = f4 * eV / (1e-10) ** 4   # Quartic constant (N/m^3)

# Reduced mass of H2 (kg)
m_H = 1.007825 * u
mu = (m_H * m_H) / (m_H + m_H)

# Angular harmonic frequency (rad/s)
omega = np.sqrt(k2 / mu)

# Harmonic frequency omega_e in spectroscopic wavenumber units (cm^-1)
omega_e = (omega / (2.0 * np.pi * c)) * 0.01

# Reduced Planck constant
hbar = h / (2.0 * np.pi)

# Canonical VPT2 equation for anharmonicity (in Joules)
term_cubic = (5.0 * (k3 ** 2)) / (3.0 * k2)
we_xe_joules = (hbar ** 2 / (8.0 * (mu ** 2) * (omega ** 2))) * (term_cubic - k4)

# Convert Anharmonicity constant to standard wavenumbers (cm^-1)
we_xe = (we_xe_joules / (h * c)) * 0.01

# =========================================================================
# 4. PRINT REPORT COMPARISON
# =========================================================================
exp_re, exp_we, exp_wexe = 0.741, 4401.0, 121.3
err_re = ((r_e - exp_re) / exp_re) * 100
err_we = ((omega_e - exp_we) / exp_we) * 100
err_wexe = ((we_xe - exp_wexe) / exp_wexe) * 100

print("=" * 66)
print(f"{'Spectroscopic Property':<25} | {'Computed':<10} | {'Exp.':<8} | {'Error (%)':<9}")
print("=" * 66)
print(f"{'Bond Length (r_e)':<25} | {r_e:<10.3f} | {exp_re:<8.3f} | {err_re:<+8.1f}%")
print(f"{'Harmonic Freq. (omega_e)':<25} | {omega_e:<10.1f} | {exp_we:<8.1f} | {err_we:<+8.1f}%")
print(f"{'Anharmonicity (omega_e_xe)':<25} | {we_xe:<10.1f} | {exp_wexe:<8.1f} | {err_wexe:<+8.1f}%")
print("=" * 66)

# =========================================================================
# 5. GRAPH THE POTENTIAL CURVE
# =========================================================================
plt.figure(figsize=(9, 5.5))

# Plot calculated points from EMT
plt.scatter(distances_wide, energies_wide, color='darkorange', s=25, label='EMT Data Points', zorder=3)

# Plot continuous polynomial representation for visual context
r_fit = np.linspace(0.4, 2.5, 300)
e_fit = np.polyval(poly_coarse, r_fit)
plt.plot(r_fit, e_fit, color='teal', linestyle='-', linewidth=2, label='Potential Curve Fit')

# Mark equilibrium bond length vertical line
plt.axvline(x=r_e, color='red', linestyle=':', label=f'Equilibrium $r_e$ ({r_e:.3f} Å)')

# Formatting details
plt.title(r'$H_2$ Mapped Potential Well & Spectroscopic Verification', fontsize=12)
plt.xlabel('Interatomic Distance $R$ (Å)', fontsize=11)
plt.ylabel('Potential Energy $E$ (eV)', fontsize=11)
plt.xlim(0.35, 2.55)
plt.ylim(min(energies_wide) - 0.2, max(energies_wide) + 0.5)
plt.grid(True, alpha=0.3)
plt.legend(loc='upper right')

plt.tight_layout()
plt.show()

