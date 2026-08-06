import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.calculators.morse import MorsePotential
from ase.optimize import BFGS

# 1. Initialize and optimize H2 to find simulation equilibrium
d_initial = 0.8
h2 = Atoms('H2', positions=[(0, 0, 0), (0, 0, d_initial)])

epsilon_param = 4.7446
sigma_param = 0.7416
rho_param = 1.44

morse_h2 = MorsePotential(epsilon=epsilon_param, sigma=sigma_param, rho=rho_param)
h2.calc = morse_h2

dyn = BFGS(h2, logfile=None)
dyn.run(fmax=0.01)

r_equilibrium = h2.get_distance(0, 1)
e_equilibrium = h2.get_potential_energy()

# 2. Generate the dissociation curve data points
# Scan H-H distances from a compressed 0.4 Å out to a separated 4.0 Å
distances = np.linspace(0.4, 4.0, 100)
energies = []

for r in distances:
    h2.set_distance(0, 1, r)
    energy = h2.get_potential_energy()
    energies.append(energy)

# 3. Reference Experimental Values
exp_re = 0.7414  # Å
exp_de = 4.7446  # eV
# Experimental minimum energy location on this absolute ASE energy scale (where infinity = 0)
exp_e_min = -exp_de 

# 4. Generate Plot
plt.figure(figsize=(8, 5))
plt.plot(distances, energies, label='Morse Potential Curve', color='blue', linewidth=2)
plt.axhline(0, color='gray', linestyle='--', alpha=0.7, label='Dissociated Limit (0 eV)')

# Mark simulation and experimental minima
plt.scatter(r_equilibrium, e_equilibrium, color='red', zorder=5, 
            label=f'Sim Equilibrium ({r_equilibrium:.3f} Å, {e_equilibrium:.3f} eV)')
plt.scatter(exp_re, exp_e_min, color='green', marker='x', s=100, zorder=5, 
            label=f'Exp Equilibrium ({exp_re:.3f} Å, {exp_e_min:.3f} eV)')

# Format chart elements
plt.title('$H_2$ Molecule Dissociation Curve (Morse Potential)', fontsize=12)
plt.xlabel('H-H Distance ($Å$)', fontsize=11)
plt.ylabel('Potential Energy ($eV$)', fontsize=11)
plt.ylim(-5.5, 5.0)
plt.grid(True, alpha=0.3)
plt.legend(loc='lower right')
plt.tight_layout()

# Display graph windows
plt.show()

