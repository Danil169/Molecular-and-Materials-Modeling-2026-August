from ase import Atoms
from ase.calculators.morse import MorsePotential
from ase.optimize import BFGS

# 1. Initialize H2 molecule near expected bond length
d_initial = 0.8
h2 = Atoms('H2', positions=[(0, 0, 0), (0, 0, d_initial)])

# Morse parameters parameterized to match the electronic potential depth
epsilon_param = 4.7446
sigma_param = 0.7416
rho_param = 1.44

morse_h2 = MorsePotential(epsilon=epsilon_param, sigma=sigma_param, rho=rho_param)
h2.calc = morse_h2

# 2. Relax the structure to its equilibrium minimum
dyn = BFGS(h2, logfile=None)
dyn.run(fmax=0.01)

e_equilibrium = h2.get_potential_energy()
r_equilibrium = h2.get_distance(0, 1)

# 3. Simulate fully broken bond state (atoms moved to 100 Angstroms away)
h2.set_distance(0, 1, 100.0)
e_dissociated = h2.get_potential_energy()

# Calculate simulation dissociation energy (D_e)
calculated_de = e_dissociated - e_equilibrium

# 4. Reference Experimental Values
exp_re = 0.7414  # Angstroms
exp_de = 4.7446  # eV (electronic well depth before Zero-Point Energy)

# 5. Printout Results with Comparison
print("=" * 55)
print(f"{'Metric':<25} | {'Simulation':<12} | {'Experiment':<10}")
print("=" * 55)
print(f"{'Equilibrium Distance (r_e)':<25} | {r_equilibrium:.4f} Å    | {exp_re:.4f} Å")
print(f"{'Dissociation Energy (D_e)':<25} | {calculated_de:.4f} eV   | {exp_de:.4f} eV")
print("-" * 55)
print(f"{'Error in Bond Length':<25} | {abs(r_equilibrium - exp_re):.4f} Å")
print(f"{'Error in Dissociation Energy':<25} | {abs(calculated_de - exp_de):.4f} eV")
print("=" * 55)

