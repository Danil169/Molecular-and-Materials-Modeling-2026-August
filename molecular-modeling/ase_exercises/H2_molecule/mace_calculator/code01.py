import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
import warnings
warnings.filterwarnings('ignore')

# Import mace_models instead of mace.calculators
try:
    import mace_models
    print("✓ mace-models imported successfully")
except ImportError:
    print("✗ Please install: pip install mace-models")
    exit(1)

# Setup constants
CELL_SIZE = 10.0  # Angstroms
EXPERIMENTAL_DE = 4.52  # eV
EXPERIMENTAL_RE = 0.7414  # Å

print("\n" + "="*60)
print("H2 Potential Energy Curve with MACE Machine Learning Potential")
print("="*60)

# Load model using mace_models
print("\n1. Loading MACE model...")
# Try different model options
model_options = ["MACE-MP-0_small", "MACE-OFF23_small", "MACE-MP-0b_small"]
model = None
model_name_used = None

for model_name in model_options:
    try:
        print(f"   Trying {model_name}...")
        model = mace_models.load(model_name)
        model_name_used = model_name
        print(f"   ✓ Loaded {model_name}")
        break
    except Exception as e:
        print(f"   ✗ {model_name} failed: {e}")

if model is None:
    print("\n   Loading default model...")
    model = mace_models.load()
    model_name_used = "default"

# Get ASE calculator
ase_calculator = model.get_calculator(device='cpu', dtype="float64")

# 2. Calculate energy of a single isolated Hydrogen atom
print("\n2. Calculating isolated H atom energy...")
isolated_atom = Atoms('H', 
                      positions=[(CELL_SIZE/2, CELL_SIZE/2, CELL_SIZE/2)],
                      cell=[CELL_SIZE, CELL_SIZE, CELL_SIZE],
                      pbc=True)
isolated_atom.calc = ase_calculator
e_single_atom = isolated_atom.get_potential_energy()
e_two_isolated_atoms = 2 * e_single_atom
print(f"   Energy of single H atom: {e_single_atom:.8f} eV")
print(f"   Energy of two isolated H atoms: {e_two_isolated_atoms:.8f} eV")

# 3. Define range of interatomic distances
distances = np.linspace(0.5, 3.5, 45)
energies = []

# 4. Loop over distances
print("\n3. Calculating potential energy curve...")
print("-" * 50)

for i, r in enumerate(distances):
    box_center = CELL_SIZE / 2
    positions = [
        (box_center - r/2, box_center, box_center),
        (box_center + r/2, box_center, box_center)
    ]
    
    molecule = Atoms('H2', 
                    positions=positions,
                    cell=[CELL_SIZE, CELL_SIZE, CELL_SIZE],
                    pbc=True)
    molecule.calc = ase_calculator
    energy = molecule.get_potential_energy()
    energies.append(energy)
    
    if (i + 1) % 10 == 0 or i == 0 or i == len(distances) - 1:
        print(f"   r = {r:.3f} Å | E = {energy:.6f} eV")

energies = np.array(energies)
print("-" * 50)

# 5. Find equilibrium bond length
print("\n4. Finding equilibrium bond length...")
min_idx = np.argmin(energies)

# Quadratic fit
if 1 <= min_idx <= len(distances) - 2:
    x_fit = distances[min_idx-1:min_idx+2]
    y_fit = energies[min_idx-1:min_idx+2]
    coeffs = np.polyfit(x_fit, y_fit, 2)
    r_eq = -coeffs[1] / (2 * coeffs[0])
    e_min = np.polyval(coeffs, r_eq)
    print(f"   Quadratic fit: r_e = {r_eq:.4f} Å, E_min = {e_min:.6f} eV")
else:
    r_eq = distances[min_idx]
    e_min = energies[min_idx]
    print(f"   Discrete minimum: r_e = {r_eq:.4f} Å, E_min = {e_min:.6f} eV")

# 6. Compute Dissociation Energy
simulated_de = e_two_isolated_atoms - e_min

print(f"\n5. Dissociation Energy Calculation:")
print(f"   Simulated De = {simulated_de:.4f} eV")
print(f"   Experimental De = {EXPERIMENTAL_DE:.4f} eV")
print(f"   Error = {abs(simulated_de - EXPERIMENTAL_DE):.4f} eV")

# 7. Save data
output_filename = f"H2_potential_curve_{model_name_used}.txt"
data_to_save = np.column_stack((distances, energies))

header_text = (
    f"H2 Potential Energy Curve (MACE)\n"
    f"Model: {model_name_used}\n"
    f"Equilibrium Bond Length: {r_eq:.5f} A\n"
    f"Experimental r_e: {EXPERIMENTAL_RE:.4f} A\n"
    f"Simulated De: {simulated_de:.5f} eV\n"
    f"Experimental De: {EXPERIMENTAL_DE:.2f} eV\n"
    f"{'Distance(A)':<15} {'Energy(eV)':<20}"
)

np.savetxt(output_filename, data_to_save, fmt="%-15.6f %-20.8f", 
           header=header_text, comments='')
print(f"\n6. Data saved to {output_filename}")

# 8. Plot results
plt.figure(figsize=(11, 8))
plt.plot(distances, energies, 'b-o', markersize=5, lw=1.8, label='MACE Data', alpha=0.8)
plt.axvline(r_eq, color='r', linestyle='--', lw=2, label=f'r$_e$ (MACE) = {r_eq:.3f} Å')
plt.axvline(EXPERIMENTAL_RE, color='orange', linestyle=':', lw=2, label=f'r$_e$ (Exp) = {EXPERIMENTAL_RE:.3f} Å')
plt.axhline(e_min, color='g', linestyle='--', lw=1.5, label=f'E$_min$ = {e_min:.2f} eV')
plt.axhline(e_two_isolated_atoms, color='m', linestyle=':', lw=2, label=f'2H isolated = {e_two_isolated_atoms:.2f} eV')

plt.title(r'Potential Energy Curve of $\mathrm{H_2}$: MACE Machine Learning Potential', fontsize=13)
plt.xlabel(r'Interatomic Distance $r$ (Å)', fontsize=12)
plt.ylabel('Potential Energy $E$ (eV)', fontsize=12)

textstr = f'Model: {model_name_used}\nDe error: {abs(simulated_de - EXPERIMENTAL_DE):.3f} eV\nr$_e$ error: {abs(r_eq - EXPERIMENTAL_RE):.3f} Å'
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

plt.legend(loc='best', fontsize=10)
plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()
plt.savefig(f'H2_potential_curve_{model_name_used}.png', dpi=150)
plt.show()

print("\n" + "="*60)
print("✓ Calculation completed successfully!")
print("="*60)
