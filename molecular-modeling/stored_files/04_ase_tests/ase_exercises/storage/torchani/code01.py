import torch
import torchani
import ase
import ase.optimize
from ase.build import molecule
from ase.visualize import view
from ase.io.trajectory import Trajectory
import numpy as np

from torchani.grad import energies_forces_and_hessians, vibrational_analysis
from torchani.utils import get_atomic_masses

# 1. Setup
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = torchani.models.ANI1x(periodic_table_index=True).to(device).double()

# 2. Optimize
atoms = molecule('H2O')
atoms.calc = model.ase()
print("Optimizing geometry...")
opt = ase.optimize.BFGS(atoms)
opt.run(fmax=1e-6)

# 3. Tensors
species = torch.tensor(np.array([atoms.get_atomic_numbers()]), device=device)
coordinates = torch.tensor(np.array([atoms.get_positions()]), 
                           requires_grad=True, device=device, dtype=torch.double)

# 4. Analysis
masses = get_atomic_masses(species)
energies, forces, hessian = energies_forces_and_hessians(model, species, coordinates)
hessian = hessian.view(-1, 9, 9)

# freq (cm^-1), modes (batch, atoms, degrees_of_freedom)
# For H2O: (1, 3, 9). Dimension 1 is atoms, Dimension 2 is modes.
freq, modes, fconstants, rmasses = vibrational_analysis(masses, hessian, mode_kind="mdu")

# 5. Extract and Sort
all_freqs = freq.detach().cpu().numpy()
vib_indices = np.argsort(all_freqs)[::-1][:3]
comp_freqs = all_freqs[vib_indices]

exp_freqs = np.array([3756.0, 3657.0, 1595.0])
mode_names = ["Asymmetric Stretch", "Symmetric Stretch", "Bending"]

print("\n--- Vibrational Frequencies (cm^-1) ---")
for i in range(3):
    calc = comp_freqs[i]
    exp = exp_freqs[i]
    error = abs(calc - exp) / exp * 100
    print(f"{mode_names[i]:20}: Calc={calc:8.2f} | Exp={exp:8.2f} | Error={error:.2f}%")

# 6. Visualization
print("\nGenerating trajectories for visualization...")
for i, idx in enumerate(vib_indices):
    traj_name = f'water_mode_{i}.traj'
    
    # NEW INDEXING: modes is (batch, atoms, 3*atoms)
    # We take all atoms (:) for the specific mode index (idx)
    # This gives us a flattened vector that we reshape to (3 atoms, 3 coordinates)
    mode_vec_flat = modes[0, :, idx].detach().cpu().numpy()
    # If the flat vector is size 9 (3 atoms * 3 coords), reshape it:
    mode_vec = mode_vec_flat.reshape(3, 3) 
    
    with Trajectory(traj_name, 'w') as traj:
        for phase in np.linspace(0, 2 * np.pi, 20):
            displacement = mode_vec * np.sin(phase) * 0.4 
            frame = atoms.copy()
            frame.set_positions(atoms.get_positions() + displacement)
            traj.write(frame)
    print(f"Saved {mode_names[i]} to {traj_name}")

# Open the Asymmetric Stretch
print("\nOpening ASE GUI. Go to 'Tools -> Movie' and hit Play.")
view(ase.io.read('water_mode_0.traj', index=':'))

