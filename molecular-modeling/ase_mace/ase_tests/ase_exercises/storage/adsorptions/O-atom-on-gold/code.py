from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.calculators.emt import EMT
from ase.optimize import BFGS

# 1. Calculate energy of an isolated Oxygen atom (Reference)
adsorbate = Atoms('O')
adsorbate.calc = EMT()
e_o = adsorbate.get_potential_energy()

# 2. Calculate energy of the clean Gold surface (Reference)
# fcc111 creates a (111) surface; vacuum=10.0 adds space above
slab = fcc111('Au', size=(2, 2, 3), vacuum=10.0)
slab.calc = EMT()
e_slab = slab.get_potential_energy()

# 3. Add the Oxygen atom to the 'fcc' hollow site
# 'fcc' is a specific hollow site on the Au(111) surface
add_adsorbate(slab, 'O', height=1.5, position='fcc')

# 4. Relax the structure (Optimization)
# This moves the Oxygen atom until the forces are nearly zero
# fmax=0.05 is the force threshold for convergence
dyn = BFGS(slab, trajectory='opt.traj')
dyn.run(fmax=0.05)

# 5. Final Energy Analysis
e_total = slab.get_potential_energy()
e_bind = e_total - (e_slab + e_o)

print(f"Clean Slab Energy: {e_slab:.3f} eV")
print(f"Isolated O Energy: {e_o:.3f} eV")
print(f"Total Relaxed System Energy: {e_total:.3f} eV")
print("-" * 35)
print(f"Adsorption (Binding) Energy: {e_bind:.3f} eV")

