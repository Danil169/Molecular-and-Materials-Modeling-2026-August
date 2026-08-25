from ase import Atoms
from ase.optimize import BFGS
from chgnet.model.dynamics import CHGNetCalculator  # Required for CHGNet

# 1. Initialize the CHGNet calculator
# This loads the pretrained weights (default version 0.3.0)
calc = CHGNetCalculator()

# 2. Set up single atom energy
# Neural networks often need a 'cell' even for single atoms
atom = Atoms('N', pbc=True, cell=[10, 10, 10], calculator=calc)
e_atom = atom.get_potential_energy()

# 3. Set up N2 molecule with experimental bond length d = 1.1 Å
d = 1.1 
molecule = Atoms('2N', 
                 positions=[(0., 0., 0.), (0., 0., d)], 
                 pbc=True, 
                 cell=[10, 10, 10], 
                 calculator=calc)

print('\n Running ASE-CHGNet calculations :')
opt = BFGS(molecule)
print(f'\n Running geometry optimization with initial d(N-N) = {d}')
opt.run(fmax=0.01)

# 4. Results
d_opt = molecule.get_distance(0, 1)
e_molecule = molecule.get_potential_energy()
e_atomization = 2 * e_atom - e_molecule

print(f'\n Optimized d(N-N): {d_opt:.3f} Ang (Experiment: 1.098 Ang)')
print(f' Nitrogen atom energy: {e_atom:.2f} eV')
print(f' Nitrogen molecule energy: {e_molecule:.2f} eV')
print(f' Atomization energy: {e_atomization:.2f} eV (Experiment: ~9.76 eV)')

