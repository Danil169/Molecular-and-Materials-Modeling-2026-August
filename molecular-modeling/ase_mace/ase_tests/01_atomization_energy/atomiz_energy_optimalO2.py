#
#
#
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS


atom = Atoms('O', calculator=EMT())
e_atom = atom.get_potential_energy()

d = 1.2 # this is the approx O2 experimental bond length (in Angs)
molecule = Atoms('2O', [(0., 0., 0.), (0., 0., d)], calculator=EMT() )

opt = BFGS(molecule)
print('\n running geometry optimization of the molecule with the initial distance d(O-O)=',d)
opt.run(fmax=0.01)

# print out optimal internuclear distance
print('d(O-O)optimiz=',molecule.get_distance(0,1),' Ang (exp. 1.21 Ang)')

e_molecule = molecule.get_potential_energy()

e_atomization = 2 * e_atom - e_molecule

print('\n Oxygen atom energy: %5.2f eV' % e_atom)
print('Oxygen molecule energy: %5.2f eV' % e_molecule)
print('Atomization energy: %5.2f eV' % e_atomization,  ' experiment cca 5.17 eV ')

