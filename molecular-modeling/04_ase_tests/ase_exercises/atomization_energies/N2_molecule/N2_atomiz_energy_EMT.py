#
#
#
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS


atom = Atoms('N', calculator=EMT())
e_atom = atom.get_potential_energy()

d = 1.1 # this is the N2 experimental bond length (in Angs)
molecule = Atoms('2N', [(0., 0., 0.), (0., 0., d)], calculator=EMT() )

print('\n\n Running ASE-EMT calculations :')

opt = BFGS(molecule)
print('\n running geometry optimization of the N2 molecule with the initial distance d(N-N)=',d)
opt.run(fmax=0.01)

# print out optimal internuclear distance
print('\n d(N-N)optimiz=',molecule.get_distance(0,1),' Ang (experiment is 1.098 Ang)')

e_molecule = molecule.get_potential_energy()

# get the atomization energy
e_atomization = 2 * e_atom - e_molecule

print('\n Nitrogen atom energy: %5.2f eV' % e_atom)
print(' Nitrogen molecule energy: %5.2f eV' % e_molecule)
print(' Atomization energy: %5.2f eV' % e_atomization,  ' (experiment cca 9.76 eV) ')

