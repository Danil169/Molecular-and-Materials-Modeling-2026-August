from ase import Atoms
from ase.build import molecule
from ase.io import write
from ase.calculators.lammpsrun import LAMMPS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

# --- Build initial system ---
methane = molecule('CH4')
oxygen1 = molecule('O2')
oxygen2 = molecule('O2')
system = methane + oxygen1 + oxygen2
system.center(vacuum=5.0)

# --- Define ReaxFF calculator ---
calc = LAMMPS(
    command='/usr/bin/lmp',
    pair_style='reax/c NULL',
    pair_coeff=['* * CHO.ff C H O'],
    fix=['qeq all reax/c 1 0.0 10.0 1e-6 reax/c'],
    files=['CHO.ff']
)

system.calc = calc

# --- Initialize velocities ---
MaxwellBoltzmannDistribution(system, temperature_K=2500)

# --- Run MD ---
dyn = VelocityVerlet(system, timestep=0.25 * units.fs)

def print_energy(a=system):
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print(f'Epot = {epot:.3f} eV/atom  Ekin = {ekin:.3f} eV/atom')

dyn.attach(print_energy, interval=10)

print("Starting combustion simulation...")
dyn.run(5000)

# --- Save trajectory ---
write('combustion.traj', system)

