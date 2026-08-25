from ase.build import molecule
from ase.calculators.lammpslib import LAMMPSlib
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

# Build methane + oxygen system
methane = molecule('CH4')
oxygen1 = molecule('O2')
oxygen2 = molecule('O2')
atoms = methane + oxygen1 + oxygen2
atoms.center(vacuum=5.0)

# ReaxFF commands (LAMMPS 2024 uses 'reaxff')
lmpcmds = [
    "pair_style reaxff NULL",
    "pair_coeff * * CHO.ff C H O",
    "fix qeq all qeq/reaxff 1 0.0 10.0 1e-6 reaxff",
    "fix species all reaxff/species 100 species.out"
]

calc = LAMMPSlib(
    lmpcmds=lmpcmds,
    atom_types={'C': 1, 'H': 2, 'O': 3},
    atom_style='charge',   # <-- set atom_style here, not in lmpcmds
    keep_alive=False,
    log_file='lammps.log'
)

atoms.calc = calc

# Initialize velocities at combustion temperature
MaxwellBoltzmannDistribution(atoms, temperature_K=2500)

# Run short MD
dyn = VelocityVerlet(atoms, timestep=0.25 * units.fs)

def print_energy(a=atoms):
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print(f"Epot = {epot:.3f} eV/atom  Ekin = {ekin:.3f} eV/atom")

dyn.attach(print_energy, interval=10)

print("Starting short combustion test...")
dyn.run(200)  # run 200 steps
print("Done. Check species.out for product tracking.")

