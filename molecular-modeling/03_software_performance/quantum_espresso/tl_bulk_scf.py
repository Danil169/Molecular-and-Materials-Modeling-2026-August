from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
import subprocess
import os
import sys

def run_pw_parallel(calc, atoms, nproc=4, mpirun_cmd='mpirun'):
    """
    Run the Quantum ESPRESSO pw.x calculation in parallel.
    """
    # First, write the input file
    calc.write_input(atoms)
    
    # Get the input file path
    input_file = os.path.join(calc.directory, 'pw.in')
    output_file = os.path.join(calc.directory, 'pw.out')
    
    # Construct the parallel command
    cmd = f"{mpirun_cmd} -np {nproc} pw.x -i {input_file} > {output_file}"
    
    print(f"Running: {cmd}")
    print(f"Using {nproc} processors in parallel")
    
    # Execute the command
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=False)
        print(f"Calculation completed successfully. Output written to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error: Calculation failed with exit code {e.returncode}")
        print(f"Check {output_file} for error messages")
        sys.exit(1)

# ============================================
# MAIN SCRIPT
# ============================================

# 1. Set up the Thallium crystal structure
atoms = bulk('Tl', crystalstructure='hcp', a=3.46, c=5.52)

# 2. Define the pseudopotentials
# Make sure this file exists in your pseudopotential directory
pseudopotentials = {
    'Tl': 'Tl.upf'
}

# 3. Set input parameters for a two-component calculation
input_data = {
    'system': {
        'ecutwfc': 60.0,
        'ecutrho': 480.0,
        'noncolin': True,
        'lspinorb': True,
        'occupations': 'smearing',
        'smearing': 'cold',
        'degauss': 0.01,
    },
    'electrons': {
        'mixing_beta': 0.7,
        'conv_thr': 1.0e-8,
    },
}

# 4. Configure the Espresso calculator with the correct path
# Using the path you found with 'which pw.x'
pw_executable = '/home/milias/miniconda3/envs/molmatmodel/bin/pw.x'

# Create an EspressoProfile with the executable path
profile = EspressoProfile(command=pw_executable)

# Create the Espresso calculator
calc = Espresso(
    profile=profile,           # <-- This provides the configuration
    directory='qe_tl_scf',
    input_data=input_data,
    pseudopotentials=pseudopotentials,
    kpts=(8, 8, 8),
)

# 5. Run the calculation in parallel
# Adjust the number of processors as needed
NPROC = 4  # You can increase this if you have more cores available
MPI_CMD = 'mpirun'  # Use 'mpiexec' if 'mpirun' doesn't work

run_pw_parallel(calc, atoms, nproc=NPROC, mpirun_cmd=MPI_CMD)

print("Calculation finished successfully!")
