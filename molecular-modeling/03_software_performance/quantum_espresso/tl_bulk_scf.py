from ase.build import bulk
from ase.calculators.espresso import Espresso
import subprocess
import os
import sys

def run_pw_parallel(calc, atoms, nproc=4, mpirun_cmd='mpirun'):
    """
    Run the Quantum ESPRESSO pw.x calculation in parallel.
    
    Parameters:
    - calc: Espresso calculator object
    - atoms: Atoms object
    - nproc: Number of processors/cores to use
    - mpirun_cmd: MPI launcher command (default: 'mpirun')
    """
    # First, write the input file
    calc.write_input(atoms)
    
    # Get the input file path
    input_file = os.path.join(calc.directory, 'pw.in')
    output_file = os.path.join(calc.directory, 'pw.out')
    
    # Construct the parallel command
    # Option 1: Using mpirun with -np flag
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
pseudopotentials = {
    'Tl': 'Tl.rel-pbe-n-kjpaw_psl.1.0.0.UPF'
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

# 4. Create the Espresso calculator
calc = Espresso(
    directory='qe_tl_scf',
    input_data=input_data,
    pseudopotentials=pseudopotentials,
    kpts=(8, 8, 8),
)

# 5. Run the calculation in parallel
# Adjust the number of processors as needed
NPROC = 8  # Change this to your desired number of processors
MPI_CMD = 'mpirun'  # Use 'mpiexec' on some systems

run_pw_parallel(calc, atoms, nproc=NPROC, mpirun_cmd=MPI_CMD)

print("Calculation finished successfully!")
