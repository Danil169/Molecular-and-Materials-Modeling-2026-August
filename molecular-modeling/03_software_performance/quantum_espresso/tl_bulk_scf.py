from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
import subprocess
import os
import sys

print("📢 Setting OMP_NUM_THREADS=1 to disable OpenMP parallelization.")
os.environ["OMP_NUM_THREADS"] = "1"

def run_pw_calculation(calc, atoms, nproc=4, mpirun_cmd='mpirun'):
    """
    Run the Quantum ESPRESSO pw.x calculation in parallel.
    """
    print("Writing input files...")
    # The properties argument is required - specify what we want to calculate
    calc.write_inputfiles(atoms, properties=['energy', 'forces', 'stress'])
    
    # Check what files were created
    calc_dir = calc.directory
    print(f"Files in {calc_dir}:")
    if os.path.exists(calc_dir):
        for f in os.listdir(calc_dir):
            print(f"  - {f}")
    else:
        print(f"Error: Directory {calc_dir} was not created")
        sys.exit(1)
    
    # Find the input file - look for multiple possible extensions
    input_files = [f for f in os.listdir(calc_dir) if f.endswith('.in') or f.endswith('.pwi')]
    if not input_files:
        print("Error: No input file found in", calc_dir)
        sys.exit(1)
    
    input_file = os.path.join(calc_dir, input_files[0])
    output_file = os.path.join(calc_dir, 'pw.out')
    
    print(f"Using input file: {input_file}")
    
    # Check if pw.x exists
    pw_path = calc.profile.command.split()[0]  # Get the executable path


    
    # Run the calculation in parallel
    if nproc > 1:
        cmd = f"{mpirun_cmd} -np {nproc} {pw_path} -i {input_file} > {output_file}"
        print(f"Running in parallel with {nproc} processors")
    else:
        cmd = f"{pw_path} -i {input_file} > {output_file}"
        print("Running in serial")
    
    print(f"Command: {cmd}")
    
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

print("Setting up Thallium calculation with spin-orbit coupling...")

# 1. Set up the Thallium crystal structure
atoms = bulk('Tl', crystalstructure='hcp', a=3.46, c=5.52)
print(f"Created bulk Thallium with {len(atoms)} atoms")

# 2. Configure the Espresso calculator
pw_executable = '/home/milias/miniconda3/envs/molmatmodel/bin/pw.x'
pseudo_dir = os.getcwd()

# Create profile
profile = EspressoProfile(
    command=pw_executable,
    pseudo_dir=pseudo_dir
)

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
    profile=profile,
    directory='qe_tl_scf',
    input_data=input_data,
    pseudopotentials={'Tl': 'Tl.upf'},
    kpts=(8, 8, 8),
)

# 5. Run the calculation in parallel
NPROC = 4  # Number of processors - adjust based on your system
MPI_CMD = 'mpirun'  # Use 'mpiexec' if 'mpirun' doesn't work

run_pw_calculation(calc, atoms, nproc=NPROC, mpirun_cmd=MPI_CMD)

print("\n=== CALCULATION FINISHED ===")
print(f"Check results in: {calc.directory}")
