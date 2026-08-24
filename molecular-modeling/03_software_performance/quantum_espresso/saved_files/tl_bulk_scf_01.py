from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
import subprocess
import os
import sys
import platform
import re

# ============================================
# SYSTEM INFORMATION
# ============================================

print("=" * 60)
print("SYSTEM INFORMATION")
print("=" * 60)

# CPU Information
print(f"Platform: {platform.platform()}")
print(f"Processor: {platform.processor()}")

try:
    # Get CPU info on Linux
    with open('/proc/cpuinfo', 'r') as f:
        cpuinfo = f.read()
    
    # Get CPU model name
    model_match = re.search(r'model name\s+:\s+(.+)', cpuinfo)
    if model_match:
        print(f"CPU Model: {model_match.group(1)}")
    
    # Get number of physical cores
    cores = re.findall(r'processor\s+:\s+\d+', cpuinfo)
    total_cores = len(cores)
    print(f"Total CPU cores: {total_cores}")
    
    # Get number of physical cores (unique core ids)
    core_ids = re.findall(r'core id\s+:\s+(\d+)', cpuinfo)
    unique_cores = len(set(core_ids))
    print(f"Physical cores: {unique_cores}")
    
    # Check for hyper-threading
    siblings = re.findall(r'siblings\s+:\s+(\d+)', cpuinfo)
    if siblings:
        print(f"Threads per core: {int(siblings[0]) // unique_cores if unique_cores > 0 else 'N/A'}")
    
except Exception as e:
    print(f"Could not read CPU info: {e}")

# Memory Information
try:
    with open('/proc/meminfo', 'r') as f:
        meminfo = f.read()
    mem_match = re.search(r'MemTotal:\s+(\d+)', meminfo)
    if mem_match:
        mem_kb = int(mem_match.group(1))
        mem_gb = mem_kb / (1024 * 1024)
        print(f"Total Memory: {mem_gb:.1f} GB")
except Exception as e:
    print(f"Could not read memory info: {e}")

# Parallelization settings
print(f"\nOMP_NUM_THREADS: {os.environ.get('OMP_NUM_THREADS', 'Not set')}")
print("=" * 60)
print()

# ============================================
# MAIN SCRIPT
# ============================================

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
        print(f"Running in parallel with {nproc} MPI processes")
    else:
        cmd = f"{pw_path} -i {input_file} > {output_file}"
        print("Running in serial")
    
    print(f"Command: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=False)
        print(f"Calculation completed successfully. Output written to {output_file}")
        
        # Extract wall time from output file
        extract_wall_time(output_file)
        
    except subprocess.CalledProcessError as e:
        print(f"Error: Calculation failed with exit code {e.returncode}")
        print(f"Check {output_file} for error messages")
        sys.exit(1)

def extract_wall_time(output_file):
    """
    Extract wall time and CPU time from the pw.out file.
    """
    print("\n" + "=" * 60)
    print("PERFORMANCE SUMMARY")
    print("=" * 60)
    
    try:
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Look for the timing information
        # Pattern: "PWSCF        :     18.32s CPU     19.57s WALL"
        pattern = r'PWSCF\s*:\s*([\d.]+)s\s+CPU\s+([\d.]+)s\s+WALL'
        match = re.search(pattern, content)
        
        if match:
            cpu_time = float(match.group(1))
            wall_time = float(match.group(2))
            
            print(f"CPU Time:  {cpu_time:.2f} s")
            print(f"Wall Time: {wall_time:.2f} s")
            print(f"Efficiency: {cpu_time/wall_time:.2f}x" if wall_time > 0 else "Efficiency: N/A")
            
            # Check if calculation completed successfully
            if "JOB DONE" in content:
                print("✅ Calculation completed successfully (JOB DONE)")
            else:
                print("⚠️  Warning: 'JOB DONE' not found in output")
                
        else:
            print("Could not find timing information in output file.")
            print("Look for 'PWSCF' lines in the output.")
            
        # Also print the end time if available
        time_pattern = r'This run was terminated on:\s*([\d:]+)\s+(\d+\w+\d+)'
        time_match = re.search(time_pattern, content)
        if time_match:
            print(f"Run terminated at: {time_match.group(1)} {time_match.group(2)}")
            
    except FileNotFoundError:
        print(f"Output file {output_file} not found.")
    except Exception as e:
        print(f"Error reading output file: {e}")
    
    print("=" * 60)

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

# Print parallelization settings before running
print("\n" + "=" * 60)
print("PARALLELIZATION SETTINGS")
print("=" * 60)
print(f"MPI processes: {NPROC}")
print(f"OMP_NUM_THREADS: {os.environ.get('OMP_NUM_THREADS', 'Not set')}")
print("=" * 60)
print()

run_pw_calculation(calc, atoms, nproc=NPROC, mpirun_cmd=MPI_CMD)

print("\n=== CALCULATION FINISHED ===")
print(f"Check results in: {calc.directory}")
