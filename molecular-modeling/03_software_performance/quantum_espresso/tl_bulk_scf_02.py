from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
import subprocess
import os
import sys
import platform
import re
import time
import json
from datetime import datetime
import itertools
from pathlib import Path

# ============================================
# SYSTEM INFORMATION
# ============================================

def get_total_cores():
    """
    Get the total number of available CPU cores.
    Returns the number of logical cores (threads) available.
    """
    try:
        # Method 1: Use /proc/cpuinfo (Linux)
        with open('/proc/cpuinfo', 'r') as f:
            cpuinfo = f.read()
        cores = re.findall(r'processor\s+:\s+\d+', cpuinfo)
        if cores:
            return len(cores)
    except:
        pass
    
    try:
        # Method 2: Use os.sched_getaffinity (Linux)
        return len(os.sched_getaffinity(0))
    except:
        pass
    
    try:
        # Method 3: Use multiprocessing
        import multiprocessing
        return multiprocessing.cpu_count()
    except:
        pass
    
    # Fallback
    print("⚠️  Could not detect number of cores. Using default value 1.")
    return 1

def print_system_info():
    """Print detailed system information."""
    print("=" * 80)
    print("SYSTEM INFORMATION")
    print("=" * 80)
    
    # CPU Information
    print(f"Platform: {platform.platform()}")
    print(f"Processor: {platform.processor()}")
    
    try:
        with open('/proc/cpuinfo', 'r') as f:
            cpuinfo = f.read()
        
        model_match = re.search(r'model name\s+:\s+(.+)', cpuinfo)
        if model_match:
            print(f"CPU Model: {model_match.group(1)}")
        
        cores = re.findall(r'processor\s+:\s+\d+', cpuinfo)
        total_cores = len(cores)
        print(f"Total CPU cores (logical): {total_cores}")
        
        core_ids = re.findall(r'core id\s+:\s+(\d+)', cpuinfo)
        if core_ids:
            unique_cores = len(set(core_ids))
            print(f"Physical cores: {unique_cores}")
            
            siblings = re.findall(r'siblings\s+:\s+(\d+)', cpuinfo)
            if siblings:
                print(f"Threads per core: {int(siblings[0]) // unique_cores if unique_cores > 0 else 'N/A'}")
        
        freq_match = re.search(r'cpu MHz\s+:\s+([\d.]+)', cpuinfo)
        if freq_match:
            print(f"CPU Frequency: {float(freq_match.group(1)):.0f} MHz")
        
        cache_match = re.search(r'L3\s+cache\s*:\s+(\d+)', cpuinfo)
        if cache_match:
            cache_kb = int(cache_match.group(1))
            cache_mb = cache_kb / 1024
            print(f"L3 Cache: {cache_mb:.1f} MB")
            
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
    
    print("=" * 80)
    print()

# ============================================
# READ PW.X PATH
# ============================================

def read_pw_path():
    """Read pw.x executable path from pw_exec_path.txt."""
    path_file = 'pw_exec_path.txt'
    
    if not os.path.exists(path_file):
        print(f"⚠️  {path_file} not found. Creating with default path...")
        default_path = '/home/milias/miniconda3/envs/molmatmodel/bin/pw.x'
        with open(path_file, 'w') as f:
            f.write(default_path)
        print(f"Created {path_file} with: {default_path}")
        return default_path
    
    try:
        with open(path_file, 'r') as f:
            path = f.read().strip()
        if os.path.exists(path):
            print(f"✓ Found pw.x at: {path}")
            return path
        else:
            print(f"⚠️  pw.x not found at: {path}")
            print("Please update pw_exec_path.txt with the correct path.")
            sys.exit(1)
    except Exception as e:
        print(f"Error reading {path_file}: {e}")
        sys.exit(1)

# ============================================
# PERFORMANCE BENCHMARK
# ============================================

def run_benchmark(profile, input_data, pseudopotentials, atoms, nproc, omp_threads, mpirun_cmd='mpirun'):
    """
    Run a single benchmark with specified MPI and OpenMP settings.
    """
    # Set OMP_NUM_THREADS
    os.environ["OMP_NUM_THREADS"] = str(omp_threads)
    
    print(f"\n{'='*60}")
    print(f"BENCHMARK: MPI={nproc}, OMP={omp_threads}")
    print(f"{'='*60}")
    
    # Create a unique directory for this run
    run_dir = f"qe_benchmark_mpi{nproc}_omp{omp_threads}"
    
    # Create the directory if it doesn't exist
    Path(run_dir).mkdir(exist_ok=True)
    
    # Create a new calculator for this run
    calc = Espresso(
        profile=profile,
        directory=run_dir,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kpts=(8, 8, 8),
    )
    
    start_time = time.time()
    
    try:
        # Write input files
        calc.write_inputfiles(atoms, properties=['energy', 'forces', 'stress'])
        
        # Find input file
        calc_dir = calc.directory
        input_files = [f for f in os.listdir(calc_dir) if f.endswith('.in') or f.endswith('.pwi')]
        if not input_files:
            print(f"Error: No input file found in {calc_dir}")
            return None
        
        input_file = os.path.join(calc_dir, input_files[0])
        output_file = os.path.join(calc_dir, 'pw.out')
        
        # Build command
        pw_path = profile.command.split()[0]
        
        if nproc > 1:
            cmd = f"{mpirun_cmd} -np {nproc} {pw_path} -i {input_file} > {output_file}"
        else:
            cmd = f"{pw_path} -i {input_file} > {output_file}"
        
        print(f"Command: {cmd}")
        
        # Run calculation
        subprocess.run(cmd, shell=True, check=True, capture_output=False)
        
        elapsed_time = time.time() - start_time
        
        # Extract performance data
        performance = extract_performance(output_file, elapsed_time, nproc, omp_threads)
        
        if performance and performance.get('status') == 'completed':
            print(f"✅ Benchmark completed: MPI={nproc}, OMP={omp_threads}")
            print(f"   Wall Time: {performance['wall_time']:.2f}s")
            print(f"   Efficiency: {performance['efficiency']:.2f}x")
            return performance
        else:
            print(f"❌ Failed to extract performance data")
            return None
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Calculation failed with exit code {e.returncode}")
        return None
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return None

def extract_performance(output_file, total_elapsed, nproc, omp_threads):
    """
    Extract performance metrics from pw.out file.
    Handles both MPI-only and hybrid (MPI+OpenMP) output formats.
    """
    performance = {
        'timestamp': datetime.now().isoformat(),
        'nproc': nproc,
        'omp_threads': omp_threads,
        'total_elapsed': total_elapsed,
        'status': 'unknown'
    }
    
    try:
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Try multiple patterns for PWSCF timing
        # Pattern 1: Standard format with seconds
        # "PWSCF        :     18.32s CPU     19.57s WALL"
        pattern1 = r'PWSCF\s*:\s*([\d.]+)s\s+CPU\s+([\d.]+)s\s+WALL'
        match = re.search(pattern1, content)
        
        if not match:
            # Pattern 2: Format with minutes and seconds
            # "PWSCF        :   3m 6.06s CPU   1m 7.80s WALL"
            pattern2 = r'PWSCF\s*:\s*(?:(\d+)m\s*)?([\d.]+)s\s+CPU\s+(?:(\d+)m\s*)?([\d.]+)s\s+WALL'
            match = re.search(pattern2, content)
            
            if match:
                # Parse minutes and seconds for CPU
                if match.group(1):  # CPU minutes
                    cpu_min = int(match.group(1))
                    cpu_sec = float(match.group(2))
                    cpu_time = cpu_min * 60 + cpu_sec
                else:
                    cpu_time = float(match.group(2))
                
                # Parse minutes and seconds for WALL
                if match.group(3):  # WALL minutes
                    wall_min = int(match.group(3))
                    wall_sec = float(match.group(4))
                    wall_time = wall_min * 60 + wall_sec
                else:
                    wall_time = float(match.group(4))
                
                performance['cpu_time'] = cpu_time
                performance['wall_time'] = wall_time
                performance['efficiency'] = cpu_time / wall_time if wall_time > 0 else 0
                performance['overhead'] = total_elapsed - wall_time
                
                if "JOB DONE" in content:
                    performance['status'] = 'completed'
                else:
                    performance['status'] = 'incomplete'
        
        if match:
            # We already have the data from pattern2
            pass
        elif 'match' in locals() and match:
            # Pattern1 matched
            cpu_time = float(match.group(1))
            wall_time = float(match.group(2))
            performance['cpu_time'] = cpu_time
            performance['wall_time'] = wall_time
            performance['efficiency'] = cpu_time / wall_time if wall_time > 0 else 0
            performance['overhead'] = total_elapsed - wall_time
            
            if "JOB DONE" in content:
                performance['status'] = 'completed'
            else:
                performance['status'] = 'incomplete'
        else:
            # Try a more general pattern
            pattern3 = r'PWSCF\s*:\s*([\d.]+)\s*s\s+CPU\s+([\d.]+)\s*s\s+WALL'
            match = re.search(pattern3, content)
            if match:
                cpu_time = float(match.group(1))
                wall_time = float(match.group(2))
                performance['cpu_time'] = cpu_time
                performance['wall_time'] = wall_time
                performance['efficiency'] = cpu_time / wall_time if wall_time > 0 else 0
                performance['overhead'] = total_elapsed - wall_time
                
                if "JOB DONE" in content:
                    performance['status'] = 'completed'
                else:
                    performance['status'] = 'incomplete'
            else:
                performance['status'] = 'no_timing'
                print(f"Warning: Could not find PWSCF timing in {output_file}")
        
        # Extract k-points
        kpoints_pattern = r'number of k points=\s*(\d+)'
        kpoints_match = re.search(kpoints_pattern, content)
        if kpoints_match:
            performance['kpoints'] = int(kpoints_match.group(1))
        
        # Extract bands
        bands_pattern = r'number of Kohn-Sham states=\s*(\d+)'
        bands_match = re.search(bands_pattern, content)
        if bands_match:
            performance['bands'] = int(bands_match.group(1))
        
        # Extract MPI/OpenMP info
        mpi_omp_pattern = r'Number of MPI processes:\s*(\d+).*?Threads/MPI process:\s*(\d+)'
        mpi_omp_match = re.search(mpi_omp_pattern, content, re.DOTALL)
        if mpi_omp_match:
            performance['actual_mpi'] = int(mpi_omp_match.group(1))
            performance['actual_omp'] = int(mpi_omp_match.group(2))
            
    except Exception as e:
        print(f"Error reading output file: {e}")
        performance['status'] = 'error'
    
    return performance

def save_results(results, filename='benchmark_results.json'):
    """Save benchmark results to JSON file."""
    # Load existing results if any
    existing = []
    if os.path.exists(filename):
        try:
            with open(filename, 'r') as f:
                existing = json.load(f)
        except:
            existing = []
    
    # Add new results
    existing.extend(results)
    
    # Save
    with open(filename, 'w') as f:
        json.dump(existing, f, indent=2)
    
    print(f"\n✓ Results saved to {filename}")

def print_summary(results):
    """Print a summary of all benchmark results."""
    if not results:
        print("\n❌ No results to summarize.")
        return
    
    print("\n" + "=" * 80)
    print("BENCHMARK SUMMARY")
    print("=" * 80)
    print(f"{'MPI':>6} {'OMP':>6} {'Wall Time':>12} {'Efficiency':>12} {'Status':>12}")
    print("-" * 80)
    
    for r in results:
        if r and r.get('status') == 'completed':
            wall = r.get('wall_time', 0)
            eff = r.get('efficiency', 0)
            print(f"{r['nproc']:>6} {r['omp_threads']:>6} {wall:>12.2f}s {eff:>12.2f}x {r['status']:>12}")
        elif r:
            print(f"{r.get('nproc', 0):>6} {r.get('omp_threads', 0):>6} {'N/A':>12} {'N/A':>12} {r.get('status', 'failed'):>12}")
    
    print("=" * 80)
    
    # Find best performance (fastest wall time)
    completed = [r for r in results if r and r.get('status') == 'completed']
    if completed:
        best = min(completed, key=lambda x: x['wall_time'])
        print(f"\n🏆 Fastest calculation: MPI={best['nproc']}, OMP={best['omp_threads']}")
        print(f"   Wall Time: {best['wall_time']:.2f}s")
        print(f"   Efficiency: {best['efficiency']:.2f}x")
        
        # Find best efficiency
        best_eff = max(completed, key=lambda x: x['efficiency'])
        if best_eff['nproc'] != best['nproc'] or best_eff['omp_threads'] != best['omp_threads']:
            print(f"\n🏆 Best efficiency: MPI={best_eff['nproc']}, OMP={best_eff['omp_threads']}")
            print(f"   Efficiency: {best_eff['efficiency']:.2f}x")
            print(f"   Wall Time: {best_eff['wall_time']:.2f}s")
        
        # Calculate speedup relative to serial
        serial = [r for r in completed if r['nproc'] == 1 and r['omp_threads'] == 1]
        if serial:
            serial_time = serial[0]['wall_time']
            print(f"\n📊 Speedup relative to serial (1 MPI, 1 OMP): {serial_time:.2f}s")
            for r in sorted(completed, key=lambda x: x['wall_time']):
                speedup = serial_time / r['wall_time']
                print(f"   MPI={r['nproc']:>2}, OMP={r['omp_threads']:>2}: {speedup:>6.2f}x speedup, {r['wall_time']:>8.2f}s")

# ============================================
# MAIN BENCHMARK SCRIPT
# ============================================

def main():
    """Main benchmarking function."""
    
    # Print system information
    print_system_info()
    
    # Read pw.x path
    pw_executable = read_pw_path()
    
    # Set up the Thallium crystal structure
    print("\nSetting up Thallium calculation with spin-orbit coupling...")
    atoms = bulk('Tl', crystalstructure='hcp', a=3.46, c=5.52)
    print(f"Created bulk Thallium with {len(atoms)} atoms")
    
    # Configure calculator
    pseudo_dir = os.getcwd()
    profile = EspressoProfile(
        command=pw_executable,
        pseudo_dir=pseudo_dir
    )
    
    # Input parameters
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
    
    pseudopotentials = {'Tl': 'Tl.upf'}
    
    # ============================================
    # DEFINE BENCHMARK PARAMETERS
    # ============================================
    
    # Get total available cores dynamically from the CPU
    total_cores = get_total_cores()
    
    # Define MPI process counts to test
    # Test powers of 2 up to total_cores
    mpi_values = []
    i = 1
    while i <= total_cores:
        mpi_values.append(i)
        i *= 2
    
    # Also include total_cores if it's not a power of 2
    if total_cores not in mpi_values:
        mpi_values.append(total_cores)
    
    mpi_values = sorted(mpi_values)
    
    # Define OpenMP thread counts to test
    omp_values = [1, 2, 4]
    # Filter to keep MPI * OMP <= total_cores
    combinations = []
    for mpi, omp in itertools.product(mpi_values, omp_values):
        if mpi * omp <= total_cores:
            combinations.append((mpi, omp))
    
    # Remove duplicates and sort
    combinations = list(set(combinations))
    combinations.sort()
    
    print(f"\n{'='*80}")
    print("BENCHMARK CONFIGURATION")
    print(f"{'='*80}")
    print(f"Total available cores: {total_cores}")
    print(f"Number of benchmarks to run: {len(combinations)}")
    print(f"MPI values: {mpi_values}")
    print(f"OMP values: {omp_values}")
    print(f"Combinations: {combinations}")
    print(f"{'='*80}")
    
    # ============================================
    # RUN BENCHMARKS
    # ============================================
    
    results = []
    total = len(combinations)
    
    for idx, (nproc, omp_threads) in enumerate(combinations, 1):
        print(f"\n📊 Benchmark {idx}/{total}")
        result = run_benchmark(profile, input_data, pseudopotentials, atoms, nproc, omp_threads)
        if result:
            results.append(result)
    
    # ============================================
    # SAVE AND SUMMARIZE RESULTS
    # ============================================
    
    if results:
        save_results(results)
        print_summary(results)
    else:
        print("\n❌ No benchmark results were collected.")
    
    print("\n=== BENCHMARK FINISHED ===")

if __name__ == "__main__":
    main()
