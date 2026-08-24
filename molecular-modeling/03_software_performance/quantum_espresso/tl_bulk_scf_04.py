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
import shutil

# ============================================
# SYSTEM INFORMATION
# ============================================

def get_total_cores(use_physical=False):
    """
    Get the total number of available CPU cores.
    
    Args:
        use_physical: If True, returns physical cores. If False, returns logical cores.
    """
    try:
        # Method 1: Use /proc/cpuinfo (Linux)
        with open('/proc/cpuinfo', 'r') as f:
            cpuinfo = f.read()
        
        if use_physical:
            # Get physical cores (unique core ids)
            core_ids = re.findall(r'core id\s+:\s+(\d+)', cpuinfo)
            if core_ids:
                return len(set(core_ids))
        
        # Get logical cores (processors)
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

def get_physical_cores():
    """Get the number of physical cores."""
    return get_total_cores(use_physical=True)

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
        
        # Get logical and physical cores
        cores = re.findall(r'processor\s+:\s+\d+', cpuinfo)
        total_logical = len(cores)
        print(f"Total CPU cores (logical): {total_logical}")
        
        core_ids = re.findall(r'core id\s+:\s+(\d+)', cpuinfo)
        if core_ids:
            total_physical = len(set(core_ids))
            print(f"Physical cores: {total_physical}")
            
            siblings = re.findall(r'siblings\s+:\s+(\d+)', cpuinfo)
            if siblings:
                print(f"Threads per core: {int(siblings[0]) // total_physical if total_physical > 0 else 'N/A'}")
        
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
# READ CONFIGURATION FILE
# ============================================

def read_config():
    """
    Read configuration from pw_configure.txt.
    
    Format:
    Line 1: Full path to pw.x executable
    Line 2: MPI configuration (optional)
            - "MPI all" - use all available MPI values
            - "MPI 2 4 6" - use only specified MPI values
            - If missing, defaults to "MPI all"
    Line 3: OMP configuration (optional)
            - "OMP all" - use all available OMP values
            - "OMP 1 2 4" - use only specified OMP values
            - If missing, defaults to "OMP all"
    """
    config_file = 'pw_configure.txt'
    
    # Default configuration
    config = {
        'pw_path': None,
        'mpi_values': None,  # None means "all"
        'omp_values': None,  # None means "all"
    }
    
    if not os.path.exists(config_file):
        print(f"⚠️  {config_file} not found. Creating with default configuration...")
        # Try to find pw.x
        pw_path = find_pw_path_auto()
        if pw_path:
            with open(config_file, 'w') as f:
                f.write(f"{pw_path}\n")
                f.write("MPI all\n")
                f.write("OMP all\n")
            print(f"Created {config_file} with default configuration.")
            config['pw_path'] = pw_path
            return config
        else:
            print(f"Created {config_file} with placeholder. Please edit it with your pw.x path.")
            with open(config_file, 'w') as f:
                f.write("/path/to/pw.x\n")
                f.write("MPI all\n")
                f.write("OMP all\n")
            print(f"Created {config_file} with placeholder. Please edit it with your pw.x path.")
            return config
    
    try:
        with open(config_file, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        
        # Line 1: pw.x path
        if len(lines) >= 1:
            config['pw_path'] = lines[0]
        
        # Line 2: MPI configuration
        if len(lines) >= 2:
            parts = lines[1].split()
            if parts[0].upper() == 'MPI':
                if len(parts) > 1 and parts[1].upper() != 'ALL':
                    config['mpi_values'] = [int(x) for x in parts[1:] if x.isdigit()]
                else:
                    config['mpi_values'] = None  # Use all
            else:
                print(f"⚠️  Invalid MPI line: {lines[1]}. Using default (all).")
        
        # Line 3: OMP configuration
        if len(lines) >= 3:
            parts = lines[2].split()
            if parts[0].upper() == 'OMP':
                if len(parts) > 1 and parts[1].upper() != 'ALL':
                    config['omp_values'] = [int(x) for x in parts[1:] if x.isdigit()]
                else:
                    config['omp_values'] = None  # Use all
            else:
                print(f"⚠️  Invalid OMP line: {lines[2]}. Using default (all).")
        
        # Validate pw.x path
        if config['pw_path'] and os.path.exists(config['pw_path']) and os.access(config['pw_path'], os.X_OK):
            print(f"✓ Found pw.x at: {config['pw_path']} (from {config_file})")
        else:
            print(f"⚠️  pw.x not found at: {config['pw_path']}")
            print("Trying to find pw.x automatically...")
            pw_path = find_pw_path_auto()
            if pw_path:
                config['pw_path'] = pw_path
                print(f"✓ Found pw.x at: {pw_path}")
                # Update config file
                with open(config_file, 'r') as f:
                    lines = f.readlines()
                lines[0] = f"{pw_path}\n"
                with open(config_file, 'w') as f:
                    f.writelines(lines)
            else:
                print(f"❌ Could not find pw.x. Please update {config_file} with the correct path.")
                return config
        
        return config
        
    except Exception as e:
        print(f"Error reading {config_file}: {e}")
        return config

def find_pw_path_auto():
    """Auto-detect pw.x path."""
    # Check if pw.x is in PATH
    pw_in_path = shutil.which('pw.x')
    if pw_in_path:
        return pw_in_path
    
    # Check common locations
    common_locations = [
        '/usr/local/bin/pw.x',
        '/usr/bin/pw.x',
        '/opt/quantum_espresso/bin/pw.x',
        os.path.expanduser('~/qe/bin/pw.x'),
        os.path.expanduser('~/quantum_espresso/bin/pw.x'),
        os.path.expanduser('~/miniconda3/envs/mace_env/bin/pw.x'),
        os.path.expanduser('~/miniconda3/envs/molmatmodel/bin/pw.x'),
        os.path.expanduser('~/miniconda3/envs/venv/bin/pw.x'),
    ]
    
    for loc in common_locations:
        if os.path.exists(loc) and os.access(loc, os.X_OK):
            return loc
    
    return None

def get_mpi_values_from_config(config, physical_cores, logical_cores):
    """
    Get MPI values based on configuration.
    """
    if config.get('mpi_values') is not None:
        # Use specified values, but filter to ensure they're valid
        mpi_values = [v for v in config['mpi_values'] if v <= logical_cores and v > 0]
        if not mpi_values:
            print("⚠️  No valid MPI values specified. Using default.")
            return None
        return sorted(mpi_values)
    return None  # Use all

def get_omp_values_from_config(config, logical_cores):
    """
    Get OMP values based on configuration.
    """
    if config.get('omp_values') is not None:
        # Use specified values, but filter to ensure they're valid
        omp_values = [v for v in config['omp_values'] if v <= logical_cores and v > 0]
        if not omp_values:
            print("⚠️  No valid OMP values specified. Using default.")
            return None
        return sorted(omp_values)
    return None  # Use all

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
        
        # Add --oversubscribe flag if nproc > physical cores
        oversubscribe = ""
        physical_cores = get_physical_cores()
        if nproc > physical_cores:
            oversubscribe = "--oversubscribe"
            print(f"⚠️  Using {oversubscribe} because MPI={nproc} > physical cores={physical_cores}")
        
        if nproc > 1:
            cmd = f"{mpirun_cmd} {oversubscribe} -np {nproc} {pw_path} -i {input_file} > {output_file}"
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
            if performance:
                print(f"   Status: {performance.get('status')}")
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
        
        # Check if JOB DONE is present
        if "JOB DONE" not in content:
            performance['status'] = 'incomplete'
            print(f"Warning: JOB DONE not found in {output_file}")
            return performance
        
        # Multiple patterns for PWSCF timing
        patterns = [
            r'PWSCF\s*:\s*([\d.]+)s\s+CPU\s+([\d.]+)s\s+WALL',
            r'PWSCF\s*:\s*(?:(\d+)m\s*)?([\d.]+)s\s+CPU\s+(?:(\d+)m\s*)?([\d.]+)s\s+WALL',
            r'PWSCF\s*:\s*([\d.]+)\s*s\s+CPU\s+([\d.]+)\s*s\s+WALL',
            r'PWSCF\s+:\s+([\d.]+)s\s+CPU\s+([\d.]+)s\s+WALL',
        ]
        
        match = None
        for pattern in patterns:
            match = re.search(pattern, content)
            if match:
                break
        
        if match:
            if len(match.groups()) == 4:
                # Pattern with minutes and seconds
                if match.group(1):  # CPU minutes
                    cpu_min = int(match.group(1))
                    cpu_sec = float(match.group(2))
                    cpu_time = cpu_min * 60 + cpu_sec
                else:
                    cpu_time = float(match.group(2))
                
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
                performance['status'] = 'completed'
                
            elif len(match.groups()) == 2:
                # Pattern with just seconds
                cpu_time = float(match.group(1))
                wall_time = float(match.group(2))
                performance['cpu_time'] = cpu_time
                performance['wall_time'] = wall_time
                performance['efficiency'] = cpu_time / wall_time if wall_time > 0 else 0
                performance['overhead'] = total_elapsed - wall_time
                performance['status'] = 'completed'
            else:
                performance['status'] = 'no_timing'
        else:
            # Fallback: look for any line with CPU and WALL
            lines = content.split('\n')
            for line in lines:
                if 'CPU' in line and 'WALL' in line and 's' in line and 'PWSCF' in line:
                    numbers = re.findall(r'([\d.]+)s', line)
                    if len(numbers) >= 2:
                        performance['cpu_time'] = float(numbers[0])
                        performance['wall_time'] = float(numbers[1])
                        performance['efficiency'] = float(numbers[0]) / float(numbers[1]) if float(numbers[1]) > 0 else 0
                        performance['overhead'] = total_elapsed - float(numbers[1])
                        performance['status'] = 'completed'
                        break
            
            if performance['status'] != 'completed':
                performance['status'] = 'no_timing'
        
        # Extract additional info
        kpoints_pattern = r'number of k points=\s*(\d+)'
        kpoints_match = re.search(kpoints_pattern, content)
        if kpoints_match:
            performance['kpoints'] = int(kpoints_match.group(1))
        
        bands_pattern = r'number of Kohn-Sham states=\s*(\d+)'
        bands_match = re.search(bands_pattern, content)
        if bands_match:
            performance['bands'] = int(bands_match.group(1))
        
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
    existing = []
    if os.path.exists(filename):
        try:
            with open(filename, 'r') as f:
                existing = json.load(f)
        except:
            existing = []
    
    existing.extend(results)
    
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
    
    completed = [r for r in results if r and r.get('status') == 'completed']
    if completed:
        best = min(completed, key=lambda x: x['wall_time'])
        print(f"\n🏆 Fastest calculation: MPI={best['nproc']}, OMP={best['omp_threads']}")
        print(f"   Wall Time: {best['wall_time']:.2f}s")
        print(f"   Efficiency: {best['efficiency']:.2f}x")
        
        best_eff = max(completed, key=lambda x: x['efficiency'])
        if best_eff['nproc'] != best['nproc'] or best_eff['omp_threads'] != best['omp_threads']:
            print(f"\n🏆 Best efficiency: MPI={best_eff['nproc']}, OMP={best_eff['omp_threads']}")
            print(f"   Efficiency: {best_eff['efficiency']:.2f}x")
            print(f"   Wall Time: {best_eff['wall_time']:.2f}s")
        
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
    
    # Read configuration
    config = read_config()
    
    if not config.get('pw_path') or not os.path.exists(config['pw_path']):
        print("❌ No valid pw.x path found. Please update pw_configure.txt")
        sys.exit(1)
    
    pw_executable = config['pw_path']
    
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
           # 'ecutwfc': 60.0,
            'ecutwfc': 40.0,
           # 'ecutrho': 480.0,
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
    
    # Get physical and logical cores
    physical_cores = get_physical_cores()
    logical_cores = get_total_cores(use_physical=False)
    
    print(f"\n{'='*80}")
    print(f"CPU CORE INFORMATION")
    print(f"{'='*80}")
    print(f"Physical cores: {physical_cores}")
    print(f"Logical cores (with hyper-threading): {logical_cores}")
    print(f"{'='*80}")
    
    # Determine MPI values
    mpi_values = get_mpi_values_from_config(config, physical_cores, logical_cores)
    if mpi_values is None:
        # Use default: powers of 2 up to physical cores
        mpi_values = []
        i = 1
        while i <= physical_cores:
            mpi_values.append(i)
            i *= 2
        if physical_cores not in mpi_values:
            mpi_values.append(physical_cores)
        mpi_values = sorted(mpi_values)
    else:
        # Filter to ensure MPI values don't exceed logical cores
        mpi_values = [v for v in mpi_values if v <= logical_cores]
        if not mpi_values:
            print("❌ No valid MPI values after filtering.")
            sys.exit(1)
    
    # Determine OMP values
    omp_values = get_omp_values_from_config(config, logical_cores)
    if omp_values is None:
        # Use default: [1, 2, 4]
        omp_values = [1, 2, 4]
    else:
        # Filter to ensure OMP values don't exceed logical cores
        omp_values = [v for v in omp_values if v <= logical_cores]
        if not omp_values:
            print("❌ No valid OMP values after filtering.")
            sys.exit(1)
    
    # Create combinations where MPI * OMP <= logical_cores
    combinations = []
    for mpi, omp in itertools.product(mpi_values, omp_values):
        if mpi * omp <= logical_cores:
            combinations.append((mpi, omp))
    
    # Remove duplicates and sort
    combinations = list(set(combinations))
    combinations.sort()
    
    print(f"\n{'='*80}")
    print("BENCHMARK CONFIGURATION")
    print(f"{'='*80}")
    print(f"Physical cores: {physical_cores}")
    print(f"Logical cores: {logical_cores}")
    print(f"Number of benchmarks to run: {len(combinations)}")
    print(f"MPI values: {mpi_values}")
    print(f"OMP values: {omp_values}")
    print(f"Combinations: {combinations}")
    if any(mpi > physical_cores for mpi, _ in combinations):
        print(f"Note: Some combinations will use --oversubscribe")
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
