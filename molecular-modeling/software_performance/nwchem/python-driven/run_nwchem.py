#!/usr/bin/env python3
"""
Python script to run NWChem with multiple NPROC values from config.txt
Input file is kept untouched - modified copies are used for each run
All outputs are redirected to a dedicated directory
Extracts and summarizes wall times from NWChem output
"""

import os
import re
import sys
import subprocess
import argparse
from pathlib import Path
from datetime import datetime


class NWChemRunner:
    def __init__(self, config_file='config.txt'):
        self.config = {}
        self.input_file = None
        self.nproc_values = []
        self.nwchem_executable = 'nwchem'
        self.mpirun_executable = 'mpirun'
        self.keep_modified = False
        self.output_dir = 'results'
        self.silent = False
        self.wall_times = {}
        
        if os.path.exists(config_file):
            self.load_config(config_file)
        else:
            print(f"Warning: Configuration file '{config_file}' not found")
    
    def load_config(self, config_file):
        """Load configuration from simple key=value format"""
        try:
            with open(config_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue
                    
                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        self.config[key] = value
            
            # Parse configuration
            self.input_file = self.config.get('input')
            
            # Parse NPROC values (comma-separated)
            nproc_str = self.config.get('nproc', '4')
            self.nproc_values = [int(x.strip()) for x in nproc_str.split(',')]
            
            # Optional settings
            self.nwchem_executable = self.config.get('executable', 'nwchem')
            self.mpirun_executable = self.config.get('mpirun', 'mpirun')
            self.output_dir = self.config.get('output_dir', 'results')
            
            # Parse keep_modified (boolean)
            keep = self.config.get('keep_modified', 'false')
            self.keep_modified = keep.lower() in ['true', 'yes', '1']
            
            # Parse silent mode
            silent = self.config.get('silent', 'false')
            self.silent = silent.lower() in ['true', 'yes', '1']
            
            print(f"Configuration loaded from {config_file}")
            
        except Exception as e:
            print(f"Error loading configuration: {e}")
            sys.exit(1)
    
    def create_modified_input(self, input_file, output_file, nproc):
        """
        Create a modified copy of the input file without the mpirun command
        The original input file remains untouched
        """
        try:
            with open(input_file, 'r') as f:
                content = f.read()
            
            # Remove any mpirun lines from the content
            lines = content.split('\n')
            cleaned_lines = []
            for line in lines:
                # Skip lines that contain mpirun commands
                if re.search(r'mpirun\s+-np\s+\d+\s+nwchem', line):
                    continue
                cleaned_lines.append(line)
            
            # Add the correct mpirun command at the end
            modified_content = '\n'.join(cleaned_lines)
            if not modified_content.endswith('\n'):
                modified_content += '\n'
            modified_content += f'\n# mpirun -np {nproc} nwchem {os.path.basename(input_file)}\n'
            
            # Write the modified content
            with open(output_file, 'w') as f:
                f.write(modified_content)
            
            return True
        
        except Exception as e:
            print(f"Error creating modified file: {e}")
            return False
    
    def extract_wall_time(self, output_file):
        """Extract wall time from NWChem output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Look for "Total times" line
            pattern = r'Total times\s+cpu:\s+([\d.]+)s\s+wall:\s+([\d.]+)s'
            match = re.search(pattern, content)
            
            if match:
                cpu_time = float(match.group(1))
                wall_time = float(match.group(2))
                return cpu_time, wall_time
            
            # Alternative pattern if not found
            pattern2 = r'Task\s+times\s+cpu:\s+([\d.]+)s\s+wall:\s+([\d.]+)s'
            match2 = re.search(pattern2, content)
            
            if match2:
                cpu_time = float(match2.group(1))
                wall_time = float(match2.group(2))
                return cpu_time, wall_time
            
            return None, None
            
        except Exception as e:
            print(f"Error extracting wall time from {output_file}: {e}")
            return None, None
    
    def run_nwchem(self, input_file, nproc):
        """Run NWChem with mpirun for a specific NPROC value"""
        if not os.path.exists(input_file):
            print(f"Error: Input file '{input_file}' not found")
            return 1
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Generate output file name
        base_name = Path(input_file).stem
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_filename = f"{base_name}_nproc{nproc}_{timestamp}.out"
        output_file = os.path.join(self.output_dir, out_filename)
        
        # Build the command - use the modified input file
        cmd = [
            self.mpirun_executable,
            '-np', str(nproc),
            self.nwchem_executable,
            input_file
        ]
        
        if not self.silent:
            print(f"\n{'='*60}")
            print(f"Running with NPROC={nproc}")
            print(f"Command: {' '.join(cmd)}")
            print(f"Output: {output_file}")
            print(f"{'='*60}\n")
        
        try:
            with open(output_file, 'w') as f:
                # Write header
                f.write(f"NWChem Run\n")
                f.write(f"Date: {datetime.now()}\n")
                f.write(f"Command: {' '.join(cmd)}\n")
                f.write(f"NPROC: {nproc}\n")
                f.write(f"{'='*60}\n\n")
                
                # Run the process - redirect all output to file
                process = subprocess.Popen(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True
                )
                
                # Wait for process to complete
                return_code = process.wait()
                
                if return_code == 0:
                    # Extract wall time from output file
                    cpu_time, wall_time = self.extract_wall_time(output_file)
                    if cpu_time is not None and wall_time is not None:
                        self.wall_times[nproc] = {
                            'cpu': cpu_time,
                            'wall': wall_time,
                            'output_file': output_file,
                            'status': 'SUCCESS'
                        }
                        if not self.silent:
                            print(f"✓ NWChem completed successfully with NPROC={nproc} (wall: {wall_time:.2f}s)")
                    else:
                        self.wall_times[nproc] = {
                            'cpu': None,
                            'wall': None,
                            'output_file': output_file,
                            'status': 'SUCCESS'
                        }
                        if not self.silent:
                            print(f"✓ NWChem completed successfully with NPROC={nproc} (wall time not found)")
                else:
                    self.wall_times[nproc] = {
                        'cpu': None,
                        'wall': None,
                        'output_file': output_file,
                        'status': f'FAILED (code {return_code})'
                    }
                    if not self.silent:
                        print(f"✗ NWChem exited with error code: {return_code} for NPROC={nproc}")
                
                return return_code
        
        except FileNotFoundError as e:
            print(f"Error: Executable not found - {e}")
            print(f"Make sure {self.nwchem_executable} and {self.mpirun_executable} are in your PATH")
            return 1
        except Exception as e:
            print(f"Error running NWChem: {e}")
            return 1
    
    def print_summary(self):
        """Print summary of wall times"""
        if not self.wall_times:
            print("\nNo timing data available")
            return
        
        print(f"\n{'='*60}")
        print("PERFORMANCE SUMMARY")
        print(f"{'='*60}")
        print(f"{'NPROC':<10} {'Wall Time (s)':<15} {'Speedup':<12} {'Efficiency':<12} {'Status':<15}")
        print(f"{'-'*60}")
        
        # Sort by NPROC
        sorted_nproc = sorted(self.wall_times.keys())
        
        # Find base time (first NPROC value)
        base_time = None
        for nproc in sorted_nproc:
            if self.wall_times[nproc]['wall'] is not None:
                base_time = self.wall_times[nproc]['wall']
                break
        
        for nproc in sorted_nproc:
            data = self.wall_times[nproc]
            wall_time = data['wall']
            status = data['status']
            
            if wall_time is not None:
                speedup = base_time / wall_time if base_time else 1.0
                efficiency = (speedup / nproc) * 100 if nproc > 0 else 0
                print(f"{nproc:<10} {wall_time:<15.2f} {speedup:<12.2f} {efficiency:<11.1f}% {status:<15}")
            else:
                print(f"{nproc:<10} {'N/A':<15} {'N/A':<12} {'N/A':<12} {status:<15}")
        
        print(f"{'='*60}")
        
        # Print best speedup
        valid_times = [(n, data['wall']) for n, data in self.wall_times.items() if data['wall'] is not None]
        if len(valid_times) > 1:
            fastest = min(valid_times, key=lambda x: x[1])
            print(f"\nFastest run: NPROC={fastest[0]} with wall time={fastest[1]:.2f}s")
            
            if base_time and fastest[0] != sorted_nproc[0]:
                speedup = base_time / fastest[1]
                print(f"Speedup over NPROC={sorted_nproc[0]}: {speedup:.2f}x")
    
    def run_all(self):
        """Run NWChem with all NPROC values"""
        if not self.input_file:
            print("Error: No input file specified in configuration")
            return 1
        
        if not os.path.exists(self.input_file):
            print(f"Error: Input file '{self.input_file}' not found")
            return 1
        
        if not self.nproc_values:
            print("Error: No NPROC values specified in configuration")
            return 1
        
        if not self.silent:
            print(f"\n{'='*60}")
            print(f"NWChem Runner")
            print(f"{'='*60}")
            print(f"Input file: {self.input_file}")
            print(f"NPROC values: {self.nproc_values}")
            print(f"Output directory: {self.output_dir}")
            print(f"Original input file will NOT be modified")
            print(f"{'='*60}\n")
        
        results = {}
        modified_files = []
        
        try:
            for nproc in self.nproc_values:
                # Create modified input file (does not modify original)
                base = Path(self.input_file).stem
                modified_file = f"{base}_nproc{nproc}.nw"
                
                if not self.silent:
                    print(f"Creating modified input: {modified_file}")
                
                if not self.create_modified_input(self.input_file, modified_file, nproc):
                    print(f"Failed to create modified input file for NPROC={nproc}")
                    continue
                
                modified_files.append(modified_file)
                
                # Run NWChem
                return_code = self.run_nwchem(modified_file, nproc)
                results[nproc] = {
                    'return_code': return_code,
                    'modified_file': modified_file
                }
            
            # Print timing summary
            self.print_summary()
            
            # Clean up modified files
            if not self.keep_modified:
                if not self.silent:
                    print(f"\nCleaning up temporary files...")
                for file in modified_files:
                    if os.path.exists(file):
                        os.remove(file)
                        if not self.silent:
                            print(f"  Removed: {file}")
            
            # Check if any runs failed
            failed = any(r['return_code'] != 0 for r in results.values())
            return 1 if failed else 0
            
        except KeyboardInterrupt:
            print("\n\nInterrupted by user")
            return 1
        except Exception as e:
            print(f"Error during execution: {e}")
            return 1


def create_sample_config():
    """Create a sample config.txt file"""
    sample = """# NWChem Configuration File
# This file contains settings for running NWChem with multiple NPROC values

# Input file (required) - this file will NOT be modified
input = ch3_zora_b3lyp_prop.nw

# NPROC values (comma-separated)
nproc = 2,4,6

# Executables (optional)
executable = nwchem
mpirun = mpirun

# Output directory (optional) - all output files go here
output_dir = results

# Keep modified input files (true/false)
keep_modified = false

# Silent mode - suppress console output (true/false)
# When true, only errors and summary are printed
silent = false
"""
    
    with open('config.txt', 'w') as f:
        f.write(sample)
    print("Sample configuration file created: config.txt")


def main():
    parser = argparse.ArgumentParser(
        description='Run NWChem with multiple NPROC values from config.txt'
    )
    parser.add_argument(
        '-c', '--config',
        default='config.txt',
        help='Configuration file (default: config.txt)'
    )
    parser.add_argument(
        '--generate-config',
        action='store_true',
        help='Generate a sample config.txt file'
    )
    parser.add_argument(
        '--silent',
        action='store_true',
        help='Suppress console output (overrides config file)'
    )
    
    args = parser.parse_args()
    
    # Generate sample config if requested
    if args.generate_config:
        create_sample_config()
        return 0
    
    # Create runner instance
    runner = NWChemRunner(args.config)
    
    # Override silent mode if specified
    if args.silent:
        runner.silent = True
    
    # Run all configurations
    return runner.run_all()


if __name__ == "__main__":
    sys.exit(main())
