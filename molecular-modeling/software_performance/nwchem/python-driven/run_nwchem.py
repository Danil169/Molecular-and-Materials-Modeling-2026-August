#!/usr/bin/env python3
"""
Python script to run NWChem with multiple NPROC values from config.txt
Input file is kept untouched - modified copies are used for each run
All outputs are redirected to a dedicated directory
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
                    if not self.silent:
                        print(f"✓ NWChem completed successfully with NPROC={nproc}")
                else:
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
            
            # Print summary
            if not self.silent:
                print(f"\n{'='*60}")
                print("Summary of runs:")
                print(f"{'='*60}")
                for nproc, result in results.items():
                    status = "✓ SUCCESS" if result['return_code'] == 0 else "✗ FAILED"
                    print(f"NPROC={nproc}: {status} (return code: {result['return_code']})")
            
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
