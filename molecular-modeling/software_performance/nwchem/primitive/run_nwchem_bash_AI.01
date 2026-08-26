#!/bin/bash

# ============================================================
# NWChem Job Submission Script
# Usage: ./script.sh [NTHREADS] [PROJECT_NAME]
# ============================================================

set -e  # Exit on error

# Function to print colored output
print_status() {
    echo -e "\n\033[1;34m>>> $1 \033[0m"
}

# Function to check errors
check_error() {
    if [ $? -ne 0 ]; then
        echo -e "\033[1;31mERROR: $1\033[0m"
        exit 1
    fi
}

# Display start time
echo -e "\n=========================================="
echo "NWChem job started: $(date)"
echo "Hostname: $(hostname)"
echo "User: $(whoami)"
echo "=========================================="

# Get CPU information
if command -v nproc &> /dev/null; then
    MAXTHREADS=$(nproc)
elif [ -f /proc/cpuinfo ]; then
    MAXTHREADS=$(grep -c processor /proc/cpuinfo)
else
    MAXTHREADS=2
    echo "WARNING: Cannot detect CPU cores, defaulting to $MAXTHREADS"
fi

echo "CPU: $(cat /proc/cpuinfo | grep 'model name' | uniq | cut -d: -f2 | xargs)"
echo "Available cores: $MAXTHREADS"

# Parse command line arguments
NTHREADS=2
proj="ch3_zora_b3lyp_prop"  # Default project name

if [ $# -ge 1 ]; then
    if [[ $1 =~ ^[0-9]+$ ]]; then
        if [ $1 -ge 1 ] && [ $1 -le $MAXTHREADS ]; then
            NTHREADS=$1
            echo "Using $NTHREADS threads"
        else
            echo "WARNING: Thread count $1 out of range (1-$MAXTHREADS). Using default: $NTHREADS"
        fi
    else
        echo "WARNING: First argument not a number. Using default threads: $NTHREADS"
        if [ $# -ge 2 ]; then
            proj=$2
        fi
    fi
fi

# If second argument is provided, use it as project name
if [ $# -ge 2 ] && [[ ! $1 =~ ^[0-9]+$ ]]; then
    proj=$2
elif [ $# -ge 1 ] && [[ $1 =~ ^[0-9]+$ ]] && [ $# -ge 2 ]; then
    proj=$2
elif [ $# -ge 1 ] && [[ ! $1 =~ ^[0-9]+$ ]]; then
    proj=$1
fi

echo "Project name: $proj"

# Check NWChem executable
export NWCHEM_EXECUTABLE=${NWCHEM_EXECUTABLE:-nwchem}
if ! command -v $NWCHEM_EXECUTABLE &> /dev/null; then
    echo "ERROR: NWChem executable '$NWCHEM_EXECUTABLE' not found in PATH"
    echo "Set NWCHEM_EXECUTABLE environment variable or install NWChem"
    exit 1
fi
echo "NWChem executable: $(which $NWCHEM_EXECUTABLE)"

# Check input file
THISDIR=$PWD
INPUT_FILE="$THISDIR/$proj.nw"
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file '$INPUT_FILE' not found!"
    echo "Current directory: $THISDIR"
    echo "Available .nw files:"
    ls -la *.nw 2>/dev/null || echo "  No .nw files found"
    exit 1
fi
echo "Input file: $INPUT_FILE"

# Create scratch directory
randstr=$(openssl rand -hex 6 2>/dev/null || date +%s%N | md5sum | head -c 12)
export TMPDIR="/tmp/${USER}/nwchem_run.${randstr}"
echo "Scratch directory: $TMPDIR"

# Check /tmp space
echo "Available space in /tmp:"
df -h /tmp

# Create and enter scratch directory
mkdir -p "$TMPDIR"
check_error "Failed to create scratch directory"
cd "$TMPDIR"
echo "Working in: $(pwd)"

# Check MPI
if command -v mpirun &> /dev/null; then
    echo "MPI: $(which mpirun)"
    mpirun --version | head -1
else
    echo "WARNING: mpirun not found! Using direct execution."
fi

# Launch NWChem
print_status "Launching NWChem with $NTHREADS thread(s)"
OUTPUT_FILE="${THISDIR}/${proj}_${NTHREADS}cores_$(date +%Y%m%d_%H%M%S).out"

echo "Output file: $OUTPUT_FILE"
echo "----------------------------------------"

if command -v mpirun &> /dev/null; then
    mpirun -np $NTHREADS $NWCHEM_EXECUTABLE "$INPUT_FILE" > "$OUTPUT_FILE" 2>&1
else
    $NWCHEM_EXECUTABLE "$INPUT_FILE" > "$OUTPUT_FILE" 2>&1
fi

# Check execution status
if [ $? -eq 0 ]; then
    echo "✅ NWChem completed successfully"
else
    echo "❌ NWChem failed. Check output file: $OUTPUT_FILE"
fi

# Show files in scratch
echo -e "\nFiles in scratch directory:"
ls -lh "$TMPDIR" | head -10
if [ $(ls -1 "$TMPDIR" | wc -l) -gt 10 ]; then
    echo "... and $(($(ls -1 "$TMPDIR" | wc -l) - 10)) more files"
fi
echo "Scratch space used: $(du -sh "$TMPDIR" | cut -f1)"

# Cleanup
print_status "Cleaning up scratch directory"
cd ..
rm -rf "$TMPDIR"
check_error "Failed to remove scratch directory"

# Return to original directory
cd "$THISDIR"
echo "Returned to: $(pwd)"

# Display completion
echo -e "\n=========================================="
echo "NWChem job finished: $(date)"
echo "=========================================="
echo "Output saved to: $OUTPUT_FILE"

# Show output file size
if [ -f "$OUTPUT_FILE" ]; then
    echo "Output size: $(du -h "$OUTPUT_FILE" | cut -f1)"
    echo -e "\nLast 5 lines of output:"
    tail -5 "$OUTPUT_FILE"
fi

exit 0
