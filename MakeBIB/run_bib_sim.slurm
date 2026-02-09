#!/bin/bash
#SBATCH --job-name=bib_sim
#SBATCH --account=m5197
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=200
#SBATCH --mem=450GB
#SBATCH --time=10:00:00
#SBATCH --output=bib_sim_%j.out
#SBATCH --error=bib_sim_%j.err

# Get particle type from command line argument
PARTICLE_TYPE=${1:-"BOTH"}
echo "Processing particle type: $PARTICLE_TYPE"

# Load environment and run simulation inside container
echo "Loading environment..."
shifter --image=docker:ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main -- /bin/bash <<EOF
set -e

# Load environment
source /opt/setup_mucoll.sh
source /global/cfs/cdirs/atlas/sferrar2/TrackMuC/mucoll-benchmarks/k4MuCPlayground/setup_digireco.sh /global/cfs/cdirs/atlas/sferrar2/TrackMuC/mucoll-benchmarks/ MAIA_v0

# Set working directory
cd /global/cfs/cdirs/atlas/sferrar2/BIB

# Create output directories
mkdir -p MUMINUS MUPLUS

echo "Starting BIB simulation processing..."
echo "Job started at: \$(date)"
echo "Using $SLURM_CPUS_PER_TASK CPU cores"

# Function to process a single file
process_file() {
    local particle=\$1
    local file_num=\$2
    local input_dir="/global/cfs/projectdirs/atlas/arastogi/MuonCollider/data/FLUKA/FLUKA_to_SLCIO/\${particle}"
    local output_dir="\${particle}"
    
    local input_file="\${input_dir}/summary\${file_num}_DET_IP.slcio"
    local output_file="\${output_dir}/summary\${file_num}_DET_IP.edm4hep.root"
    
    # Check if input file exists
    if [[ ! -f "\$input_file" ]]; then
        echo "WARNING: Input file not found: \$input_file"
        return 1
    fi
    
    # Skip if output already exists
    if [[ -f "\$output_file" ]]; then
        echo "SKIP: Output already exists: \$output_file"
        return 0
    fi
    
    # Run ddsim
    echo "Processing \${particle} file \${file_num}..."
    ddsim \\
        --steeringFile ../TrackMuC/mucoll-benchmarks/simulation/steer_baseline.py \\
        --numberOfEvents 1 \\
        --inputFiles "\$input_file" \\
        --outputFile "\$output_file"
    
    local exit_code=\$?
    if [[ \$exit_code -eq 0 ]]; then
        echo "SUCCESS: \${particle} file \${file_num} completed"
    else
        echo "ERROR: \${particle} file \${file_num} failed with exit code \$exit_code"
    fi
    
    return \$exit_code
}

# Export function for parallel execution
export -f process_file

# Create job list
echo "Generating job list..."
job_list=\$(mktemp)

# Add jobs based on particle type
if [[ "$PARTICLE_TYPE" == "MUMINUS" || "$PARTICLE_TYPE" == "BOTH" ]]; then
    echo "Adding MUMINUS jobs..."
    for i in \$(seq 1 6666); do
        echo "MUMINUS \$i" >> "\$job_list"
    done
fi

if [[ "$PARTICLE_TYPE" == "MUPLUS" || "$PARTICLE_TYPE" == "BOTH" ]]; then
    echo "Adding MUPLUS jobs..."
    for i in \$(seq 1 6666); do
        echo "MUPLUS \$i" >> "\$job_list"
    done
fi

echo "Total jobs to process: \$(wc -l < "\$job_list")"

# Run jobs in parallel using xargs
echo "Starting parallel execution with $SLURM_CPUS_PER_TASK processes..."
cat "\$job_list" | xargs -n 2 -P $SLURM_CPUS_PER_TASK -I {} bash -c 'process_file {}'

# Clean up
rm -f "\$job_list"

echo "Job completed at: \$(date)"

# Summary statistics
echo ""
echo "=== JOB SUMMARY ==="
echo "MUMINUS files processed: \$(ls MUMINUS/*.root 2>/dev/null | wc -l)"
echo "MUPLUS files processed: \$(ls MUPLUS/*.root 2>/dev/null | wc -l)"
echo "Total output files: \$(ls MUMINUS/*.root MUPLUS/*.root 2>/dev/null | wc -l)"

# Check for any failed jobs by looking for missing outputs
echo ""
echo "=== MISSING FILES CHECK ==="
missing_muminus=0
missing_muplus=0

for i in \$(seq 1 6666); do
    if [[ ! -f "MUMINUS/summary\${i}_DET_IP.edm4hep.root" ]]; then
        ((missing_muminus++))
    fi
    if [[ ! -f "MUPLUS/summary\${i}_DET_IP.edm4hep.root" ]]; then
        ((missing_muplus++))
    fi
done

echo "Missing MUMINUS files: \$missing_muminus"
echo "Missing MUPLUS files: \$missing_muplus"
echo "Total missing files: \$((missing_muminus + missing_muplus))"

if [[ \$((missing_muminus + missing_muplus)) -gt 0 ]]; then
    echo ""
    echo "=== RERUN SCRIPT FOR MISSING FILES ==="
    echo "To reprocess missing files, you can rerun this job - it will skip existing outputs."
fi
EOF