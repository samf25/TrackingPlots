#!/bin/bash

# Show current queue status
echo "Current queue status for your jobs:"
squeue -u $USER

echo ""
read -p "Submit the BIB simulation jobs? (y/N): " -n 1 -r
echo

if [[ $REPLY =~ ^[Yy]$ ]]; then
    cd /global/cfs/cdirs/atlas/sferrar2/BIB
    
    # Submit MUMINUS job
    echo "Submitting MUMINUS job..."
    muminus_job_id=$(sbatch run_bib_sim.slurm MUMINUS | grep -o '[0-9]*')
    echo "MUMINUS job submitted with ID: $muminus_job_id"
    
    # Submit MUPLUS job
    echo "Submitting MUPLUS job..."
    muplus_job_id=$(sbatch run_bib_sim.slurm MUPLUS | grep -o '[0-9]*')
    echo "MUPLUS job submitted with ID: $muplus_job_id"
else
    echo "Job submission cancelled."
fi