#!/bin/bash
#SBATCH --job-name=MolSim             
#SBATCH --output=MolSim_output.log   
#SBATCH --error=MolSim_error.log      
#SBATCH --get-user-env=L              
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny        
#SBATCH --cpus-per-task=112             
#SBATCH --export=NONE                
#SBATCH --time=01:00:00 

# Load required modules
module load slurm_setup
module load gcc/13.2.0                # GCC module for compilation
module load cmake/3.26.3              # Ensure the correct version of CMake
module load xerces-c                  # Xerces-C for XML processing

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to the build directory
cd ~/MolSim/build

# Run the program
./MolSim -e 0.5 -i ../input/xml_file/sheet4-task2-larger-example.xml -l off -n

