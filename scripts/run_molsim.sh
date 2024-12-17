#!/bin/bash
#SBATCH --job-name=MolSim             # Job name
#SBATCH --output=MolSim_output.log    # Standard output log
#SBATCH --error=MolSim_error.log      # Error log
#SBATCH --time=01:00:00               # Time limit (hh:mm:ss)
#SBATCH --partition=cm4_tiny          # Partition to use
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=4             # Number of CPU cores per task
#SBATCH --mem=8G                      # Memory per node
#SBATCH --export=NONE                 # Ensure clean environment
#SBATCH --get-user-env=L              # Load user environment

# Load required modules
module load gcc/13.2.0                # GCC module for compilation
module load cmake/3.26.3              # Ensure the correct version of CMake
module load xerces-c                  # Xerces-C for XML processing

# Navigate to the build directory
cd ~/MolSim/build

# Run the program
./MolSim -e 5 -d 0.0002 -i ../input/collision-of-two-bodies.txt -o ../output -l warn -n

