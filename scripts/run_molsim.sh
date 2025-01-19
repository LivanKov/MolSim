#!/bin/bash
#SBATCH --job-name=MolSim             
#SBATCH --output=MolSim_output.log   
#SBATCH --error=MolSim_error.log      
#SBATCH --get-user-env=L              
#SBATCH --clusters=serial
#SBATCH --partition=serial_std        
#SBATCH --mem=8G 
#SBATCH --cpus-per-task=1             
#SBATCH --export=NONE                
#SBATCH --time=01:00:00 

# Load required modules
module load gcc/13.2.0                # GCC module for compilation
module load cmake/3.26.3              # Ensure the correct version of CMake
module load xerces-c                  # Xerces-C for XML processing

# Navigate to the build directory
cd ~/MolSim/build

# Run the program
./MolSim -e 0.5 -i ../input/xml_file/sheet4-task2-larger-example.xml -l off -n

