#!/bin/bash

#SBATCH -J /work/06037/tg853369/stampede2/feasst/tutorial/9_co2/4_surface/TRAPPE_molecules/alkenes/ethylene
#SBATCH -o feasst_ar-co2.o%j    # Name of stdout output file
#SBATCH -e feasst_ar-co2.e%j    # Name of stderr error file
#SBATCH -p normal                # Use the KNL Partition
#SBATCH -N 1                    # Total # of nodes (must be 1 for OpenMP)
#SBATCH --ntasks-per-node=64
#SBATCH -t 48:00:00             # Run time (hh:mm:ss)
#SBATCH --mail-user=chr218@lehigh.edu
#SBATCH --mail-type=all         # Send email at begin and end of job

# Speak to me...
pwd
date

# Stampede2 Config:
#  KNL Nodes: 68 cores @ 4 hardware threads/core -> 272 cores total
#    1.4GHz; 96GB DDR4
#  SKX Nodes: 48 cores in two CPUs @ 2 hardware threads/core -> 96 cores total
#    2.1GHz nominal (1.4-3.7GHz depending on load/instruction)

# Set thread count (have to find a way to do this elegantly)
export OMP_NUM_THREADS=64   #run with 48 cores

# Launch FEASST
$HOME/WORK/feasst/tools/run.sh test.py

# Finish
date
