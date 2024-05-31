#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --job-name=FabF_WT_dimer_C4_1
#SBATCH --output=FabF_WT_C4_dimer_1.%j.out
#SBATCH --mail-type=ALL

# Load gromacs
module purge
ml gcc openmpi gromacs

# Make sure you are in the right directory
cd /projects/anth4580/KS_Mutants/FabF_WT_C4

# Run production run
gmx mdrun -deffnm md_0_50
