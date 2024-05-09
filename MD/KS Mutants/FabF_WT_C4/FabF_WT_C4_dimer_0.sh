#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --job-name=FabF_WT_dimer_C4_0
#SBATCH --output=FabF_WT_C4_dimer_0.%j.out
#SBATCH --mail-type=ALL

# Load gromacs
module purge
ml gcc openmpi gromacs

# Make sure you are in the right directory
cd /projects/anth4580/KS_Mutants/FabF_WT_C4

# Set box size
gmx editconf -f FabF_WT_dimer_C4.gro -o box.gro -bt dodecahedron -d 1.0

# Fill box with water
gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro

# Make ions tpr file
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

# Add ions
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# Make energy minimization tpr file
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

# Run energy minimization
gmx mdrun -deffnm em

# Preprocessing for NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

# Run NVT equilibration
gmx mdrun -deffnm nvt

# Preprocessing for NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr

# Run NPT equilibration
gmx mdrun -deffnm npt

# Preprocessing for production run
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_50.tpr

# Run production run
gmx mdrun -deffnm md_0_50
