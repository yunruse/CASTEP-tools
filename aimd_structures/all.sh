#!/bin/bash -l

# Wallclock (h:mm:ss)
#$ -l h_rt=48:00:00

# Budget
#$ -P Free
#$ -A UKCP_ED_P

# Request nodes with RAM and TMPDIR
#$ -pe mpi 24
#$ -l mem=1G
#$ -l tmpfs=15G

#$ -N $NAME
#$ -cwd

# Load modules
module remove compilers/intel/2018/update3
module remove mpi/intel/2018/update3/intel

module load compilers/intel/2019/update4
module load mpi/intel/2019/update4/intel
module load castep/19.1.1/intel-2019

# Run CASTEP on given cell (and print time to .e file)
time gerun castep.mpi $NAME
