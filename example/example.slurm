#!/bin/bash
#
#SBATCH --job-name=example
#SBATCH --output=example.log
#
#SBATCH --ntasks=90
#SBATCH --cpus-per-task=1
# don't use this SBATCH --time=10:00
#SBATCH --partition=my-mep

# This is needed to run intel MPI from slurm
# see https://software.intel.com/en-us/articles/how-to-use-slurm-pmi-with-the-intel-mpi-library-for-linux
# and https://slurm.schedmd.com/mpi_guide.html#intel_srun
# (note, this is now set in .bashrc, so don't need it here, but leaving to document)
# export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# These variables allow system load to show a decrease as processes finish their jobs
export I_MPI_THREAD_YIELD=3
export I_MPI_THREAD_SLEEP=100

# these directory references likely need to be changed for your system, but this assumes
# everything starts in src directory, which is likely wrong.
cd ./example
srun  ../src/AutoHetDet_0D_Opt example

