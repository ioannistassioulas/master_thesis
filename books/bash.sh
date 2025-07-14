#!/bin/sh
#
#SBATCH --job-name="job_name"
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G
#SBATCH --account=research-<faculty>-<department>

module load 2023r1
module load intel/oneapi-all

export I_MPI_PMI_LIBRARY=/cm/shared/apps/slurm/current/lib64/libpmi2.so

srun ./executable