#!/bin/sh
#
#SBATCH --job-name="TFIM TEST JOB"
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=education-as-msc-ap

module load python/3.12.3

srun python3 script.py
