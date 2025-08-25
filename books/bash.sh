#!/bin/sh
#
#SBATCH --job-name="TFIM TEST JOB"
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=education-as-msc-ap

module load 2024r1
module load python/3.10.12

source ~/venvs/dmrg_example/bin/activate

SIM_NAME="Sim: h=0,2;k=-2,2;j=1"
SCRATCH_DIR="/scratch/${USER}/${SIM_NAME}_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_DIR"

cd "$SCRATCH_DIR" || exit 1
cp "$SLURM_SUBMIT_DIR/script.py" .

#SBATCH --output=/scratch/%u/${SIM_NAME}_%j/slurm_%j.out
#SBATCH --error=/scratch/%u/${SIM_NAME}_%j/slurm_%j.err

srun python3 script.py 2 2 11
