#!/bin/sh
#
#SBATCH --job-name="MPO_First_Round"
#SBATCH --partition=compute
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=education-as-msc-ap
#SBATCH --output=/scratch/%u/${SIM_NAME}_%j/slurm_%j.out
#SBATCH --error=/scratch/%u/${SIM_NAME}_%j/slurm_%j.err

echo "Job started at $(date)"
echo "Job directory: $SLURM_SUBMIT_DIR"
echo "Scratch directory: $SCRATCH_DIR"

module load 2024r1
module load python/3.10.12

source ~/venvs/dmrg_example/bin/activate
echo "Using python from: $(which python)"

SIM_NAME="$(date '+%Y-%m-%d')::hx=(0, 2);k =(-2, 2);J = 1"
SCRATCH_DIR="/scratch/${USER}/${SIM_NAME}_${SLURM_JOB_ID}"
mkdir -p "$SCRATCH_DIR"

cd "$SCRATCH_DIR" || exit 1
echo "Successfully entered scratch"
echo "Listing files in SLURM_SUBMIT_DIR:"
ls "$SLURM_SUBMIT_DIR"
cp "$SLURM_SUBMIT_DIR/script.py" . || { echo "Failed to copy script.py"; exit 1; }

echo "Running Python script with arguments 2 2 11"
srun python3 script.py 2 1 11

echo "Job completed at $(date)"
