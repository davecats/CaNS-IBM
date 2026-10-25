#!/usr/bin/env bash
# Job script used for the reference runs on HoreKA (KIT): one node, 4 x A100, one MPI rank per GPU.
# usage: CANS_DIR=<CaNS-IBM checkout> sbatch -A <project> -p accelerated -t 48:00:00 -J <name> run_horeka.sh <run dir>
# The run directory must contain input.nml (and the restart files fld_*.bin in data/ for restart = T).
#SBATCH -N 1 --ntasks-per-node=4 --gres=gpu:4 --cpus-per-task=16
#SBATCH --distribution=block:block
#SBATCH -o %x-%j.out
set -euo pipefail
CANS_DIR=${CANS_DIR:?set CANS_DIR to the CaNS-IBM checkout (built with FCOMP=NVIDIA GPU=1)}
module purge > /dev/null 2>&1 || true
module load toolkit/nvidia-hpc-sdk/24.9
source $NVHPC_ROOT/comm_libs/12.6/hpcx/hpcx-2.20/hpcx-init.sh; hpcx_load
export LD_LIBRARY_PATH=$CANS_DIR/dependencies/cuDecomp/build/lib:$LD_LIBRARY_PATH
D=$1; cd $D; mkdir -p data
echo "job $SLURM_JOB_ID on $SLURM_JOB_NODELIST start $(date)" >> jobs.txt
cp input.nml input-$SLURM_JOB_ID.nml
mpirun -n 4 --map-by slot --bind-to none $CANS_DIR/build/cans > log-$SLURM_JOB_ID.txt 2>&1 || echo "CANS rc=$?"
echo "job $SLURM_JOB_ID end $(date)" >> jobs.txt
