#!/bin/sh
#
#SBATCH --job-name="FPV"
#SBATCH --partition=rome
#SBATCH --time=12:00:00
#SBATCH -n 32
#SBATCH -o stdout-batch/slurm-%j-%4t-%n.out
#SBATCH -e stdout-batch/slurm-%j-%4t-%n.err
#SBATCH --mem=56G
##SBATCH --exclusive

# source ./compile/modules_snellius.sh

INITIAL_CASE=1
FINAL_CASE=6
for i in $(seq $INITIAL_CASE $FINAL_CASE)
do
    echo "case: $i"
    export CASE_ID=$i
    # mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_serial.jl")' &
    srun -N1 -n1 -c1 --mem-per-cpu 8000MB --exact \
	julia --project=. -O3 --check-bounds=no \
	./scripts/fpv/fpvRieke/paper3_initialize_freq_time_hourly_batch.jl &
done
wait


