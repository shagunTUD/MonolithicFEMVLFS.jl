#!/bin/sh
#
#SBATCH --job-name="fpvVLFS"
#SBATCH --partition=compute
#SBATCH --time=96:00:00
#SBATCH -n 18
#SBATCH -o stdout-batch/slurm-%j-%4t-%n.out
#SBATCH -e stdout-batch/slurm-%j-%4t-%n.err
#SBATCH --mem=80G
#SBATCH --account=research-ceg-he

# source ./compile/modules.sh

INITIAL_CASE=1
FINAL_CASE=18
for i in $(seq $INITIAL_CASE $FINAL_CASE)
do
    echo "case: $i"
    export CASE_ID=$i
    # mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_serial.jl")' &
    srun --exclusive -N1 -n1 -c1 --mem-per-cpu 4000MB --exact julia --project=. -O3 --check-bounds=no ./scripts/fpv/fpvData2/paper2_initialize_freq_time_hourly_batch.jl &
done
wait
