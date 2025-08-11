#!/bin/sh
#
#SBATCH --job-name="FPV"
#SBATCH --partition=rome
#SBATCH --time=12:00:00
#SBATCH -n 1
#SBATCH -o stdout-batch/slurm-%j-%4t-%n.out
#SBATCH -e stdout-batch/slurm-%j-%4t-%n.err
#SBATCH --mem=28G
##SBATCH --exclusive

# source ./compile/modules_snellius.sh

# mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_serial.jl")' &
srun -N1 -n1 -c1 --mem-per-cpu 8000MB --exact \
julia --project=. -O3 --check-bounds=no \
./scripts/fpv/fpvRieke/paper_initialize_sim_empty.jl

