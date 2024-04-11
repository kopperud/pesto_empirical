#!/usr/bin/env sh
#SBATCH --job-name=rate_shift_type
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=logs/rate_shift_type.log
#SBATCH --error=logs/rate_shift_type.err
#SBATCH --qos=low_prio_res
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=100
#SBATCH --partition=krypton

#module load R/4.2.3 gnu openblas
module load R/4.3.2 gnu openblas

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"

julia --threads ${SLURM_CPUS_PER_TASK} scripts/06_rate_shift_type_inference.jl > logs/rate_shift_type_inference.txt
