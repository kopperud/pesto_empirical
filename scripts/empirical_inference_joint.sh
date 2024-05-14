#!/usr/bin/env sh
#SBATCH --job-name=empirical_joint
#SBATCH --mail-type=END
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem-per-cpu=4GB
#SBATCH --output=logs/empirical_joint.log
#SBATCH --error=logs/empirical_joint.err
#SBATCH --qos=low_prio_res
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --partition=krypton

#module load R/4.2.3 gnu openblas
#module load openblas
module load gnu/7
#module load R/4.3.2
module load R/4.2.3

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"
echo ${SLURM_CPUS_PER_TASK} > output/ntasks.txt

julia --threads ${SLURM_CPUS_PER_TASK} scripts/empirical_inference_joint.jl > output/screen.txt

