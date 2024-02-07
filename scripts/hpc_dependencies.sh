#!/usr/bin/env sh
#SBATCH --job-name=install_deps
#SBATCH --mail-type=END
#SBATCH --mail-user=b.kopperud@lmu.de
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=logs/dependencies.log
#SBATCH --error=logs/dependencies.err
#SBATCH --qos=high_prio
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --partition=krypton

#module load R/4.2.3 gnu openblas
module load R/4.3.2 gnu openblas

export R_HOME="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R"
export LD_LIBRARY_PATH="/opt/cres/lib/hpc/gcc7/R/4.2.3/lib64/R/lib"
julia scripts/setup/dependencies.jl > logs/jldep.txt
#Rscript scripts/setup/dependencies.R > logs/Rdep.txt
#julia --threads 8 script.jl > output/screen.txt

