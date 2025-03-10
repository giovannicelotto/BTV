#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G                         
#SBATCH --partition=short              # Specify your cluster partition
#SBATCH --time=1:00:00  
#SBATCH --output=/t3home/gcelotto/slurm/output/BTV.out  # Output file for stdout
#SBATCH --error=/t3home/gcelotto/slurm/output/BTV.out    # Output file for stderr
#SBATCH --dependency=singleton

source_dir="/scratch"
fileName="$1"
fileNumber="$2" 
flatPath="$3"
python /work/gcelotto/BTV/scripts/tuplizer/ntupleLinker_SVD.py -fileName $fileName -fileNumber $fileNumber -mE 500


xrdcp -f -N $source_dir/"TTToH_"$fileNumber.root root://t3dcachedb.psi.ch:1094//$flatPath