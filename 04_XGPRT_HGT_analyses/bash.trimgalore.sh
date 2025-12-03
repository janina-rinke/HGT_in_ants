#!/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=normal,express,requeue,r0bornb
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name=Readtrimming
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lfranke1@uni-muenster.de
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

# Load module
module load palma/2022a foss/2022a Trim_Galore/0.6.10

# Set directories
#cd /home/l/lfranke1/projects/LGTRNAseq/0.data/raw_reads
OUTDIR="/home/l/lfranke1/projects/LGTRNAseq/2.pipeline/trimmed_reads"

# Loop over all FASTQ files
echo "Trimming $file..."
trim_galore -o "$OUTDIR" "$file"

