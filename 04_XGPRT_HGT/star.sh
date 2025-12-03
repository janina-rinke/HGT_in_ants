#!/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=normal,express,requeue,r0bornb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=36
#SBATCH --mem=20G
#SBATCH --job-name=aligning
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lfranke1@uni-muenster.de
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

# Load module
module load palma/2024a foss/2024a STAR

# Set directories
OUTDIR="/home/l/lfranke1/projects/LGTRNAseq/2.pipeline/STAR/alignments"
# Extract just the filename without path or extension
FILENAME=$(basename "$file" .fq.gz)

# Loop over all FASTQ files
echo "alinging $file..."

  STAR \
  --runThreadN 36 \
  --genomeDir /home/l/lfranke1/projects/LGTRNAseq/2.pipeline/STAR/genome_index \
  --readFilesIn "$file" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUTDIR/${FILENAME}_" \
  --outSAMtype BAM SortedByCoordinate