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
module load palma/2024a  GCC/13.3.0

# Set directories
OUTDIR="/scratch/tmp/lfranke1/projects/LGTRNAseq/2.pipeline/STAR/alignments"
# Extract just the filename without path or extension
FILENAME=$(basename "$file" .sortedByCoord.out.bam)

#one for Cobs
featureCounts \
  -T 36 \
  -t exon \
  -g gene \
  -a /home/l/lfranke1/projects/LGTRNAseq/0.data/Cobs3.1_BR_gene_annotation.gff \
  -o /home/l/lfranke1/projects/LGTRNAseq/2.pipeline/Read_counting/gene_counts_Cobs_$FILENAME.txt \
  "$file"

  #one for Westeberhardia
  featureCounts \
  -T 36 \
  -t CDS \
  -g Parent \
  -a /home/l/lfranke1/projects/LGTRNAseq/0.data/Wobs_gene_annotation.gff \
  -o /home/l/lfranke1/projects/LGTRNAseq/2.pipeline/Read_counting/gene_counts_Wobs_$FILENAME.txt \
  "$file"