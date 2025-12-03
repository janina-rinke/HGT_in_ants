#!/bin/bash
# trimming reads using Trim galore
module load palma/2024a foss/2024a parallel/20240722
ls /home/l/lfranke1/projects/LGTRNAseq/0.data/raw_reads/*.fastq.gz | parallel --dryrun 'sbatch --export=ALL,file="{}" --job-name {/.} --time 00:50:00 1.code/bash.trimgalore.sh' > 1.code/run.trimgalore.sh
bash 1.code/run.trimgalore.sh

# repare genome index (run once):
mkdir -p /home/l/lfranke1/projects/LGTRNAseq/2.pipeline/STAR/genome_index
STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /home/l/lfranke1/projects/LGTRNAseq/2.pipeline/STAR/genome_index \
  --genomeFastaFiles /cloud/wwu1/r_gadau/data/schradel/projects/LGT/RNAseq/Cobs3.1_BR_Wobs_v01_genome.fna \
  --sjdbGTFfile /cloud/wwu1/r_gadau/data/schradel/projects/LGT/RNAseq/Cobs3.1_BR_Wobs_v01_gene_annotation.gff \
  --sjdbOverhang 99

# map trimmed reads against genome index
module load palma/2024a foss/2024a parallel/20240722
ls /home/l/lfranke1/projects/LGTRNAseq/2.pipeline/trimmed_reads/*.fq.gz | parallel --dryrun 'sbatch --export=ALL,file="{}" --job-name {/.} --time 00:50:00 1.code/star.sh' > 1.code/run.star.sh
bash 1.code/run.star.sh

# quantify counts per gene
module load palma/2024a  GCC/13.3.0 parallel/20240722
ls /home/l/lfranke1/projects/LGTRNAseq/2.pipeline/STAR/alignments/*.sortedByCoord.out.bam | parallel --dryrun 'sbatch --export=ALL,file="{}" --job-name {/.} --time 00:50:00 1.code/featurecounts.sh' > 1.code/run.featurecounts.sh
cd /home/l/lfranke1/projects/LGTRNAseq
bash 1.code/run.featurecounts.sh