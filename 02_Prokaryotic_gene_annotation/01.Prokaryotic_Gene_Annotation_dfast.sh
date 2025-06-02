#!/usr/bin/bash


# Prepare HGT candidates file to run prokaryotic gene annotation with DFAST
cat GAGA.LGTs.allcoordinates.tsv | parallel --colsep '\t' "samtools faidx {1}*.fasta {2}:{3}-{4} > ./2_analysis/candidates.fasta/{1}.{2}.{3}-{4}.fa"

# Prokaryotic gene annotation with DFAST
## DFAST installation
conda create --name dfast python=3.8

conda activate dfast

conda install -c bioconda dfast=1.2.14

## Download databases for DFAST
dfast_file_downloader.py --protein dfast
dfast_file_downloader.py --cdd Cog --hmm TIGR

## Run DFAST
dfast --genome all.candidates.fa  --force --minimum_length 100 --metagenome \
-o ./2_analysis/gene_annotation/dfast

# Run DFAST as batch job for all HGTs
find . -name "GAGA*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh
find . -name "NCBI*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh
find . -name "OUT*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh



# Prokaryotic gene annotation with Prodigal
## Load module for Prodigal
source /usr/share/modules/init/bash  # enables the module package
module use /global/projects/programs/modules/
module avail   # list all available modules
module load annotation/prodigal/2.6.3 

## Run Prodigal
prodigal -p meta -i all.candidates.fa -a protein.translations.faa \
-d annotated_genes_prodigal.cds -w prodigal.statistics -f gff -o \
./2_analysis/gene_annotation/prodigal/prodigal_annotated_genes.gff