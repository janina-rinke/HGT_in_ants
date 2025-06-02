#$ -S /bin/bash
#$ -N batchdfastjob
#$ -cwd
#$ -w e
#$ -V
#$ -pe smp 20
#$ -l h_vmem=6G

conda activate dfast

echo "Running on file: $file"

dfast -g $file --force --metagenome --cpu 20 --debug --use_original_name t --minimum_length 100 --database ./databases/dfast/uniprot_bacteria-0.9.ref -o ./2_analysis/gene_annotation/dfast/$file --config custom_config.py