#$ -S /bin/bash
 #$ -N plusLGTs
 #$ -cwd
 #$ -w e
 #$ -V
 #$ -pe smp 10
 #$ -l h_vmem=2G

 conda activate dfast

 echo "Running on file: $file"

 dfast -g $file --force --cpu 10 --debug --use_original_name t --minimum_length 100 --database ~/relevant.bacteria.uni90.ref -o ./2_analysis/gene_annotation/blast_output_dfast/$file --config /home/j/j_rink02/anaconda3/envs/dfast/bin/custom_config.py