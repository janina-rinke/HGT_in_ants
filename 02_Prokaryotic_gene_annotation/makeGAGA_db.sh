#$ -S /bin/bash
#$ -N makeGAGAdb
#$ -cwd
#$ -pe smp 20
#$ -l h_vmem=3G


makeblastdb -dbtype nucl -in ./0_data/assemblies_copy/GAGA_genomes_database