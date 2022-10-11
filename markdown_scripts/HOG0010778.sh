#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Job name (comment out the next line to get the name of the script used as the job name)
##PBS -N test
### Output files (comment out the next 2 lines to get the job name used instead)
##PBS -e test.err
##PBS -o test.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=1:thinnode
### Memory
#PBS -l mem=8gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 20:00:00, 20 hours)
#PBS -l walltime=20:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
WORKDIR=/global/scratch2/j_rink02/master/lgt/2_analysis/selection_analysis/Etherases
echo Working directory is $WORKDIR, going to this script directory
cd $WORKDIR

# Load all required modules for the job

#module load tools
#module load ngs
#module load anaconda3/4.4.0
#module load fasttree/2.1.11
#module load perl/5.24.0
#module load trimal/1.4.1
module load msa/muscle/3.8.31
module load phylo/iqtree/2.1.3
#module load gcc/11.1.0
module load msa/prank/170427
#module load openmpi/gcc
#module load usr/bin/hyphy

#prank -d=../HOG0010778.cds.fasta -o=../PRANK/HOG0010778.cds.aln -codon -F
prank -d=etherases.fa -o=etherases.aln -codon -F
echo "checking if alignment file exists..."

#FILE=../PRANK/HOG0010778.cds.aln.best.fas
FILE=etherases.aln.best.fas
if [ -f "$FILE" ]; then
	echo "$FILE exists."
else
	echo "$FILE does not exist. Running again"
	prank -d=../HOG0010778.cds.fasta -o=../PRANK/HOG0010778.cds.aln -codon -F
fi
if [ ! -f "$FILE" ]; then
	echo "$FILE not exists." >> ../errors.txt
	exit
fi

# Check alignment quality using zorro
#/global/scratch2/j_rink02/master/programs/zorro_linux_x86_64 -sample etherases.aln.best.fas > etherases.aln.zorro.txt
#FILE1=etherases.aln.zorro.txt
#if [ -f "$FILE1" ]; then
#	echo "$FILE1 exists."
#else
#	echo "$FILE1 does not exist. Running again"
#	/home/projects/ku_00039/people/joeviz/programs/zorro_linux_x86_64 -sample ../PRANK/HOG0010778.cds.aln.best.fas > ../ZORRO/HOG0010778.cds.aln.zorro.txt
#fi
#if [ ! -f "$FILE1" ]; then
#	echo "$FILE1 not exists." >> ../errors.txt
#	exit
#fi

hyphy CPU=1 gard --alignment etherases.aln.best.fas --mode Faster --output etherases_gard.json --output-lf etherases.best-gard > etherases_gard.output
perl /home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/analyze_gard_v2_singleOG.pl etherases.aln.best.fas etherases_gard.output ../HYPHY/alignments/HOG0010778
for gard_part in ../HYPHY/alignments/HOG0010778*aln; do
iqtree2 -s "$gard_part" -B 1000 -T 1 -m MFP
echo "checking if treefile file exists..."
FILE="${gard_part}.treefile"
if [ -f "$FILE" ]; then
	echo "$FILE exists."
else
	echo "$FILE does not exist. Running again"
	iqtree2 -s "$gard_part" -B 1000 -T 1 -m MFP
fi
if [ ! -f "$FILE" ]; then
	echo "$FILE not exists." >> ../errors.txt
	exit
fi

hyphy CPU=1 absrel --alignment "$gard_part" --tree "${gard_part}.treefile" --output "${gard_part}_ABSREL.json" > "${gard_part}_ABSREL.output" 
python3 /home/projects/ku_00039/people/joeviz/scripts_ignasi/significant_aBSREL_2_table.py "${gard_part}_ABSREL.json"  "${gard_part}_ABSREL_significant_table.txt" 0.05 0.00001
python3 /home/projects/ku_00039/people/joeviz/scripts_ignasi/aBSREL_2_table.py "${gard_part}_ABSREL.json" "${gard_part}_ABSREL_table.txt"
python3 /home/projects/ku_00039/people/joeviz/scripts_ignasi/label_relax.py "${gard_part}.treefile" /home/projects/ku_00039/people/joeviz/selection_pipeline/Species_selection_Worker_ovary_sterility_clades.txt "${gard_part}.treefile.labeled"
hyphy CPU=1 relax --alignment "$gard_part" --tree "${gard_part}.treefile.labeled" --output "${gard_part}_RELAX.json" --test test --reference reference > "${gard_part}_RELAX.output" 
done
