#!/bin/bash

echo "LGT_pipeline.sh" $(date)

# before starting
source /usr/share/modules/init/bash  # enables the module package
module use /global/projects/programs/modules/
module load seq-search/diamond/2.0.11

### change working directory and output directory
WD=/global/scratch2/s_tret01/research_project_genomics
OUTPUT=$WD/2_pipeline/LGT_pipeline/store

### blast LGT candidate against bacterial database
# define your input file (LGT candidate protein sequence)
INPUT=$WD/0_data/external/GAGA-0363.Scaffold36.521520-521743.fa/protein.faa
diamond blastp -d /global/scratch2/databases/diamond/uniprot90_bacteria_20210707/uniprot_bacteria-0.9.dmnd -q $INPUT --outfmt 6 qseqid sseqid slen full_sseq evalue bitscore score pident stitle > $OUTPUT/blast_search.txt

### Get fasta file with sequences and headers for blast search
awk 'BEGIN{OFS="\t";FS="\t"} {print ">"$9,$3,"\n"$4}' $OUTPUT/blast_search.txt > $OUTPUT/blast_search.fa

### get top 20 blast hits
head -40 $OUTPUT/blast_search.fa > $OUTPUT/top_20_hits.fa

### combine LGT fasta file with top 20 hits
cat $INPUT $OUTPUT/top_20_hits.fa > $OUTPUT/LGT_top20.fa

# creates multiple sequence alignments of prot sequences
prank -d=$OUTPUT/LGT_top20.fa -o=$OUTPUT/msa.fa

# create phylogeny with IQtree
iqtree -s $OUTPUT/msa.fa.best.fas -nt 12 -b 100
