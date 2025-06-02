#!/usr/bin/bash

# Extract the best uniprot hit from DFAST for each HGT candidate
ls ./|egrep "^GAGA|^OUT|^NCBI" |parallel \
"awk -v var=$1  '/##FASTA/ {exit} {if (/^[^##]/) {print var"\t"$0}}' $1/genome.gff | \
perl -pe 's/(.*?)\t(.*?)\t.*note=(.*?\])\,.*/$1\t$2\t$3/g' > LGTs_uniProt_besthit.tsv"


# Extract start- and stop codons from HGT candidates
for i in * ; do seqkit fx2tab $i/cds.fna |awk -F '\t' '{print $1,substr($2,1,3),substr($2,length($2)-2,length($2))}' \
> $i.start.stop.codons.tsv ; done

# Prints out the file directory names as lines
printf "%s\n" *.start.stop.codons.tsv | xargs -n1 -d $'\n' bash -c 'xargs -n1 -d $'\''\n'\'' printf "%s,%s\n" "$1" <"$1"' -- > all.start.stop.codons.tsv


# Count number of predicted CDS for every HGT candidate
grep -c CDS */genome.gff