# reannotation with DFAST

## setup environment
```bash
cd ~/data/Projects/LGT/reannotation.dfast
conda activate dfast
```

# define locus
```bash
#define locus
locus="GAGA-0221.Scaffold3.8820002-8820310"
# define start, stop, scaffold, gagaID
## see https://stackoverflow.com/questions/10520623/how-to-split-one-string-into-multiple-variables-in-bash-shell

# this defines variables $GAGAid, $scaffold, $region, $start, and $stop for a given locus
IFS='.' read GAGAid scaffold region <<< "$locus"
IFS='-' read start stop <<< "$region"

# define size of the extension around LGT locus to use for the annotation
shiftSize="5000"
```

```bash
# the genome.gff_mod.gff is not required for this script to work. 
# retrieved genome.gff_mod.gff from Janina
#ls genome.gff_mod.gff
ln -s /global/scratch2/j_rink02/master/lgt/0_data/assemblies/GAGA-0221_curated_nextpolish_final_dupsrm_filt.fasta .

# this defines the fasta file to use as reference genome
genome=$(ls ${GAGAid}*.fasta)

# index fasta
samtools faidx ${genome}

# create "genome file" for bedtools
cut -f 1,2 ${genome}.fai > ${genome}.genome

```


## extract fasta +- 1000bp around locus
```bash
# create bedfile extended by 5000 bp around the locus
echo -e ${scaffold}"\t"${start}"\t"${stop}"\t"${locus} |bedtools slop -b ${shiftSize} -g ${genome}.genome > ${locus}.${shiftSize}.bed

# extract fasta
bedtools getfasta -bed ${locus}.${shiftSize}.bed -fi ${genome} -fo ${locus}.${shiftSize}.fa

# run dfast
#conda activate dfast
#dfast -g ${locus}.1kb.fa --force --metagenome --cpu 2 --debug --use_original_name t --minimum_length 100 --database /global/scratch2/databases/dfast/uniprot_bacteria-0.9.ref -o ${locus}.${shiftSize}.fa.out --config /global/scratch2/j_rink02/master/lgt/0_data/custom_config.py
dfast -g ${locus}.${shiftSize}.fa --force --metagenome --cpu 20 --debug --use_original_name t --minimum_length 100 -o ${locus} --config /global/scratch2/j_rink02/master/lgt/0_data/custom_config.py
```

## update coordinates to genome-level
```bash

#remove fasta from gff (bedtools doesn't like it)
#change first field to keep only the scaffold
# shift by original start coordinate
#shift by -n kb

sed '/^##FASTA$/,$d'  ${locus}/genome.gff| \
perl -pe 's/^('${scaffold}').*?\t/$1\t/g' | \
bedtools shift -i - -g ${genome}.genome -s ${start} -header | \
bedtools shift -i - -g ${genome}.genome -s -${shiftSize} -header \
> ${locus}/genome_reannotate_mod.gff

## check if extraction is right
#bedtools getfasta  -fi ${genome} -bed ${locus}/genome_reannotate_mod.gff -s

```

## intersect with `bedtools`

```bash
# this command retrieves all the entries from the new "expanded" gff that overlap with an entry in the old gff.
#bedtools intersect -b genome.gff_mod.gff -a ${locus}/genome_reannotate_mod.gff  > ${locus}/genome_reannotate_mod_intersect.gff

# better solution:
# find all the entries in the "expanded" gff that overlap with the LGT locus.
echo -e ${scaffold}"\t"${start}"\t"${stop}"\t"${locus}|bedtools intersect -b stdin -a ${locus}/genome_reannotate_mod.gff  > ${locus}/genome_reannotate_mod_intersect.gff
```

