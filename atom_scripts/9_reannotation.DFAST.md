# Reannotation with DFAST

This script only needs a GAGA genome as a FASTA file (path: `/global/scratch2/j_rink02/master/lgt/0_data/assemblies`) and the extended LGT locus by user-defined size.

### 1) Setup environment
```bash
cd /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/loci
```

### 2) Extend the locus by 1000 bp on each side.

Submit the job as a script called `rerun.dfast.sh` to GridEngine and loop over each row in the file `locus.txt`.

```bash
#$ -S /bin/bash
#$ -N rerunDFAST
#$ -cwd
#$ -w e
#$ -V
#$ -pe smp 30
#$ -l h_vmem=6G
#$ -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/loci/tmp/locus.out
#$ -e /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/loci/tmp/locus.err
#$ -wd /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/loci

# define size of the extension around LGT locus to use for the annotation
shiftSize="1000"

# define locus
locus=${GAGAid}.${scaffold}.${start}.${stop}

# create a softlink in the current working directory that links to the fasta index for the current GAGAid
ln -s /global/scratch2/j_rink02/master/lgt/0_data/assemblies/${GAGAid}*.fasta .

# create a softlink in the current working directory that links to the fasta index for the current GAGAid
ln -s /global/scratch2/j_rink02/master/lgt/0_data/assemblies/${GAGAid}*fasta.fai .

# define ${genome} to point to the (soft-link) genome fasta in the current directory
genome=$(ls ${GAGAid}*.fasta)

# create "genome file" for bedtools
cut -f 1,2 ${genome}.fai > ${genome}.genome

### Extract fasta +- 1000 bp around the locus
bedtools getfasta -bed ${locus}.${shiftSize}.bed -fi ${genome} -fo ${locus}.${shiftSize}.fa

# create bedfile extended by 1000 bp around the locus
echo -e ${scaffold}"\t"${start}"\t"${stop}"\t"${locus} | bedtools slop -b ${shiftSize} -g ${genome}.genome > ${locus}.${shiftSize}.bed


### Re-run dfast with the extended locus
conda activate dfast

source /usr/share/modules/init/bash  # enables the module package
module use /global/projects/programs/modules/
module load seq-search/mmseqs/sse2-13-45111

dfast -g ${locus}.${shiftSize}.fa --force --metagenome --cpu 30 --debug --use_original_name t --minimum_length 100 --database /global/scratch2/databases/dfast/uniprot_bacteria-0.9.ref -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out --config /global/scratch2/j_rink02/master/lgt/0_data/custom_config.py

### Update the coordinates to genome level

#remove fasta from gff (bedtools doesn't like it)
#change first field to keep only the scaffold
# shift by original start coordinate
#shift by -n kb

sed '/^##FASTA$/,$d' /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out/genome.gff| \
perl -pe 's/^('${scaffold}').*?\t/$1\t/g' | \
bedtools shift -i - -g ${genome}.genome -s ${start} -header | \
bedtools shift -i - -g ${genome}.genome -s -${shiftSize} -header \
> /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out/genome_reannotate_mod.gff

## check if extraction is right
#bedtools getfasta  -fi ${genome} -bed ${locus}/genome_reannotate_mod.gff -s

### Intersect with `bedtools`
# At the end, you will have a new gff ${locus}.1000.out/genome_reannotate_mod_intersect.gff, which has all predicted CDS that overlap with the LGT locus.

# find all the entries in the "expanded" gff that overlap with the LGT locus.
echo -e ${scaffold}"\t"${start}"\t"${stop}"\t"${locus} | bedtools intersect -b stdin -a /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out/genome_reannotate_mod.gff > /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out/genome_reannotate_mod_intersect.gff
```

Submit the script `rerun.dfast.sh` to GridEngine:
```bash
cat locus.txt | parallel --col-sep "\t" --dryrun --max-args 1 qsub -v GAGAid={1} -v scaffold={2} -v start={3} -v stop={4} rerun.dfast.sh
```
