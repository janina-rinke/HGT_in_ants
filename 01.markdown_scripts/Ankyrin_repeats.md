# Analysis of Ankyrin Repeat Proteins in GAGA ants

### 1 Plotting the distribution of ankyrin repeats across the GAGA phylogeny
HGT candidates were checked for ankyrin repeat domains at the following website:
```bash
http://smart.embl-heidelberg.de/
```

In case of detection of ANK domains, these proteins were incorporated in the ANK repeat analysis.

After that, all proteins [CDS] were extracted and plotted onto the GAGA phylogeny with the following R-Script:
```bash
./02.r_scripts/11_ANKs_analysis.Rmd
```

## 2.1 Metagenomic analysis based on dfast gene annotation (re-annotation)
### Use ankyrin repeats coordinates from `./02.r_scripts/11_ANKs_analysis.Rmd` to extract ANK nucleotide sequences from genome assemblies

Use the following file to extract FASTA sequences:
```bash
/global/scratch2/j_rink02/master/lgt/0_data/candidatefiles
```

Structure of the file:
```bash
GAGA-0335       Scaffold10      5599332 5600854 A0A060Q434      OTU domain-containing protein   Wolbachia endosymbiont of Cimex lectularius
GAGA-0513       Scaffold27      2132477 2133991 A0A060Q434      OTU domain-containing protein   Wolbachia endosymbiont of Cimex lectularius
GAGA-0028       Scaffold283     301904  302239  A0A176Q8J2      Uncharacterized ANK     Wolbachia endosymbiont of Laodelphax striatellus
GAGA-0114       Scaffold15      8938304 8939053 A0A176Q8J2      Uncharacterized ANK     Wolbachia endosymbiont of Laodelphax striatellus
GAGA-0463       Scaffold11      1830581 1833058 A0A176Q8J2      Uncharacterized ANK     Wolbachia endosymbiont of Laodelphax striatellus
```

Find the corresponding genome sequences:
```bash
cat /global/scratch2/j_rink02/master/lgt/2_analysis/ANKs/GAGA.ANKs.tsv | parallel --colsep '\t' "samtools faidx /global/scratch2/j_rink02/master/lgt/0_data/assemblies/{1}*.fasta {2}:{3}-{4} > /global/scratch2/j_rink02/master/lgt/2_analysis/ANKs/{1}.{2}.{3}-{4}.{5}.{6}.{7}.fa"
```

### Concatenate all ANK fasta candidate files to obtain one FASTA file for all ANK loci
Keep the filenames as sequence headers
```bash
# Using sed: 's/text_to_replace/new_text/g'
for file in *.fa; do sed -i "s/^>.*/>${file%.*}/g" "$file"; done

### The above sed command will replace every line starting (^) with ">" with the filename until ".fa" is matched.
# -i: --in-place: File will be changed
# ^.*: every line starting with ">" and everything that comes afterwards will be substituted
# >${file%.*}: replace with ">" and filename

# EXPLANATION:
# s= substitute
# g= global, the pattern will be replaced in every single occurrence
```

Concatenate the files
```bash
cat *.fa >> ANKs.DFAST.annotation.fa
```

## 2.2 Metagenomic analysis based on GAGA gene annotation 
```bash
basedir=/global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/results
echo ${basedir}

# For one locus 
## GAGA-0513.Scaffold27.2132518-2132991
locus=GAGA-0513.Scaffold27.2132518-2132991
cd $basedir/${locus}.1000.out
start=$(cat genome_reannotate_mod_intersect.gff |cut -f 4)
stop=$(cat genome_reannotate_mod_intersect.gff |cut -f 5)
scf=$(cat genome_reannotate_mod_intersect.gff |cut -f 1)
id=$(echo $locus|cut -f 1 -d ".")

# Grep scaffold from GAGA annotation files 
grep "$scf" /global/homes/jg/schradel/projects/GAGA/annotations/${id}_final_annotation_repfilt_addreannot_noparpse_representative.gff3 |

# awk -v: var=val assigns a value to a variable which will be used by awk
# -F: Field-separator
awk -v start=$start -v stop=$stop -F"\t" '{if ($4>start-5000 && $5<stop+5000) print $0}'

geneID=$(grep "$scf" /global/homes/jg/schradel/projects/GAGA/annotations/${id}_final_annotation_repfilt_addreannot_noparpse_representative.gff3 |awk -v start=$start -v stop=$stop -F"\t" '{if ($4>start-5000 && $5<stop+5000) print $0}'|awk '{if ($3=="gene") print $9}'|perl -pe 's/.*ID=(.*?)\;.*/$1/g')

grep "$geneID" /global/homes/jg/schradel/projects/GAGA/annotations/${id}_final_annotation_repfilt_addreannot_noparpse_representative.pep.fasta -A 1
```

### 2.2.1 For all ANK loci:
Write a script called `ANKs.GAGA.annotation.sh` to intersect ANK loci with GAGA annotation

```bash
#!/bin/bash
#SBATCH -N intersect.ANKs
#SBATCH -cwd
#SBATCH --mem 200M
#SBATCH -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci/tmp/ANK.loci.out
#SBATCH -e /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci/tmp/ANK.loci.err
#SBATCH -wd /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci

basedir=/global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation

# define locus
locus=${GAGAid}.${scaffold}.${start}-${stop}


# locus=GAGA-0020.Scaffold17.156258-159230
# GAGAid=GAGA-0020
# scaffold=Scaffold17
# start=156258
# stop=159230


# Keep only Ankyrin.loci.1000.out files 
echo -e ${scaffold}"\t"${start}"\t"${stop}"\t"${locus} | cp -r ${basedir}/reannotation.dfast/results/${locus}.1000.out /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci

# Intersect dfast.extended reannotation with GAGA annotation
# At the end, you will have a new gff ${locus}.1000.out/GAGA_dfast_ANKs_intersect.gff, which has all predicted CDS that overlap with the ANK reannotation locus.

# find all the entries in the "genome_reannotate_mod_intersect" gff that overlap with the respective GAGA annotations.
# -a: File A to intersect with B
# -b: File B to be intersected with A
bedtools intersect -wb \
-a /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci/${locus}.1000.out/genome_reannotate_mod_intersect.gff \
-b /global/homes/jg/schradel/projects/GAGA/annotations/${GAGAid}_final_annotation_repfilt_addreannot_noparpse_representative.gff3 > /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci/${locus}.1000.out/GAGA_dfast_ANKs_intersect.gff

```

Submit bash script:
```bash
#parallel --dryrun --max-args 1 qsub -v GAGAid=GAGA-0020 -v scaffold=Scaffold17 -v start=156258 -v stop=159230 ANKs.GAGA.annotation.sh
cat ANKs.loci.txt | parallel --col-sep "\t" --dryrun --max-args 1 qsub -v GAGAid={1} -v scaffold={2} -v start={3} -v stop={4} ANKs.GAGA.annotation.sh
```

#### Quality check of all files:

Check if all files have completed the dfast re-annotation:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/genome_reannotate_mod_intersect.gff ]] && echo "$dir"; done
```

Grep fasta sequence from newly intersected files:
```bash

```