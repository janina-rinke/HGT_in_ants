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
/global/scratch2/j_rink02/master/lgt/2_analysis/ANKs/GAGA.ANKs.tsv
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
```

### 2.2.1 For one locus 
```bash
locus=GAGA-0335.Scaffold10.5598799-5599854;
GAGAid=GAGA-0335; 
cp -r /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out \
/global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci/;
bedtools intersect -wa \
-a /global/homes/jg/schradel/projects/GAGA/annotations/GAGA-0335_final_annotation_repfilt_addreannot_noparpse_representative.gff3 \
-b /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci/${locus}.1000.out/genome_reannotate_mod_intersect.gff \
> /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci/${locus}.1000.out/GAGA_dfast_ANKs_intersect.gff;
id=$(cat /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci/${locus}.1000.out/GAGA_dfast_ANKs_intersect.gff|awk '{if ($3=="mRNA") print $0}'|perl -pe 's/.*ID=(.*?)\;.*/$1/g');
# id Mbic_g01422
echo ${id} > ids.txt;
seqkit grep -nrif ids.txt /global/homes/jg/schradel/projects/GAGA/annotations/${GAGAid}_final_annotation_repfilt_addreannot_noparpse_representative.cds.fasta > ${locus}.out.fa
```

### 2.2.2 For all ANK loci:

```bash
basedir=/global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation
echo ${basedir}

# ANKs.loci.txt file structure:
# GAGA-0028       Scaffold283     299552  303004

# Parallel --dryrun
cat ANKs.loci.txt | parallel --col-sep "\t" --dryrun \
"locus={1}.{2}.{3}-{4}; 
cp -r ${basedir}/reannotation.dfast/results/{1}.{2}.{3}-{4}.1000.out /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci;
bedtools intersect -wa \
-a /global/homes/jg/schradel/projects/GAGA/annotations/{1}_final_annotation_repfilt_addreannot_noparpse_representative.gff3 \
-b ${basedir}/ANKs/loci/{1}.{2}.{3}-{4}.1000.out/genome_reannotate_mod_intersect.gff \
> ${basedir}/ANKs/loci/{1}.{2}.{3}-{4}.1000.out/GAGA_dfast_ANKs_intersect.gff;
cat ${basedir}/ANKs/loci/{1}.{2}.{3}-{4}.1000.out/GAGA_dfast_ANKs_intersect.gff|awk '{if (\$3==\"mRNA\") print \$0}'|perl -pe 's/.*ID=(.*?)\;.*/\$1/g' >> ids.txt;
seqkit grep -nrif ids.txt /global/homes/jg/schradel/projects/GAGA/annotations/{1}_final_annotation_repfilt_addreannot_noparpse_representative.cds.fasta > {1}.{2}.{3}-{4}.out.fa"


# Use Parallel on all ankyrin files
cat ANKs.loci.txt | parallel --col-sep "\t" "locus={1}.{2}.{3}-{4}; 
cp -r ${basedir}/reannotation.dfast/results/{1}.{2}.{3}-{4}.1000.out /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/ANKs/loci;
bedtools intersect -wa \
-a /global/homes/jg/schradel/projects/GAGA/annotations/{1}_final_annotation_repfilt_addreannot_noparpse_representative.gff3 \
-b ${basedir}/ANKs/loci/{1}.{2}.{3}-{4}.1000.out/genome_reannotate_mod_intersect.gff \
> ${basedir}/ANKs/loci/{1}.{2}.{3}-{4}.1000.out/GAGA_dfast_ANKs_intersect.gff;
cat ${basedir}/ANKs/loci/{1}.{2}.{3}-{4}.1000.out/GAGA_dfast_ANKs_intersect.gff|awk '{if (\$3==\"mRNA\") print \$0}'|perl -pe 's/.*ID=(.*?)\;.*/\$1/g' >> ids.txt;
seqkit grep -nrif ids.txt /global/homes/jg/schradel/projects/GAGA/annotations/{1}_final_annotation_repfilt_addreannot_noparpse_representative.cds.fasta > {1}.{2}.{3}-{4}.out.fa"
```

#### Make file names as fasta headers and concatenate all files 
```bash
# For one file, keeping the GAGA gene id (e.g.Ankyrin_Obir_g10246_i1)
for file in NCBI-0001.NC_039516.1.740424-747898.out.fa; do sed -i "s/^>/>${file%out.*}/g" "$file"; done

# For all files
for file in *.out.fa; do sed -i "s/^>/>${file%out.*}/g" "$file"; done

# Concatenate all changed files 
cat *out.fa >> ANKs.GAGA.annotation.fa
```


#### Quality check of all files:

```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/genome_reannotate_mod_intersect.gff ]] && echo "$dir"; done
```

Grep fasta sequence from newly intersected files:
```bash

```