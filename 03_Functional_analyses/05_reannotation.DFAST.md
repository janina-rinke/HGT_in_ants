# Reannotation with DFAST

This script only needs a GAGA genome as a FASTA file (path: `/global/scratch2/j_rink02/master/lgt/0_data/assemblies`) and the extended HGT locus by user-defined size.

### 1) Extend the locus by 1000 bp on each side.

Submit the job as a script called `rerun.dfast.sh` to GridEngine and loop over each row in the file `locus.txt`.

```bash
#$ -S /bin/bash
#$ -N rerunDFAST
#$ -cwd
#$ -w e
#$ -V
#$ -pe smp 30
#$ -l h_vmem=6G

# define size of the extension around LGT locus to use for the annotation
shiftSize="1000"

# define locus
locus=${GAGAid}.${scaffold}.${start}-${stop}

# create a softlink in the current working directory that links to the fasta index for the current GAGAid
ln -s ./0_data/assemblies/${GAGAid}*.fasta .

# create a softlink in the current working directory that links to the fasta index for the current GAGAid
ln -s ./0_data/assemblies/${GAGAid}*fasta.fai .

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

dfast -g ${locus}.${shiftSize}.fa --force --metagenome --cpu 30 --debug --use_original_name t --minimum_length 100 --database ./databases/dfast/uniprot_bacteria-0.9.ref -o ./2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out --config ./0_data/custom_config.py

### Update the coordinates to genome level

#remove fasta from gff (bedtools doesn't like it)
#change first field to keep only the scaffold
# shift by original start coordinate
#shift by -n kb

sed '/^##FASTA$/,$d' ./2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out/genome.gff| \
perl -pe 's/^('${scaffold}').*?\t/$1\t/g' | \
bedtools shift -i - -g ${genome}.genome -s ${start} -header | \
bedtools shift -i - -g ${genome}.genome -s -${shiftSize} -header \
> ./2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out/genome_reannotate_mod.gff

## check if extraction is right
#bedtools getfasta  -fi ${genome} -bed ${locus}/genome_reannotate_mod.gff -s

### Intersect with `bedtools`
# At the end, you will have a new gff ${locus}.1000.out/genome_reannotate_mod_intersect.gff, which has all predicted CDS that overlap with the LGT locus.

# find all the entries in the "expanded" gff that overlap with the LGT locus.
echo -e ${scaffold}"\t"${start}"\t"${stop}"\t"${locus} | bedtools intersect -wa -b stdin -a ./2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out/genome_reannotate_mod.gff > ./2_analysis/gene_annotation/reannotation.dfast/results/${locus}.1000.out/genome_reannotate_mod_intersect.gff
```

Submit the script `rerun.dfast.sh` to GridEngine:
```bash
cat locus.txt | parallel --col-sep "\t" --dryrun --max-args 1 qsub -v GAGAid={1} -v scaffold={2} -v start={3} -v stop={4} rerun.dfast.sh
```

In the first version of the script, bedtools intersect did not actualize the coordinates for the re-annotated CDS  `genome_reannotate_mod_intersect.gff` correctly. After completing all dfast jobs, the `bedtools intersect` command was run again.

```bash
cd ./2_analysis/gene_annotation/reannotation.dfast/results

cat locus2.txt | parallel --col-sep "\t" --max-args 1 "echo -e '{2}\t{3}\t{4}\t{5}' | bedtools intersect -wa -b stdin -a ./2_analysis/gene_annotation/reannotation.dfast/results/{5}.1000.out/genome_reannotate_mod.gff > ./2_analysis/gene_annotation/reannotation.dfast/results/{5}.1000.out/genome_reannotate_mod_intersect.gff"
```

Quality check of all files:

Check if all files have completed the dfast re-annotation:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/genome_reannotate_mod_intersect.gff ]] && echo "$dir"; done
```

No files are appearing here, which means that all candidates have completed the pipeline.

However some files are empty (`find */*gff -empty`):
```bash
GAGA-0025.Scaffold15.4178206-4178515.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0085.Scaffold5.1024757-1024905.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0087.Scaffold383.14647-14826.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0099.Scaffold4.5794430-5794610.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0200.Scaffold27.1009029-1849104.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0221.Scaffold11.416033-416312.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0222.Scaffold130.329306-329454.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0222.Scaffold85.858212-858361.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0224.Scaffold24.2173785-2173936.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0229.Scaffold18.2465400-2465653.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0245.Chr7.18043257-18043592.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0246.Scaffold7.9510351-9510585.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0256.Scaffold51.1261353-1261551.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0288.Scaffold75.646615-647586.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0302.Scaffold187.20714-20890.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0302.Scaffold406.14472-14650.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0328.Scaffold1.9848465-9848878.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0328.Scaffold17.636160-636366.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0333.Scaffold8.1681954-1682103.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0350.Scaffold9.243113-243382.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0353.Scaffold133.191513-191764.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0353.Scaffold33.2046976-2047194.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0354.Scaffold108.264230-264804.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0363.Scaffold15.2721767-2721972.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0363.Scaffold2.13656233-13656485.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0364.Scaffold43.1414194-1414359.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0365.Scaffold3.164904-165353.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0365.Scaffold3.225903-226228.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0382.Scaffold12.6451061-6451382.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0401.Scaffold102.52501-52756.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0405.Scaffold7.485313-485463.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0502.Scaffold24.3179491-3179882.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0505.Scaffold84.55509-55735.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0512.Scaffold26.1740048-1740197.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0512.Scaffold65.432822-433169.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0513.Scaffold118.10998-11334.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0524.Scaffold92.200610-200775.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0527.Scaffold16.2537796-2538281.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0527.Scaffold7.11475564-11475713.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0528.Scaffold116.534097-534334.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0532.Scaffold73.609296-610495.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0578.Scaffold21.4767189-4768777.1000.out/genome_reannotate_mod_intersect.gff
GAGA-0579.Scaffold13.389277-389689.1000.out/genome_reannotate_mod_intersect.gff
OUT-0001.CM020815.1.6441834-6442103.1000.out/genome_reannotate_mod_intersect.gff
OUT-0002.Scaffold14.831336-831605.1000.out/genome_reannotate_mod_intersect.gff
```
The whole folder `reannotate.dfast` was used to run the RNAseq mapping again. The reannotated RNAseqmapping folder is called  `run_final_keepingbam_pergagaid_latest`.

##### Re-Run commands to update summary table

### 2) Re-calculate start and stop codons from the dfast reannotation.

The first codon (`start codon`) should usually be ATG (theoretically: GTG or rarely TTG, see https://en.wikipedia.org/wiki/Start_codon).
The last codon (`stop-codon`) should be TAA , TAG , or TGA.

```bash
seqkit fx2tab cds.fna |awk -F '\t' '{print $1,substr($2,1,3),substr($2,length($2)-2,length($2))}'
```

For all candidates:
```bash
for i in * ; do seqkit fx2tab $i/cds.fna |awk -F '\t' '{print $1,substr($2,1,3),substr($2,length($2)-2,length($2))}' > $i.start.stop.codons.tsv ; done

#prints out the file directory names as lines, does not work with MacOS
printf "%s\n" *.start.stop.codons.tsv | xargs -n1 -d $'\n' bash -c 'xargs -n1 -d $'\''\n'\'' printf "%s,%s\n" "$1" <"$1"' -- > all.reannotated.start.stop.codons.tsv

#printf: used for format and print data, %s presents a string place holder for the space separated string that follows and \n represents a line feed. In this case, "%s" will then be substituted for the input of *.start.stop.codons.tsv followed by a line feed.
```

Output of the file `all.reannotated.start.stop.codons.tsv`:
```bash
GAGA-0020.Scaffold17.156258-159230.fa	MGA_1 LOCUS_10 Ankyrin_repeat_protein ATG TAA
GAGA-0020.Scaffold17.156258-159230.fa	MGA_2 LOCUS_20 hypothetical_protein ATG TGT
GAGA-0020.Scaffold267.55721-56766.fa	MGA_1 LOCUS_10 hypothetical_protein ATG TAA
GAGA-0020.Scaffold31.542876-547130.fa	MGA_1 LOCUS_10 hypothetical_protein ATG GTG
GAGA-0020.Scaffold31.613867-616833.fa	MGA_1 LOCUS_10 Ankyrin_repeat_protein GAA CCT
GAGA-0020.Scaffold31.648280-651311.fa	MGA_1 LOCUS_10 hypothetical_protein ATG TGA
GAGA-0020.Scaffold31.648280-651311.fa	MGA_2 LOCUS_20 hypothetical_protein ATG TAA
```

### 3) Merge all reannotated bam files and extract unique read counts

Merge all reannotatedRNAseqbamfiles:
```bash
cd ./RNAseqmapping/reannotation.dfast_RNAseq_mapping
for i in */*.1000.out; do samtools merge  $i/mergedRNAseq.bam $i/*LGTregion.bam; done
```

Check if any files are empty:
```bash
find */*/mergedRNAseq.bam -empty
```
No file appears here, so all candidates do have a `mergedRNAseq.bam` file.

Extract only the unique read counts out of every `mergedRNAseq.bam` file. `Sambamba` works on Linux Computers only.
```bash
# For one candidate only:
sambamba view -h -F [NH]==1 mergedRNAseq.bam > uniquely_mapped.LGTregion.bam

# For all candidates:
cd ./RNAseqmapping/reannotation.dfast_RNAseq_mapping

for i in */*; do sambamba view -h -F [NH]==1 $i/mergedRNAseq.bam > $i/uniquely_mapped.LGTregion.bam; done
# -F: filter
# -h: header
```

To check if all read counts have been extracted correctly, go to any directory of your choice and compare the `mergedRNAseq.bam` file and `uniquely_mapped.LGTregion.bam`:
```bash
samtools view mergedRNAseq.bam | awk '{print $1}' | sort | uniq | wc -l

samtools view uniquely_mapped.LGTregion.bam | awk '{print $1}' | sort | uniq | wc -l
```

The `uniquely_mapped.LGTregion.bam` file still needs to be sorted, otherwise the file cannot be read in R and it cannot be indexed.
```bash
for i in */*; do samtools sort -o $i/uniquely_mapped.LGTregion.sorted.bam $i/uniquely_mapped.LGTregion.bam; done
```
