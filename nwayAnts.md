# 2-way analysis ants

This is the documentation for the analysis of 3 ant genomes with the nway program available at http://retrogenomics2.uni-muenster.de/tools/nway/index.hbi?



## Define working environment

All work for this project is done in the folder `/global/scratch2/j_rink02/hiwi_projects/nway3ants`

```bash
# define $base variable as basedir to work from
base=/global/scratch2/j_rink02/hiwi_projects/nway3ants

# create folder
mkdir ${base}
cd ${base}

mkdir ${base}/0_data #folder to store all input data (e.g. genomes, etc.)
mkdir ${base}/1_code #folder for all code (e.g. this script)
mkdir ${base}/2_analysis #intermediate results and output of intermediate steps
mkdir ${base}/3_output #final output that will be used for reporting, exporting, or in general outside processing

```

## Download genomes from NCBI

- 3 "closely related" species:
    - Target: Monomorium pharaonis https://www.ncbi.nlm.nih.gov/genome/?term=txid307658[orgn]
    - Solenopsis invicta https://www.ncbi.nlm.nih.gov/genome/?term=txid13686[orgn]
    - Acromyrmex echinatior https://www.ncbi.nlm.nih.gov/genome/?term=Acromyrmex+echinatior

Phylogeny: ((Mpha,Sinv),Aech);

Genomic data was downloaded from NCBI's ftp server at
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/373/865/GCF_013373865.1_ASM1337386v2/

Data for the three species can be found on the ftp server by:
- searching for the species
- clicking on the assembly
- use the best assembly (usually the newest one)
- RefSeq accession looks like this: `GCF_016802725.1 (latest)`
- Click on accession to get to index of all files for the species

```bash Download data
cd ${base}/0_data/
## Monomorium
# download genomic data from ftp server
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/373/865/GCF_013373865.1_ASM1337386v2/GCF_013373865.1_ASM1337386v2_genomic.fna.gz

## Solenopsis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/802/725/GCF_016802725.1_UNIL_Sinv_3.0/GCF_016802725.1_UNIL_Sinv_3.0_genomic.fna.gz

## Acromyrmex
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/204/515/GCF_000204515.1_Aech_3.9/GCF_000204515.1_Aech_3.9_genomic.fna.gz

# unzip all files
gzip -d file.gz
```

We used the `*_genomic.fna.gz` file, as this is the FASTA format of genomic sequence(s) in the assembly.
Further information on all files can be found in the README file:

```bash
*_genomic.fna.gz file
      FASTA format of the genomic sequence(s) in the assembly. Repetitive
      sequences in eukaryotes are masked to lower-case (see below).
      The FASTA title is formatted as sequence accession.version plus
      description. The genomic.fna.gz file includes all top-level sequences in
      the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds,
      unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds
      that are part of the chromosomes are not included because they are
      redundant with the chromosome sequences; sequences for these placed
      scaffolds are provided under the assembly_structure directory.
```

## Check if genomes are soft masked

By counting the number of lower-case letters, we determined whether fasta sequences are softmasked.
We also manually screened the fasta sequences to make sure repetitive stretches were soft-masked.

```bash Check for soft masking
cd ${base}/0_data/

## count number of lower case letters in a file
grep -c "a" GCF_013373865.1_ASM1337386v2_genomic.fna
grep -c "t" GCF_013373865.1_ASM1337386v2_genomic.fna
grep -c "c" GCF_013373865.1_ASM1337386v2_genomic.fna
grep -c "g" GCF_013373865.1_ASM1337386v2_genomic.fna
```

Check manually whether genomes are soft-masked:
```
>NC_050467.1 Monomorium pharaonis isolate MP-MQ-018 chromosome 1, ASM1337386v2, whole genome shotgun sequence
TCGTCAATTAAATGACAAATATAATAGattgacagaaaaaaaaaacacttgtGAAAAACTTAAGGCAAAAGTGTGGTGAA
GAAGATAAAGAGTAGCATTTCGTTCTTGAAAggatatataataacaaaaagatCACATAGCATCAGACATGTAGCACGTG
AATTTAACCTGTTTTAATTTACCAGAGAACCAGGGCGATTAGAGTAAGTCAAccatgtaatttataaaatgacaaaaaga

>NC_052664.1 Solenopsis invicta isolate M01_SB chromosome 1, UNIL_Sinv_3.0, whole genome shotgun sequence
AAAGCGTTGACAAGAAGCATAAAATACATAGATCGTTTATAGAGATGATACTATCATAATACtacagagctaatactaga
gataataatacagagctaatactagagataagaatacagagctaatactagagatgataatacagagctaatactagaga
tcataatacagagctaatactagagataagaatacagagctaatactagagatgataatacagagctaatactagagatC
ATAATACAGAGTTATgcgggctgcttgcccgtacggtcggggagtcgcgggacgcggcgggattggctaacgcggtggcG
TGAGCGTGTGCGCGAGCGCGAGTCGTTGCCGATGCGGTCgtcgccgcgattcgcccgagcgcgccccgtgcaagtggagc

>NW_011623521.1 Acromyrmex echinatior unplaced genomic scaffold, Aech_3.9 C1841224, whole genome shotgun sequence
GTTTATTCAAACTTAGATACAGATGAAGAATAAGAATGTTCATatgcaagaaatatttttttcttattattatatgtaca
atacatTCTTTATGCTCTTACTTTCAAATGTTTCTAATTAAATGAGAGAAAGCCTATATTCGATCTCGTTTCTAATGTAT
ATTGTTTGTATTTCCTTTTACATATCAGATCTCTTGAATC
>NW_011623522.1 Acromyrmex echinatior unplaced genomic scaffold, Aech_3.9 C1841250, whole genome shotgun sequence
atattttattataatttatgttatttatttagttttatatttatttaatttttatatttatgtgcaAGGATACTAAatga
ttcaatattattataaatgtgggctataattattattcaataatgttagtgatgttgaaataatagattttgtgtttagt
```

## Annotate introns in gffs

We used genometools' (http://genometools.org/) `gt gff3` to annotate introns in the gffs.

```bash Annotate introns

## download gff files from NCBI for all three species and annotate introns
# Monomorium
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/373/865/GCF_013373865.1_ASM1337386v2/GCF_013373865.1_ASM1337386v2_genomic.gff.gz

/global/projects/programs/source/genometools-1.5.10/bin/gt gff3 -retainids -addintrons GCF_013373865.1_ASM1337386v2_genomic.gff > /global/scratch2/j_rink02/hiwi_projects/nway3ants/2_analysis/Mphar_introns_gff3


# Solenopsis
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/802/725/GCF_016802725.1_UNIL_Sinv_3.0/GCF_016802725.1_UNIL_Sinv_3.0_genomic.gff.gz

/global/projects/programs/source/genometools-1.5.10/bin/gt gff3 -retainids -addintrons GCF_016802725.1_UNIL_Sinv_3.0_genomic.gff > /global/scratch2/j_rink02/hiwi_projects/nway3ants/2_analysis/Sinvi_introns_gff3


# Acromyrmex
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/204/515/GCF_000204515.1_Aech_3.9/GCF_000204515.1_Aech_3.9_genomic.gff.gz

/global/projects/programs/source/genometools-1.5.10/bin/gt gff3 -retainids -addintrons GCF_000204515.1_Aech_3.9_genomic.gff > /global/scratch2/j_rink02/hiwi_projects/nway3ants/2_analysis/Aech_gff3


```
## Keep only intron coordinates

To keep only the intron coordinates in the gffs, we used `awk` in bash.

```bash get intron coordinates
## Monomorium
cat Mphar_gff3 | awk ' $3=="intron" { print $1"\t"$4"\t"$5}' > Mphar.intron.coordinates.tsv
# Nr of intron coordinates: 237633

## Solenopsis
cat Sinvi_gff3 | awk ' $3=="intron" { print $1"\t"$4"\t"$5}' > Sinvi.intron.coordinates.tsv
# Nr of intron coordinates: 252999

## Acromyrmex
cat Aech_gff3 | awk ' $3=="intron" { print $1"\t"$4"\t"$5}' > Aech.intron.coordinates.tsv
# Nr of intron coordinates: 152183
```

## Run 2-ways
Next, the following 2-ways were run in the browser with the `*_genomic.fna` files of all three ants:

1) Mphar (target) - Sinvi (query)
Name: `Mphar_GCF_013373865.1_vs_UNIL_Sinv_3.0`
ID: `twoway 1614612353336`

2) Mphar (target) - Aech (query)
Name: `Mphar_GCF_013373865.1_vs_Aech_3.9`
ID: `twoway 1614612772558`

3) Sinvi(target) - Aech (query)
Name: `UNIL_Sinv_3.0_vs_Aech_3.9`
ID: `twoway 1614613162852`

4) Sinvi(target) - Mphar (query)
Name:`UNIL_Sinv_3.0_vs_Mphar_ASM`
ID:`twoway 1614613580306`

Last saved & changed on 01/03/21 
