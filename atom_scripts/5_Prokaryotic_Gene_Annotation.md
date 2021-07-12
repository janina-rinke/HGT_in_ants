## Prokaryotic gene annotation with `dfast` and `prodigal`

##### Janina Rinke, 02.07.2021
dfast: https://dfast.ddbj.nig.ac.jp/


## 1. Prokaryotic gene annotation with `dfast`
DFAST is a flexible and customizable pipeline for prokaryotic genome annotation.
#### 1.1 Run `dfast` on the website

File with LGT candidates:
```bash
/global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/all.candidates.fa
```
Result of `dfast`run on the website:
```bash
Total Sequence Length (bp):	3843330
Number of Sequences:	407
Longest Sequences (bp):	840076
N50 (bp):	38848
Gap Ratio (%):	0.094527
GCcontent (%):	33.3
Number of CDSs:	2021
Average Protein Length:	256.1
Coding Ratio (%):	40.4
Number of rRNAs:	6
Number of tRNAs:	52
Number of CRISPRs:	2
```

######Remove the *Blochmannia* LGT from the file, to exclude this single large region:

Use `samtools`:
1. Create an index of the fasta file
```bash
samtools faidx all.candidates.fa
```
2. Filter the fasta file with the help of the index:
```bash
samtools faidx -o all.candidates.nB.fa all.candidates.fa ids…
```

####1.2 Run `dfast` on the cluster
The dfast command-line tool is also referred to as dfast-core to differentiate it from its online version.
##### Install `dfast` with `conda`:
```bash
# Create a conda environment with python
conda create --name dfast python=3.8

# Activate the new environment
conda activate dfast

# Install dfast and provide newest version
conda install -c bioconda dfast=1.2.14
```

##### Download the databases for dfast:
```bash
dfast_file_downloader.py --protein dfast
dfast_file_downloader.py --cdd Cog --hmm TIGR
```
Files are written into `/global/homes/jg/j_rink02/anaconda3/envs/dfast/opt/dfast-1.2.14/db/protein`
How to install databases for `dfast`: https://github.com/nigyta/dfast_core/#installation

##### Run `dfast`:
```bash
dfast --genome all.candidates.fa  --force --minimum_length 100 --metagenome -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast
```

##### Run `dfast` against the UniProt bacterial database only:
```bash
# activate environment first
conda activate dfast

dfast --genome all.candidates.fa --force --use_original_name t --minimum_length 100 --metagenome --database /global/scratch2/databases/dfast/uniprot_bacteria-0.9.ref -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast/all.candidates.uniprot.bacteria  
```
-------------------------------------------------------------------------
## 2. Prokaryotic gene annotation with `prodigal`

####2.1 Run `prodigal` on the cluster
Use the option `-p meta` to annotate sequences with less than 20 kb length.

##### Load the module for `prodigal`:
```bash
source /usr/share/modules/init/bash  # enables the module package
module use /global/projects/programs/modules/
module avail   # list all available modules
module load annotation/prodigal/2.6.3
```  

Run the gene annotation with:
```bash
prodigal -p meta -i all.candidates.fa -a protein.translations.faa -d annotated_genes_prodigal.cds -w prodigal.statistics -f gff -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/prodigal/prodigal_annotated_genes.gff
```

By default, `prodigal`produces one output file, which consists of gene coordinates and some metadata associated with each gene.
To produce four more output files:
`-a`option: protein translations
`-d`option: nucleotide sequences
`-s`option: complete listing of all start/stop pairs with score information
`-w`option: summary of statistical information about the genome or metagenome
`-f`option: specify another output type e.g. `gff` format

Produced output files:
```bash
annotated_genes_prodigal.cds  prodigal_annotated_genes.gff  protein.translations.faa
annotated_genes_prodigal.gbk  prodigal.statistics
```
Fields in this header are as follows:
seqnum: An ordinal ID for this sequence, beginning at 1.
seqlen: Number of bases in the sequence.
seqhdr: The entire FASTA header line.
version: Version of Prodigal used to analyze this sequence.
run_type: "Ab initio" for normal mode, "Anonymous" for anonymous mode.
model (Anonymous mode only): Information about the preset training file used to analyze the sequence.
gc_cont: % GC content of the sequence.
transl_table: The genetic code used to analyze the sequence.
uses_sd: Set to 1 if Prodigal used its default RBS finder, 0 if it scanned for other motifs.

Number of predicted LGTs which have an annotation with `prodigal` and `dfast`:
```bash
# prodigal
cat protein.translations.faa | grep ">" | wc -l
#2182

#dfast
cat protein.faa | grep ">" | wc -l
#1980
```
