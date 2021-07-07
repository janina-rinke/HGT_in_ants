## Prokaryotic gene annotation with `dfast`, `prodigal` and `metageneannotator`

##### Janina Rinke, 02.07.2021

metageneannotator: http://metagene.nig.ac.jp/
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
prodigal -p meta -i all.candidates.fa -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/annotated_genes_prodigal.gbk
```
Output:
```bash
DEFINITION  seqnum=1;seqlen=2530;seqhdr="GAGA-0020.Scaffold107.144838-147367";version=Prodigal.v2.6.3;run_type=Metagenomic;model="39|Rickettsia_conorii_Malish_7|B|32.4|11|1";gc_cont=32.40;transl_table=11;uses_sd=1
FEATURES             Location/Qualifiers
     CDS             complement(<2..>2530)
                     /note="ID=1_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.361;conf=99.99;score=348.99;cscore=347.38;sscore=1.61;rscore=0.00;uscore=0.00;tscore=1.61;"
//
DEFINITION  seqnum=2;seqlen=168;seqhdr="GAGA-0020.Scaffold1.5322877-5323044";version=Prodigal.v2.6.3;run_type=Metagenomic;model="0|Mycoplasma_bovis_PG45|B|29.3|4|1";gc_cont=29.30;transl_table=4;uses_sd=1
FEATURES             Location/Qualifiers
     CDS             complement(38..>166)
                     /note="ID=2_1;partial=01;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.295;conf=98.17;score=17.33;cscore=14.11;sscore=3.22;rscore=0.00;uscore=0.00;tscore=3.22;"
//
```

####2.2 Compare prokaryotic genes (predicted with `prodigal`) to a prokaryotic database with `mmseqs`
