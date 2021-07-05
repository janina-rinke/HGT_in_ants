# Blast of LGT candidates with MMseqs2
The file with all 410 LGT candidates can be found here:
```bash
/global/scratch2/j_rink02/master/lgt/0_data/candidatefiles/blastCandidates.tsv
```

Structure of the file:
```bash
GAGA-0020	Scaffold1	5322877	5323044
GAGA-0020	Scaffold107	144838	147367
GAGA-0020	Scaffold17	156258	159230
GAGA-0020	Scaffold267	55721	56766
GAGA-0020	Scaffold31	542876	547130
GAGA-0020	Scaffold31	613867	651311
GAGA-0020	Scaffold31	704726	708354
GAGA-0020	Scaffold38	658653	660450          
```
## Databases to use:
- NR                  	Aminoacid 	       -	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- NT                  	Nucleotide	       -	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- UniRef90            	Aminoacid 	     yes	https://www.uniprot.org/help/uniref

## Workflow:

1. in bash
   - loop over the rows. Start with :
```
GAGA-0020	Scaffold1	5322877	5323044
```

2. Find the corresponding genome fasta file:
```
id=GAGA-0003
/0_data/GAGAgenomes/$id.fa
```
3. extract the correct region from the genome fasta ($id.fa).
​
All together with the `parallel` command:
```bash
cat /global/scratch2/j_rink02/master/lgt/0_data/blastCandidates.tsv | parallel --colsep '\t' "samtools faidx /global/scratch2/j_rink02/master/lgt/0_data/assemblies/{1}*.fasta {2}:{3}-{4} > /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/{1}.{2}.{3}-{4}.fa"
```

What the above command does:
1. give out the file `blastCandidates.tsv`
2. recognize the file and that it is tab-separated, parallel assigns automatically values to each tab separated column {1},{2} etc...
3. use samtools to extract the gaga id {1}, scaffold{2} and start-end {3}-{4}
4. give as output one file per candidate which is named after id, scaffold, start-end.fa

The files in the directory `candidates.fasta` now look like the following:
```bash
GAGA-0020.Scaffold17.156258-159230.fa			 GAGA-0363.Scaffold10.3446272-3446687.fa
GAGA-0020.Scaffold267.55721-56766.fa			 GAGA-0363.Scaffold10.7932710-7932924.fa
GAGA-0020.Scaffold31.542876-547130.fa			 GAGA-0363.Scaffold15.2721767-2721972.fa
GAGA-0020.Scaffold31.613867-651311.fa			 GAGA-0363.Scaffold15.6376837-6377052.fa
GAGA-0020.Scaffold31.704726-708354.fa			 GAGA-0363.Scaffold1.837744-837938.fa
GAGA-0020.Scaffold38.658653-660450.fa			 GAGA-0363.Scaffold2.13656233-13656485.fa
```

-------------------------------------------------------------------------------------------------------
#### Load the module for `mmseqs2`
```bash
source /usr/share/modules/init/bash #enables the module
module use /global/projects/programs/modules/

module load seq-search/mmseqs/sse2-13-45111 #to load the module

# To check which modules are loaded
module list
```
Output of `module list`if the module is loaded:
```
Currently Loaded Modulefiles:
 1) seq-search/mmseqs/sse2-13-45111
 ```

 `module whatis seq-search/mmseqs/sse2-13-45111`:
 Output:
 ```
 seq-search/mmseqs/sse2-13-45111: seq-research/mmseqs/sse2-13-45111 - sets the Environment for mmseqs-sse2-13-45111
               provides:
               - mmseqs

               databases can be found here: /global/scratch2/databases/mmseq/
```

`which mmseqs`
Output, if module is loaded right:
```
/global/projects/programs/source/mmseqs/mmseqs2-sse2-13-45111/bin//mmseqs
```
-------------------------------------------------------------------------------------------------------

### Run mmseqs2 to find the best blast hit
MMseqs2 is an ultra fast and sensitive sequence search
It is a software suite which can cluster huge protein and nucleotide sequence sets. It can run 10000 faster than BLAST.


```bash
# Databases for mmseqs can be found here:
/global/scratch2/databases/mmseq
```
Structure of command for mmseqs modules:
```bash  
mmseqs module input_db output_db args [options]
```

#### 1. Pre-format the databases

######QUERY DB: `all.candidates.fa`
######TARGET DB: `uniRef90`, `NT`, `NR`
```bash
# Query database
mmseqs createdb /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/all.candidates.fa all.candidates.queryDB

# Target database
mmseqs createdb /global/scratch2/databases/mmseq/uniRef90/uniRef90-2021_06_22 uniRef90-2021_06_22
mmseqs createdb /global/scratch2/databases/mmseq/nr/nr-2021_06_22 nr
```  
#### 2. Start the search with mmseqs

#####2.1 Try to run `mmseqs`with one file only to see if it works:
```bash
# Prepare your file of choice as a query database
mmseqs createdb /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/GAGA-0020.Scaffold107.144838-147367.fa GAGA-0020.Scaffold107.144838-147367.fa.queryDB

# start a Screen session to not loose any output
screen -R mmseqs2

# Run mmseqs with your prepared query database against the targetdb `UniRef90`
mmseqs easy-search -s 1 --split-memory-limit 25G  GAGA-0020.Scaffold107.144838-147367.fa /global/scratch2/databases/mmseq/uniRef90/uniRef90-2021_06_22 GAGA-0020.Scaffold107.144838-147367.bls tmp
```   

##### 2.2 Submit the job to the cluster with Gridengine

###### 2.2.1 Search against the `UniRef90` database:
```bash
#$ -S /bin/bash
#$ -N mmseqsjob
#$ -cwd
#$ -pe smp 62
#$ -l h_vmem=2G
#$ -o /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/all.candidates.fasta.out
#$ -e /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/all.candidates.fasta.err
#$ -wd /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta

source /usr/share/modules/init/bash  # enables the module package
module use /global/projects/programs/modules/
module load seq-search/mmseqs/sse2-13-45111

mmseqs easy-search -s 1 --split-memory-limit 25G all.candidates.fa /global/scratch2/databases/mmseq/uniRef90/uniRef90-2021_06_22 all.candidates.bls tmp
```
 `-pe smp` are the cores, while `-l option defines the RAM you'll need. So if you need 124 GB RAM on a server (e.g. from ebbsrv03) you need to ask for all 62 kernels and 2G (2 x 62 = 124 GB RAM).

Submit the job with:
```bash
qsub -l hostname=ebbsrv03 mmseqs_script.sh
```

###### 2.2.2 Search against the `NR` database:
```bash
#$ -S /bin/bash
#$ -N mmseqsjob
#$ -cwd
#$ -pe smp 62
#$ -l h_vmem=2G
#$ -o /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/all.candidates.fasta.nr.out
#$ -e /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/all.candidates.fasta.nr.err
#$ -wd /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta

source /usr/share/modules/init/bash  # enables the module package
module use /global/projects/programs/modules/
module load seq-search/mmseqs/sse2-13-45111

mmseqs easy-search -s 1 --split-memory-limit 25G all.candidates.fa /global/scratch2/databases/mmseq/nr/nr-2021_06_22 all.candidates.nr.bls tmp
```
Submit the job with:
```bash
qsub -l hostname=ebbsrv03 mmseqs_script_nr.sh
```
-----------------------------------------------------------------------
Expected output:
```bash
GAGA-0003.Scaffold1:3000-5000 200 400 Wolbachia pbla 6949218 6811581 1e-160
GAGA-0003.Scaffold1:3000-5000 300 320 ecoli 1949218 1811581 1e-10
```

Then select best hit based on evalue/bitscore
   `cat $id.$scf.$start.$end.bls | sort -k4,4 -nr -k5,5`
