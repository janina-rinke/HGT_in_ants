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
3. use samtools to extract the gaga id {1}, scaffold{2} and start-end {3}-{4} from every genome fasta file in the directory containing the assemblies
4. give as output one file per candidate which is named after id, scaffold, start-end.fa and save it in the directory `candidates.fasta`

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
It is a software suite which can cluster huge protein and nucleotide sequence sets. It can run 10000 faster than BLAST. MMseqs2 provides easy workflows to cluster, search and assign taxonomy.


```bash
# Databases for mmseqs can be found here:
/global/scratch2/databases/mmseq
```
Structure of command for mmseqs modules:
```bash  
mmseqs esay-search input_db output_db args [options]
```

The `easy-search` searches directly with a FASTA file against either another FASTA file or an already existing MMseqs2 database.

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
qsub mmseqs_script.sh
```

First ten lines of the output file `all.candidates.bls`:
```bash
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_UPI001250B95F	0.933	135	3	0	140	6	19	63	3.569E-18	87
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_A0A1A9VKG7	0.933	135	3	0	140	6	681	7253.569E-18	87
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_A0A178GX68	0.911	135	4	0	140	6	19	63	1.270E-17	86
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_A0A178GYM1	0.911	135	4	0	140	6	19	63	1.270E-17	86
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_A0A176Q5F8	0.911	135	4	0	140	6	19	63	1.745E-17	85
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_A0A7J6YK53	0.911	135	4	0	140	6	19	63	1.745E-17	85
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_C0R374	0.911	135	4	0	140	6	19	63	2.397E-17	85
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_UPI00157A6D4B	0.888	135	5	0	140	6	19	63	3.293E-17	85
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_A0A059PCR4	0.888	135	5	0	140	6	19	63	3.293E-17	85
GAGA-0350.Scaffold21.1074942-1075284	UniRef90_A0A178GWB3	0.888	135	5	0	140	6	19	63	3.293E-17	85
```
The output is formatted as a tab-separated list with 12 columns:
(1) query id
(2) target id
(3) sequence identity
(4) alignment length
(5) number of mismatches
(6) number of gap openings
(7-8) domain start and end in query
(9-10) domain start and end in target
(11) E-value
(12) bit score  

###### 2.2.2 Search against the `NR` database:
```bash
#$ -S /bin/bash
#$ -N mmseqsjob
#$ -cwd
#$ -pe smp 50
#$ -l h_vmem=6G
#$ -o /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/all.candidates.fasta.nr.out
#$ -e /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/all.candidates.fasta.nr.err
#$ -wd /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta

source /usr/share/modules/init/bash  # enables the module package
module use /global/projects/programs/modules/
module load seq-search/mmseqs/sse2-13-45111

mmseqs easy-search -s 1 --split-memory-limit 300G all.candidates.fa /global/scratch2/databases/mmseq/nr/nr-2021_06_22 all.candidates.nr.bls tmp
```
Submit the job with:
```bash
qsub mmseqs_script_nr.sh
```

The job will only run on servers which provide enough RAM, so you do not need to specify a server with `hostname`!

Output of `all.candidates.nr.bls`:
```bash
GAGA-0511.Scaffold5.3127011-3129187	WP_065094361.1	0.908	1470	45	0	2106	637	1	490	8.552E-278	861
GAGA-0511.Scaffold5.3127011-3129187	WP_182319239.1	0.902	1470	48	0	2106	637	1	490	3.740E-276	857
GAGA-0511.Scaffold5.3127011-3129187	WP_182365617.1	0.902	1470	48	0	2106	637	1	490	9.619E-276	856
GAGA-0511.Scaffold5.3127011-3129187	WP_182183394.1	0.846	1470	75	0	2106	637	1	490	1.205E-255	797
GAGA-0511.Scaffold5.3127011-3129187	MBR9983826.1	0.840	1470	78	0	2106	637	1	490	2.261E-255	797
GAGA-0511.Scaffold5.3127011-3129187	WP_149168725.1	0.844	1470	76	0	2106	637	1	490	3.097E-255	796
GAGA-0511.Scaffold5.3127011-3129187	WP_077188281.1	0.844	1470	76	0	2106	637	1	490	5.812E-255	796
GAGA-0511.Scaffold5.3127011-3129187	WP_006014162.1	0.842	1470	77	0	2106	637	1	490	1.091E-254	795
GAGA-0511.Scaffold5.3127011-3129187	WP_077190377.1	0.842	1470	77	0	2106	637	1	490	3.842E-254	793
GAGA-0511.Scaffold5.3127011-3129187	WP_141457075.1	0.840	1470	78	0	2106	637	1	490	7.210E-254	792
 ```
-----------------------------------------------------------------------
#### 3. Process the results of mmseqs

##### 3.1 Find the best blast hit for every candidate in the `NR`file

```bash
cat all.candidates.nr.bls |sort -k1,1 -k11,11g -k12,12gr | sort -u -k1,1 --merge > all.candidates.nr.bestHit.bls
```
`-k` option sorts the output based on a certain column.
We will sort first based on e-value (column 12) and after that based on bitscore (column 11). Now the hit which is at the top should always have the highest bitscore.

###### 3.2 Find the taxonomy for the best hits with `eutils`
Loop over each hit in the file `all.candidates.nr.bestHit.bls` with `parallel`. Restrict the jobs to `-j 3` to not crash the NCBI server. 

```bash
proteinQuery=WP_065094361.1
esearch -db protein -query ${proteinQuery}|elink -target taxonomy|efetch -mode xml|xtract -pattern Lineage -element Lineage
```

Output:
```bash
cellular organisms; Bacteria; Proteobacteria; Alphaproteobacteria; Rickettsiales; Anaplasmataceae; Wolbachieae; Wolbachia
```
