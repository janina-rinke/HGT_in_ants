# Blast of good candidate sequences
File can be found here:
```bash
/Users/Janina/sciebo/Master.LGTs/blastCandidates.tsv
```

Structure of the file:
```bash
GAGAid        scaffold         start    stop
GAGA-0003     Scaffold1        3000     5000            
GAGA-0010     Scaffold5        7000     9000           
```

## Workflow:

1. in bash
   - loop over the rows. Start with :
   ```
   GAGA-0003     Scaffold1        3000     5000
   ```

2. Find the corresponding genome fasta file:
    ```
    id=GAGA-0003
    /data/GAGAgenomes/$id.fa
    ```
3. extract the correct region from the genome fasta ($id.fa)
​
​All together with the `parallel` command:
```bash
cat blastCandidates.tsv| parallel --colsep '\t' "samtools faidx /global/scratch2/j_rink02/master/lgt/0_data/assemblies/{1}*.fasta {2}:{3}-{4} > candidates.fasta/{1}.{2}.{3}-{4}.fa"
```

What the above command does:
1. give out the file `blastCandidates.tsv`
2. recognize the file and that it is tab-separated, parallel assigns automatically values to each tab separated column {1},{2} etc...
3. use samtools to extract the gaga id {1}, scaffold{2} and start-end {3}-{4}
4. give as output one file per candidate which is named after id, scaffold, start-end.fa


-------------------------------------------------------------------------------------------------------
#### Load the module for `mmseqs2`
```bash
source /usr/share/modules/init/bash #enables the module
module use /global/projects/programs/modules/

module load seq-search/mmseqs/avx2-13-45111 #to load the module

# To check which modules are loaded
module list
```
Output of `module list`if the module is loaded:
```
Currently Loaded Modulefiles:
 1) seq-search/mmseqs/avx2-13-45111
 ```

 `module whatis seq-search/mmseqs/avx2-13-45111`:
 Output:
 ```
 seq-search/mmseqs/avx2-13-45111: seq-research/mmseqs/avx2-13-45111 - sets the Environment for mmseqs-avx2-13-45111
               provides:
               - mmseqs

               databases can be found here: /global/scratch2/databases/mmseq/
```

`which mmseqs`
Output, if module is loaded right:
```
/global/projects/programs/source/mmseqs/mmseqs2-avx2-13-45111/bin//mmseqs
```
-------------------------------------------------------------------------------------------------------

### Run mmseqs2 to find the best blast hit
Ultra fast and sensitive sequence search
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
```  
#### 2. Start the search with mmseqs
```bash
mmseqs search all.candidates.fa /global/scratch2/databases/mmseq/uniRef90/uniRef90-2021_06_22 all.candidates.bls
```

Expected output:
```bash
GAGA-0003.Scaffold1:3000-5000 200 400 Wolbachia pbla 6949218 6811581 1e-160
GAGA-0003.Scaffold1:3000-5000 300 320 ecoli 1949218 1811581 1e-10
```
​
5. select best hit based on evalue/bitscore
   `cat $id.$scf.$start.$end.bls | sort -k4,4 -nr -k5,5`
​
# databases to use:
- NR                  	Aminoacid 	       -	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- NT                  	Nucleotide	       -	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- UniRef90            	Aminoacid 	     yes	https://www.uniprot.org/help/uniref
