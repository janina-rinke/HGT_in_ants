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

##### Output files from first `dfast-core`run:
`cds.fna`
```bash
>MGA_1|LOCUS_00010 flavin prenyltransferase UbiX
GTGTTAACGTGGGGTGTAAAAGAAATAATTATGCAGCAACAATTACGTCTAGTAGTGGGAATATCAGGGGCATCTGGTGGTATTTATGGAATTCGTGCATTAACTATTTTAAAAAAATGTTCTATGAATATAGAATCTCATCTTATCGTTACTAGAAATGCGCTGATTACATTGCAACAAGAACTAAAGATGAGCAAACAAGCTGTATATGAATTAGCAGATGTTGTACATTTTCCGCAGGACATCGGAGCATCAGTTGCTAGTGGATCATATCCAACTGTAGGTATGTTAATTGCTCCATGCTCTATTAGAACTATGTCAGAAATTAGTTCAGGAATGACATCGTCATTGATTGGTCGTGTTGCAGATGTAATGTTGAAAGAAAAACGTAGACTAGTGTTAATGATTAGAGAAACTCCATTGCATCTGGGACATTTACGTACAATGGTTAAATTAACAGAATTTGGAGCGGTAATTATGCCGCCTGTTCCTGCATATTATACACATCCAAAAACTATAGATGACGTTATTTGCCACACTGTGGCAAGAGCTCTAAATTTATTTGGAATTGATACCAATATTTCTTCTACTTGGTTAGGAATACGGCATCAAAATCATGATAATCTTATGAAACTAAATTGA
>MGA_2|LOCUS_00020 hypothetical protein
ATGTTAAATAATTATTACAATTTAGACAATATTTTTAATCCTGAAACATCATTACAAAAATGTAACATACATTACCCAAAAAGTGGTCATCAACTATATATCCAAACTAAATTTCCACACTTTTTCGAATTTACTGGAAAAATAAAATTGGAACAATTTATTTATGACAAAAAATACAAAAAAATTTTTAATAAAAAAATTAGTGACTATTATTTATCGTTAGATCTGAATTATCAACCCATTCCCATGCTAGGTTTTAGCATAAATAATATTTTTGTAAACAAGCAATATAATAATACTATTTGTCAAGTATTAATAGATTACCAATTCGGTACTCCTATTATGGAACAAATTAATTATATAAATAAAAAAAATAAATCAATACTAAATAACCTTGATACTGTAATACAACCGTTTATACCTACTATAATTCAATATCGCGATTACATTCCAATTAACGACCATAATCATCTTCCATCCTTGCGGTATACTCAAAAAATAATAGGTTATCCTGGCGAAATTAAAATAATTAAAGTCGATGACAATAATGATAAATCTATGAAATGGGATTTTGAATCATTGCAAAATCACGGAGGTAATATTGTTGCAATAACAAATAATACTTATGCGCTCTATTTCCCTAACTATCCTATTAAACAAGAAGACAATATTTTTGTTTCATACATTACTAATAATGCAACACAAATCCCTCAGGAACAGAAAAAACAAAATATACATATTGTAATAAAAGATTTTTCACAAAAAAAATTATTAAATACGAATACACGAAAAAATAATGCTTCCATTATTAACATTAATAATGATTCCGGAATTAACGTTACAGAAAATACAGCTGAACATACATTAATTTATCATGACGATGATCAGCGTAGCATAAAGAATGATACATCGGTCTCAAATGTCATTACGCAGAATAATATTGATACTCTAACAGAAAACGCGCAGAATGCTGCAATAGAAGATGATTTCACATTTTCAGCTCCTCCTCCCCCTCCAATTCCTTTTCCTTTTTTAGAAAAAGATTCGTCGTCATTTATAGCATCATCAACACTGAACGATACATTTCTGTTGTCCCGATCCGTGTCTCATGAAGATCAACAGACTCATTCTGAAAAAAATAAAGATAATAAAGGAATACACGATGCTTTTTCTGAATATTCTAAAATTGAATTAAATCAAAGTGATGACTTATCGCATCGCTTATCTACACATAGGCAATCTAAGTTTGCTTCAGTAGGCACTAAAGAACATATAAATAAACTTGAAAATACCATAAGAGAACGAAGAAAGCTAAACACCTCAGTGATATGGAATAAATGTTTGTCAAACTAA
```
`genbank.tbl`
```bash
>Feature sequence001
96	737	gene
			locus_tag	LOCUS_00010
			gene	ubiX
96	737	CDS
			product	flavin prenyltransferase UbiX
			inference	COORDINATES:ab initio prediction:MetaGeneAnnotator
			inference	similar to AA sequence:UniProtKB:Q92XP7
			protein_id	gnl|my_center|LOCUS_00010
2300	948	gene
```
`genbank.fsa`
```bash
>sequence001 [organism=] [gcode=11]
GTATATGTAAGATTTTTTAATTCATGATAGCATTCAGAAATATTTTTAGTGTTATGTATA
GTGTTTATTATTAAAAGTCATAAGTATATTTATTTGTGTTAACGTGGGGTGTAAAAGAAA
TAATTATGCAGCAACAATTACGTCTAGTAGTGGGAATATCAGGGGCATCTGGTGGTATTT
ATGGAATTCGTGCATTAACTATTTTAAAAAAATGTTCTATGAATATAGAATCTCATCTTA
TCGTTACTAGAAATGCGCTGATTACATTGCAACAAGAACTAAAGATGAGCAAACAAGCTG
TATATGAATTAGCAGATGTTGTACATTTTCCGCAGGACATCGGAGCATCAGTTGCTAGTG
GATCATATCCAACTGTAGGTATGTTAATTGCTCCATGCTCTATTAGAACTATGTCAGAAA
TTAGTTCAGGAATGACATCGTCATTGATTGGTCGTGTTGCAGATGTAATGTTGAAAGAAA
AACGTAGACTAGTGTTAATGATTAGAGAAACTCCATTGCATCTGGGACATTTACGTACAA
```
`genome.gff`
```bash
##gff-version 3
sequence001	MetaGeneAnnotator	CDS	96	737	.	+	0	ID=MGA_1;product=flavin prenyltransferase UbiX;inference=COORDINATES:ab initio prediction:MetaGeneAnnotator,similar to AA sequence:UniProtKB:Q92XP7;translation=MLTWGVKEIIMQQQLRLVVGISGASGGIYGIRALTILKKCSMNIESHLIVTRNALITLQQELKMSKQAVYELADVVHFPQDIGASVASGSYPTVGMLIAPCSIRTMSEISSGMTSSLIGRVADVMLKEKRRLVLMIRETPLHLGHLRTMVKLTEFGAVIMPPVPAYYTHPKTIDDVICHTVARALNLFGIDTNISSTWLGIRHQNHDNLMKLN;locus_tag=LOCUS_00010;gene=ubiX;EC_number=2.5.1.129;note=Q92XP7 flavin prenyltransferase UbiX (Sinorhizobium meliloti 1021) [pid:57.8%25%2C q_cov:87.8%25%2C s_cov:93.9%25%2C Eval:1.1e-57],COG:COG0163:UbiX 3-polyprenyl-4-hydroxybenzoate decarboxylase [Category:H%2C Aligned:16-202%2C Eval:1.2e-92%2C score:266.8],MGA_1
sequence001	MetaGeneAnnotator	CDS	948	2300	.	-	0	ID=MGA_2;product=hypothetical protein;inference=COORDINATES:ab initio prediction:MetaGeneAnnotator;translation=MLNNYYNLDNIFNPETSLQKCNIHYPKSGHQLYIQTKFPHFFEFTGKIKLEQFIYDKKYKKIFNKKISDYYLSLDLNYQPIPMLGFSINNIFVNKQYNNTICQVLIDYQFGTPIMEQINYINKKNKSILNNLDTVIQPFIPTIIQYRDYIPINDHNHLPSLRYTQKIIGYPGEIKIIKVDDNNDKSMKWDFESLQNHGGNIVAITNNTYALYFPNYPIKQEDNIFVSYITNNATQIPQEQKKQNIHIVIKDFSQKKLLNTNTRKNNASIININNDSGINVTENTAEHTLIYHDDDQRSIKNDTSVSNVITQNNIDTLTENAQNAAIEDDFTFSAPPPPPIPFPFLEKDSSSFIASSTLNDTFLLSRSVSHEDQQTHSEKNKDNKGIHDAFSEYSKIELNQSDDLSHRLSTHRQSKFASVGTKEHINKLENTIRERRKLNTSVIWNKCLSN;locus_tag=LOCUS_00020;note=MGA_2
```

`rna.fna`
```bash
>Aragorn_1|LOCUS_t00010 tRNA-Val
TGCGTTCTTAGCTCAGTTGGTTAGAGCACTACCTTGACATGGTAGAGGTCGATGGTTCAAGTCCATTAGAACGCAT
>Aragorn_2|LOCUS_t00020 tRNA-Met
TAGCTACGTAGCTCAGTGGGTTAGAGCGCAGCACTCATAATGCTGAGGTCACAGGTTCAAGTCCCGTCGTAGCTAT
>Aragorn_3|LOCUS_t00030 tRNA-Leu
```
`protein.faa`
```bash
>MGA_1|LOCUS_00010 flavin prenyltransferase UbiX
MLTWGVKEIIMQQQLRLVVGISGASGGIYGIRALTILKKCSMNIESHLIVTRNALITLQQELKMSKQAVYELADVVHFPQDIGASVASGSYPTVGMLIAPCSIRTMSEISSGMTSSLIGRVADVMLKEKRRLVLMIRETPLHLGHLRTMVKLTEFGAVIMPPVPAYYTHPKTIDDVICHTVARALNLFGIDTNISSTWLGIRHQNHDNLMKLN
>MGA_2|LOCUS_00020 hypothetical protein
MLNNYYNLDNIFNPETSLQKCNIHYPKSGHQLYIQTKFPHFFEFTGKIKLEQFIYDKKYKKIFNKKISDYYLSLDLNYQPIPMLGFSINNIFVNKQYNNTICQVLIDYQFGTPIMEQINYINKKNKSILNNLDTVIQPFIPTIIQYRDYIPINDHNHLPSLRYTQKIIGYPGEIKIIKVDDNNDKSMKWDFESLQNHGGNIVAITNNTYALYFPNYPIKQEDNIFVSYITNNATQIPQEQKKQNIHIVIKDFSQKKLLNTNTRKNNASIININNDSGINVTENTAEHTLIYHDDDQRSIKNDTSVSNVITQNNIDTLTENAQNAAIEDDFTFSAPPPPPIPFPFLEKDSSSFIASSTLNDTFLLSRSVSHEDQQTHSEKNKDNKGIHDAFSEYSKIELNQSDDLSHRLSTHRQSKFASVGTKEHINKLENTIRERRKLNTSVIWNKCLSN
>MGA_3|LOCUS_00030 hypothetical protein
MPTFIILGIIFQIFTYNEIAWCSTLKGCIKSNVSNNIFQNDLYQQEMKLYTHDHIHHTLNFYPYTTNKLRAHAYNYRSPFSSTYKSKMQLQNDSINMFHSFHTQHSKHKTNLSFMQLGIHNLLSEQIFNFGGGKRHLTNDKCAIGYNTFYHCPISNQSSQPYSINVGVEYWLYRVSGKF
```

`statistics.txt`
```bash
Total Sequence Length (bp)	3843330
Number of Sequences	407
Longest Sequences (bp)	840076
N50 (bp)	38848
Gap Ratio (%)	0.094527
GCcontent (%)	33.3
Number of CDSs	1980
Average Protein Length	258.4
Coding Ratio (%)	39.9
Number of rRNAs	7
```

##### Run `dfast` against the UniProt bacterial database only:

For one candidate only to test if script is working:

`nano dfast_bacteria_onecandidate.sh`:

```bash
#$ -S /bin/bash
#$ -N dfastjob
#$ -cwd
#$ -pe smp 50
#$ -l h_vmem=1G
#$ -o /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/dfast.bacteria.onecandidate.out
#$ -e /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/dfast.bacteria.onecandidate.err
#$ -wd /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta

conda activate dfast

dfast -g GAGA-0020.Scaffold107.144838-147367.fa --force --cpu 50 --debug --use_original_name t --minimum_length 100 --database /global/scratch2/databases/dfast/uniprot_bacteria-0.9.ref -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast/one.candidate.uniprot.bacteria --config /home/j/j_rink02/anaconda3/envs/dfast/bin/custom_config.py
```
For another candidate with slightly different options:
```bash
#$ -S /bin/bash
#$ -N dfastjob
#$ -cwd
#$ -pe smp 20
#$ -l h_vmem=6G
#$ -o /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/dfast.bacteria.GAGA-0024.Scaffold1.16539177-16545088.out
#$ -e /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/dfast.bacteria.GAGA-0024.Scaffold1.16539177-16545088.err
#$ -wd /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta

conda activate dfast

dfast -g GAGA-0024.Scaffold1.16539177-16545088.fa --force --metagenome --cpu 20 --debug --use_original_name t --minimum_length 100 --database /global/scratch2/databases/dfast/uniprot_bacteria-0.9.ref -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast/GAGA-0024.Scaffold1.16539177-16545088.uniprot.bacteria --config /home/j/j_rink02/anaconda3/envs/dfast/bin/custom_config.py
```

For all candidates together:

Write a GridEngine Script to run the process:

nano `batch_dfast_job.sh`:

```bash
#$ -S /bin/bash
#$ -N batchdfastjob
#$ -cwd
#$ -w e
#$ -V
#$ -pe smp 20
#$ -l h_vmem=6G
#$ -o /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/$file.out
#$ -e /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/$file.err
#$ -wd /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/by.genome

conda activate dfast

echo "Running on file: $file"

dfast -g $file --force --metagenome --cpu 20 --debug --use_original_name t --minimum_length 100 --database /global/scratch2/databases/dfast/uniprot_bacteria-0.9.ref -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast/$file --config /home/j/j_rink02/anaconda3/envs/dfast/bin/custom_config.py
```

Submit the script and parse bash variables to the script using `-v file="file1.fa"`

To submit jobs more easily, use `parallel`:
```bash
find . -name "GAGA*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh

find . -name "NCBI*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh

find . -name "OUT*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh
```

After running the script, not all candidates contain all files from the complete pipeline (Possible bug in GridEngine). Files, which do not contain all resulting dfast files need to be identified and the `dfast`pipeline needs to be re-run for them.

####1.3 Debugging of incomplete `dfast` candidates

Find all directories which do not contain a `cds.fna` file and save them into a list:
```bash
cd /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast

find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/cds.fna ]] && echo "$dir"; done > ../dfast2repeat.lst
```

Re-run the `dfast_batch_job.sh` script on all genomes in the `candidates.fasta` folder based on the elements given by the `dfast2repeat.lst`:

```bash
cd /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta

cat ../gene_annotation/dfast2repeat.lst| cut -f 2 -d "/"| parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh -o /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/%.out -e /global/scratch2/j_rink02/master/lgt/2_analysis/candidates.fasta/tmp/%.err
```

Control, that all candidates have a completed pipeline:
```bash
cd /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast

find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/cds.fna ]] && echo "$dir"; done | wc -l
```
If all candidates have successfully completed the pipeline, no more candidate directories should be showing up here.

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

## 3. Prokaryotic gene annotation with `kraken`

Kraken is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic studies.


`Input:` Fasta files of all GAGA-genomes with evaluated good candidates.
Output was produced by Ding He with the `kraken` tool.

 All candidates were re-evaluated for being good candidates based on the prokaryotic gene annotation from Kraken.

The LGT in GAGA-0515 is a good example how the annotation for a good
LGT candidate should look like:


##### Output Files:
```bash
GAGA-0515.k2_bracken.report
GAGA-0515.k2.class.bracken
GAGA-0515.k2.krona
GAGA-0515.k2.out
GAGA-0515.k2.report
GAGA-0515.k2.species.bracken
GAGA-0515.krona.html
```

Output of `GAGA-0515.k2_bracken.report`:
```bash
0.00	0	0	U	0	unclassified
100.00	1	0	R	1	root
100.00	1	0	R1	131567	  cellular organisms
100.00	1	0	D	2	    Bacteria
100.00	1	0	P	1224	      Proteobacteria
100.00	1	0	C	1236	        Gammaproteobacteria
100.00	1	0	O	91347	          Enterobacterales
100.00	1	0	F	1903411	            Yersiniaceae
100.00	1	1	G	613	              Serratia
```

Output of `GAGA-0515.k2.class.bracken`:
```bash
name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
Gammaproteobacteria	1236	C	1	0	1	1.00000
```

The `GAGA-0515.krona.html` file (file:///Users/Janina/Downloads/taxonomy_classification/GAGA-0515.krona.html) shows the taxonomic assignment graphically.

Based on this good example, all candidates were evaluated.
