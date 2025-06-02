# 1 Prepare fasta files for all HGT candidates to use in a prokaryotic gene annotation with DFAST

Structure of the file with coordinates of all 497 high-quality HGT candidates `GAGA.LGTs.allcoordinates.tsv`:
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
id=GAGA-0020
/0_data/GAGAgenomes/$id.fa
```
3. extract the correct region from the genome fasta (`$id.fa`).
​
All together with the `parallel` command:
```bash
cat GAGA.LGTs.allcoordinates.tsv | parallel --colsep '\t' "samtools faidx {1}*.fasta {2}:{3}-{4} > ./2_analysis/candidates.fasta/{1}.{2}.{3}-{4}.fa"
```

What the above command does:
1. give out the file `GAGA.LGTs.allcoordinates.tsv`
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

All files in the `candidates.fasta` directory are ready to be used for the prokaryotic gene annotation.

# 2 Prokaryotic gene annotation with `dfast` and `prodigal`

dfast: https://dfast.ddbj.nig.ac.jp/

## 2.1 Prokaryotic gene annotation with `dfast`
DFAST is a flexible and customizable pipeline for prokaryotic genome annotation.
#### 2.1.1 Run `dfast` on the website

File with LGT candidates:
```bash
all.candidates.fa
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

#### 2.1.2 Run `dfast` on the cluster
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

How to install databases for `dfast`: https://github.com/nigyta/dfast_core/#installation

##### Run `dfast`:
```bash
dfast --genome all.candidates.fa  --force --minimum_length 100 --metagenome -o ./2_analysis/gene_annotation/dfast
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

conda activate dfast

dfast -g GAGA-0020.Scaffold107.144838-147367.fa --force --cpu 50 --debug --use_original_name t --minimum_length 100 --database ./databases/dfast/uniprot_bacteria-0.9.ref -o ./2_analysis/gene_annotation/dfast/one.candidate.uniprot.bacteria --config custom_config.py
```
For another candidate with slightly different options:
```bash
#$ -S /bin/bash
#$ -N dfastjob
#$ -cwd
#$ -pe smp 20
#$ -l h_vmem=6G

conda activate dfast

dfast -g GAGA-0024.Scaffold1.16539177-16545088.fa --force --metagenome --cpu 20 --debug --use_original_name t --minimum_length 100 --database ./databases/dfast/uniprot_bacteria-0.9.ref -o ./2_analysis/gene_annotation/dfast/GAGA-0024.Scaffold1.16539177-16545088.uniprot.bacteria --config custom_config.py
```

For all candidates together:

Write a Script to run the process:

nano `batch_dfast_job.sh`:

```bash
#$ -S /bin/bash
#$ -N batchdfastjob
#$ -cwd
#$ -w e
#$ -V
#$ -pe smp 20
#$ -l h_vmem=6G

conda activate dfast

echo "Running on file: $file"

dfast -g $file --force --metagenome --cpu 20 --debug --use_original_name t --minimum_length 100 --database ./databases/dfast/uniprot_bacteria-0.9.ref -o ./2_analysis/gene_annotation/dfast/$file --config custom_config.py
```

Submit the script and parse bash variables to the script using `-v file="file1.fa"`

To submit jobs more easily, use `parallel`:
```bash
find . -name "GAGA*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh

find . -name "NCBI*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh

find . -name "OUT*" | parallel -I% --max-args 1 qsub -v file="%" batch_dfast_job.sh
```

Control, that all candidates have a completed pipeline:
```bash
cd ./2_analysis/gene_annotation/dfast

find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/cds.fna ]] && echo "$dir"; done | wc -l
```
If all candidates have successfully completed the pipeline, no more candidate directories should be showing up here.

## 2.2 Prokaryotic gene annotation with `prodigal`

### 2.2.1 Run `prodigal` on the cluster
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
prodigal -p meta -i all.candidates.fa -a protein.translations.faa -d annotated_genes_prodigal.cds -w prodigal.statistics -f gff -o \
./2_analysis/gene_annotation/prodigal/prodigal_annotated_genes.gff
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

Number of predicted HGTs which have an annotation with `prodigal` and `dfast`:
```bash
# prodigal
cat protein.translations.faa | grep ">" | wc -l

#dfast
cat protein.faa | grep ">" | wc -l
```

## 2.3 Prokaryotic gene annotation with `kraken`

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
