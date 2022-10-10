# Selection analysis of selected LGTs using HyPhy and RELAX

Here, we will test for signatures of selection acting at different branches among the gene trees of LGTs.

1. Identify the 20 best BLAST hits in bacterial genomes
2. Make an alignment: Retrieve and align the coding sequence of genes of interest
   -> Using mafft-linsi or PRANK (msa)
3. Check alignment quality using zorro and visualize the alignment quality with seaview
   
4. Clean the alignments and compute a solid gene tree with IQtree
   
5. Test the alignment for evidence of recombination
   
6. Test for evidence of positive selection using HyPhy
   
7. Identify branches that show significant signatures of positive selection

After that, we will test for selection in other ants with relevant LGTs.


### Calculate selection for Etherase LGTs

Compare the etherases from all ant species and the closest bacterial hit plus two outgroups and test for signatures of positive selection across different branches in the phylogeny. 

Our working hypothesis is that etherases were adapted in the respective ants and allowed them to use the bacterially-derived gene.

1. Make an alignment of the coding sequence (nucleotide sequence) of the ant and bacterial lysozyme genes.

Extract all necessary CDS for the ant etherases (8 ant species with etherases):

Coordinates of etherases:
```bash
GAGA-0221	Scaffold3	8819326	8820234
NCBI-0005	NW_020229769.1	105240	104332	
GAGA-0361	Scaffold4	3723826	3724734	
GAGA-0362	Scaffold15	5724715	5725623	
GAGA-0396	Scaffold28	2026308	2027216
GAGA-0200	Scaffold111	363210	364118
GAGA-0374	Scaffold20	2858966	2859874
GAGA-0360	Scaffold16	5909052	5909960
```

```bash
samtools faidx genome Scaffold:start-stop > /global/scratch2/j_rink02/master/lgt/2_analysis/selection_analysis/Etherases/etherases.fa
```

Extract CDS for the bacterial etherase from closest bacterial hit and from two outgroups from NCBI BLAST:

```bash
#closest bacterial hit: Spiroplasma sp.
elink -db protein -id "APE13601.1" -target nuccore | efilter -molecule mrna| efetch -format fasta_cds_na >> /global/scratch2/j_rink02/master/lgt/2_analysis/selection_analysis/Etherases/etherases.fa

#Outgroup 1: Paraclostridium bifermentans
elink -db protein -id "UAG17450.1" -target nuccore | efilter -molecule mrna| efetch -format fasta_cds_na >> /global/scratch2/j_rink02/master/lgt/2_analysis/selection_analysis/Etherases/etherases.fa

#Outgroup 2:
elink -db protein -id "QQK08127.1" -target nuccore | efilter -molecule mrna| efetch -format fasta_cds_na >> /global/scratch2/j_rink02/master/lgt/2_analysis/selection_analysis/Etherases/etherases.fa

```