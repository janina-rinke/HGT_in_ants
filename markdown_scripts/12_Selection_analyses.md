# Selection analysis of selected LGTs using HyPhy and RELAX

Here, we will test for signatures of selection acting at different branches among the gene trees of LGTs.


### Calculate selection for Etherase LGTs

Compare the etherases from all ant species and the closest bacterial hit plus two outgroups and test for signatures of positive selection across different branches in the phylogeny. 

Our working hypothesis is that etherases were adapted in the respective ants and allowed them to use the bacterially-derived gene.

### 1. Make a FASTA file of the coding sequence (nucleotide sequence) of the ant etherase genes.

#### 1.1 Extract all necessary CDS for the ant etherases (8 ant species with etherases):

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
samtools faidx genome Scaffold:start-stop > /global/scratch2/j_rink02/master/lgt/2_analysis/selection_analysis/Etherases/etherases.cds.fa
```

After extracting all sequences, the FASTA file with the nucleotide sequences should look like this:
!You need to make sure that all sequences are in the same orientation!

```bash
>Cfel
TTATTTACTATTAATAATATCAATTATTCGATTGCTATGTTGTGCTAATAATTTTTGACT
TGTTTTAAAATCAACATTCTTTAATAACATAATAATTGCATGTTTACAAGAATAATTACA
TGCATTTAAAGCTTTCTCAATTTCTTCTTGTGAACAATCTGTTACTTCACTAACAATTTT
AGCTGTTCTAACTTTTAATTTTTCATTTGTTGCCACCAAATCAACCATTAAATTTTGATA
AACCTTTCCTGATTTAATCATTAATGTTGTCGATATCATATTACATACTAATTTTGTTGC
TGTCCCACTTTTCATTCTTGTACTTCCTGTGACAACTTCTGACCCTGTTTCAATTATAAT
TGTTTCATCAGCAATTTCTGACATTTTTGAATTTTTTGTCATACATAATCCAATTGCTAA
ACCATTAATTGATTTACAATATTTTAAGCCAGTAAGAACATATGGTGTTCTTCCTGATGC
AGCAATTCCAATCACAGTATCTAATGAATTTAAGTTACGATCTTTTAAATCTTTAATAGC
TAAATTACTATTATCTTCCGCACCTTCAATTGGCATTCTAATTGCAATATCACCACCTGC
TATTAATCCCAGAATACGATCATTGCCAACATTATATGTGGGTAACATCTCTGAAGCATC
TAACACACCAAGTCTACCTGATGTTCCTGCTCCAATATAAATTAATCTCCCACCTTTTGT
AAAACGATCAAAAATTAAATCAATCACCTTGGCAATGATTGGAATTTTATCTTTAATTGC
TCCTGCTATTTTAGCATCTTCATTATTGATTATTGTTAAAATTCCAGTTGTATTTTCTTT
ATCAATTGATAAACTATTTTGATTACGTTGTTCTGTATCAATCTTGGTTAGATCAATTTT
CATGTTCAT
>Cflo
TTA....
```
#### 1.2 Remove stop codons and anything that is not a multiple of three
To obtain dNdS ratios and for prank and hyphy to work correctly, we need to transform the FASTA CDS sequences into a codon-fasta file with a custom perl script.
```bash
perl replace_stop_codons.pl etherases.cds.fa etherases_codoncorrect.cds.fa 
```

The corrected fasta cds file now looks like this:
```bash
>Colo
TTATTTACTATTAATAATATCAATTATTCGATTGTTATGTTGTACNNNNNNTTTTTGACTTGTTTTAAAATCAATATTCTTNNNNNNCATAACAATTGCATGTTTACAAGAATAATTACATGCATTNNNAGCTTTCTCAGCTTCTTCTTGNNNACAATTTGTTACTTCACTAACAATTTTAGCTGTTCTAACCCTNNNTTTTTCATTTGTTGCCACCAAATCAATCATNNNATTCTGATAAACCTTTCCNNNTTTAATCATNNNTGTTGTNNNTATCATATTACATACNNNTTTTGTTGCTGTTCCACTTTTCATTCTTGTACTTCCTGTGATAACTTCNNNCCCTGTTTCAATTATAATTGTTTCATCAGCAATTTCNNNCATTTTNNNATTTTTTGTCATACANNNTCCAACTGCNNNACCATTAATNNNTTTACAATATTTNNNACCAGTAAGAACATATGGTGTTCTTCCNNNTGCAGCAATTCCAATGACAGTATCNNNNNNATTNNNATTGTGATCTTTNNNATCTTTAATAGCNNNATTATTATTATCTTCTGCACCTTCAATTGGTATTCTAATCGCAATATCACCACCTGCTATNNNTCCCAGAATACGATCACTGCCAACATTATATGTTGGNNNCATCTCNNNAGCATCNNNCACACCAAGTCTACCNNNTGTTCCTGCTCCAATATAAATNNNCCTCCCGCCTTTTGTAAAACGATCAAAGATNNNATCAATCACTTTGGCAATGATNNNGATTTTATCTTTAATTGCTCCTGCTATTTTAGCATCTTCATTATTGATTATTCTNNNAATTTCAATTGTATTTTCTTTATCAATNNNNNNACTATTTTGATTACGTTGTTCTGTCTCAATTTTGGTNNNATCAATTTTCATGTTCAT
>Cfed
TTATTTACTATTAATAATATCAATTATTCGATTGCTATGTTGTGCNNNNNNTTTTTGACTTGTTTTAAAATCAACATTCTTNNNNNNCATAATAATTGCATGTTTACAAGAATAATTACATGCATTNNNAGCTTTCTCAATTTCTTCTTGNNNACAATCTGTTATTTCACTAACAATTTTAGCTGTTCTAACTTTNNNTTTTTCATTTGTTGCCACCAAATCAACCATNNNATTTTGATAAACCTTTCCNNNTTTAATCATNNNTGTTGTCGATATCATATTACATACNNNTTTTGTTGCTGTCCCACTTTTCATTCTTGTACTTCCTGTGACAACTTCNNNCCCTGTTTCAATTATAATTGTTTCATCAGTAATTTCNNNCATTTTNNNATTTTTTGTCATACANNNTCCAATCGTNNNACCATTAATNNNTTTACAATATTTNNNGCCAGTAAGAACATATGGTGTTCTTCCNNNTGCAGCAATTCCAATCACAGTATCNNNNNNATTNNNGTTATGATCTTTNNNATCTTTAATAGCNNNATTACTATTATCTTCCGCACCTTCAATTGGCATTCTAATTGCAATATCACCACCTGCTATNNNTCCCAGAATACGATCATTACCAACATTATATGTGGGNNNCATCTCNNNAGCATCNNNCACACCAAGTCTACCNNNTGTTCCTGCTCCAATATAAATNNNTCTCCCACCTTTCATAAAACGATCAAAAATNNNATCAATTACCTTGGCAATGATTGGAATTTTATCTTTAATTGCTCCTGCTATTTTAGCATCTTCATTATTGATTATTGTNNNAATTCCAGTTGTATTTTCTTTATCAATNNNNNNACTATTTTGATTACGTTGTTCTGTATCAATCTTGGTNNNATCAATTTTCATGTTCAT
```


### 2 Make an alignment from the CDS FASTA file and compute a gene tree
For this, we will use a custom script called `dNdS.sh` which takes our corrected CDS FASTA file as input, makes the alignment and computes a gene tree (`IQtree`) using the alignment

```bash
bash dNdS.sh
```

Check all files which have been produced by the bash script:
```bash
dNdS.sh							etherases_codoncorrect.cds.aln.best.fas.iqtree
etherases.cds.fa					etherases_codoncorrect.cds.aln.best.fas.log
etherases.pep.fa					etherases_codoncorrect.cds.aln.best.fas.mldist
etherases_codoncorrect.cds.aln.best.fas			etherases_codoncorrect.cds.aln.best.fas.model.gz
etherases_codoncorrect.cds.aln.best.fas.bionj		etherases_codoncorrect.cds.aln.best.fas.splits.nex
etherases_codoncorrect.cds.aln.best.fas.ckp.gz		etherases_codoncorrect.cds.aln.best.fas.treefile
etherases_codoncorrect.cds.aln.best.fas.contree		etherases_codoncorrect.cds.fa
```

Now we have a phylogenetic tree and an alignment and we can start running some selection analyses.

### 3 Test for recombination using `hyphy SBP`

Provide the following options to the menu prompts:
1. Enter `14` to test an alignment for recombination.
2. Enter `3` to search an alignment for a single breakpoint using the SPB method.
3. Enter `2` as we are using codon data.
4. Enter the absolute path to your alignment file: `/path/to/data/files/etherases_codoncorrect.cds.aln.best.fas` .
5. Enter `1` to select the universal codon table.
6. Enter `1` to use only (c)AIC to measure goodness of fit and skip the statistical inference using KH
resampling.
7. Enter `MG94CUSTOMCF3X4` . This selects the Muse-Gaut 94 model crossed with a nucleotide substituion model, using a corrected CF3x4 codon frequency estimator.
8. Enter `2` to select that model parameters should be shared by all branches, branch lengths are estimated independently.
9. Enter `012345` to set the nucleotide substitution model to GTR (general time reversible) with 5 independent rate categories.
  
10. Enter a file where the results should be stored, e.g. `/path/to/data/files/sbpResults/runRecombination.SBP` .

Output from `sbp`:
```bash
Best supported breakpoint is located at position 138
AIC = 2926.383928344459 : an improvement of 7.562666873956005 AIC points

AIC-c

Best supported breakpoint is located at position 138
AIC = 2934.206150566681 : an improvement of 2.42595701922528 AIC points

BIC

There seems to be NO recombination in this alignment
```

We do find evidence for recombination at position 138 of the alignment.

Checking the alternative phylogenies from file sbpResults_cAIC.splits:
```bash
0-138
((Colo,Pila),(((Cfel,Cfed),Cspec),Cjap),(Cflo,Csin))
139-908
((((Colo,Cjap),Pila),(Cfel,Cflo)),(Cfed,Csin),Cspec)
```

##### Estimate a single alignment-wide ω
1. Enter `5` to select Basic Analyses.
2. Enter `1` to select "Analyze codon data with a variety of standard models using a given tree".
3. Enter `1` to select the universal genetic code.
4. Enter the absolute path to your alignment file: `/path/to/data/files/etherases_codoncorrect.cds.aln.best.fas` .
5. The screen then shows a variety of model options to fit to this dataset. We will use HyPhy’s preferred model MG94xREV. Select this option by entering `MG94CUSTOMCF3X4` .
6. Enter `2` , which tells HyPhy that the estimated ω should be global, i.e. shared across the entire alignment tree.
7. Enter the string `012345` , which tells HyPhy to use the GTR mutation model in the MG94 model.
8. Enter the absolute path to your IQtree file:
`/path/to/data/files/etherases.aln.best.fas.treefile` .
9. Enter `1` to tell HyPhy to use maximum likelihood to estimate and optimize branch lengths on the provided tree.


```bash
Current Max: -1445.2395     (92 % done) LF Evals/Sec: 2860    CPU Load: 4.92    
______________RESULTS______________
Log Likelihood = -1445.23953088664;
Global Parameters:
CT=1.243129585435661;
CG=0.2405341567135409;
GT=0.02473892222772022;
AC=0.02357639595615491;
R=4.063556035280313;
AT=0.03930980742771926;

Tree givenTree=((((Colo:0.03941250355226409,Cjap:0.009184983743003177)Node4:0.003708049796125132,Pila:0.03221354672533978)Node3:0.004287686466103441,Cflo:0.007726856025186305)Node2:0,Csin:0.01031822637827303,(Cfel:0.002531559596540717,(Cfed:0.008980788797870651,Cspec:0.007714654868827347)Node12:0)Node10:0.001270054233219844);
```

The alignment wide ω is: `R=4.064`.

### Run BUSTED to test for positive selection in a "branch-site manner"

###### Full tree analysis
```bash
 hyphy2.5 busted --alignment etherases.aln.best.fas --tree etherases.aln.best.fas.treefile --output BUSTED.full.json
 ```

###### Performing the full (dN/dS > 1 allowed) branch-site model fit

```bash
>Save BUSTED model fit to this file ['/dev/null' to skip] (`/global/scratch2/j_rink02/master/lgt/2_analysis/selection_analysis/Etherases/`) BUSTED3.full.json
* Log(L) = -1440.85, AIC-c =  2956.88 (37 estimated parameters)
* For *test* branches, the following rate distribution for branch-site combinations was inferred

|          Selection mode           |     dN/dS     |Proportion, %|               Notes               |
|-----------------------------------|---------------|-------------|-----------------------------------|
|        Negative selection         |     0.999     |    0.000    |       Not supported by data       |
|         Neutral evolution         |     1.000     |   41.607    |       Collapsed rate class        |
|      Diversifying selection       |     6.473     |   58.393    |                                   |

## Branch-site unrestricted statistical test of episodic diversification [BUSTED]
Likelihood ratio test for episodic diversifying positive selection, **p =   0.0000**.
```

A proportion of sites (~58.4 %) is evolving with ω > 1 in a subset of branches. The reported p-value returned p = 0.
Probably, positive selection has left a signature somewhere in this gene during the evolution of the different species in the phylogeny.

### Run aBSREL to test for signatures of positive selection 
aBSREL = adaptive branch site random effects likelihood

aBSREL will successively scan the phylogeny for all the branches where selection may have operated.

##### aBSREL: Full tree analysis 
```bash
### Testing selected branches for selection

|              Branch               |  Rates   |     Max. dN/dS     |      Test LRT      |Uncorrected p-value |
|-----------------------------------|----------|--------------------|--------------------|--------------------|
|               Colo                |     1    |   2.48 (100.00%)   |        4.68        |       0.03495      |
|               Pila                |     1    |   6.74 (100.00%)   |       13.17        |       0.00047      |
|               Csin                |     1    |   1.52 (100.00%)   |        0.33        |       0.36692      |
|               Cjap                |     1    |  >1000 (100.00%)   |        7.69        |       0.00752      |
|               Cfed                |     2    |   51.48 ( 2.78%)   |        5.61        |       0.02170      |
|               Cspec               |     1    |  >1000 (100.00%)   |        7.74        |       0.00734      |
|               Cflo                |     1    |   4.52 (100.00%)   |        2.47        |       0.11051      |
|               Node3               |     1    |  >1000 (100.00%)   |        4.13        |       0.04637      |
|               Node2               |     1    |  >1000 (100.00%)   |        3.41        |       0.06758      |
|               Cfel                |     1    |  >1000 (100.00%)   |        2.58        |       0.10424      |
|               Node5               |     1    |  >1000 (100.00%)   |        1.29        |       0.20826      |
----
### Adaptive branch site random effects likelihood test 
Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p =   0.0500_ found **1** branches under selection among **11** tested.

* Pila, p-value =  0.00519
```
