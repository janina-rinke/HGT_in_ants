# HGT in Ants

## Introduction to the Code Repositories

This repository holds all the code to evaluate HGTs (synonym: LGTs) within ant genomes (subproject of the Global Ant Genomics Alliance).
In here, you can find five different folders storing all the data to analyse horizontal gene transfers from bacteria to ants.

## Information about Code & Scripts

Local and remote copies of files and folders: 

```bash
# mac
macbase=/Users/Janina/sciebo/MASTER/Masterarbeit/lgt
# cluster
cluster=/global/scratch2/j_rink02/master/HGT
```
### 01. Markdown Scripts (`.md`)
This folder holds all markdown scripts for different analyses within the Linux environment.

#### MMseqs2
`mmseqs2` was run to annotate putative prokaryotic HGT sequences in ant genomes.
```bash
./markdown_scripts/01_Blast_mmseqs.md
```

#### DFAST

DFAST was run as another prokaryotic gene annotation tool.
```bash
./markdown_scripts/02_Prokaryotic_Gene_Annotation_dfast.md
```

#### Calculation of reads overlapping HGT boundaries

I used bedtools and samtools to calculate the reads overlapping the start, end and total sequence length of HGT candidates.
```bash
./markdown_scripts/03_countBoundaryReads.md
```

Read counts only for remaining HQ HGT candidates
```bash
./markdown_scripts/04_ReadCounts_FilteredCandidates.md
```

#### Check completeness of HGT HGT candidates
Input files are files obtained from `dfast`.
```bash
./markdown_scripts/05_HGT_completeness_check.md
```

#### Re-annotation of incomplete HGT candidates
All candidates were extended by 1000 bp in 5' and 3' direction and re-annotated with `dfast`.
```bash
./markdown_scripts/06_reannotation.DFAST.md
```

#### Check for additional candidates within all GAGA genomes
We used this script to identify putative HGTs that were previously filtered out due to our strict filtering criteria. All resulting candidates were intersected with initially predicted HGT candidates by the automated HGT finder pipeline.
```bash
./markdown_scripts/07_Blast_GAGA_genomes.md
```

#### Running gene trees on HGTs of interest
We ran gene trees on the clade-specific prokaryotic protein HGT sequences (Lysozymes, MurNAc etherases, CFA synthases) and their closest blast homologs to conduct phylogenetic analyses.
```bash
./markdown_scripts/08_Run_GeneTrees.md
```

#### Checking for selection among clade-specific HGTs
We checked for selection using HyPhy on the clade-specific HGTs (Lysozymes, MurNAc etherases, CFA synthases).
```bash
./markdown_scripts/09_Selection_analyses.md
```

### 02. Bash Scripts (`.sh`)
Contains scripts written in bash which can directly be submitted to any cluster.

Blast any HGT candidate against a bacterial database and obtain the 20 best blast hits to create a phylogeny
```bash
./bash_scripts/HGT_pipeline.sh
./bash_scripts/proteinBlast.sh
```
### 03. Perl Scripts (`.pl`)
Code to analyse RNAseq reads obtained from GAGA and conduct the RNAseq mapping for all high-quality HGT candidates to investigate gene expression.

#### Modify all gff files obtained from dfast to obtain coordinates in reference to the genomes
```bash
./perl_scripts/1_Modify_gffs.pl
```

#### Run RNAseq mapping with StringTie
```bash
./perl_scripts/2_run_rna_mapping_wholegenome.pl
```
### 04. Python Scripts (`.py`)
This folder contains a script to plot a pie chart of all bacterial donors of HGTs.

The script was used to create a donut plot of the bacterial donors of HGTs, which is visualized in Figure 2. 
```bash
./python_scripts/1_Bacterial_distribution_piechart.ipynb
```
### 05. R Scripts (`.Rmd`)
In this folder, the main analyses for HGT candidates can be found. All scripts have been written in RStudio 4.3.1.

#### GAGA HGT finder analysis
Identification of HGTs with protein-coding hits and calculation of overview plots for all HGT candidates detected by the automated HGT finder pipeline. These overview plots were used as a starting point for all downstream analyses and to determine HGT quality, as well as filters for the identification of true HGTs.
```bash
./r_scripts/01_analyseLGTs.Rmd
```

#### HGT filtering of predicted candidates with pre-defined thresholds
This script contains all filtering steps which were used to exclude false-positives and filter the predicted HGTs. Filters were determined based on manual assessment of HGT candidates from seven randomly chosen GAGA genomes, as well as plotted parameter distributions of all predicted HGTs.
```bash
./r_scripts/02_LGTfiltering.Rmd
```

#### Prepare candidates for BLAST searches
```bash
./r_scripts/03_Preparation_blast_candidates.Rmd
```
#### Visualization plots of HGT candidates
```bash
./r_scripts/04_Avg_length_plot.Rmd
./r_scripts/05_Distribution_plots_avg_read_length.Rmd
```

#### Visualize HGT candidates and plot GAGA phylogeny with HGTs
This script was used to plot all HGTs detected by the automatic HGT finder pipeline, which is visualized in Figure 1. 
```bash
./r_scripts/06_map2phylogeny.Rmd
```

#### Obtain all information from dfast gffs for candidates
```bash
./r_scripts/07_screenDfastGFFs.Rmd
```

#### Phylogeny creator for single HGT candidates with information from UniProt
```bash
./r_scripts/08_Gene_trees_HGTs.Rmd
```

#### Plot RNAseq coverage for all candidates using data from RNAseq mapping
We used unique read counts as input and added up all read counts for all life stages within each respective HGT candidate.

```bash
./r_scripts/09_plotCoverage.Rmd
```

#### Calculate the expression and completeness of Candidates
Here, we obtained a summary CDS-level table with information for all candidates, as well as a locus-level table with summarized information for all high-quality HGT candidates.
```bash
./r_scripts/10_calculateExpression.Rmd
```

#### Retrieve additional information from UniProt
We used the package *uniprotR* to obtain additional information for all HGT candidates and integrated this information into our summarized tables.
```bash
./r_scripts/11_retrieveUniprotInfo.Rmd
```

####Comparison of original annotation and re-reannotation
We visualized this comparison as a CDS length extension plot.
```bash
./r_scripts/12_ReannotationPlots.Rmd
```

#### Synteny of clade-specific HGTs
```bash
./r_scripts/13_Synteny.Rmd
```

#### Detail analyses of unique and species-specific HGTs and creation of plots for global HGTs
```bash
./r_scripts/14_Unique_HGTs.Rmd 
```
