# Horizontal Gene Acquisition in Ants

## Introduction to the Code Repositories

This repository holds all the code to evaluate horizontally acquired genes (HGAs) deriving from bacteria in 163 ant genomes, previously identified by an automated Horizontal Gene Transfer (HGT) finder pipeline. All the code for the automated detection of HGTs/LGTs in ant genomes can be found [HERE](https://github.com/dinhe878/GAGA-Metagenome-LGT).


This project is a part of the Global Ant Genomics Alliance (GAGA).

## Information about Code & Scripts

### 01. Markdown Scripts (`.md`)
This folder holds all markdown scripts for different analyses to be conducted within a Linux bash coding environment.

#### Calculation of reads overlapping HGT boundaries
`bedtools` and `samtools` were used to map raw sequencing reads to the respective genomes and to calculate reads overlapping the start, end, and total sequence length of HGT candidates.
```bash
#Calculate extension based on average read length distribution, count one read value per boundary
./markdown_scripts/01_countBoundaryReads.md
# Read counts only for remaining HQ HGT candidates, count 3 values: start, end, complete HGT + extension
./markdown_scripts/02_count3ValueBoundaryReads_HQ_LGTs.md
```

#### Prokaryotic gene annotation 
`prodigal`, `kraken2`, and `dfast` were run to annotate putative prokaryotic HGT sequences in ant genomes. For all downstream analyses, only `dfast` results were used.

```bash
./markdown_scripts/03_Prokaryotic_Gene_Annotation_dfast.md
```

#### Check completeness of HGT candidates
Input files are files obtained from `dfast`.
```bash
./markdown_scripts/04_HGT_completeness_check.md
```

#### Re-annotation of incomplete HGT candidates
All candidates were extended by 1000 bp in 5' and 3' direction and re-annotated with `dfast` to identify HGTs which were annotated incompletely.
```bash
./markdown_scripts/05_reannotation.DFAST.md
```

#### Check for additional candidates within all GAGA genomes
We used this script to identify putative HGTs that were previously filtered out due to our strict filtering criteria. All resulting candidates were intersected with initially predicted HGT candidates by the automated HGT finder pipeline.
```bash
./markdown_scripts/06_Blast_GAGA_genomes.md
```

#### Running gene trees on HGTs of interest
We ran gene trees on the clade-specific prokaryotic protein HGT sequences (Lysozymes, MurNAc etherases, CFA synthases) and their closest blast homologs to conduct phylogenetic analyses.
```bash
./markdown_scripts/07_Run_GeneTrees.md
```

#### Checking for selection among clade-specific HGTs
We checked for selection using `HyPhy` on the clade-specific HGTs (Lysozymes, MurNAc etherases, CFA synthases).
```bash
./markdown_scripts/08_Selection_analyses.md
```
### 02. R Scripts (`.Rmd`)
All scripts have been written in RStudio 4.3.1.

#### GAGA HGT finder analysis
Identification of HGTs with protein-coding hits and calculation of overview plots for all HGT candidates detected by the automated HGT finder pipeline (`https://github.com/dinhe878/GAGA-Metagenome-LGT`). These overview plots were used as a starting point for all downstream analyses and to determine HGT quality, as well as filters, for the identification of true HGTs.
```bash
./r_scripts/01_analyseLGTs.Rmd
```

#### HGT filtering of predicted candidates with pre-defined thresholds
This script contains all filtering steps which were used to exclude false-positives and filter the predicted HGTs. Filters were determined based on manual assessment of HGT candidates from seven randomly chosen GAGA genomes, as well as plotted parameter distributions of all predicted HGTs.
```bash
./r_scripts/02_LGTfiltering.Rmd
```

#### Visualization of average read length in PacBio and stLFR genomes
We calculated average read lengths for both PacBio and stLFR genomes separately. This was used as a basis to determine HGT boundary expansions in the `./markdown_scripts/01_countBoundaryReads.md` script to calculate read overlaps accordingly. 
```bash
./r_scripts/03_Distribution_plots_avg_read_length.Rmd
```

#### Calculate expression of HGTs and candidate completeness
Here, we obtained a summary CDS-level table, as well as a locus-level table with summarized information, including expression based on RNAseq data, for all high-quality HGT candidates.
```bash
./r_scripts/04_calculateExpression.Rmd
```

#### Plot coverage of gene expression for all HGTs
We used unique read counts as input and added up all read counts for all life stages within each respective HGT candidate. RNAseq coverage was visualized for each HGT candidate.

```bash
./r_scripts/05_plotRNAseqCoverage.Rmd
```

#### Obtain all information from dfast gffs for candidates
Query coverage and subject coverage were added as additional information for every HGT candidate to estimate completeness and functionality of the HGT.
```bash
./r_scripts/06_screenDfastGFFs.Rmd
```

#### Visualize HGT candidates and plot GAGA phylogeny with HGTs
This script was used to plot all HGTs detected by the automatic HGT finder pipeline, which is visualized in **Figure 1**. 
```bash
./r_scripts/07_map2phylogeny.Rmd
```

#### Retrieve additional information from UniProt
We used the package *uniprotR* to obtain additional information for all HGT candidates and integrated this information into our summarized tables for functional analyses and as a basis for in-depth analyses of clade-specific LGTs.
```bash
./r_scripts/08_retrieveUniprotInfo.Rmd
```

#### Synteny of clade-specific HGTs
All clade-specific HGT regions were extended by 40 kb up- and downstream to infer synteny of genes.
```bash
./r_scripts/09_Synteny.Rmd
```

#### Detail analyses of unique and species-specific HGTs and creation of plots for global HGTs
This script was used to investigate unique and species-specific HGTs.
```bash
./r_scripts/10_Unique_HGTs.Rmd 
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
Script used to create a pie chart visualization of the bacterial donors of HGTs. 
```bash
./python_scripts/1_Bacterial_distribution_piechart.ipynb
```