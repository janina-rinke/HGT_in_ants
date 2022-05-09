# HGT in Ants

**Introduction to the Code Repositories**
This repository holds all the code to evaluate LGTs within ant genomes (subproject of the Global Ant Genomics Alliance).
In here, you can find four different folders storing all the data to analyse horizontal gene transfers from bacteria to ants.

## 1. Atom Scripts
This folder holds all initial workflows and Markdown scripts for different analyses within the Linux environment, as well as useful information about data and file types.

## 2. bash Scripts
This folder holds scripts written in bash which can directly be submitted to the GridEngine Cluster.

## 3. Perl Scripts
This folder holds the code to analyse RNAseq reads and conduct the RNAseq mapping for all high-quality LGT candidates.

## 4. Python Scripts
This folder contains a script to plot a pie chart of all bacterial donors of LGTs.

## 5. R Scripts
This folder contains scripts that have been written in RStudio 4.3.1.

**Information about code & scripts**

Local and remote copies of files and folders

```bash
# mac
macbase=/Users/Janina/sciebo/Master.LGTs
```

### 1.1 Atom scripts
General & Useful Information about data and file types
```bash
./atom_scripts/1_Data_and_file_types.md
```
Initially used commands to work with PacBio reads within the targeted ant genomes.
```bash
./atom_scripts/2_LGT_Radanalysis.md
```

Initial starting point to calculate coverage of LGT Candidates
```bash
./atom_scripts/3_Calculation_Coverage.md
```

#### MMseqs2
`mmseqs2` was run to annotate putative prokaryotic LGT sequences in ant genomes.
```bash
./atom_scripts/4_Blast_mmseqs.md
```

#### DFAST

DFAST was run as another prokaryotic gene annotation tool.
```bash
./atom_scripts/5_Prokaryotic_Gene_Annotation_dfast.md
```

#### Calculation of reads overlapping LGT boundaries

I used bedtools and samtools to calculate the reads overlapping the start, end and total sequence length of LGT candidates.
```bash
./atom_scripts/6_countBoundaryReads.md
```

Read counts only for remaining HQ LGT candidates
```bash
./atom_scripts/7_ReadCounts_FilteredCandidates.md
```

#### Check completeness of HGT LGT candidates
Input files are files obtained from `dfast`.
```bash
./atom_scripts/8_LGT_completeness_check.md
```

#### Re-annotation of incomplete LGT candidates
All candidates were extended by 1000 bp in 5' and 3' direction and re-annotated with `dfast`.
```bash
./atom_scripts/9_reannotation.DFAST.md
```

#### Check for additional candidates within all GAGA genomes
We used this script to identify putative LGTs that were previously filtered out due to our strict filtering criteria. All resulting candidates were intersected with initially predicted LGT candidates by the automated LGT finder pipeline.
```bash
./atom_scripts/10_Blast_GAGA_genomes.md
```

#### Running gene trees on LGTs of interest
We ran gene trees on prokaryotic protein LGT sequences and their closest blast homologs to conduct phylogenetic analyses.
```bash
./atom_scripts/11_Run_GeneTrees.md
```

### 2.1 Bash scripts
Blast any LGT candidate against a bacterial database and obtain the 20 best blast hits to create a phylogeny
```bash
./bash_scripts/LGT_pipeline.sh
./bash_scripts/proteinBlast.sh
```

### 3.1 Perl scripts

#### Modify all gff files obtained from dfast to obtain coordinates in reference to the genomes
```bash
./perl_scripts/1_Modify_gffs.pl
```

#### Run RNAseq mapping with StringTie
```bash
./perl_scripts/2_run_rna_mapping_wholegenome.pl
```

### 4.1 Python scripts
Create a donut chart for bacterial origin of LGTs
```bash
./python_scripts/1_pie_chart.py
```

### 5.1 R scripts
In this folder, the main analyses for LGT candidates can be found. All scripts have been written in RStudio 4.3.1.

#### GAGA LGT finder analysis
```bash
./r_scripts/0_analyseLGTs.Rmd
```

#### Script to plot LGT candidates from automated prediction
```bash
./r_scripts/1_Candidate_plot.Rmd
```

#### LGT filtering of predicted candidates with pre-defined thresholds
```bash
./r_scripts/2_LGTfiltering.Rmd
```

#### Prepare candidates for BLAST searches
```bash
./r_scripts/3_Preparation_blast_candidates.Rmd
```
#### Visualization plots of LGT candidates
```bash
./r_scripts/4_Avg_length_plot.Rmd
./r_scripts/5_Distribution_plots_avg_read_length.Rmd
```

#### Visualize LGT candidates and plot GAGA phylogeny with LGTs
```bash
./r_scripts/6_map2phylogeny.Rmd
```

#### Obtain all information from dfast gffs for candidates
```bash
./r_scripts/7_screenDfastGFFs.Rmd
```

#### Phylogeny creator for single LGT candidates with information from UniProt
```bash
./r_scripts/8_Gene_trees_LGTs.Rmd
```

#### Plot RNAseq coverage for all candidates using data from RNAseq mapping
We used unique read counts as input and added up all read counts for all life stages within each respective HGT candidate.

```bash
./r_scripts/9_plotCoverage.Rmd
```

#### Calculate the expression and completeness of Candidates
Here, we obtained a summary CDS-level table with information for all candidates, as well as a locus-level table with summarized information for all high-quality HGT candidates.
```bash
./r_scripts/10_calculateExpression.Rmd
```

#### Retrieve additional information from UniProt
We used the package *uniprotR* to obtain additional information for all LGT candidates and integrated this information into our summarized tables.
```bash
./r_scripts/11_retrieveUniprotInfo.Rmd
```

####Comparison of original annotation and re-reannotation
We visualized this comparison as a CDS length extension plot.
```bash
./r_scripts/12_ReannotationPlots.Rmd
```

#### Synteny of clade-specific LGTs
```bash
./r_scripts/13_Synteny.Rmd
```

#### Detail analyses of unique and species-specific LGTs and creation of plots for global LGTs
```bash
./r_scripts/14_Unique_LGTs.Rmd 
```
