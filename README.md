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

### 1.1 Atom Scripts
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

### 2.1 bash Scripts
