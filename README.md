# Horizontal Gene Transfer in Ants

## Introduction to the Code Repositories

This repository holds all the code to evaluate horizontally transferred genes (HGTs) from prokaryotic donors in 163 ant genomes, previously identified by an automated Horizontal Gene Transfer (HGT) finder pipeline. All the code for the previous, automated detection of potential HGTs/LGTs in ant genomes can be found [HERE](https://github.com/dinhe878/GAGA-Metagenome-LGT).


## 1 HGT detection
The detection of potential HGTs with protein-coding hits in the 163 ant genomes of interest obtained by our bioinformatic HGT finder pipeline is explained in detail at (`https://github.com/dinhe878/GAGA-Metagenome-LGT`). 

This HGT finder pipeline yielded 13,664 initial HGT candidates for downstream analyses. We produced overview plots for each candidate to aid the determination of HGT quality and processing of HGT candidates in downstream analyses:

[*01.analyseLGTs.Rmd*](01_HGT_detection/01.analyseLGTs.Rmd)



#### HGT filtering of predicted candidates with pre-defined thresholds
This script contains all filtering steps which were used to exclude false-positives and filter the predicted HGTs. Filters were determined based on manual assessment of HGT candidates from seven randomly chosen GAGA genomes, as well as plotted parameter distributions of all predicted HGTs.
```bash
./r_scripts/02_LGTfiltering.Rmd
```
## 2 Validation and quality assessment of HGT candidates


### Calculation of reads overlapping HGT boundaries

[*01.countBoundaryReads.sh*](01_HGT_validation/01.countBoundaryReads.sh)

We used `bedtools` and `samtools` to map raw PacBio sequencing reads (from `bam` files) to the respective genomes and to calculate reads overlapping the start, end, and total sequence length of HGT candidates. All coordinates were stored in bed files (`LGTs.candidateloci.loose.bed`). HGT boundaries were expanded based on the average read length distribution of the respective genome. Based on plots produced by this script, we decided to expand HGT boundaries in PacBio genomes by 1000 bp each side and in stLFR genomes by 25 bp each side. All reads overlapping the HGT boundaries were counted and stored in a bed file with start+stop coordinates of HGTs as separate entries. To extract reads overlapping with the expanded HGT boundary, `bedtools intersect` was used and the option `-F 1` was added to require that the entire stretch defined in the bed file is covered by a given read. We further filtered reads that did not map well to the HGT boundary region (e.g. they could also map somewhere else) and removed multimapping reads.

After filtering, we counted boundary reads again only for the 497 remaining high-quality HGT candidates using the summary script [countReads_HQ_HGTcandidates.sh](01_HGT_validation/countReads_HQ_HGTcandidates.sh). Finally, all final read count files were merged together.

### Prokaryotic gene annotation 
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

#### Predict gene models for selected HGTs
We reconstructed gene models for specific HGTs, such as Lysozymes using `StringTie`.
```bash
./markdown_scripts/09_run.stringtie.md
```

#### Additional analysis for Ankyrin repeat HGAs in ants
This script deals with ankyrin repeat HGAs detected in ants and extracts ANK nucleotide sequences. 
```bash
./markdown_scripts/10_Ankyrin_repeats.md
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

#### Detail analyses of other HGTs
This script was used to investigate HGTs other than ankyrin repeats and clade-specific HGAs.
```bash
./r_scripts/10_Other_LGTs.Rmd 
```

#### Detail analyses of Ankyrin repeat HGAs
This script was used to investigate in detail the ankyrin-repeat HGAs in ants including plotting of ankyrin repeats onto the GAGA phylogeny.
```bash
./r_scripts/11_ANKs_analysis.Rmd 
```

#### Plot a subset of branches of the GAGA phylogeny
This script was used to plot only a subset of relevant branches of the inferred GAGA ant phylogeny.
```bash
./r_scripts/12_plotGAGA.tree.Rmd 
```

#### Overview of species-specific LGTs
This script produces some summarized overview plots related to putatively uniquely-occurring HGAs in ants.

```bash
./r_scripts/13_Unique_LGTs.Rmd 
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