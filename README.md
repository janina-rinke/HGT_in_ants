# Horizontal Gene Transfer in Ants

This repository holds all the code to evaluate horizontally transferred genes (HGTs) from prokaryotic donors in 163 ant genomes, previously identified by an automated Horizontal Gene Transfer (HGT) finder pipeline. All the code for the previous, automated detection of potential HGTs/LGTs in ant genomes can be found [HERE](https://github.com/dinhe878/GAGA-Metagenome-LGT).


## 1 Detection, validation and quality assessment of HGT candidates
The detection of potential HGTs with protein-coding hits in the 163 ant genomes of interest obtained by our bioinformatic HGT finder pipeline is explained in detail [HERE](https://github.com/dinhe878/GAGA-Metagenome-LGT).

This HGT finder pipeline yielded 13,664 initial HGT candidates for downstream analyses. We produced overview plots for each candidate to aid the examination of HGT quality and processing of HGT candidates in downstream analyses:

[*01.analyseLGTs.Rmd*](01_HGT_detection_and_validation/01.analyseLGTs.Rmd)


After that, HGT filtering steps were employed based on the previously produced and examined HGT candidate overview plots. Filters were further determined based on manual assessment of HGT candidates from seven randomly chosen GAGA genomes, as well as plotted parameter distributions of all predicted HGTs:

[*02.LGTfiltering.Rmd*](01_HGT_detection_and_validation/02.LGTfiltering.Rmd)

The following filter thresholds were applied based on a step-by-step filtering process:

0) Unfiltered: 13664 candidates, 163 genomes
1) Filtered by e-value (>1e-25, removed 7616 candidates)
2) Filtered by ct4 (ct4 > 0.25, removed 2.646 candidates)
3) Filtered by ce (ce > 1.5, removed 485 candidates)
4) Filtered by BitDiffSum (BitDiffSum > 150, removed 4610 candidates)
5) Filtered by candidate length (length > 100, removed 2209 candidates)
6) Filtering out candidates within first 1000 bp of scaffold (removed 182 candidates)
7) Filtering out candidates within last 1000 bp of scaffold (removed 92 candidates)
8) Filtering out candidates with 0/1 read at boundaries (reads_start & reads_end > 1, removed 251 candidates)
9) Merging candidates in 20 kb range: 1148 candidate regions.


### Calculation of reads overlapping HGT boundaries

We used `bedtools` and `samtools` to map raw PacBio sequencing reads (from `bam` files) to the respective genomes and to calculate reads overlapping the start, end, and total sequence length of HGT candidates. All coordinates were stored in bed files (`LGTs.candidateloci.loose.bed`). HGT boundaries were expanded based on the average read length distribution of the respective genomes, calculated and plotted with:

[*03.Distribution_plots_avg_read_lengths.Rmd*](01_HGT_detection_and_validation/03.Distribution_plots_avg_read_lengths.Rmd)


Based on plots produced by this script, we decided to expand HGT boundaries in PacBio genomes by 1000 bp each side and in stLFR genomes by 25 bp each side. All reads overlapping the HGT boundaries were counted and stored in a bed file with start+stop coordinates of HGTs as separate entries. To extract reads overlapping with the expanded HGT boundary, `bedtools intersect` was used and the option `-F 1` was added to require that the entire stretch defined in the bed file is covered by a given read. We further filtered reads that did not map well to the HGT boundary region (e.g. they could also map somewhere else) and removed multimapping reads:

[*04.countBoundaryReads.sh*](01_HGT_detection_and_validation/04.countBoundaryReads.sh)

After filtering, we counted boundary reads again only for the 497 remaining high-quality HGT candidates using the summary script [countReads_HQ_HGTcandidates.sh](01_HGT_detection_and_validation/countReads_HQ_HGTcandidates.sh). Finally, all final read count files were merged together.

## 2 Prokaryotic HGT annotation 
`prodigal`, `kraken2`, and `dfast` were run searching against `nr` (https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA), `nt` https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA), and `UniRef90` (https://www.uniprot.org/help/uniref) databases downloaded from NCBI as well as against `COG` and `tigrfam` databases, to annotate the remaining high-quality HGT candidate sequences in the ant genomes. For all downstream analyses following prokaryotic gene annotation, only `dfast` results were used.

[*01.Prokaryotic_Gene_Annotation_dfast.sh*](02_Prokaryotic_gene_annotation/01.Prokaryotic_Gene_Annotation_dfast.sh)

### Check HGT completeness
All annotated HGT candidates were checked for their completeness using the obtained DFAST results. We extracted for every HGT candidate the best uniprot hit from DFAST. Further, we extracted the start and stop codon of each HGT candidate to check for complete ORFs using `seqkit` and counted the number of predicted CDS per HGT:

[*02.HGT_completeness_check.sh*](02_Prokaryotic_gene_annotation/02.HGT_completeness_check.sh)

To obtain all information from the DFAST gffs for our HGT candidates, query coverage and subject coverage were added as additional information for every HGT candidate to estimate completeness and functionality of the HGT:

[*03.screenDfastGFFs.Rmd*](02_Prokaryotic_gene_annotation/03.screenDfastGFFs.Rmd)

### Re-annotation of incomplete HGT candidates
All candidates were extended by 1000 bp in 5' and 3' direction and re-annotated with DFAST to identify HGTs which were previously annotated incompletely. We further modified `gff` files obtained from DFAST to obtain coordinates in reference to the genomes using the script `Modify_gffs.pl`.

[*04.reannotation.DFAST.sh*](02_Prokaryotic_gene_annotation/04.reannotation.DFAST.sh)

### Check for additional HGT candidates within all GAGA genomes
Finally, we used the annotated HGT genes to identify putative additional HGT candidates that were previously filtered out due to our strict filtering criteria. All resulting candidates were intersected with initially predicted HGT candidates by the automated HGT finder pipeline.

[*05.Blast_GAGA_genomes.sh*](02_Prokaryotic_gene_annotation/05.Blast_GAGA_genomes.sh)


### Integrate information from Uniprot
We used the package *uniprotR* to obtain additional information for all HGT candidates and integrated this information into our summarized tables for functional analyses and as a basis for in-depth analyses of clade-specific HGTs.

[*06.retrieveUniproInfo.Rmd*](02_Prokaryotic_gene_annotation/06.retrieveUniproInfo.Rmd)



## 3 Functional and Comparative analyses
We mapped paired and unpaired RNAseq reads to the genomes using STAR. RNAseq data was made available and collected by the GAGA project. Information about available RNAseq data for each investigated species, including the sampling of different ant castes and developmental stages, can be obtained from Tab. S1A in Vizueta et al (2025). RNAseq data was available for 130 of the 163 studied ant species and subsequently mapped to the genomes to assess gene expression of the HGT loci: 

[*01.run_rna_mapping_wholegenome.pl*](03_Functional_Comparative_analyses/01.run_rna_mapping_wholegenome.pl)

### Gene expression analysis of HGTs
We obtained a summary CDS-level table, as well as a locus-level table with summarized information as described above. Subsequently, we investigated HGT expression for all high-quality HGT candidates. Here, we calculated the sum of unique read counts per CDS from bam files for each HGT region:

[*02.calculateExpression.Rmd*](03_Functional_Comparative_analyses/02.calculateExpression.Rmd)

To plot the coverage of gene expression for all HGTs, we used unique read counts as input and added up all read counts for all life stages within each respective HGT candidate. RNAseq coverage was visualized for each HGT candidate in R:

[*03.plotRNAseqCoverage.Rmd*](03_Functional_Comparative_analyses/03.plotRNAseqCoverage.Rmd)

We further reconstructed gene models for specific HGTs, such as Lysozymes using `StringTie`:

[*04.run.stringtie.sh*](03_Functional_Comparative_analyses/04.run.stringtie.sh)


### Synteny of clade-specific HGTs
To infer synteny of genes, all clade-specific HGT regions were extended by 40 kb up- and downstream and all ant genes and protein sequences within this flanking region were extracted. `Minimap2` was then used to conduct an all-vs-all alignment after which OrthoFinder was used to determine orthogroups across species. The extent of synteny was subsequently plotted with the R package `gggenomes`.

[*05.Synteny.Rmd*](03_Functional_Comparative_analyses/05.Synteny.Rmd)


### Phylogenetic analyses
We ran gene trees on the clade-specific prokaryotic protein HGT sequences (Lysozymes, MurNAc etherases, CFA synthases) and their closest blast homologs to conduct phylogenetic analyses. We retrieved the five most similar homologs for each candidate via BLASTp against the NCBI non-redundant (nr) protein database. Multiple sequence alignments were performed using MAFFT with default settings. Phylogenetic trees were then constructed using maximum likelihood inference in IQTREE2, with node support assessed through 100 bootstrap replicates. Substitution model selection was performed automatically using ModelFinder, which is integrated within IQ-TREE and selects the best-fitting model based on statistical criteria such as the Bayesian Information Criterion (BIC):

[*06.Run_GeneTrees.sh*](03_Functional_Comparative_analyses/06.Run_GeneTrees.sh)

## 4 XGPRT HGT functional analyses
We analysed in detail one HGT case from an enterobacterial donor into the genome of the ant *Cardiocondyla obscurior*. Using previously published RNAseq datasets, we tested for HGT co-expression between the *Cardiocondyla*-associated intracellular endosymbiont *Candidatus Westeberhardia cardiocondylae* and the ant host genome.

### XGPRT Gene expression analysis
We used the following reference genomes and gene annotations for the analysis:

- *Cardiocondyla obscurior* reference genome [GCF_019399895.1.fna](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_019399895.1/)
- *Cardiocondyla obscurior* gene annotation [GCF_019399895.1.gff](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_019399895.1/)
- *Westeberhardia* reference genome [GCF_001242845.1.fna](https://www.ebi.ac.uk/ena/browser/view/LN774881)
- *Westeberhardia* gene annotation [GCF_019399895.1.gff](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_001242845.1/)

The XGPRT HGT locus was inserted into the gene annotation manually with the CDS sequence being defined based on annotation [KAL0100757.1](https://www.ncbi.nlm.nih.gov/protein/KAL0100757.1/). Both host and symbiont reference genomes and gene annotations were then concatenated for subsequent mapping.


RNAseq Short-read files were obtained from NCBI´s main raw data repository:
- [Queen larvae](https://www.ncbi.nlm.nih.gov/sra/SRX692538)
- [Winged male larvae](https://www.ncbi.nlm.nih.gov/sra/SRX879678)
- [Ergatoid male larvae](https://www.ncbi.nlm.nih.gov/sra/SRX879676)
- [Worker larvae](https://www.ncbi.nlm.nih.gov/sra/SRX879674)


We trimmed raw RNAseq reads using Trim Galore, mapped them against the concatenated reference genomes using STAR, and counted mapped reads using feature counts:

[*01.RNAseq_read_processing.sh*](04_XGPRT_HGT_analyses/01.RNAseq_read_processing.sh)

Read counts were converted to log2 counts per million (library size normalization). We performed a pearson correlation of the XGPRT HGT against each of the annotated host and endosymbiont genes, respectively, followed by Bonferroni adjustments of the p-values. After that, we specifically tested for a gene expression correlation between the XGPRT HGT and the *Westeberhardia* genes.

[*02.Gene_expression_normalization_and_correlation.Rmd*](04_XGPRT_HGT_analyses/02.Gene_expression_normalization_and_correlation.Rmd)

[*03.Westeberhardia_correlations.Rmd*](04_XGPRT_HGT_analyses/03.Westeberhardia_correlations.Rmd)

To prepare for a weighted gene correlation network analysis (WGCNA) we created an adjacency network, using a softPower = 13 and type = unassigned. This is used to compute a TOM (Topological Overlap Matrix). A subset of this TOM which includes XGPRT and all 82 genes that significantly correlate with the XGPRT HGT in their expression patterns was then used to generate input data creating a weighted gene correlation network in Cytoscape.

[*04.WGCNA.Rmd*](04_XGPRT_HGT_analyses/03.WGCNA.Rmd)