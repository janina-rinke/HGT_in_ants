

# 4. RNA seq analysis
## 4.1 Preparing Reference Genomes and Annotations

Both host and symbiont reference genomes and gene annotations were concatenated for mapping.
### 4.1.1 Reference genomes & gene annotations

We used the following reference genomes and gene annotations for the analysis:

- Cardiocondyla obscurior BR reference genome [GCF_019399895.1.fna](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_019399895.1/)
- Cardiocondyla obscurior BR gene annotation [GCF_019399895.1.gff](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_019399895.1/)
- Westeberhardia reference genome [GCF_001242845.1.fna](https://www.ebi.ac.uk/ena/browser/view/LN774881)
- Westeberhardia gene annotation [GCF_019399895.1.gff](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_001242845.1/)

### 4.1.2 Inserting XGPRT_LGT locus
The XGPRT locus was inserted manually. The CDS sequence was defined based on annotation [KAL0100757.1](https://www.ncbi.nlm.nih.gov/protein/KAL0100757.1/)

[01.XGPRT_insertion.txt](1.code/01.XGPRT_insertion.txt)

### 4.1.3 Concatenating genomes and annotations

Both host and symbiont reference genomes and gene annotations were concatenated for mapping:

[02.Prepareing_concatenated_genomes_and_annotations.sh](1.Prepareing_concatenated_genomes_and_annotations01.Prepareing_concatenated_genomes_and_annotations.sh)

## 4.2. Collecting RNA-seq Data
Short-read files were linked from NCBI´s main raw data repository:
- [queen larvae](https://www.ncbi.nlm.nih.gov/sra/SRX692538)
- [winged male larvae](https://www.ncbi.nlm.nih.gov/sra/SRX879678)
- [ergatoid male larvae](https://www.ncbi.nlm.nih.gov/sra/SRX879676)
- [worker larvae](https://www.ncbi.nlm.nih.gov/sra/SRX879674)

## 4.3. Sequence trimming, mapping and read counting
We trimmed raw RNAseq reads using Trim Galore, mapped them againste the concatenated reference genomes using STAR and counted gene reads using feature counts

[03.RNAseq_read_processing.sh](1.code/03.RNAseq_read_processing.sh)

## 4.4 Count normalization and gene expression normalization

Read counts were converted to log2 counts per million (library size normalization). We performed a pearson correlation of XGPRT against each host and endosymbiont genes, followed by Bonferroni adjustments of the p-values

[04.gene_expression_normalization_and_correlation.Rmd](1.code/04.gene_expression_normalization_and_correlation.Rmd)

## 4.5 WGCNA

To prepare for a weighted gene correlation network analysis we create an adjacency network, using a softPower = 13 and tpye = unsigned. This is used to compute a TOM (Topological Overlap Matrix). A subset of this TOM which includes XGPRT and all 82 genes that significantly correlate with XGPRT in their expression patterns is used to generate input data to create a weighted gene correlation network in Cytoscape.

[05.WGCNA.Rmd](1.code/05.WGCNA.Rmd)