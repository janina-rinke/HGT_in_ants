# Gene Trees for in-depth analysis of HGT candidates


## Protein gene trees

### 1. Prepare a `fasta` file with all protein sequences of the clade-specific HGT detected in the respective ant genomes and the best NCBI BLAST results.

The candidates have been checked for completeness before.

Lyzozyme query sequence to search for homologous proteins in other species with NCBI blastp (nr database):
```bash
# Protein blast with default settings
>GAGA-0098.Harpagoxenus_sublaevis_BLAST_query_sequence
MMENVIVDLSHWNQSVNFKLAKEDGIIGVLHKATQGLKYVDPTYAERRKAAEEENLMWGAYHFGVGENGTDQAEHFLKTVGNVSKVLLALDIEENRNGKNITPKQAEEFVNRVHEATGRFPLIYGSPYFLKDFATPILTKSLLWMAKWDVRPILPKGWKKWVLWQYTDGRKGPEPHSVRGIGPCDRDKFNGTLEELKDFWISKPEDT
```
For lyzozymes, we chose the ten best bacterial hits (all *Wolbachia*), one hit against a putative HGT transferred into the genome of *Drosophila ananassae* and three bacterial species, other than *Wolbachia*, as outgroups.

Prepare a `fasta` file called `<HGT>_all_seqs_proteins.fa` with all relevant HGT protein sequences and all bacterial sequences which should be aligned to the HGT.

### 1.1 Run MAFFT alignment on `<HGT>_all_seqs_proteins.fa`:
```bash
mafft <HGT>_all_seqs_proteins.fa > <HGT>_all_seqs_proteins.aln
```

### 1.2 Run IQtree
```bash
/global/projects/programs/source/iqtree-2.1.3-Linux/bin/iqtree2 -s <HGT>_all_seqs_proteins.aln -nt 12 -b 100

# -s ALIGNMENT input file
# -nt: Number of iterations
#- b: Bootstrap
```
After running IQtree, you will get a resulting `<HGT>_all_seqs_proteins.aln.treefile`. This file can be used to visualize the three with common tree viewing programs such as iTOL etc.
