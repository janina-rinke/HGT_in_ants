# Gene Trees for LGT candidates
### Janina Rinke, 01.03.2022


## CDS gene trees

### 1. Prepare CDS gene tree with all detected Lyzozyme candidates and best BLAST bacteria results

The candidates have been checked for completeness before. Only the longest CDS was considered for the gene tree.


Prepare FASTA file called `Lyzozyme_all_seqs.fa` with all the bacterial sequences you want to align to the Lyzozyme LGTs.

### 1.1 Run MAFFT alignment on `Lyzozymes_all_seqs.fa`:
```bash
mafft Lyzozymes_all_seqs.fa > Lyzozymes_all_seqs.aln
```

### 1.2 Run IQtree
```bash
/global/projects/programs/source/iqtree-2.1.3-Linux/bin/iqtree2 -s Lyzozyme_all_seqs.aln -nt 12 -b 100

# -s ALIGNMENT input file
# -nt: Number of iterationsls
#- b: Bootstrap
```

### 2. Prepare protein gene tree with Lyzozymes
