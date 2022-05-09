# Gene Trees for in-depth analysis of LGT candidates


## Protein gene trees

### 1. Prepare fasta file with all detected Lyzozyme candidates and best BLAST bacteria results

The candidates have been checked for completeness before.


Prepare FASTA file called `HGTcandidate_all_seqs.fa` with all relevant LGT protein sequences and all the bacterial sequences you want to align to the LGT(s).

### 1.1 Run MAFFT alignment on `HGTcandidate_all_seqs.fa`:
```bash
mafft HGTcandidate_all_seqs.fa > HGTcandidate_all_seqs.aln
```

### 1.2 Run IQtree
```bash
/global/projects/programs/source/iqtree-2.1.3-Linux/bin/iqtree2 -s HGTcandidate_all_seqs.aln -nt 12 -b 100

# -s ALIGNMENT input file
# -nt: Number of iterationsls
#- b: Bootstrap
```
After running IQtree, you will get a resulting `HGTcandidate_all_seqs.aln.treefile`. This file can be used to visualize the three with common tree viewing programs such as iTOL etc.
