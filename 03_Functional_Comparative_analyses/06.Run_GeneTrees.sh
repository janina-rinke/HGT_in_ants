#!/usr/bin/bash

# Gene Trees for in-depth analysis of HGT candidates

## Alignment with MAFFT
mafft <HGT>_all_seqs_proteins.fa > <HGT>_all_seqs_proteins.aln

## Run IQtree2
/global/projects/programs/source/iqtree-2.1.3-Linux/bin/iqtree2 -s <HGT>_all_seqs_proteins.aln -nt 12 -b 100