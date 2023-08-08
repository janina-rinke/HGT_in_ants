#!/bin/bash

cd /global/scratch2/s_tret01/research_project_genomics

### Identify proteins where GAGA-0200 LGT has more genes than all other Blochmannia species

## with a local BLAST
#esearch -db protein -query "WP_015344634.1 OR WP_015344631.1 OR WP_011282911.1 OR WP_015344625.1 OR WP_015344617.1 OR WP_015344616.1 OR WP_015344615.1 OR WP_015344614.1 OR WP_015344604.1 OR WP_015344602.1 OR WP_015344600.1 OR WP_015344598.1 OR WP_015344508.1 OR WP_015344506.1 OR WP_015344504.1 OR WP_011282775.1 OR WP_011282773.1 OR WP_015344500.1 OR WP_011282770.1 OR WP_015344498.1 OR WP_015344497.1 OR WP_015344496.1 OR WP_015344495.1 OR WP_015344494.1 OR WP_015344493.1 OR WP_011282897.1 OR WP_015344499.1 OR WP_015344492.1 OR WP_015344343.1 OR WP_015344635.1 OR WP_171805071.1 OR WP_015344632.1 OR WP_015344629.1 OR WP_015344603.1 OR WP_041567607.1 OR WP_015344599.1 OR WP_015344597.1 OR WP_015344628.1 OR WP_015344596.1 OR WP_015344507.1 OR WP_015344503.1 OR WP_015344627.1 OR WP_015344505.1 OR WP_015344501.1 OR WP_015344502.1 OR WP_015344626.1 OR WP_015344630.1 OR WP_015344610.1 OR WP_015344509.1" | efetch -db protein -format docsum | xtract -pattern DocumentSummary -element Caption Title > TableOfGenes.tsv

#esearch -db protein -query "WP_015344634.1" | efetch -db protein -format docsum | xtract -pattern DocumentSummary -element Caption Title Organism

## with the annotation file



### Identify GAGA-0200 LGT unique proteins

#cat ./2_pipeline/phylogenetic_analysis_v2/store/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups_UnassignedGenes.tsv | grep -o "LOCUS.*" | tr -d '	' > ./2_pipeline/proteinBlast/store/uniqueProteinList.txt
#seqkit grep -nrif ./2_pipeline/proteinBlast/store/uniqueProteinList.txt ./0_data/external/protein.fasta > ./2_pipeline/proteinBlast/temporary/uniqueProteins.fasta
#blastp -db /global/databases/swissprot/2020_11/uniprot_sprot.fasta -query ./2_pipeline/proteinBlast/temporary/uniqueProteins.fasta -outfmt 7 -evalue 0.0000001 -out ./2_pipeline/proteinBlast/store/uniqueProteins.tsv
diamond blastp --fast -d /global/scratch2/databases/diamond/uniprot90_bacteria_20210707/uniprot_bacteria-0.9.dmnd -q ./2_pipeline/proteinBlast/temporary/uniqueProteins.fasta --outfmt 6 qseqid sseqid full_sseq evalue bitscore score pident stitle > ./2_pipeline/proteinBlast/store/uniqueProteins_bacdb.txt

#esearch -db protein -query ./2_pipeline/proteinBlast/temporary/uniqueProteins.fasta | efetch -db protein -format docsum #| xtract -pattern DocumentSummary -element Caption Title Organism

