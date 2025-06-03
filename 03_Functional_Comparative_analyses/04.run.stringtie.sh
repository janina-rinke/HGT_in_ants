#!/usr/bin/bash

# Run StringTie using RNAseq data to predict gene models for HGT candidates

basedir=./RNAseqmapping/reannotation.dfast_RNAseq_mapping

## Example: GAGA-0528
species=GAGA-0528
HGT_locus=GAGA-0528.Scaffold97.419863-420872.1000.out
echo $basedir

cd $basedir


# Merge all bam files in a new file called <genome>.RNAseq.bam
samtools merge ${basedir}/${species}/${HGT_locus}/${species}.RNAseq.bam ./RNAseq/run_final_keepingbam_pergagaid/GAGA-0099/GAGA-0099.Scaffold20.2658162-2659328.fa/*.bam

# run stringtie on GAGA-0099.RNAseq.bam
/Users/Janina/software/stringtie-2.2.1.OSX_x86_64/stringtie ${basedir}/${species}/${HGT_locus}/mergedRNAseq.bam > ${HGT_locus}.RNAseq.stringtie.gff