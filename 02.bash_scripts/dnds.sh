#!/bin/bash

cd /global/scratch2/s_tret01/research_project_genomics/

#while read p; do
#	echo "$p"
#	B_C_nipp=$(sed -n '1{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	B_chro=$(sed -n '2{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	B_flor=$(sed -n '3{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	B_nipp=$(sed -n '4{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	B_obli=$(sed -n '5{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	B_penn=$(sed -n '6{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	B_turn=$(sed -n '7{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	B_vafe=$(sed -n '8{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	GAGA_0200=$(sed -n '9{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#	W_card=$(sed -n '10{p;q}' ./2_pipeline/phylogenetic_analysis_v2/temporary/$p.headerlist.fa)
#echo "($W_card,(($B_C_nipp,$B_obli),((((($B_chro,$B_penn),$GAGA_0200),$B_nipp),($B_flor,$B_vafe)),$B_turn)));" > ./2_pipeline/dnds/store/Species_Trees/$p.sp.nwk
#done <./2_pipeline/phylogenetic_analysis_v2/store/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups_SingleCopyOrthologues.txt

#date
#echo "currently running orthogroup 038"
#HYPHYMPI aBSREL "./2_pipeline/phylogenetic_analysis_v2/store/cds_alignments/OG0000038.cds.aln.fa" "./2_pipeline/dnds/store/Species_Trees/OG0000038.sp.nwk" > ./2_pipeline/dnds/store/aBSREL_results/OG0000038.aBSREL.results


#while read p; do
#	date
#	echo "currently running orthogroup $p"
#	HYPHYMPI aBSREL "./2_pipeline/phylogenetic_analysis_v2/store/cds_alignments03/$p.cds.aln.fa" "./2_pipeline/dnds/store/Species_Trees/$p.sp.nwk" > ./2_pipeline/dnds/store/aBSREL_results/$p.aBSREL.results
#done <./2_pipeline/phylogenetic_analysis_v2/store/OrthoFinder/Results_Jul13/Orthogroups/Orthogroups_SingleCopyOrthologues.txt


