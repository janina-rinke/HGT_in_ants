# BLAST sequence against GAGA genomes

## Prepare the GAGA genomes as a database
for f in * ; do mv -n "$f" "${f%%\_*}"; done

## Add GAGAid to FASTA headers
for f in *; do sed -i '' -e "s/^>/>${f%\_}_/g" "${f}"; done

## Concatenate GAGA genomes into one large genome database
cat * >> GAGA_genomes_database

## Local BLAST against GAGA genomes database
blastn -query Sequence.txt \
-db ../assemblies_copy/GAGA_genomes_database \
-outfmt "6 std qlen" -out blast_Sequence_GAGAgenomes.out



# Intersect newly detected sequences with previously obtained HGT loci by automatic pipeline (Example: GAGA-0328)
bedtools intersect -wa -a LGTs.candidateloci.GAGAid.bed \
-b blast_Lyzozyme_GAGA-0328.bed > results.blast.GAGA-0328.bed

# Extract Fasta sequence for resulting overlapping coordinates
bedtools getfasta -fi ./0_data/assemblies_copy/GAGA-0328 -bed results.blast.GAGA-0328.bed -fo results.blast.GAGA-0328.fa



# Change all bed files fitting to blast bed format
for i in */*; do cp $i/LGTs.candidateloci.loose.bed $i/LGTs.candidateloci2.loose.bed; done
for i in */*; do cp $i/LGTs.nA.candidateloci.loose.bed $i/LGTs.nA.candidateloci2.loose.bed; done

# Insert GAGAid in all files at beginning of each line
for i in GAGA-*; do echo $i ; cat $i/results/LGTs.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.candidateloci.GAGAid.bed; done

for i in OUT-*; do echo $i ; cat $i/results/LGTs.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.candidateloci.GAGAid.bed; done
for i in OUT-*; do echo $i ; cat $i/results/LGTs.nA.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.nA.candidateloci.GAGAid.bed; done

for i in NCBI-*; do echo $i ; cat $i/results/LGTs.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.candidateloci.GAGAid.bed; done
for i in NCBI-*; do echo $i ; cat $i/results/LGTs.nA.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.nA.candidateloci.GAGAid.bed; done



# Run bedtools intersect on all files
for i in */results; do bedtools intersect -wa -a $i/LGTs.candidateloci.GAGAid.bed -b ./0_data/local_blast/blast_Lyzozyme_GAGA-0328.bed >> ./0_data/local_blast/results_Lyzozymes.bed; done
for i in */results; do bedtools intersect -wa -a $i/LGTs.nA.candidateloci.GAGAid.bed -b ./0_data/local_blast/blast_Lyzozyme_GAGA-0328.bed >> ./0_data/local_blast/results_Lyzozymes_nA.bed; done

## For the Etherases:
for i in */results; do bedtools intersect -wa -a $i/LGTs.candidateloci.GAGAid.bed -b ./0_data/local_blast/blast_Etherase_GAGAgenomes.bed >> ./0_data/local_blast/results_Etherase.bed; done
for i in */results; do bedtools intersect -wa -a $i/LGTs.nA.candidateloci.GAGAid.bed -b ./0_data/local_blast/blast_Etherase_GAGAgenomes.bed >> ./0_data/local_blast/results_Etherase_nA.bed; done

## Extract FASTA sequence for resulting coordinates (Example: Lysozymes)
for i in *; do bedtools getfasta -fi ./0_data/assemblies_copy/$i -bed ./0_data/local_blast/results_Lyzozymes.bed -fo ./0_data/local_blast/dfast_candidates/$i.Lyzozyme.fa; done
for i in *; do bedtools getfasta -fi ./0_data/assemblies_copy/$i -bed ./0_data/local_blast/results_Lyzozymes_nA.bed -fo ./0_data/local_blast/dfast_candidates/$i.Lyzozyme.nA.fa; done

## Extract FASTA sequence for resulting coordinates (Example: Etherases)
for i in *; do bedtools getfasta -fi ./0_data/assemblies_copy/$i -bed ./0_data/local_blast/results_Etherase.bed -fo ./0_data/local_blast/dfast_candidates/$i.Etherase.fa; done
for i in *; do bedtools getfasta -fi ./0_data/assemblies_copy/$i -bed ./0_data/local_blast/results_Etherase_nA.bed -fo ./0_data/local_blast/dfast_candidates/$i.Etherase.nA.fa; done


awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' GAGA*
awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' NCBI*
awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' OUT*

# Remove _ and replace with .
for i in * ; do mv ${i} `echo ${i} | sed 's/_/\./'`; done

# Remove : and replace with .
for i in * ; do mv ${i} `echo ${i} | sed 's/\:/\./'`; done

# Remove duplicates
for i in *; do cat $i | seqkit rmdup -s -o $i; done



# Build a database for all bacteria of interest to run dfast for the additional group-specific candidates:
## Extract relevant bacteria
seqkit grep -n -r -p "Wolbachia|Sodalis|Serratia|Yersinia|Rahnella|Escherichia|Spiroplasma|Mycoplasma" ./databases/dfast/uniprot_bacteria-0.9.fasta > ~/relevant.bacteria.uni90.fa

## Build a dfast database
source /usr/share/modules/init/bash
module use /global/projects/programs/modules
module load annotation/dfast/1.2.12

/global/projects/programs/source/dfast/dfast_core-1.2.12/scripts/reference_util.py fasta2dfast -i relevant.bacteria.uni90.fa -o relevant.bacteria.uni90.ref
/global/projects/programs/source/dfast/dfast_core-1.2.12/scripts/reference_util.py formatdb -i relevant.bacteria.uni90.ref


## Submit the jobs:
find . -name "GAGA*" | parallel -I% --max-args 1 qsub -v file="%" run_dfast_BLAST_candidates.sh
find . -name "NCBI*" | parallel -I% --max-args 1 qsub -v file="%" run_dfast_BLAST_candidates.sh
find . -name "OUT*" | parallel -I% --max-args 1 qsub -v file="%" run_dfast_BLAST_candidates.sh