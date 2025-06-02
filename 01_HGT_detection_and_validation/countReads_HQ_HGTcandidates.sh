#!/usr/bin/bash

# To run this scipt on a bed file with high quality HGT candidates only:
## cat GAGA.LGTs.allcoordinates.tsv | parallel --colsep "\t" bash countallReads.sh {1} {2} {3} {4}

GAGAid=$1
scf=$2
start=$3
end=$4

# Step 1: create bed file with start, stop and complete coordinates
echo -e "$scf\t$start\t$end" | parallel --colsep "\t"  echo -e '{1}"\t"{2}"\t"{2}"\t"{1}"-"{2}":"{3}.start"\n"{1}"\t"{3}"\t"{3}"\t"{1}"-"{2}":"{3}.end"\n"{1}"\t"{2}"\t"{3}"\t"{1}"-"{2}":"{3}.complete' > $GAGAid.$scf.$start.$end.LGTboundaries.bed

# Step 2: Expand the required overlap by 1000 bp for PacBio genomes
bedtools slop -i $GAGAid.$scf.$start.$end.LGTboundaries.bed -g ../../0_data/$GAGAid/results/genome.file -b 1000 > $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.up.down.bed

# Step 3: Extract reads completely overlapping the +- 1000 bp of the boundary
bedtools intersect -F 1 -abam ../../0_data/$GAGAid/results/merged.candidateloci.loose.bam -b $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.up.down.bed > $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.PacBio.overlap.bam

# Step 4: Filter reads based on mapping quality and remove duplicated reads with awk command
samtools view -F 1024 -F 256 -F 2048 $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.PacBio.overlap.bam -h | awk '!visited[$0]++|| $1 ~ /^@/' | samtools view -bS - > $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.bam

# Step 5: Count the number of reads per boundary
bedtools coverage -f 1 -b $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.bam -a $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.up.down.bed -counts > $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.bed

# Step 6: Write the result in a tsv file
cat $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.bed | paste - - - > $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.tsv

#Step 7: Prepare tsv file to combine read counts for all candidates
paste <(echo -e $GAGAid'\t'$scf'\t'$start'\t'$end) <(awk '{print $5"\t"$10"\t"$15}' $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.tsv) > $GAGAid.$scf.$start.$end.final.read.counts.tsv