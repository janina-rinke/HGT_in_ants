#!/usr/bin/bash


# Calculate the number of reads overlapping HGT boundaries and create a bed file with start+stop coordinates of HGTs as separate entries

# for all HGT candidates
find */results/LGTs.candidateloci.loose.bed | parallel -I% --max-args 1 qsub -v file="%" \
makeHGTboundaryfile.sh -o ./tmp/$file.out -e ./tmp/$file.err


# for all no ant HGT candidates
find */results/LGTs.nAo.candidateloci.loose.bed | parallel -I% --max-args 1 qsub -v file="%" \
makeHGTboundaryfile.sh -o ./tmp/$file.out -e ./tmp/$file.err



# Check if all genomes have this file:
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed ]] && echo "$dir"; done

find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/LGTs.nAo.candidateloci.loose.bed.LGTboundaries.bed ]] && echo "$dir"; done


# Merge the bam files for nA and other
for i in * ; do samtools merge $i/results/merged.candidateloci.loose.bam \
$i/results/LGTs.nAo.candidateloci.loose.PacBio.bam $i/results/LGTs.candidateloci.loose.PacBio.bam; done

# Check if all genomes have the file
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/merged.candidateloci.loose.bam ]] && echo "$dir"; done


# Calculate read overlap for HGT boundaries with bedtools intersect
for i in * ; do bedtools intersect -abam $i/results/merged.candidateloci.loose.bam \
-b $i/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed \
> $i/results/$i.LGTboundaries.PacBio.overlap.bam; done

# index the bam file
for i in * ; do samtools index $i/results/$i.LGTboundaries.PacBio.overlap.bam; done

# Check if all genomes have the file
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.PacBio.overlap.bam ]] && echo "$dir"; done

# For all nA genomes
for i in * ; do bedtools intersect -abam $i/results/merged.candidateloci.loose.bam \
-b $i/results/LGTs.nAo.candidateloci.loose.bed.LGTboundaries.bed \
> $i/results/$i.nAo.LGTboundaries.PacBio.overlap.bam; done

# index the bam file
for i in * ; do samtools index $i/results/$i.nAo.LGTboundaries.PacBio.overlap.bam; done

find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.nAo.LGTboundaries.PacBio.overlap.bam ]] && echo "$dir"; done

# Expand the required overlap by 1000 bp for PacBio genome assemblies
for i in *; do bedtools slop -i $i/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed \
-g $i/results/genome.file -b 1000 > $i/results/$i.LGTboundaries.1000bp.up.down.bed; done

for i in *; do bedtools slop -i $i/results/LGTs.nAo.candidateloci.loose.bed.LGTboundaries.bed \
-g $i/results/genome.file -b 1000 > $i/results/$i.nAo.LGTboundaries.1000bp.up.down.bed; done

# Extract reads completely overlapping the +-1000 bp of the boundary
for i in *; do bedtools intersect -F 1 -abam $i/results/merged.candidateloci.loose.bam \
-b $i/results/$i.LGTboundaries.1000bp.up.down.bed > $i/results/$i.LGTboundaries.1000bp.PacBio.overlap.bam; done

for i in *; do bedtools intersect -F 1 -abam $i/results/merged.candidateloci.loose.bam \
-b $i/results/$i.nAo.LGTboundaries.1000bp.up.down.bed > $i/results/$i.nAo.LGTboundaries.1000bp.PacBio.overlap.bam; done



# Expand the required overlap by 25 bp for stLFR genome assemblies
for i in *; do bedtools slop -i $i/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed \
-g $i/results/genome.file -b 25 > $i/results/$i.LGTboundaries.25bp.up.down.bed; done

for i in *; do bedtools slop -i $i/results/LGTs.nAo.candidateloci.loose.bed.LGTboundaries.bed \
-g $i/results/genome.file -b 25 > $i/results/$i.nAo.LGTboundaries.25bp.up.down.bed; done

# Extract reads completely overlapping the +-25 bp of the boundary
for i in *; do bedtools intersect -F 1 -abam $i/results/merged.candidateloci.loose.bam \
-b $i/results/$i.LGTboundaries.25bp.up.down.bed > $i/results/$i.LGTboundaries.25bp.PacBio.overlap.bam; done

for i in *; do bedtools intersect -F 1 -abam $i/results/merged.candidateloci.loose.bam \
-b $i/results/$i.nAo.LGTboundaries.25bp.up.down.bed > $i/results/$i.nAo.LGTboundaries.25bp.PacBio.overlap.bam; done



# Filter reads based on mapping quality

# For PacBio genomes
for i in *; do samtools view -F 1024 -F 256 -F 2048 \
$i/results/$i.LGTboundaries.1000bp.PacBio.overlap.bam -h |awk '!visited[$0]++|| $1 ~ /^@/' | samtools view -bS - > $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bam; done

for i in *; do samtools view -F 1024 -F 256 -F 2048 \
$i/results/$i.nAo.LGTboundaries.1000bp.PacBio.overlap.bam -h |awk '!visited[$0]++|| $1 ~ /^@/' | samtools view -bS - > $i/results/$i.nAo.LGTboundaries.1000bp.good.PacBio.overlap.bam; done

# For stLFR genomes
for i in *; do samtools view -F 1024 -F 256 -F 2048 \
$i/results/$i.LGTboundaries.25bp.PacBio.overlap.bam -h |awk '!visited[$0]++|| $1 ~ /^@/' | samtools view -bS - > $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bam; done

for i in *; do samtools view -F 1024 -F 256 -F 2048 \
$i/results/$i.nAo.LGTboundaries.25bp.PacBio.overlap.bam -h |awk '!visited[$0]++|| $1 ~ /^@/' | samtools view -bS - > $i/results/$i.nAo.LGTboundaries.25bp.good.stLFR.overlap.bam; done



# Count the number of reads per boundary

# For PacBio genomes
for i in *; do bedtools coverage -f 1 \
-b $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bam \
-a $i/results/$i.LGTboundaries.1000bp.up.down.bed  \
-counts > $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bed; done

for i in *; do cat $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bed | paste - - > \
$i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.tsv; done

for i in *; do bedtools coverage -f 1 \
-b $i/results/$i.nAo.LGTboundaries.1000bp.good.PacBio.overlap.bam \
-a $i/results/$i.nAo.LGTboundaries.1000bp.up.down.bed  \
-counts > $i/results/$i.nAo.LGTboundaries.1000bp.good.PacBio.overlap.bed; done

for i in *; do cat $i/results/$i.nAo.LGTboundaries.1000bp.good.PacBio.overlap.bed | paste - - > \
$i/results/$i.nAo.LGTboundaries.1000bp.good.PacBio.overlap.tsv; done

# For stLFR genomes
for i in *; do bedtools coverage -f 1 \
-b $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bam \
-a $i/results/$i.LGTboundaries.25bp.up.down.bed \
-counts > $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bed; done

for i in *; do cat $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bed | paste - - > \
$i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.tsv; done

for i in *; do bedtools coverage -f 1 \
-b $i/results/$i.nAo.LGTboundaries.25bp.good.stLFR.overlap.bam \
-a $i/results/$i.nAo.LGTboundaries.25bp.up.down.bed \
-counts > $i/results/$i.nAo.LGTboundaries.25bp.good.stLFR.overlap.bed; done

for i in *; do cat $i/results/$i.nAo.LGTboundaries.25bp.good.stLFR.overlap.bed | paste - - > \
$i/results/$i.nAo.LGTboundaries.25bp.good.stLFR.overlap.tsv; done