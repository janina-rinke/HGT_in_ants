# Quick commands to calculate the number of reads overlapping LGT boundaries

```bash
# change to where you want to work
#cd /global/homes/jg/schradel/data/GAGA/GAGA-0200-LGT/results/boundaryReads
```

To start, we need 2 files: One `bam` file holding the PacBio reads and one `bed` file containing the LGT start and stop coordinates.

We start with a file containing the LGT candidates you want to screen (I used `LGTs.candidateloci.bed`, see if that makes sense) to create a `bed` file with the start and stop positions as separate entries.

```bash
cat ../../data/GAGA-0200/LGTs.candidateloci.bed|parallel --colsep "\t"  echo -e '{1}"\t"{2}"\t"{2}"\t"{1}"-"{2}":"{3}.start"\n"{1}"\t"{3}"\t"{3}"\t"{1}"-"{2}":"{3}.end' >  GAGA-0200.LGTboundaries.bed
```

Next, we can extract all the PacBio reads overlapping these boundaries
```bash
## merge the two bam files for nAo and other
#samtools merge <outfile.bam> <infile1.bam> <infile2.bam>
samtools merge merged.candidateloci.loose.bam ../../data/GAGA-0200/LGTs.nAo.candidateloci.loose.PacBio.bam ../../data/GAGA-0200/LGTs.candidateloci.loose.PacBio.bam


# extract reads overlapping the boundaries
bedtools intersect -abam merged.candidateloci.loose.bam -b GAGA-0200.LGTboundaries.bed > GAGA-0200.LGTboundaries.PacBio.overlap.bam
# index the bam file
samtools index GAGA-0200.LGTboundaries.PacBio.overlap.bam
```

## Slop

Now, if we want to expand the required overlap, and want e.g. that reads overlap the LGT start position by 1000 bp on both sides, we can do that with `bedtools slop` and `bedtools intersect`

First, we `slop` the boundaries
```bash
# create a genome file containing the length of all scaffolds in the assembly (required by bedtools)
samtools faidx ../../data/GAGA-0200_SLR-superscaffolder_final_dupsrm_filt.fasta
cut -f 1-2 ../../data/GAGA-0200_SLR-superscaffolder_final_dupsrm_filt.fasta.fai > GAGA-0200.genome

# then slop the bed file by 50 bp.
bedtools slop -i GAGA-0200.LGTboundaries.bed -g GAGA-0200.genome -b 50 > GAGA-0200.LGTboundaries.50bp.up.down.bed
```

## intersect
Now we can `intersect` the reads with the expanded LGT boundary. Here, we add the `-F 1` option, to require that the entire stretch defined in the bed file is covered by a given read.

```bash
# extract reads completely overlapping the +-50 bp of the boundary
bedtools intersect -F 1 -abam merged.candidateloci.loose.bam -b GAGA-0200.LGTboundaries.50bp.up.down.bed > GAGA-0200.LGTboundaries.50bp.PacBio.overlap.bam
```

## Filter reads based on mapping quality
Finally, we filter reads that are not mapped very well to the region (e.g. they could also map somewhere else). Also, we remove reads that are in the file twice (because we merged the nAo and the other bam file above)
```bash
# remove not primary alignments and supplementary alignments from bam file. Google "sam flags explained" for details
## I also remove duplicated reads (due to the merging of nAo Bam and other bam above), with the awk command
samtools view -F 1024 -F 256 -F 2048 GAGA-0200.LGTboundaries.50bp.PacBio.overlap.bam -h |awk '!visited[$0]++|| $1 ~ /^@/' |samtools view -bS - > GAGA-0200.LGTboundaries.50bp.good.PacBio.overlap.bam
```

## Stats
```bash
# count how many reads you get per boundary
bedtools coverage -f 1 -b GAGA-0200.LGTboundaries.50bp.good.PacBio.overlap.bam -a GAGA-0200.LGTboundaries.50bp.up.down.bed  -counts > GAGA-0200.LGTboundaries.50bp.good.PacBio.overlap.bed
cat GAGA-0200.LGTboundaries.50bp.good.PacBio.overlap.bed|paste - - > GAGA-0200.LGTboundaries.50bp.good.PacBio.overlap.tsv
# calculate average read length (not sure this is useful :) )
bedtools bamtobed -bed12 -i GAGA-0200.LGTboundaries.50bp.good.PacBio.overlap.bam|cut -f 11|awk '{ total += $1; count++ } END { print total/count }'
```

The file `GAGA-0200.LGTboundaries.50bp.good.PacBio.overlap.tsv` contains now one row for each candidate and gives the following columns

```
col1: Scaffold
col2: start region - 50 bp
col3: start region + 50 bp
col4: "complete candidate region.start"
col5: number of reads overlapping start +- 50 bp
col6: Scaffold
col7: end region - 50 bp
col8: end region + 50 bp
col9: "complete candidate region.end"
col20: number of reads overlapping end +- 50 bp
```
