Calculate reads for remaining filtered candidates.

In this script, we want to get reads only for the filtered, remaining 497 high-quality HGA candidates.
We would like to get reads, which
A) overlap the start of the HGT by +- 1000 bp on each side (or 25 bp for stLFR assembled genomes)
B) overlap the end of the HGT by +- 1000bp on each side (or 25bp for stLFR assembled genomes)
C) completely overlap the whole HGT from start to end +- 1000 bp.


# For one candidate:
Start with one candidate to try it out (file: `GAGA_example.tsv`):

```bash
GAGA-0513	Scaffold118	10998	11334
```

What we want to get is a bed file with the following data:

`GAGA-0513.Scaffold118.10998-11334.bed`:

```bash
Scaffold118   9998  11998  start
Scaffold118  10334  12334  end
Scaffold118   9998  12334  complete
```

Then we need to find the corresponding bam file for intersecting.

###1. Write a bash script to create the files:

`nano countReads.sh`

```bash
#!/usr/bin/bash

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
cat $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.bed | paste - - > $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.tsv

#Step 7: Prepare tsv file to combine read counts for all candidates
paste <(echo -e $GAGAid'\t'$scf'\t'$start'\t'$end) <(awk '{print $5"\t"$10"\t"$15}' $GAGAid.$scf.$start.$end.LGTboundaries.1000bp.good.PacBio.overlap.tsv) > $GAGAid.LGTboundaries.finalcounts.tsv
```

To run the script for one candidate:
```bash
cat GAGA_example.tsv | parallel --colsep "\t" bash countReads.sh {1} {2} {3} {4}
```

#### Output file step 1: Create bed file
The resulting bed file `GAGA-0513.Scaffold118.10998.11334.LGTboundaries.bed` now contains the start and stop coordinates as separate entries:
```bash
Scaffold118     10998   10998   Scaffold118-10998:11334.start
Scaffold118     11334   11334   Scaffold118-10998:11334.end
Scaffold118     10998   11334   Scaffold118-10998:11334.complete
```

#### Output file step 2: Expand the required overlap of the candidates by 1000 bp
(File: `GAGA-0513.LGTboundaries.1000bp.up.down.bed`):
```bash
Scaffold118     9998    11998   Scaffold118-10998:11334.start
Scaffold118     10334   12334   Scaffold118-10998:11334.end
Scaffold118     9998    12334   Scaffold118-10998:11334.complete
```

#### Output file step 5: Count the number of reads per boundary:
(File: `GAGA-0513.LGTboundaries.1000bp.up.down.bed`):
```bash
Scaffold118     9998    11998   Scaffold118-10998:11334.start   54
Scaffold118     10334   12334   Scaffold118-10998:11334.end     55
Scaffold118     9998    12334   Scaffold118-10998:11334.complete 53
```

#### Output file step 7: All results in a final tsv file:
(File: `GAGA-0513.LGTboundaries.finalcounts.tsv`)
```bash
GAGA-0513       Scaffold118     10998   11334   54      55      53
```

# For all candidates:

The remaining 497 candidates were all assembled by PacBio sequencing, thus the candidates did not need to be separated into stLFR assembled genomes and PacBio assembled genomes. All candidates were expanded by 1000 bp, as previously calculated as a suitable expansion of the LGT boundary.

nano `countallReads.sh`
```bash
#!/usr/bin/bash

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
```

To run the script for all candidates:
```bash
cat GAGA.LGTs.allcoordinates.tsv | parallel --colsep "\t" bash countallReads.sh {1} {2} {3} {4}
```

### Merge all final read count files together
```bash
cat *.final.read.counts.tsv > Candidates.final.read.counts.tsv
```

### Add headers to the final read count file
```bash
echo -e "GAGAid\tcand.scaffold\tcand.start\tcand.end\tstart_overlap\tend_overlap\tcomplete_overlap" | cat - Candidates.final.read.counts.tsv > Candidates.final.read.counts.headers.tsv
```

At the end, the output file `Candidates.final.read.counts.headers.tsv` looks like this:
```bash
GAGAid  cand.scaffold   cand.start      cand.end        start_overlap   end_overlap     complete_overlap
GAGA-0020       Scaffold107     144838  147367  73      75      67
GAGA-0020       Scaffold1       5322877 5323044 64      64      64
GAGA-0020       Scaffold17      156258  159230  69      59      56
GAGA-0020       Scaffold267     55721   56766   42      44      41
GAGA-0020       Scaffold31      542876  547130  50      30      26
GAGA-0020       Scaffold31      613867  616833  60      67      49
GAGA-0020       Scaffold31      648280  651311  37      50      29
GAGA-0020       Scaffold31      704726  708354  49      55      46
GAGA-0020       Scaffold38      658653  660450  43      34      32
GAGA-0020       Scaffold42      25702   28653   75      71      59
GAGA-0020       Scaffold43      15725   21639   92      92      62
```
