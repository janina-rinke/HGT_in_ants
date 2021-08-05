# Calculate the number of reads overlapping LGT boundaries

```bash
# change to where you want to work
#cd /global/scratch2/j_rink02/master/lgt/0_data
```

To start, we need 2 files: One `bam` file holding the PacBio reads and one `bed` file containing the LGT start and stop coordinates.

We start with a file containing the LGT candidates you want to screen to create a `bed` file with the start and stop positions as separate entries.

### 1. Use the file `LGTs.candidateloci.bed` to create a bed file with start+stop coordinates of the LGT as separate entries.

```bash
cat GAGA-0515/results/LGTs.candidateloci.bed | parallel --colsep "\t"  echo -e '{1}"\t"{2}"\t"{2}"\t"{1}"-"{2}":"{3}.start"\n"{1}"\t"{3}"\t"{3}"\t"{1}"-"{2}":"{3}.end' > GAGA-0515/results/GAGA-0515.LGTboundaries.bed
```

Output of the above file:
```bash
Scaffold10      88994   88994   Scaffold10-88994:89221.start
Scaffold10      89221   89221   Scaffold10-88994:89221.end
Scaffold10      90067   90067   Scaffold10-90067:91033.start
Scaffold10      91033   91033   Scaffold10-90067:91033.end
Scaffold10      94503   94503   Scaffold10-94503:94704.start
Scaffold10      94704   94704   Scaffold10-94503:94704.end
Scaffold14      30015   30015   Scaffold14-30015:30210.start
Scaffold14      30210   30210   Scaffold14-30015:30210.end
Scaffold10      95463   95463   Scaffold10-95463:95627.start
Scaffold10      95627   95627   Scaffold10-95463:95627.end
Scaffold14      41250   41250   Scaffold14-41250:42034.start
Scaffold14      42034   42034   Scaffold14-41250:42034.end
Scaffold10      96605   96605   Scaffold10-96605:97289.start
Scaffold10      97289   97289   Scaffold10-96605:97289.end
Scaffold50      177554  177554  Scaffold50-177554:178683.start
Scaffold50      178683  178683  Scaffold50-177554:178683.end
Scaffold14      57096   57096   Scaffold14-57096:57916.start
Scaffold14      57916   57916   Scaffold14-57096:57916.end
Scaffold14      58799   58799   Scaffold14-58799:59000.start
Scaffold14      59000   59000   Scaffold14-58799:59000.end
Scaffold150     31481   31481   Scaffold150-31481:33376.start
Scaffold150     33376   33376   Scaffold150-31481:33376.end
Scaffold96      1       1       Scaffold96-1:22253.start
Scaffold96      22253   22253   Scaffold96-1:22253.end
Scaffold8       1760778 1760778 Scaffold8-1760778:1761241.start
Scaffold8       1761241 1761241 Scaffold8-1760778:1761241.end
```

For all LGT candidates:
```bash
ls | parallel --dryrun cat {}/results/LGTs.candidateloci.bed "|" parallel --colsep "\t"  echo -e '{1}"\t"{2}"\t"{2}"\t"{1}"-"{2}":"{3}.start"\n"{1}"\t"{3}"\t"{3}"\t"{1}"-"{2}":"{3}.end' ">" {}/results/{}.LGTboundaries.bed
```
###### DOES NOT WORK YET.

Write a GridEngine Script to produce this file for all GAGA genomes:
`nano makeLGTboundaryfile.sh`
```bash
#$ -S /bin/bash
#$ -N LGTboundarybedfile
#$ -cwd
#$ -pe smp 10
#$ -l h_vmem=1G
#$ -o /global/scratch2/j_rink02/master/lgt/0_data/tmp/bedfile.out
#$ -e /global/scratch2/j_rink02/master/lgt/0_data/tmp/bedfile.err
#$ -wd /global/scratch2/j_rink02/master/lgt/0_data

cat $file | parallel --colsep "\t"  echo -e '{1}"\t"{2}"\t"{2}"\t"{1}"-"{2}":"{3}.start"\n"{1}"\t"{3}"\t"{3}"\t"{1}"-"{2}":"{3}.end' > $file.LGTboundaries.bed
```
Execute the script with:
`find */results/LGTs.candidateloci.bed | parallel -I% --max-args 1 qsub -v file="%"
 makeLGTboundarybedfile.sh -o ./tmp/$file.out -e ./tmp/$file.err`

Check if all genomes have this file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/L
GTs.candidateloci.bed.LGTboundaries.bed ]] && echo "$dir"; done
```
### 2. Extract all PacBio reads overlapping these boundaries.

#### 2.1 Merge the bam files for nAo and other
```bash
#samtools merge <outfile.bam> <infile1.bam> <infile2.bam>
#for one file
samtools merge GAGA-0515/results/merged.candidateloci.loose.bam GAGA-0515/results/LGTs.nAo.candidateloci.loose.PacBio.bam GAGA-0515/results/LGTs.candidateloci.loose.PacBio.bam

#for all files
for i in * ; do samtools merge $i/results/merged.candidateloci.loose.bam $i/result
s/LGTs.nAo.candidateloci.loose.PacBio.bam $i/results/LGTs.candidateloci.loose.PacBio.bam; done
```
Check if all genomes have the file `merged.candidateloci.loose.bam`:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/merged.candidateloci.loose.bam ]] && echo "$dir"; done
```

#### 2.2 Extract the reads overlapping the boundaries:
```
bedtools intersect -abam GAGA-0515/results/merged.candidateloci.loose.bam -b GAGA-0515/results/GAGA-0515.LGTboundaries.bed > GAGA-0515/results/GAGA-0515.LGTboundaries.PacBio.overlap.bam
# index the bam file
samtools index GAGA-0515/results/GAGA-0515.LGTboundaries.PacBio.overlap.bam
```

### 3. Expand the required overlap
##### Slop

Now, if we want to expand the required overlap, and want e.g. that reads overlap the LGT start position by 1000 bp on both sides, we can do that with `bedtools slop` and `bedtools intersect`

First, we `slop` the boundaries
```bash
# create a genome file containing the length of all scaffolds in the assembly (required by bedtools)
samtools faidx assemblies/GAGA-0515_SLR-superscaffolder_final_dupsrm_filt.fasta
cut -f 1-2 assemblies/GAGA-0515_SLR-superscaffolder_final_dupsrm_filt.fasta.fai > GAGA-0515.genome

# then slop the bed file by 50 bp.
bedtools slop -i GAGA-0515/results/GAGA-0515.LGTboundaries.bed -g GAGA-0515.genome -b 50 > GAGA-0515.LGTboundaries.50bp.up.down.bed
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
