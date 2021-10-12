# Calculate the number of reads overlapping LGT boundaries

```bash
# change to where you want to work
#cd /global/scratch2/j_rink02/master/lgt/0_data
```

To start, we need 2 files: One `bam` file holding the PacBio reads and one `bed` file containing the LGT start and stop coordinates (`LGTs.candidateloci.loose.bed`).

We start with the file containing the LGT candidates you want to screen to create a `bed` file with the start and stop positions as separate entries.

### 1. Use the file `LGTs.candidateloci.loose.bed` to create a bed file with start+stop coordinates of the LGT as separate entries.

#### 1.1 For one genome only, with `GAGA-0515` as example:
```bash
cd /global/scratch2/j_rink02/master/lgt/0_data

cat GAGA-0515/results/LGTs.candidateloci.loose.bed | parallel --colsep "\t"  echo -e '{1}"\t"{2}"\t"{2}"\t"{1}"-"{2}":"{3}.start"\n"{1}"\t"{3}"\t"{3}"\t"{1}"-"{2}":"{3}.end' > GAGA-0515/results/GAGA-0515.LGTboundaries.bed
```
The above command creates a new file called `GAGA-0515.LGTboundaries.bed` which now contains only the start and stop coordinates of each LGT as a separate entry.

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

#### 1.2 For all LGT candidates:

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
`find */results/LGTs.candidateloci.loose.bed | parallel -I% --max-args 1 qsub -v file="%"
 makeLGTboundaryfile.sh -o ./tmp/$file.out -e ./tmp/$file.err`

The file should be called `LGTs.candidateloci.loose.bed.LGTboundaries.bed` and should be found in every GAGA genome folder.

##### 1.2.1 For all no ant LGT candidates:
Execute the script with:
`find */results/LGTs.nAo.candidateloci.loose.bed | parallel -I% --max-args 1 qsub -v file="%"
 makeLGTboundaryfile.sh -o ./tmp/$file.out -e ./tmp/$file.err`

The file should be called `LGTs.nAo.candidateloci.loose.bed.LGTboundaries.bed` and should be found in every GAGA genome folder.

#### 1.3 Check if all genomes have this file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed ]] && echo "$dir"; done
```

##### 1.3.1 Check if all nA genomes have this file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/LGTs.nAo.candidateloci.loose.bed.LGTboundaries.bed ]] && echo "$dir"; done
```

If all GAGA genomes contain this file, no GAGA genome folder should appear here.
### 2. Extract all reads overlapping these boundaries.

#### 2.1 Merge the bam files for nAo and other
```bash
#command: samtools merge <outfile.bam> <infile1.bam> <infile2.bam>
#for one file (GAGA-0515 as example file):
samtools merge GAGA-0515/results/merged.candidateloci.loose.bam GAGA-0515/results/LGTs.nAo.candidateloci.loose.PacBio.bam GAGA-0515/results/LGTs.candidateloci.loose.PacBio.bam

#for all files:
for i in * ; do samtools merge $i/results/merged.candidateloci.loose.bam $i/result
s/LGTs.nAo.candidateloci.loose.PacBio.bam $i/results/LGTs.candidateloci.loose.PacBio.bam; done
```
Check if all genomes have the file `merged.candidateloci.loose.bam`:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/merged.candidateloci.loose.bam ]] && echo "$dir"; done
```
#### 2.2 Extract the reads overlapping the boundaries:

For one genome only, with `GAGA-0515` as example:
```bash
bedtools intersect -abam GAGA-0515/results/merged.candidateloci.loose.bam -b GAGA-0515/results/GAGA-0515.LGTboundaries.bed > GAGA-0515/results/GAGA-0515.LGTboundaries.PacBio.overlap.bam

# -abam: BAM/BED file A. Each BAM alignment in A is compared to B in search of overlaps.
# -b: One or more BAM/BED file(s).
# index the bam file
samtools index GAGA-0515/results/GAGA-0515.LGTboundaries.PacBio.overlap.bam
```

`bedtools intersect` allows one to screen for overlaps between two sets of genomic features, in this case whether any of the reads from the `merged.candidateloci.loose.bam` file overlap with the LGT boundaries from the `LGT.boundaries.bed` file. The resulting overlap between the LGT boundaries and the reads is written into a new file with the name `$GAGA-id.LGTboundaries.PacBio.overlap.bam`.

For all genomes:
```bash
# calculate overlap with bedtools intersect
for i in * ; do bedtools intersect -abam $i/results/merged.candidateloci.loose.bam -b $i/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed > $i/results/$i.LGTboundaries.PacBio.overlap.bam; done

# index the bam file
for i in * ; do samtools index $i/results/$i.LGTboundaries.PacBio.overlap.bam; done
```

Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.PacBio.overlap.bam ]] && echo "$dir"; done
```

For all nAo genomes:
```bash
# calculate overlap with bedtools intersect
for i in * ; do bedtools intersect -abam $i/results/merged.candidateloci.loose.bam -b $i/results/LGTs.nAo.candidateloci.loose.bed.LGTboundaries.bed > $i/results/$i.nAo.LGTboundaries.PacBio.overlap.bam; done

# index the bam file
for i in * ; do samtools index $i/results/$i.nAo.LGTboundaries.PacBio.overlap.bam; done
```

Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.nAo.LGTboundaries.PacBio.overlap.bam ]] && echo "$dir"; done
```

### 3. Expand the required overlap
##### Slop
`bedtools slop` will increase the size of each feature in a file by a user-defined number of bases. It will restrict resizing to the size of the chromosome (no end above chromosome size and no start below 0). In order to prevent extension beyond scaffold boundaries, `bedtools slop` requires a genome file, defining the length of each scaffold.

The required `genome.file`already exists for every genome in the respective folder (e.g. `GAGA-0515/results/genome.file`).

#### 3.1 Sort the genomes into short-read data and long-read data

```bash
mkdir PacBio
mkdir stLFR
```
Sort according to sequencing information from GAGA project.

Genomes with short-read data only (directory `stLFR`):
```bash
GAGA-0002  GAGA-0331  GAGA-0389  GAGA-0495  GAGA-0553  NCBI-0003  NCBI-0008  NCBI-0013	NCBI-0017
GAGA-0055  GAGA-0366  GAGA-0392  GAGA-0535  GAGA-0554  NCBI-0004  NCBI-0010  NCBI-0014
GAGA-0220  GAGA-0386  GAGA-0406  GAGA-0536  GAGA-0577  NCBI-0006  NCBI-0011  NCBI-0015
GAGA-0234  GAGA-0387  GAGA-0408  GAGA-0537  GAGA-0580  NCBI-0007  NCBI-0012  NCBI-0016
```

Genomes with long-read data (directory `PacBio`):
```bash
GAGA-0001  GAGA-0085  GAGA-0223  GAGA-0307  GAGA-0351  GAGA-0378  GAGA-0463  GAGA-0522	GAGA-0552
GAGA-0003  GAGA-0087  GAGA-0224  GAGA-0328  GAGA-0352  GAGA-0379  GAGA-0485  GAGA-0524	GAGA-0578
GAGA-0004  GAGA-0090  GAGA-0229  GAGA-0330  GAGA-0353  GAGA-0380  GAGA-0491  GAGA-0527	GAGA-0579
GAGA-0014  GAGA-0098  GAGA-0245  GAGA-0332  GAGA-0354  GAGA-0382  GAGA-0494  GAGA-0528	NCBI-0001
GAGA-0020  GAGA-0099  GAGA-0246  GAGA-0333  GAGA-0356  GAGA-0384  GAGA-0502  GAGA-0530	NCBI-0002
GAGA-0024  GAGA-0103  GAGA-0256  GAGA-0334  GAGA-0358  GAGA-0391  GAGA-0503  GAGA-0531	NCBI-0005
GAGA-0025  GAGA-0109  GAGA-0266  GAGA-0335  GAGA-0359  GAGA-0393  GAGA-0505  GAGA-0532	NCBI-0009
GAGA-0026  GAGA-0114  GAGA-0275  GAGA-0336  GAGA-0360  GAGA-0395  GAGA-0510  GAGA-0533	OUT-0001
GAGA-0028  GAGA-0177  GAGA-0288  GAGA-0337  GAGA-0361  GAGA-0396  GAGA-0511  GAGA-0534	OUT-0002
GAGA-0063  GAGA-0187  GAGA-0300  GAGA-0338  GAGA-0362  GAGA-0401  GAGA-0512  GAGA-0538
GAGA-0074  GAGA-0198  GAGA-0301  GAGA-0340  GAGA-0363  GAGA-0404  GAGA-0513  GAGA-0539
GAGA-0080  GAGA-0199  GAGA-0302  GAGA-0341  GAGA-0364  GAGA-0405  GAGA-0515  GAGA-0540
GAGA-0082  GAGA-0200  GAGA-0303  GAGA-0343  GAGA-0365  GAGA-0407  GAGA-0517  GAGA-0541
GAGA-0083  GAGA-0221  GAGA-0304  GAGA-0346  GAGA-0374  GAGA-0413  GAGA-0520  GAGA-0543
GAGA-0084  GAGA-0222  GAGA-0306  GAGA-0350  GAGA-0376  GAGA-0454  GAGA-0521  GAGA-0550
```

#### 3.2 Genomes with long-read PacBio data
`cd PacBio`
##### Slop
We `slop` the boundaries, so that the boundary is expanded by 1000 bp on each side.

For one candidate (with `GAGA-0515` as an example):
```bash
# Expand the LGT boundaries by 1000 bp for long-read data with bedtools slop
bedtools slop -i GAGA-0515/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed -g GAGA-0515/results/genome.file -b 1000 > GAGA-0515/results/GAGA-0515.LGTboundaries.1000bp.up.down.bed
# -i: input file
# -g: genome file
# -b: increase the BED/GFF file by the same number of bp in each direction.
```
Output of the file `GAGA-0515.LGTboundaries.1000bp.up.down.bed`:
```bash
Scaffold10      87994   89994   Scaffold10-88994:89221.start
Scaffold10      88221   90221   Scaffold10-88994:89221.end
Scaffold10      89067   91067   Scaffold10-90067:91033.start
Scaffold10      90033   92033   Scaffold10-90067:91033.end
Scaffold10      94463   96463   Scaffold10-95463:95627.start
Scaffold10      94627   96627   Scaffold10-95463:95627.end
Scaffold10      93503   95503   Scaffold10-94503:94704.start
Scaffold10      93704   95704   Scaffold10-94503:94704.end
Scaffold14      56096   58096   Scaffold14-57096:57916.start
Scaffold14      56916   58916   Scaffold14-57096:57916.end
Scaffold10      95605   97605   Scaffold10-96605:97289.start
Scaffold10      96289   98289   Scaffold10-96605:97289.end
Scaffold14      29015   31015   Scaffold14-30015:30210.start
Scaffold14      29210   31210   Scaffold14-30015:30210.end
Scaffold14      40250   42250   Scaffold14-41250:42034.start
Scaffold14      41034   43034   Scaffold14-41250:42034.end
Scaffold96      1       1001    Scaffold96-1:22253.start
Scaffold96      21253   23253   Scaffold96-1:22253.end
Scaffold14      57799   59799   Scaffold14-58799:59000.start
Scaffold14      58000   60000   Scaffold14-58799:59000.end
Scaffold150     30481   32481   Scaffold150-31481:33376.start
Scaffold150     32376   34376   Scaffold150-31481:33376.end
Scaffold50      176554  178554  Scaffold50-177554:178683.start
Scaffold50      177683  179683  Scaffold50-177554:178683.end
Scaffold8       1759778 1761778 Scaffold8-1760778:1761241.start
Scaffold8       1760241 1762241 Scaffold8-1760778:1761241.end
```

For all genomes:
```bash
# Expand the LGT boundaries by 1000 bp for long-read data with bedtools slop
for i in *; do bedtools slop -i $i/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed -g $i/results/genome.file -b 1000 > $i/results/$i.LGTboundaries.1000bp.up.down.bed; done
```

Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.1000bp.up.down.bed ]] && echo "$dir"; done
```
##### intersect
Now we can `intersect` the reads with the expanded LGT boundary. Here, we add the `-F 1` option, to require that the entire stretch defined in the bed file is covered by a given read.

```bash
# extract reads completely overlapping the +-1000 bp of the boundary
for i in *; do bedtools intersect -F 1 -abam $i/results/merged.candidateloci.loose.bam -b $i/results/$i.LGTboundaries.1000bp.up.down.bed > $i/results/$i.LGTboundaries.1000bp.PacBio.overlap.bam; done
```
Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.1000bp.PacBio.overlap.bam ]] && echo "$dir"; done
```

#### 3.3 Genomes with short-read stLFR data
`cd stLFR`
##### Slop
We `slop` the boundaries, so that the boundary is expanded by 25 bp on each side.

For all genomes:
```bash
# Expand the LGT boundaries by 25 bp for short-read data with bedtools slop
for i in *; do bedtools slop -i $i/results/LGTs.candidateloci.loose.bed.LGTboundaries.bed -g $i/results/genome.file -b 25 > $i/results/$i.LGTboundaries.25bp.up.down.bed; done
```
Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.25bp.up.down.bed ]] && echo "$dir"; done
```
##### intersect
```bash
# extract reads completely overlapping the +-25 bp of the boundary
for i in *; do bedtools intersect -F 1 -abam $i/results/merged.candidateloci.loose.bam -b $i/results/$i.LGTboundaries.25bp.up.down.bed > $i/results/$i.LGTboundaries.25bp.PacBio.overlap.bam; done
```
Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.25bp.PacBio.overlap.bam ]] && echo "$dir"; done
```
## 4. Filter reads based on mapping quality
Finally, we filter reads that are not mapped very well to the region (e.g. they could also map somewhere else). Also, we remove reads that are in the file twice (because we merged the nAo and the other bam file above)

#### 4.1 Genomes with PacBio data:
`cd PacBio`

For one genome only (with `GAGA-0200` as an example):
```bash
# remove not primary alignments and supplementary alignments from bam file. Google "sam flags explained" for details
## Additionally remove duplicated reads (due to the merging of nAo Bam and other bam above), with the awk command
samtools view -F 1024 -F 256 -F 2048 GAGA-0200/results/GAGA-0200.LGTboundaries.1000bp.PacBio.overlap.bam -h |awk '!visited[$0]++|| $1 ~ /^@/' | samtools view -bS - > GAGA-0200/results/GAGA-0200.LGTboundaries.1000bp.good.PacBio.overlap.bam

# -F 1024: read is PCR or optical duplicate
# -F 256: not primary alignment
# -F 2048: supplementary alignment
```
For all genomes with long-read PacBio data:
```bash
for i in *; do samtools view -F 1024 -F 256 -F 2048 $i/results/$i.LGTboundaries.1000bp.PacBio.overlap.bam -h |awk '!visited[$0]++|| $1 ~ /^@/' | samtools view -bS - > $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bam; done
```
Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.1000bp.good.PacBio.overlap.bam ]] && echo "$dir"; done
```
#### 4.2 Genomes with stLFR data:
`cd stLFR`

For all genomes with short-read stLFR data:
```bash
for i in *; do samtools view -F 1024 -F 256 -F 2048 $i/results/$i.LGTboundaries.25bp.PacBio.overlap.bam -h |awk '!visited[$0]++|| $1 ~ /^@/' | samtools view -bS - > $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bam; done
```
Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.25bp.good.stLFR.overlap.bam ]] && echo "$dir"; done
```
## 5. Stats

### 5.1 Count the number of reads per boundary

`bedtools coverage` computes the depth and breadth of coverage of features in a file B on features in file A. It counts the number of features that overlap an interval in file A, but also computes the fraction of bases in the interval A that were overlapped by one or more features.

#### 5.1.1 Long-read data:
`cd PacBio`
For one genome only (`GAGA-0515` as example):
```bash
# Count how many reads you get per boundary
bedtools coverage -f 1 -b GAGA-0515/results/GAGA-0515.LGTboundaries.1000bp.good.PacBio.overlap.bam -a GAGA-0515/results/GAGA-0515.LGTboundaries.1000bp.up.down.bed  -counts > GAGA-0515/results/GAGA-0515.LGTboundaries.1000bp.good.PacBio.overlap.bed
# -a: BAM/BED/GFF file A: Each feature in A is compared to B in search of overlaps
# -b: One or more BAM/BED/GFF files.
# -f: Minimum overlap required as a fraction of A
# -counts: Only report the count of overlaps
```

Output of the file `GAGA-0515.LGTboundaries.1000bp.good.PacBio.overlap.bed`:
```bash
Scaffold10      87994   89994   Scaffold10-88994:89221.start    0
Scaffold10      88221   90221   Scaffold10-88994:89221.end      2
Scaffold10      89067   91067   Scaffold10-90067:91033.start    8
Scaffold10      90033   92033   Scaffold10-90067:91033.end      2
Scaffold10      94463   96463   Scaffold10-95463:95627.start    7
Scaffold10      94627   96627   Scaffold10-95463:95627.end      2
Scaffold10      93503   95503   Scaffold10-94503:94704.start    16
Scaffold10      93704   95704   Scaffold10-94503:94704.end      14
Scaffold14      56096   58096   Scaffold14-57096:57916.start    336
Scaffold14      56916   58916   Scaffold14-57096:57916.end      313
Scaffold10      95605   97605   Scaffold10-96605:97289.start    7
Scaffold10      96289   98289   Scaffold10-96605:97289.end      16
Scaffold14      29015   31015   Scaffold14-30015:30210.start    0
Scaffold14      29210   31210   Scaffold14-30015:30210.end      0
Scaffold14      40250   42250   Scaffold14-41250:42034.start    25
Scaffold14      41034   43034   Scaffold14-41250:42034.end      0
Scaffold96      1       1001    Scaffold96-1:22253.start        0
Scaffold96      21253   23253   Scaffold96-1:22253.end  0
Scaffold14      57799   59799   Scaffold14-58799:59000.start    289
Scaffold14      58000   60000   Scaffold14-58799:59000.end      284
Scaffold150     30481   32481   Scaffold150-31481:33376.start   0
Scaffold150     32376   34376   Scaffold150-31481:33376.end     0
Scaffold50      176554  178554  Scaffold50-177554:178683.start  72
Scaffold50      177683  179683  Scaffold50-177554:178683.end    78
Scaffold8       1759778 1761778 Scaffold8-1760778:1761241.start 57
Scaffold8       1760241 1762241 Scaffold8-1760778:1761241.end   54
```

Write the result in a `tsv` file (`GAGA-0515` as an example):
```bash
cat GAGA-0515/results/GAGA-0515.LGTboundaries.1000bp.good.PacBio.overlap.bed | paste - - > GAGA-0515/results/GAGA-0515.LGTboundaries.1000bp.good.PacBio.overlap.tsv
```

Output of the `GAGA-0515.LGTboundaries.1000bp.good.PacBio.overlap.tsv`file:
```bash
Scaffold10	87994	89994	Scaffold10-88994:89221.start	0	Scaffold10	88221	90221	Scaffold10-88994:89221.end	2
Scaffold10	89067	91067	Scaffold10-90067:91033.start	8	Scaffold10	90033	92033	Scaffold10-90067:91033.end	2
Scaffold10	94463	96463	Scaffold10-95463:95627.start	7	Scaffold10	94627	96627	Scaffold10-95463:95627.end	2
Scaffold10	93503	95503	Scaffold10-94503:94704.start	16	Scaffold10	93704	95704	Scaffold10-94503:94704.end	14
Scaffold14	56096	58096	Scaffold14-57096:57916.start	336	Scaffold14	56916	58916	Scaffold14-57096:57916.end	313
Scaffold10	95605	97605	Scaffold10-96605:97289.start	7	Scaffold10	96289	98289	Scaffold10-96605:97289.end	16
Scaffold14	29015	31015	Scaffold14-30015:30210.start	0	Scaffold14	29210	31210	Scaffold14-30015:30210.end	0
Scaffold14	40250	42250	Scaffold14-41250:42034.start	25	Scaffold14	41034	43034	Scaffold14-41250:42034.end	0
Scaffold96	1	1001	Scaffold96-1:22253.start	0	Scaffold96	21253	23253	Scaffold96-1:22253.end	0
Scaffold14	57799	59799	Scaffold14-58799:59000.start	289	Scaffold14	58000	60000	Scaffold14-58799:59000.end	284
Scaffold150	30481	32481	Scaffold150-31481:33376.start	0	Scaffold150	32376	34376	Scaffold150-31481:33376.end	0
Scaffold50	176554	178554	Scaffold50-177554:178683.start	72	Scaffold50	177683	179683	Scaffold50-177554:178683.end	78
Scaffold8	1759778	1761778	Scaffold8-1760778:1761241.start	57	Scaffold8	1760241	1762241	Scaffold8-1760778:1761241.end	54
```
The file `GAGA-0515.LGTboundaries.1000bp.good.PacBio.overlap.tsv` contains now one row for each candidate and gives the following columns
```
col1: Scaffold
col2: start region - 1000 bp
col3: start region + 1000 bp
col4: "complete candidate region.start"
col5: number of reads overlapping start +- 1000 bp
col6: Scaffold
col7: end region - 1000 bp
col8: end region + 1000 bp
col9: "complete candidate region.end"
col10: number of reads overlapping end +- 1000 bp
```

For all PacBio genomes:
```bash
# Count how many reads you get per boundary
for i in *; do bedtools coverage -f 1 -b $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bam -a $i/results/$i.LGTboundaries.1000bp.up.down.bed  -counts > $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bed; done

for i in *; do cat $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bed | paste - - > $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.tsv; done
```
Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.1000bp.good.PacBio.overlap.tsv ]] && echo "$dir"; done
```
#### 5.1.2 Short-read data:
`cd stLFR`

For all short-read genomes:
```bash
# Count how many reads you get per boundary
for i in *; do bedtools coverage -f 1 -b $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bam -a $i/results/$i.LGTboundaries.25bp.up.down.bed  -counts > $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bed; done

for i in *; do cat $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bed | paste - - > $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.tsv; done
```
Check if all genomes have the file:
```bash
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/results/$dir.LGTboundaries.25bp.good.stLFR.overlap.tsv ]] && echo "$dir"; done
```


### 5.2 Calculate the average read length across all genomes to estimate LGT boundary expansion
`bedtools bamtobed`is a conversion utility that converts sequence alignments in BAM format into BED.
By default, each alignment in the BAM file is converted to a 6 column BED.

Plot the average read lengths for stLFR-assembled and PacBio-assembled genomes in R and after that, decide how large the LGT boundary expansion should be.

```bash
# calculate average read length for PacBio genomes
for i in *; do bedtools bamtobed -bed12 -i $i/results/$i.LGTboundaries.1000bp.good.PacBio.overlap.bam|cut -f 11|awk '{ total += $1; count++ } END { print total/count }'; done
# -bed12: Write a "blocked" BED format. This converts "spliced" BAM alignments to BED12.
# -i: input file in BAM format

# calculate average read length for stLFR genomes
for i in *; do bedtools bamtobed -bed12 -i $i/results/$i.LGTboundaries.25bp.good.stLFR.overlap.bam|cut -f 11|awk '{ total += $1; count++ } END { print total/count }'; done
```
