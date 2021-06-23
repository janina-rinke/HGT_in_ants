# Calculating the coverage of candidate sequences and within the genome
# Janina Rinke, 30.05.2021

- Based on mapped PacBio reads
- Summary: Mapping of reads against genome -> Divide into windows -> calculate coverage in every window

### Commands to calculate coverage

1. Mapping of PacBio reads against the genome with `minimap2`

```
minimap2 -t 40 -ax sr genome.fa longreads.fq.gz > mapping/${id}.longread.sam
```

2. Dividing the genome in 2000 bp lengthy sliding windows which overlap by 500 bp
```
bedtools makewindows -g genome.lengths.tsv -w $windowSize -s $stepSize > genome.overlappingwindows.bed
```
Output:
scaffold1 1        2000
scaffold1 1500 3500
scaffold1 3000 4500

3. Calculate the coverage of PacBio reads in every window with `bedtools coverage`
```
bedtools coverage -sorted -a genome.overlappingwindows.bed -b mapping/${id}.longread.bam > mapping/genome.overlappingwindows.cov.tsv
```

Output:
Scaffold1	0	    2000	44	2000	2000	1.0000000
Scaffold1	1500	3500	46	2000	2000	1.0000000
Scaffold1	3000	5000	58	2000	2000	1.0000000
Scaffold1	4500	6500	66	2000	2000	1.0000000
Scaffold1	6000	8000	73	2000	2000	1.0000000
Scaffold1	7500	9500	66	2000	2000	1.0000000
Scaffold1	9000	11000	74	2000	2000	1.0000000
Scaffold1	10500	12500	79	2000	2000	1.0000000
Scaffold1	12000	14000	73	2000	2000	1.0000000
Scaffold1	13500	15500	49	2000	2000	1.0000000

4. Connecting LGT candidates with the coverage
```
bedtools slop -b 20000 -i LGTs.candidateloci.loose.bed -g genome.file|bedtools intersect -b genome.overlappingwindows.cov.tsv -a - -wao |bedtools merge -i - -c 4,10,11,12,13,14 -o first,collapse,collapse,collapse,collapse,collapse|bedtools sort -i > LGTs.candidateloci.loose.coverage.bed
```
Detailed explanation of the above command:
```
# Candidates get enlarged by 20kb up- and downstream
bedtools slop -b 20000 -i LGTs.candidateloci.loose.bed -g genome.file

# Then "enlarged" LGT candidates get overlapped with coverage
bedtools intersect -b genome.overlappingwindows.cov.tsv -a - -wao
```
Hereby, several sliding windows overlap with the LGT candidate.
Output:
Scaffold  1  0  4305 LGTcandidate1 Scaffold1	0	    2000	44	2000	2000	1.0000000
Scaffold  1  0  4305 LGTcandidate1 Scaffold1	1500	3500	46	2000	2000	1.0000000
Scaffold  1  0  4305 LGTcandidate1 Scaffold1	3000	5000	58	2000	2000	1.0000000

`bedtools merge` now summarizes all lines which have the same LGT candidate. `-o`option tells how columns get summarized out of different lines. "first" only keeps the first entry (e.g. Scaffold1 only once). `collapse` adds all entries one after another separated by a comma (e.g. 0,1500,3000).

```
bedtools merge -i - -c 4,10,11,12,13,14 -o first,collapse,collapse,collapse,collapse,collapse
```
Output:
Scaffold  1  0  4305 LGTcandidate1 Scaffold1	0,1500,3000	    2000,3500,5000	44,46,58	2000,2000,2000	2000,2000,2000	1.0000000,1.0000000,1.0000000

Then, everything gets sorted by LGT candidate:
```
bedtools sort -i > LGTs.candidateloci.loose.coverage.bed
```

Explanation of columns in .tsv file:
covWindowStart: Gibt die Startkoordinaten der “sliding windows” an die mit dem (um 20kb erweiterten) LGT Kandidaten überlappen

covWindowEnd: Gibt die Endkoordinate der “sliding windows” an

cov: Gibt an wie viele gemappte PacBio reads mit dem “sliding window” überlappen.

covBases: Gibt laut bedtools manual an: The number of bases that had non-zero coverage. Also wie viele Basen in einem sliding window mindestens von einem PacBio read gecovered waren.

windowLength: Gibt and covWindowEnd-covWindowStart für jedes window (ist fast immer 2000, außer am Ende des scaffolds!)

5. Plot the coverage - calculate average coverage over whole genome
Based on file `genome.overlappingwindows.cov.tsv`
```
globalCoverage<-read.csv("genome.overlappingwindows.cov.tsv",sep="\t",F)
colnames(globalCoverage) <- c("scaffold","windowstart","windowend","coverage","covBases","windowlength","fractionCovered")
averageGlobalCoverage    <- median(globalCoverage$coverage)
```

For every sliding window, calculate coverage relative to "average coverage" and then plot it as log2:
```
p1 <-     ggplot(lgt.stackedL.good[[selection]],aes(x=covWindowCenter,y=log(cov/averageGlobalCoverage,2)))
```


-> Entweder Kandidaten, die sich um einen bestimmten Wert von der average coverage in dem sliding window unterscheiden rausschmeißen ODER wenn die coverage in dem sliding window der Kandidaten sehr low, bzw 0 ist rausschmeißen (damit werden Kandidaten, die keine reads haben, entfernt wie z.B. `GAGA-0343.Scaffold25-3332043-3332282`).

Z.B sowas filtern wie: subset(cov = 0) column `cov` und `c_overlap` in Betracht ziehen

-> Wenn column `c_overlap`= 0 heißt es dann, in diesem Bereich gibt es höchstwahrscheinlich keine mapped PacBio reads?
