# Data and file types

`.gz files` gunzip

- compression algorithm that reduces the size of a file
- often used to compress web elements for faster page loading


`.gff files` general feature format
- one line per feature
- each line contains 9 columns of Data

sample GFF output from Ensembl data dump:
```
   [1]         [2]         [3]       [4]       [5]      ....    [9]
[seqname]   [source]    [feature]  [start]    [end]     .... [attribute]

     X       Ensembl      gene     2410918   2410992    ....    hid-trf
```
file example:
`GCF_000204515.1_Aech_3.9_genomic.gff`


`gtf files` general transfer format
- identical to gff version 2
- fields are tab-separated


`SAM file` Sequence Alignment/Map Format
- text format for storing sequence data in tab delimited columns
- consists of header section (optional) and alignment section
- Header lines start with @
- Each alignment line has 11 mandatory fields for essential alignment information (e.g. mapping position, variable number of optional fields for flexible or aligner specific information)


`.bam files` BAM file
- compressed binary version of a SAM file
- used to represent aligned sequences up to 128 MB
- Header: contains information about entire file (sample name, sample length, alignment method..); Alignments in the alignment section are associated with specific information in the header section
- Alignments: contains read name, read sequence, read quality, alignment information and custom tags
- Read name: includes the chromosome, start coordinate, alignment quality, and match descriptor string
- The alignments section includes the following information for each read or pair:

`RG`: Read group, which indicates the number of reads for a specific sample
`BC`: Barcode tag, which indicates demultiplexed sample ID associated with the read
`SM`: Single-end alignment quality
`AS`: Paired-end alignment quality
`NM`: Edit distance tag, which records the Levenshtein distance between the read and the reference
`XN`: Amplicon name tag, which records the amplicon tile ID associated with the read

file example:
`LGTs.5kb.candidateregions.PacBio.bam`


`.bam.bai files ` BAM index files
- provide an index of the corresponding BAM file

Create an indexed file of a BAM file:
```bash
samtools index *.bam *.bam.bai
```

## Bedtools
`bedtools bamtobed`: converts sequence alignments in `BAM` format into `BED`
By default, each alignment in the BAM file is converted to a 6 column BED

```
bedtools bamtobed -i *.bam
bedtools bamtobed -i *.bam | grep "Scaffold8"
```
Output:
```
Scaffold8	1741882	1765016	m54220_190716_142950/68354179/833_24479	60	+
Scaffold8	1742240	1764794	m54212_190712_185746/12910675/11749_35641	60	+
Scaffold8	1744601	1761860	m54220_190716_041618/14025124/0_19153	60	+
```


`bedtools makewindows`: Makes adjacent or sliding windows across a genome or BED file.
```
bedtools makewindows -g <genome>
                     -b <bed>

                     -w <window size>
                     -s <step size>
                     -n <number of windows>


# -w: Divide to fixed-size windows. (e.g. -w 500, divides into windows of 500 bp)
# -s: Step size, i.e. how many basepairs to step before creating a new window. Used to create overlapping windows.
# -n: Divide each input interval to a fixed number of windows, with varying window sizes.
```


`bedtools intersect`: Asking whether features overlap with each other = feature intersection. Allows one to screen for overlaps between two sets of genomic features. Works with both `BED/GFF/VCF`and `BAM` files as input.

```bash
bedtools intersect -a <file>.
                   -b <file1, file2 ...>



# -a: BAM/BED/GFF/VCF file, Each feature in A is compared to B in search of overlaps.
# -b: one ore more files of the above format.
# -abam: BAM file A. Each BAM file is compared to B in search of overlaps.
# -c: For each entry in A, report the number of hits in B while restricting to -f. Reports 0 for A entries that have no overlap with B. Restricted -f, -F, -r and -s
# -f: Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1 bp)
# -F: Minimum overlap required as a fraction of B.
```

Last changed 2021/05/11
#### Added to GitHub 
