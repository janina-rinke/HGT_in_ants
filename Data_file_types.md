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
