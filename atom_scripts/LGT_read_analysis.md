LGT analysis: Extraction of border regions in PACBIO reads

## Viewing PacBio reads in the IGV web browser application

Create a fasta index file of the genome file to view the data:

`samtools` implements various utilities for post-processing alignments in SAM, BAM and CRAM formats, including indexing, variant calling and a simple alignment viewer.

```bash
# With samtools
samtools faidx GAGA-0515_nextpolish_final.fasta
```
Upload the indexed genome file in the IGV web browser
https://igv.org/app/

Create an indexed file of the BAM track file you would like to upload:
```bash
samtools index LGTs.5kb.candidateregions.PacBio.bam LGTs.5kb.candidateregions.PacBio.bai
```

You can then search for the scaffold of interest.

## Extract PacBio reads from a region of interest

To extract reads of a SAM/BAM file by reference and region:
```bash
# -b needs to be specified if you want the output to be a .bam file
# -h needs to be specified to include the header in the output
samtools view -b -h input.bam "Scaffold96:22250-22260" > output.BAM
```
This extracts reads of the border region which you have specified.

To extract multiple regions of a bam file:
```bash
samtools view -b input.bam "Scaffold96:22250-22260 Scaffold96:30000-30010 Scaffold96:40000-40010" > output.BAM
```

## Count the number of reads in a region of interest

```bash
# if you have extracted a region of interest before from bam file:
samtools view -c output.bam
# view reads in specified region in original .bam file:
samtools view LGTs.5kb.candidateregions.PacBio.bam "Scaffold96:22250-22260" | wc -l
```

To get the headers from a fasta
```bash
awk < GAGA-0515_nextpolish_final.fasta '/^>/ { print $0 }' #get all headers
awk < GAGA-0515_nextpolish_final.fasta '/^>Scaffold9/ { print $0 }' #get only headers starting with 9
```

## Identify start and stop coordinates of all reads from LGT candidates

```bash
# Shows all start and stop coordinates from this file
bedtools bamtobed -i *.bam
```
## Identify missassemblies by using bedtools
```bash
# Divide the border region on the scaffold of interest (e.g. Scaffold96) into windows, requires BED file
bedtools bamtobed -i LGTs.5kb.candidateregions.PacBio.bam > LGTs.coordinates.bed
cat LGTs.coordinates.bed | awk '$1=="Scaffold96" && $2>22000 { print }' > LGTs.Scaffold96_22000.coordinates.bed

# Windows will be created for each interval in the file, Scaffold96 is now divided into windows of 100 bp.
bedtools makewindows -b LGTs.Scaffold96_22000.coordinates.bed -w 100 > LGTs.Scaffold96_22000.windows.bed

# Find overlaps of reads between the windows of 100 bp and the read coordinates from LGTs.coordinates.bed
bedtools intersect -a LGTs.Scaffold96.windows.bed -b LGTs.coordinates.bed -f 1 -c
```

## FASTA manipulation
`seqkit`
`seqtk`
`biowak`

Last updated 2021/05/11
#### Added the file to GitHub on 11/05/2021. 
