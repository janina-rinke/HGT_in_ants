# Reannotation with DFAST

## Extend the locus by 1000 bp on each side and rerun DFAST bash script
cat locus.txt | parallel --col-sep "\t" --dryrun --max-args 1 \
qsub -v GAGAid={1} -v scaffold={2} -v start={3} -v stop={4} rerun.dfast.sh

## Re-calculate start and stop codons from the dfast reannotation.
for i in * ; do seqkit fx2tab $i/cds.fna |awk -F '\t' '{print $1,substr($2,1,3),substr($2,length($2)-2,length($2))}' \
> $i.start.stop.codons.tsv ; done

## Prints out the file directory names as lines, does not work with MacOS
printf "%s\n" *.start.stop.codons.tsv | xargs -n1 -d $'\n' bash -c 'xargs -n1 -d $'\''\n'\'' printf "%s,%s\n" "$1" <"$1"' -- \
> all.reannotated.start.stop.codons.tsv

## Merge all reannotated bam files and extract unique read counts
for i in */*.1000.out; do samtools merge  $i/mergedRNAseq.bam $i/*LGTregion.bam; done

## Extract only the unique read counts out of every `mergedRNAseq.bam` file. 
for i in */*; do sambamba view -h -F [NH]==1 $i/mergedRNAseq.bam > $i/uniquely_mapped.LGTregion.bam; done

# Sanity check
samtools view mergedRNAseq.bam | awk '{print $1}' | sort | uniq | wc -l
samtools view uniquely_mapped.LGTregion.bam | awk '{print $1}' | sort | uniq | wc -l

for i in */*; do samtools sort \
-o $i/uniquely_mapped.LGTregion.sorted.bam $i/uniquely_mapped.LGTregion.bam; done
