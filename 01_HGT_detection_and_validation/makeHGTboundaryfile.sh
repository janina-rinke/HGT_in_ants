#$ -S /bin/bash
#$ -N HGTboundarybedfile
#$ -cwd
#$ -pe smp 10
#$ -l h_vmem=1G

cat $file | parallel --colsep "\t"  echo -e '{1}"\t"{2}"\t"{2}"\t"{1}"-"{2}":"{3}.start"\n"{1}"\t"{3}"\t"{3}"\t"{1}"-"{2}":"{3}.end' \
> $file.LGTboundaries.bed