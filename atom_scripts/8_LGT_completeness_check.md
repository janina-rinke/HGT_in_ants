# Check LGT candidates for completeness
##### Janina Rinke, 29.10.2021


## 1. Extract for every candidate the best uniprot hit from `dfast`.

Working directory:
```bash
/global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast
```
`dfast` file to use: `genome.gff`.

Example file from `GAGA-0020.Scaffold107.144838-147367.fa`:

```bash
##gff-version 3
Scaffold107:144838-147367	MetaGeneAnnotator	CDS	1	2530	.	-	0	ID=MGA_1;product=hypothetical protein;inference=COORDINATES:ab initio prediction:MetaGeneAnnotator;t
ranslation=LLNRGANVDAKGIDCKTSLHIAAERGYLQIVEHLLKHGARVNSPCFCTCCENYTPLHLAVEKGSEEMTKLLLSRGANVNAKAEDGITALHIAAERGYFQIVEHLLNHEADVNSAITGTYPKGPTPLYLAVEKGSEEMTKLLLSRGANVNAKAEDGITSLHIAAERGYLQIVEHLLNHGADVNSAINWIHYTPL
HLAVEKGNEEIVKLLLSRGANIDAKNQYGKTLYHNAVASKNMKIIELFLNREANVNARNNNSTTLLHLAAKEGDEEIVKLLLSKGANVDAKAEDGITALHIAAERGYLQIVEHLLKHGAHINSTYTSIYREDYTPLHLAVQIGNKEIVKLLLSRGANVDAKGKSGNTSLHIAAEKGYLQIVEHLLNHEADVNSAITGTYPKGYT
PLHLAVEKGNEEIVKLLLSRGANVDAKGIDCKTSLHIAAKRGYSQIVEHLLKHGARVNSPCFCTCCENYTPLHLAVEKGSEEMTKLLLSRGANVNAKAEDGITALHIAAKRGYLQIVEHLLKHEARVNSAYTSTCRKGHTPLYLAVEEGNEEIVKLLLSRGANVDAKGKVGITSLHIAAEKEYLQIVKHLLKYGARVNSAYTST
CREGYTPLHLAVEKGNEEITKLLVSRGADVDATGKYGITSLHIAAEKGYLQIVEHLLNHEADVNSAITGTYPKGYTPLHLAVEKGNEEIVKLLLSRGANVDAKGIDCKTSLHIAAERGYLQIVEHLLKHGARVNSQCFCTCCENYTPLHLAVEKGSEEMTKLLLSRGANVNAKAEDGITALHIAAERGYLQIVEHLLKHGAHIN
STYTSTYREDYTPLHLATERGSKEIIKLLLIRGANVNA;locus_tag=LOCUS_10;note=uniprot_bacteria-0.9.ref:UniRef90_A0A1Q3WAM0 ANK_REP_REGION domain-containing protein (Fragment) (Candidatus Amoebophilus sp. 36-38) [pid:39.
6%25%2C q_cov:99.9%25%2C s_cov:64.6%25%2C Eval:2.6e-153%2C partial hit],WP_010957294.1 hypothetical protein (Treponema denticola ATCC 35405) [pid:28.3%25%2C q_cov:84.1%25%2C s_cov:74.6%25%2C Eval:3.1e-52%
2C partial hit],TIGR:TIGR00870%3B trp: transient-receptor-potential calcium channel protein [Name:TIGR00870%2C Eval:8.9e-117%2C score:389.0%2C bias:56.9],COG:COG0666:ANKYR Ankyrin repeat [Category:T%2C Al
igned:200-388%2C Eval:1.3e-21%2C score:91.8],COG:COG0666:ANKYR Ankyrin repeat [Category:T%2C Aligned:122-256%2C Eval:1.2e-20%2C score:89.1%2C N-term missing],COG:COG0666:ANKYR Ankyrin repeat [Category:T%2
C Aligned:566-739%2C Eval:2.1e-17%2C score:79.5],COG:COG0666:ANKYR Ankyrin repeat [Category:T%2C Aligned:427-597%2C Eval:2.4e-17%2C score:79.5],COG:COG0666:ANKYR Ankyrin repeat [Category:T%2C Aligned:675-
835%2C Eval:7.3e-16%2C score:74.9],COG:COG0666:ANKYR Ankyrin repeat [Category:T%2C Aligned:9-150%2C Eval:2.4e-13%2C score:67.5%2C N-term missing],MGA_1
```

To extract the best uniprot hit for an LGT candidate:

```bash
# go to directory of an LGT candidate within dfast directory
cd GAGA-0020.Scaffold107.144838-147367.fa/

#extract best uniprot hit from genome.gff file
id=$(pwd|perl -pe 's/.*\/(.*?)\..*/$1/g')
awk -v var=${id}  '/##FASTA/ {exit} {if (/^[^##]/) {print var"\t"$0}}' genome.gff | \
perl -pe 's/(.*?)\t(.*?)\t.*note=(.*?\])\,.*/$1\t$2\t$3/g'
```

Output of the above command:
```bash
GAGA-0020	Scaffold107:144838-147367	uniprot_bacteria-0.9.ref:UniRef90_A0A1Q3WAM0 ANK_REP_REGION domain-containing protein (Fragment) (Candidatus Amoebophilus sp. 36-38) [pid:39.6%25%2C q_cov:99.9%25%2C s_cov:64.6%25%2C Eval:2.6e-153%2C partial hit]
```

### 1.1 Extract best uniprot hit for all candidates:
```bash
cd /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast
```

`nano extractUniProt.sh`
```bash
#!/usr/bin/bash
awk -v var=$1  '/##FASTA/ {exit} {if (/^[^##]/) {print var"\t"$0}}' $1/genome.gff | \
perl -pe 's/(.*?)\t(.*?)\t.*note=(.*?\])\,.*/$1\t$2\t$3/g'
```

To execute the bash script:
```bash
bash extractUniprot.sh  OUT-0002.Scaffold1.173150-173918.fa
```

Run it in parallel on all LGTs:
```bash
ls ./|egrep "^GAGA|^OUT|^NCBI" |parallel "bash extractUniprot.sh {}" > LGTs_uniProt_besthit.tsv
```



## 2. Extract start- and stop codons from LGT candidates

The first codon (`start codon`) should usually be ATG (theoretically: GTG or rarely TTG, see https://en.wikipedia.org/wiki/Start_codon).
The last codon (`stop-codon`) should be TAA , TAG , or TGA.

```bash
seqkit fx2tab cds.fna |awk -F '\t' '{print $1,substr($2,1,3),substr($2,length($2)-2,length($2))}'
```

For all candidates:
```bash
for i in * ; do seqkit fx2tab $i/cds.fna |awk -F '\t' '{print $1,substr($2,1,3),substr($2,length($2)-2,length($2))}' > $i.start.stop.codons.tsv ; done

#prints out the file directory names as lines
printf "%s\n" *.start.stop.codons.tsv | xargs -n1 -d $'\n' bash -c 'xargs -n1 -d $'\''\n'\'' printf "%s,%s\n" "$1" <"$1"' -- > all.start.stop.codons.tsv
```

## 3. Count number of predicted CDS for every LGT candidate

```bash
cd /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/dfast/

grep -c CDS */genome.gff
```

Output:
```bash
GAGA-0020.Scaffold107.144838-147367.fa/genome.gff:1
GAGA-0020.Scaffold1.5322877-5323044.fa/genome.gff:1
GAGA-0020.Scaffold17.156258-159230.fa/genome.gff:2
GAGA-0020.Scaffold267.55721-56766.fa/genome.gff:1
```

A column with the number of CDS sequences was added to the file `GAGA.LGT.filtered.completeness.xlsx`. Additionally, the percentage of complete CDS for every LGT candidate (column: `Complete_CDS`) was calculated, based on `q_cov` and `s_cov`in the dfast gff files for every candidate. A candidate CDS was considered as complete, if both q_cov and s_cov were above 75 %. For example, if an LGT candidate had two predicted CDS sequences and one of them was complete, a value of 50% was calculated in the column `Complete_CDS`.

Moreover, the expression of every CDS from every LGT candidate was calculated and visualized with the script `9_calculateCoverage.Rmd` in R. Here, the percentage of expressed CDS per LGT region was again calculated as described for the column `Complete_CDS`. 
