# BLAST sequence against GAGA genomes

This script describes the process how LGT sequences can be blasted against all GAGA genomes to possibly identify more LGTs that were previously filtered out as being a good candidate.

In case of additional LGT sequences being identified as good LGT candidates, these candidates will be intersected with all initially predicted LGT loci of the LGT finder pipeline within the respective genome and run again with the dfast prokaryotic gene annotation.


### 1 Prepare the GAGA genomes as a database
##### 1.1 Add GAGAid to fasta headers:
Copy genomes to another folder to work on them, e.g. in `assemblies_copy` folder:
```bash
cd /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy

# Remove everything after the GAGAid
# e.g. GAGA-0001_SLR-superscaffolder_final_dupsrm_filt.fasta turns into GAGA-0001

for f in * ; do mv -n "$f" "${f%%\_*}"; done

# To check if it will work out the way it is expected:
for f in * ; do echo mv -n "$f" "${f%%\_*}"; done
```

For every Scaffold, add the GAGAid to the FASTA header (which starts with a `>`)
```bash
for f in *; do sed -i '' -e "s/^>/>${f%\_}_/g" "${f}"; done

### The above sed command will replace every line starting (^) with ">" with the filename until ".fasta" is matched.
# -i: --in-place File will be changed
# -e: expression

# EXPLANATION:
# sed 's/text_to_replace/new_text/g'

# s= substitute
# g= global, the pattern will be replaced in every single occurrence
```
Now, every file consists of the GAGAid and the scaffold name, which helps identifying the GAGA genome and location in case of a BLAST hit.

##### 1.2 Concatenate all GAGA genomes into one large genome database
Concatenate all files and make a large database to blast against.
```bash
cat * >> GAGA_genomes_database
```
### 2 Sequence search against the large GAGA genomes database
Make a local BLAST search to detect whether a query sequence exists in any of the GAGA genomes.

GridEngine Script to make `GAGA_genomes_database` a BLAST db:

`cd /global/scratch2/j_rink02/master/lgt/0_data`
```bash
#$ -S /bin/bash
#$ -N makeGAGAdb
#$ -cwd
#$ -pe smp 20
#$ -l h_vmem=3G
#$ -o /global/scratch2/j_rink02/master/lgt/0_data/local_blast/tmp/makeblastdb.out
#$ -e /global/scratch2/j_rink02/master/lgt/0_data/local_blast/tmp/makeblastdb.err
#$ -wd /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy

makeblastdb -dbtype nucl -in /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy/GAGA_genomes_database
```

Blast any sequence of interest against the GAGA genomes:
```bash
blastn -query Sequence.txt -db ../assemblies_copy/GAGA_genomes_database -outfmt "6 std qlen" -out blast_Sequence_GAGAgenomes.out
```

##### 2.1 Sequence search with the Lyzozyme sequence from GAGA-0098 (Harpagoxenus sublaevis).

Input file:`Lyzozyme_sequence_GAGA-0098.txt`

Output file: `blast_Lyzozyme_GAGA-0098.out`

```bash
GAGA-0098.LGT	GAGA-0098_Scaffold31	100.000	624	0	0	1	624	1153305	1153928	0.0	1153	624
GAGA-0098.LGT	GAGA-0099_Scaffold20	97.756	624	14	0	1	624	2658767	2658144	0.0	1075	624
GAGA-0098.LGT	GAGA-0103_Scaffold3	86.446	605	82	0	1	605	176889	177493	0.0	664	624
GAGA-0098.LGT	GAGA-0222_Scaffold8	92.447	609	46	0	1	609	5585551	5584943	0.0	870	624
GAGA-0098.LGT	GAGA-0223_Scaffold1	92.118	609	48	0	1	609	1849423	1848815	0.0	859	624
GAGA-0098.LGT	GAGA-0224_Scaffold2	91.790	609	50	0	1	609	5524378	5524986	0.0	848	624
GAGA-0098.LGT	GAGA-0256_Scaffold2	89.109	606	64	2	1	605	7916749	7917353	0.0	752	624
GAGA-0098.LGT	GAGA-0288_Scaffold13	92.118	609	48	0	1	609	4155591	4154983	0.0	859	624
GAGA-0098.LGT	GAGA-0407_Scaffold17	90.728	604	56	0	1	604	738017	738620	0.0	806	624
GAGA-0098.LGT	GAGA-0463_Scaffold5	97.596	624	15	0	1	624	10067288	10067911	0.0	1070	624
GAGA-0098.LGT	GAGA-0510_Scaffold38	92.775	609	44	0	1	609	1251324	1251932	0.0	881	624
GAGA-0098.LGT	GAGA-0511_Scaffold84	91.790	609	50	0	1	609	595715	596323	0.0	848	624
GAGA-0098.LGT	GAGA-0512_Scaffold4	91.954	609	49	0	1	609	5228196	5227588	0.0	854	624
GAGA-0098.LGT	GAGA-0513_Scaffold6	91.790	609	50	0	1	609	3023405	3022797	0.0	848	624
```
This BLAST search did not result in any other species than the ones detected before which were predicted to incorporate a Lyzozyme. However, the three Carebara spezies with a predicted Lyzozyme did not appear here.

##### 2.2 Sequence search with the Lyzozyme sequence from GAGA-0328 (Carebara)

Input file: `Lyzozyme_sequence_GAGA-0328.txt`
Output file: `blast_Lyzozyme_GAGA-0328.out`

```bash
GAGA-0328.LGT	GAGA-0328_Scaffold36	100.000	618	0	0	1	618	505676	506293	0.0	1142	618
GAGA-0328.LGT	GAGA-0331_Scaffold2506	89.286	616	57	5	1	616	5742	6348	0.0	763	618
GAGA-0328.LGT	GAGA-0378_Scaffold29	92.857	616	44	0	1	616	1056536	1057151	0.0	894	618
GAGA-0328.LGT	GAGA-0382_Scaffold1	97.735	618	14	0	1	618	29206274	29205657	0.0	1064	618
GAGA-0328.LGT	GAGA-0533_Scaffold25	83.359	637	76	11	1	616	1172729	1173356	1.33e-156	562	618
GAGA-0328.LGT	GAGA-0578_Scaffold3	90.792	619	56	1	1	618	7430847	7430229	0.0	826	618
GAGA-0328.LGT	GAGA-0579_Scaffold2	93.861	619	36	2	1	618	762894	763511	0.0	931	618
```
Here, all Carebara species occurred as hits showing that these are two independent and separate LGT events.

### 3 Intersect newly detected sequences with LGT loci
Take the `LGTs.candidateloci.loose.bed` file and intersect this file with the coordinates from the BLAST search. Use those intersected hits and re-run dfast with the coordinates of these LGT loci.
