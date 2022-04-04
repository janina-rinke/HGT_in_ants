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

##### 2.3 Sequence search with the Etherase sequence from GAGA-0221 (Camponotus fellah).

Input file:`Etherase_sequence.txt`

Output file: `blast_Etherase_GAGAgenomes.out`

```bash
GAGA-0221.LGT  GAGA-0221_Scaffold3    100.000 909    0      0      1      909    8819326 8820234 0.0    1679   909
GAGA-0221.LGT  NCBI-0005_NW_020229769.1       99.010 909    9      0      1      909    105240 104332 0.0    1629   909
GAGA-0221.LGT  GAGA-0361_Scaffold4    99.010 909    9      0      1      909    3723826 3724734 0.0    1629   909
GAGA-0221.LGT  GAGA-0362_Scaffold15   98.900 909    10     0      1      909    5724715 5725623 0.0    1624   909
GAGA-0221.LGT  GAGA-0396_Scaffold28   98.680 909    12     0      1      909    2026308 2027216 0.0    1613   909
GAGA-0221.LGT  GAGA-0200_Scaffold111  98.130 909    17     0      1      909    363210 364118 0.0    1585   909
GAGA-0221.LGT  GAGA-0374_Scaffold20   96.370 909    33     0      1      909    2858966 2859874 0.0    1496   909
GAGA-0221.LGT  GAGA-0360_Scaffold16   95.820 909    38     0      1      909    5909052 5909960 0.0    1469   909
```

In the BLAST search with GAGA-0221 Etherase sequence as input, in total eight species occurred as hits. All of these species were *Camponotus* species, as well as *Polyrhachis* and *Colobopsis*.

### 3 Intersect newly detected sequences with LGT loci
Take the `LGTs.candidateloci.loose.bed` file and intersect this file with the coordinates from the BLAST search. Use those intersected hits and re-run dfast with the coordinates of these LGT loci.

#### 3.1 For one file only (GAGA-0331 as example):

`cd /global/scratch2/j_rink02/master/lgt/0_data/GAGA-0331/results`

##### 3.1.1 Convert BLAST output file into a bed file

All sequences which are reverse complement and thus have a higher start then stop coordinate, need to be switched to be readable by `bedtools`.

```bash
cat blast_Lyzozyme_GAGA-0328.out | awk '{ if ($10 > $9) print($2"\t"$9"\t"$10); else if ($9 > $10) print($2"\t"$10"\t"$9);}' > blast_Lyzozyme_GAGA-0328.bed
```
Output file `blast_Lyzozyme_GAGA-0328.bed`:
```bash
GAGA-0328_Scaffold36	505676	506293
GAGA-0382_Scaffold1	29205657	29206274
GAGA-0579_Scaffold2	762894	763511
GAGA-0378_Scaffold29	1056536	1057151
GAGA-0578_Scaffold3	7430229	7430847
GAGA-0331_Scaffold2506	5742	6348
GAGA-0533_Scaffold25	1172729	1173356
```
##### 3.1.2 Change GAGA genome bed file to make it fitting to blast bed format
All `LGTs.candidateloci.loose.bed` files need to be changed to be in accordance with the BLAST output bed file (starting with `GAGAid_Scaffold` as ID)

`LGTs.candidateloci.loose.bed` for GAGA-0331:

```bash
Scaffold10016   10      116     CP002371.1;640122;640229;2.481E-09;72;108;76.300;Candidatus Liberibacter solanacearum CLso-ZC1  72      -1      73      0
Scaffold10911   548     642     CP003418.1;796953;797046;9.251E-12;80;94;78.700;Ignavibacterium album JCM 16511 80,94,70,96,75,95,75    -1,-1,-1,-1,-1,-1,-1    81,95,71,97,76,96,76    0,0,0,0,0,0,0
Scaffold11      9074    9257    AP012210.1;1446066;1446250;4.621E-09;71;185;37.400;Candidatus Arthromitus sp. SFB-rat-Yit       71,71   -1,-1   72,72   0,0
Scaffold1300    6482    8312    CP001391.1;741827;743353;0.000E+00;2185;1527;91.800;Wolbachia sp. wRi   2185,1176       2079,700        106,476 1517,812
Scaffold1300    10504   12833   CP001391.1;743798;745788;0.000E+00;2500;1994;87.700;Wolbachia sp. wRi   2500,827        2363,750        137,77  1993,804
Scaffold1463    18709   18952   CP001739.1;1490348;1490601;7.711E-13;84;254;66.900;Sebaldella termitidis ATCC 33386     84,91   -1,-1   85,92   0,0
Scaffold1463    19109   19246   CP003261.1;563097;563235;2.068E-10;76;139;71.900;Clostridium pasteurianum BC1   76      -1      77      0
Scaffold2142    3168    3557    CP003884.1;935415;935756;2.775E-73;285;342;78.700;Wolbachia endosymbiont of Drosophila simulans wHa     285,338 130,130 155,208 10,10
Scaffold21535   668     740     CP000738.1;2016679;2016749;6.898E-09;70;74;81.000;Sinorhizobium medicae WSM419  70      -1      71      0
Scaffold2506    5741    6367    CP003883.1;614723;615361;3.441E-116;427;639;74.800;Wolbachia endosymbiont of Drosophila simulans wNo    427,457 -1,-1   428,458 0,0
Scaffold26106   404     525     CP001391.1;738049;738171;1.515E-14;88;123;75.600;Wolbachia sp. wRi      88      -1      89      0
Scaffold2695    7513    9070    AP013028.1;967212;968194;8.170E-169;602;986;73.400;Wolbachia endosymbiont of Cimex lectularius  602,494,685,497,473,516,371,564,495,336,557,417,514,364,502,447,390,251,247,395,511,146,213,96,173,207,346,198,292,180,260,169,160      184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,79,79,-1,79,79,79,79,62,62,62,62,62 418,310,501,313,289,332,187,380,311,152,373,233,330,180,318,263,206,67,63,211,327,67,134,97,94,128,267,119,230,118,198,107,98   274,262,253,228,203,193,121,121,120,109,107,97,59,58,58,56,56,56,56,56,56,151,273,0,273,273,273,273,53,53,53,53,53
Scaffold3073    7572    8258    AM999887.1;361350;362058;4.978E-80;307;709;70.700;Wolbachia endosymbiont of Culex quinquefasciatus Pel  307,398,397,317,390,397,192,192,381,167 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1   308,399,398,318,391,398,193,193,382,168 0,0,0,0,0,0,0,0,0,0
Scaffold3478    6559    6991    CP001391.1;740953;741380;5.359E-66;260;434;73.500;Wolbachia sp. wRi     260     -1      261     0
Scaffold3478    9526    10039   CP001391.1;737722;738242;5.219E-101;376;521;76.100;Wolbachia sp. wRi    376     269     107     444
Scaffold3839    4008    4140    AM422018.1;681423;681520;1.601E-08;70;98;75.500;Candidatus Phytoplasma australiense     70,85,84,83,92,71,92,89,81,84,83,70,76,70,73,73,72,71,72,71,71  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1  71,86,85,84,93,72,93,90,82,85,84,71,77,71,74,74,73,72,73,72,72  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Scaffold3839    4283    4403    AP012292.1;2674970;2675076;1.196E-13;87;114;76.300;Selenomonas ruminantium subsp. lactilytica TAM6421   87,74,70,74,73,71,79,74,72,79,73,77,74,71,73,73,74,70   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1   88,75,71,75,74,72,80,75,73,80,74,78,75,72,74,74,75,71   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Scaffold506     33791   34276   CP001391.1;737675;738171;2.483E-94;354;497;75.800;Wolbachia sp. wRi     354     261     93      484
Scaffold506     34663   35069   CP001391.1;739581;739952;1.857E-48;202;408;71.000;Wolbachia sp. wRi     202     62      140     54
Scaffold506     36023   37255   CP001391.1;740170;741373;1.065E-231;810;1234;74.400;Wolbachia sp. wRi   810     386     424     594
Scaffold506     59622   60484   CP001391.1;740541;741380;9.171E-148;532;864;73.100;Wolbachia sp. wRi    532     187     345     375
Scaffold506     61928   62470   CP001391.1;737699;738242;4.000E-71;277;544;60.200;Wolbachia sp. wRi     277     213     64      472
Scaffold506     75276   75778   CP001391.1;737666;737893;1.790E-38;169;228;77.200;Wolbachia sp. wRi     169,358 114,274 55,84   192,501
Scaffold506     76649   78057   CP001391.1;740002;741389;2.540E-294;1018;1410;76.100;Wolbachia sp. wRi  1018    736     282     757
Scaffold6938    289     419     CT573213.2;6612046;6612174;3.847E-10;75;132;72.700;Frankia alni ACN14a  75      -1      76      0
Scaffold7395    1       402     CP001391.1;741660;742028;3.080E-120;440;403;83.300;Wolbachia sp. wRi    440     268     172     275
Scaffold8478    3272    3352    BA000021.3;254826;254904;7.165E-10;74;82;83.100;Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis 74,74   -1,-1   75,75   0,0
```

Command to change file:
```bash
sed -e "s/^/GAGA-0331_/g" LGTs.candidateloci.loose.bed > LGTs.candidateloci.GAGAid.bed
```

Output `LGTs.candidateloci.GAGAid.bed`:
```bash
GAGA-0331_Scaffold10016 10      116     CP002371.1;640122;640229;2.481E-09;72;108;76.300;Candidatus Liberibacter solanacearum CLso-ZC1  72      -1      73      0
GAGA-0331_Scaffold10911 548     642     CP003418.1;796953;797046;9.251E-12;80;94;78.700;Ignavibacterium album JCM 16511 80,94,70,96,75,95,75    -1,-1,-1,-1,-1,-1,-1    81,95,71,97,76,96,76    0,0,0,0,0,0,0
GAGA-0331_Scaffold11    9074    9257    AP012210.1;1446066;1446250;4.621E-09;71;185;37.400;Candidatus Arthromitus sp. SFB-rat-Yit       71,71   -1,-1   72,72   0,0
GAGA-0331_Scaffold1300  6482    8312    CP001391.1;741827;743353;0.000E+00;2185;1527;91.800;Wolbachia sp. wRi   2185,1176       2079,700        106,476 1517,812
GAGA-0331_Scaffold1300  10504   12833   CP001391.1;743798;745788;0.000E+00;2500;1994;87.700;Wolbachia sp. wRi   2500,827        2363,750        137,77  1993,804
GAGA-0331_Scaffold1463  18709   18952   CP001739.1;1490348;1490601;7.711E-13;84;254;66.900;Sebaldella termitidis ATCC 33386     84,91   -1,-1   85,92   0,0
GAGA-0331_Scaffold1463  19109   19246   CP003261.1;563097;563235;2.068E-10;76;139;71.900;Clostridium pasteurianum BC1   76      -1      77      0
GAGA-0331_Scaffold2142  3168    3557    CP003884.1;935415;935756;2.775E-73;285;342;78.700;Wolbachia endosymbiont of Drosophila simulans wHa     285,338 130,130 155,208 10,10
GAGA-0331_Scaffold21535 668     740     CP000738.1;2016679;2016749;6.898E-09;70;74;81.000;Sinorhizobium medicae WSM419  70      -1      71      0
GAGA-0331_Scaffold2506  5741    6367    CP003883.1;614723;615361;3.441E-116;427;639;74.800;Wolbachia endosymbiont of Drosophila simulans wNo    427,457 -1,-1   428,458 0,0
GAGA-0331_Scaffold26106 404     525     CP001391.1;738049;738171;1.515E-14;88;123;75.600;Wolbachia sp. wRi      88      -1      89      0
GAGA-0331_Scaffold2695  7513    9070    AP013028.1;967212;968194;8.170E-169;602;986;73.400;Wolbachia endosymbiont of Cimex lectularius  602,494,685,497,473,516,371,564,495,336,557,417,514,364,502,447,390,251,247,395,511,146,213,96,173,207,346,198,292,180,260,169,160      184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,79,79,-1,79,79,79,79,62,62,62,62,62 418,310,501,313,289,332,187,380,311,152,373,233,330,180,318,263,206,67,63,211,327,67,134,97,94,128,267,119,230,118,198,107,98   274,262,253,228,203,193,121,121,120,109,107,97,59,58,58,56,56,56,56,56,56,151,273,0,273,273,273,273,53,53,53,53,53
GAGA-0331_Scaffold3073  7572    8258    AM999887.1;361350;362058;4.978E-80;307;709;70.700;Wolbachia endosymbiont of Culex quinquefasciatus Pel  307,398,397,317,390,397,192,192,381,167 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1   308,399,398,318,391,398,193,193,382,168 0,0,0,0,0,0,0,0,0,0
GAGA-0331_Scaffold3478  6559    6991    CP001391.1;740953;741380;5.359E-66;260;434;73.500;Wolbachia sp. wRi     260     -1      261     0
GAGA-0331_Scaffold3478  9526    10039   CP001391.1;737722;738242;5.219E-101;376;521;76.100;Wolbachia sp. wRi    376     269     107     444
GAGA-0331_Scaffold3839  4008    4140    AM422018.1;681423;681520;1.601E-08;70;98;75.500;Candidatus Phytoplasma australiense     70,85,84,83,92,71,92,89,81,84,83,70,76,70,73,73,72,71,72,71,71  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1  71,86,85,84,93,72,93,90,82,85,84,71,77,71,74,74,73,72,73,72,72  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
GAGA-0331_Scaffold3839  4283    4403    AP012292.1;2674970;2675076;1.196E-13;87;114;76.300;Selenomonas ruminantium subsp. lactilytica TAM6421   87,74,70,74,73,71,79,74,72,79,73,77,74,71,73,73,74,70   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1   88,75,71,75,74,72,80,75,73,80,74,78,75,72,74,74,75,71   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
GAGA-0331_Scaffold506   33791   34276   CP001391.1;737675;738171;2.483E-94;354;497;75.800;Wolbachia sp. wRi     354     261     93      484
GAGA-0331_Scaffold506   34663   35069   CP001391.1;739581;739952;1.857E-48;202;408;71.000;Wolbachia sp. wRi     202     62      140     54
GAGA-0331_Scaffold506   36023   37255   CP001391.1;740170;741373;1.065E-231;810;1234;74.400;Wolbachia sp. wRi   810     386     424     594
GAGA-0331_Scaffold506   59622   60484   CP001391.1;740541;741380;9.171E-148;532;864;73.100;Wolbachia sp. wRi    532     187     345     375
GAGA-0331_Scaffold506   61928   62470   CP001391.1;737699;738242;4.000E-71;277;544;60.200;Wolbachia sp. wRi     277     213     64      472
GAGA-0331_Scaffold506   75276   75778   CP001391.1;737666;737893;1.790E-38;169;228;77.200;Wolbachia sp. wRi     169,358 114,274 55,84   192,501
GAGA-0331_Scaffold506   76649   78057   CP001391.1;740002;741389;2.540E-294;1018;1410;76.100;Wolbachia sp. wRi  1018    736     282     757
GAGA-0331_Scaffold6938  289     419     CT573213.2;6612046;6612174;3.847E-10;75;132;72.700;Frankia alni ACN14a  75      -1      76      0
GAGA-0331_Scaffold7395  1       402     CP001391.1;741660;742028;3.080E-120;440;403;83.300;Wolbachia sp. wRi    440     268     172     275
GAGA-0331_Scaffold8478  3272    3352    BA000021.3;254826;254904;7.165E-10;74;82;83.100;Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis 74,74   -1,-1   75,75   0,0
```

##### 3.1.3 Run `bedtools intersect` to find overlaps:
```bash
bedtools intersect -wa -a LGTs.candidateloci.GAGAid.bed -b blast_Lyzozyme_GAGA-0328.bed > results.blast.GAGA-0328.bed
```

Output for detected overlap in file `results.blast.GAGA-0328.bed`:
```bash
GAGA-0331_Scaffold2506	5741	6367	CP003883.1;614723;615361;3.441E-116;427;639;74.800;Wolbachia endosymbiont of Drosophila simulans wNo	427,457	-1,-1	428,458	0,0
```

##### 3.1.4 Extract Fasta sequence for resulting overlapping coordinates

```bash
bedtools getfasta -fi /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy/GAGA-0331 -bed results.blast.GAGA-0328.bed -fo results.blast.GAGA-0328.fa

# -fi: <input FASTA>
# -bed BED input file with coordinates
# -fo: Specify an output file name
```

Output file `results.blast.GAGA-0328.fa`:
```bash
>GAGA-0331_Scaffold2506:5741-6367
ATGTCAGAAAATGTAATAGTAGATCTCTCACATTGGAATAAGAATGTAGATTTTACATTAGCGAAAGAAGATGGAATATTGGGAATCATTCATAAAGCAACACAAGGATTAAAGTATGTGGATTCAACATAT
GCAGAAAGAAGAAAAGCTGCTGAAGATGAAGGTTTAATGTGGGGTGCATACCATTTTGGCATAGGGGAAAATGGTGAAGATCAAGCTGATCATTTCTTAAAAATGATTGGTGATAATACTAAAATATTACTT
GCCCTTGATATTGAAGAAAACCAAAATGGAAAAAATATGACATCGAAACAAGCAGAGAATTTCGTTAATAAAGTAACAGGACGTTTGCCTTTAATTTATGGAAATGCTTATTTTTTAAAAGATTTTGCTACA
TCAATTTTAACGAAATGTCCATTGTGGATAACAAGCTGGGAAGTACAACCAGTATTACCAAGAGGGTGGACAAAATGGATTTTATGGCAATATACTAATGGTAAAAGTGGAGAAAAGCCACATGAAGTACAA
GGAATAGGTCCATGTGATCGAAATAAATTCAATGGAACATTAGAAGAGTTAAAAGATTTTTGGATATCAAAATGTTATTAAAATTGTTTGAACAATTC
```

Now, `dfast` can be run with the extracted FASTA sequence.

#### 3.2 For all files:
##### 3.2.1 Change all bed files fitting to blast bed format

The pipeline produces two bed files: one file called `LGTs.candidateloci.loose.bed` which consists of candidates blasted against the prokaryotic database and another file `LGTs.nA.candidateloci.loose.bed` consisting of candidates without any ant genomes. We need to change both bed files and run the code against each file separately to identify additional candidates.

```bash
cd /global/scratch2/j_rink02/master/lgt/0_data

# Copy file to work on it
for i in */*; do cp $i/LGTs.candidateloci.loose.bed $i/LGTs.candidateloci2.loose.bed; done

for i in */*; do cp $i/LGTs.nA.candidateloci.loose.bed $i/LGTs.nA.candidateloci2.loose.bed; done

# Control that all directories have this file
ls */*/LGTs.candidateloci2.loose.bed
ls */*/LGTs.nA.candidateloci2.loose.bed

# Insert GAGAid in all files at beginning of each line
for i in GAGA-*; do echo $i ; cat $i/results/LGTs.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.candidateloci.GAGAid.bed; done

for i in GAGA-*; do echo $i ; cat $i/results/LGTs.nA.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.nA.candidateloci.GAGAid.bed; done
-------------------------------------------------------------------------

for i in OUT-*; do echo $i ; cat $i/results/LGTs.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.candidateloci.GAGAid.bed; done

for i in OUT-*; do echo $i ; cat $i/results/LGTs.nA.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.nA.candidateloci.GAGAid.bed; done
-------------------------------------------------------------------------

for i in NCBI-*; do echo $i ; cat $i/results/LGTs.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.candidateloci.GAGAid.bed; done

for i in NCBI-*; do echo $i ; cat $i/results/LGTs.nA.candidateloci2.loose.bed | awk -v ID=$i '{print (ID "_" $0)}' > $i/results/LGTs.nA.candidateloci.GAGAid.bed; done

# awk 'pattern{action}' file
# -v: introduces a variable
```

Example Output GAGA-0001 `LGTs.candidateloci.GAGAid.bed`:
```bash
GAGA-0001_Scaffold1     3010448 3010660 CP001131.1;4880087;4880286;1.601E-08;70;206;71.800;Anaeromyxobacter sp. K       70,79,70,84,77,88,88,70,70      -1,-1,-1,-1,-1,-1,-1,-1,-1      71,80,71,85,78,89,89,71,71      0,0,0,0,0,0,0,0,0
GAGA-0001_Scaffold10    3180417 3180639 CP000263.1;384231;384451;7.711E-13;84;224;72.100;Buchnera aphidicola BCc        84      -1      85      0
GAGA-0001_Scaffold11    5696020 5696198 CP000263.1;234106;234274;1.435E-12;83;180;71.800;Buchnera aphidicola BCc        83      -1      84      0
GAGA-0001_Scaffold114   737     1293    CP001359.1;4660106;4660668;1.996E-34;156;563;68.400;Anaeromyxobacter dehalogenans 2CP-1 156     105     51      322
GAGA-0001_Scaffold19    4418033 4418085 CP000084.1;141842;141895;4.621E-09;71;54;88.800;Candidatus Pelagibacter ubique HTCC1062 71      -1      72      0
```

Example Output NCBI-0001 `LGTs.candidateloci.GAGAid.bed`:
```bash
NCBI-0001_NC_039507.1   6928806 6929000 CP002659.1;1899187;1899369;5.162E-22;114;183;73.700;Sphaerochaeta coccoides DSM 17374   114,102,87,94,73,77,102,74,88,75,79,91,85,70,74,94,100,81,80,85,82,71,79,76,88,87,86,86,77,84,74,73,75,74,71,77,92,75,70,70,76,81,70,71,70,81,70,70,72  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1      115,103,88,95,74,78,103,75,89,76,80,92,86,71,75,95,101,82,81,86,83,72,80,77,89,88,87,87,78,85,75,74,76,75,72,78,93,76,71,71,77,82,71,72,71,82,71,71,73  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
NCBI-0001_NC_039507.1   6930093 6930292 AE010300.2;310857;311054;1.854E-14;89;201;71.400;Leptospira interrogans serovar Lai str. 56601  89      -1      90      0
NCBI-0001_NC_039507.1   11842783        11842965        CP003379.1;4871843;4871955;9.256E-12;80;114;77.300;Terriglobus roseus DSM 18391 80,84   -1,-1   81,85   0,0
NCBI-0001_NC_039518.1   1216867 1216941 CP002048.1;2119145;2119219;1.334E-09;73;75;81.300;Syntrophothermus lipocalidus DSM 12680        73,70,70,72,70  -1,-1,-1,-1,-1  74,71,71,73,71  0,0,0,0,0
NCBI-0001_NC_039518.1   8533069 8533318 CP000248.1;2062540;2062696;4.142E-13;85;157;72.700;Novosphingobium aromaticivorans DSM 12444    85,70,122,93,89,72,88,83,95,70,81       -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1        86,71,123,94,90,73,89,84,96,71,82       0,0,0,0,0,0,0,0,0,0,0
```

##### 3.2.2 Run `bedtools intersect` on all files

Input file: `blast_Lyzozyme_GAGA-0328.bed` (BLAST output)

```bash
# Run on LGTs.candidateloci.GAGAid.bed
for i in */results; do bedtools intersect -wa -a $i/LGTs.candidateloci.GAGAid.bed -b /global/scratch2/j_rink02/master/lgt/0_data/local_blast/blast_Lyzozyme_GAGA-0328.bed >> /global/scratch2/j_rink02/master/lgt/0_data/local_blast/results_Lyzozymes.bed; done

# Run on LGTs.nA.candidateloci.GAGAid.bed
for i in */results; do bedtools intersect -wa -a $i/LGTs.nA.candidateloci.GAGAid.bed -b /global/scratch2/j_rink02/master/lgt/0_data/local_blast/blast_Lyzozyme_GAGA-0328.bed >> /global/scratch2/j_rink02/master/lgt/0_data/local_blast/results_Lyzozymes_nA.bed; done

# -b: BAM/BED/GFF file to compare to
# -wa: write the original entry in A for each overlap
```

As an output file (`results_Lyzozymes.bed`), all intersected hits appear from the `LGTs.candidateloci.GAGAid.bed` and `LGTs.nA.candidateloci.GAGAid.bed` files from every directory.

Output file `results_Lyzozymes.bed`:
```bash
GAGA-0328_Scaffold36    505675  506652  CP003883.1;614723;615361;3.700E-119;437;639;75.400;Wolbachia endosymbiont of Drosophila simulans
GAGA-0331_Scaffold2506  5741    6367    CP003883.1;614723;615361;3.441E-116;427;639;74.800;Wolbachia endosymbiont of Drosophila simulans
GAGA-0378_Scaffold29    1056525 1057476 CP003883.1;614734;615371;9.234E-114;419;638;74.700;Wolbachia endosymbiont of Drosophila simulans
GAGA-0382_Scaffold1     29205332        29206276        AM999887.1;257710;257844;6.423E-14;88;135;76.600;Wolbachia endosymbiont of Culex
GAGA-0533_Scaffold25    1172728 1173697 CP003883.1;614753;615361;1.545E-83;319;624;71.200;Wolbachia endosymbiont of Drosophila simulans w
GAGA-0578_Scaffold3     7430181 7430859 CP001391.1;1355371;1356031;1.064E-136;495;661;76.900;Wolbachia sp. wRi  495,455 -1,-1   496,456 0
GAGA-0579_Scaffold2     762893  763897  CP003883.1;614723;615361;2.392E-101;378;639;73.600;Wolbachia endosymbiont of Drosophila simulans
```

Output file `results_Lyzozymes_nA.bed`:
```bash
GAGA-0328_Scaffold36    505675  506652  CP003883.1;614723;615361;3.700E-119;437;639;75.400;Wolbachia endosymbiont of Drosophila simulans wNo    437,
GAGA-0331_Scaffold2506  5741    6367    CP003883.1;614723;615361;3.441E-116;427;639;74.800;Wolbachia endosymbiont of Drosophila simulans wNo    427,
GAGA-0378_Scaffold29    1056525 1057476 CP003883.1;614734;615371;9.234E-114;419;638;74.700;Wolbachia endosymbiont of Drosophila simulans wNo    419,
GAGA-0382_Scaffold1     29205332        29206276        AM999887.1;257710;257844;6.423E-14;88;135;76.600;Wolbachia endosymbiont of Culex quinquefasc
GAGA-0533_Scaffold25    1172728 1173697 CP003883.1;614753;615361;1.545E-83;319;624;71.200;Wolbachia endosymbiont of Drosophila simulans wNo     319,
GAGA-0578_Scaffold3     7430181 7430859 CP001391.1;1355371;1356031;1.064E-136;495;661;76.900;Wolbachia sp. wRi  495,455 -1,-1   496,456 0,0
GAGA-0579_Scaffold2     762893  763897  CP003883.1;614723;615361;2.392E-101;378;639;73.600;Wolbachia endosymbiont of Drosophila simulans wNo    378,
```


#### For the Etherases:
Input file: `blast_Etherase_GAGAgenomes.bed` (BLAST output)

```bash
cd /global/scratch2/j_rink02/master/lgt/0_data

# Run on LGTs.candidateloci.GAGAid.bed
for i in */results; do bedtools intersect -wa -a $i/LGTs.candidateloci.GAGAid.bed -b /global/scratch2/j_rink02/master/lgt/0_data/local_blast/blast_Etherase_GAGAgenomes.bed >> /global/scratch2/j_rink02/master/lgt/0_data/local_blast/results_Etherase.bed; done

# Run on LGTs.nA.candidateloci.GAGAid.bed
for i in */results; do bedtools intersect -wa -a $i/LGTs.nA.candidateloci.GAGAid.bed -b /global/scratch2/j_rink02/master/lgt/0_data/local_blast/blast_Etherase_GAGAgenomes.bed >> /global/scratch2/j_rink02/master/lgt/0_data/local_blast/results_Etherase_nA.bed; done

# -b: BAM/BED/GFF file to compare to
# -wa: write the original entry in A for each overlap
```

Output file `results_Etherase.bed`:
```bash
GAGA-0221_Scaffold3    8820002 8820310 BX293980.2;283224;283430;4.979E-29;138;207;75.500;Mycoplasma mycoides subsp. mycoides SC str. PG1      138,147
GAGA-0361_Scaffold4    3723812 3724743 CP002082.1;345109;346044;9.206E-131;475;936;71.900;Spiroplasma mirum ATCC 29335 475,536,461,451,471,78,468,191,
GAGA-0362_Scaffold15   5725502 5725614 CP005077.1;787673;787786;4.142E-13;85;114;76.300;Spiroplasma chrysopicola DF-1 85     -1     86     0
```

Output file `results_Etherase_nA.bed`:
```bash
GAGA-0200_Scaffold111   363196  364682  CP002082.1;345738;346044;1.926E-41;179;307;73.400;Spiroplasma mirum ATCC 29335  179,474,532,209,213,145,568,
GAGA-0221_Scaffold3     8819312 8820310 CP002082.1;345109;346044;1.189E-132;482;936;72.000;Spiroplasma mirum ATCC 29335 482,537,464,451,483,82,474,2
GAGA-0360_Scaffold16    5909048 5910576 CP002082.1;345088;346035;2.759E-124;454;948;71.200;Spiroplasma mirum ATCC 29335 454,583,531,165,115,449,546,
GAGA-0361_Scaffold4     3723812 3724743 CP002082.1;345109;346044;9.206E-131;475;936;71.900;Spiroplasma mirum ATCC 29335 475,536,461,451,471,78,468,1
GAGA-0362_Scaffold15    5724696 5725632 CP002082.1;345109;346046;1.650E-137;498;938;72.300;Spiroplasma mirum ATCC 29335 498,551,486,469,487,87,487,2
GAGA-0374_Scaffold20    2858962 2859996 CP011856.1;383809;384320;3.573E-109;404;512;78.000;Spiroplasma eriocheiris      404,482,527,179,143,302,483,
GAGA-0396_Scaffold28    2026294 2027752 CP002082.1;345338;346044;8.915E-104;386;707;72.700;Spiroplasma mirum ATCC 29335 386,423,344,345,363,85,363,8
NCBI-0005_NW_020229769.1        103912  105255  CP001668.1;276469;277784;3.688E-136;493;1316;68.900;Mycoplasma mycoides subsp. capri str. GM12  493,
```

In this step, only three out of the originally eight detected hits appeared in the prokaryotic database with the ant genomes. This is due to the fact, that the other detected candidates have been identified within the `no.Ant database`. Here, the file `LGTs.nA.candidateloci.bed` was successfully intersected and found the other candidates.

##### 3.2.3 Extract FASTA sequence for resulting coordinates
```bash
cd /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy

# For the Lysozymes
# prokaryotic database candidates
for i in *; do bedtools getfasta -fi /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy/$i -bed /global/scratch2/j_rink02/master/lgt/0_data/local_blast/results_Lyzozymes.bed -fo /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates/$i.Lyzozyme.fa; done

# noAnt database candidates
for i in *; do bedtools getfasta -fi /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy/$i -bed /global/scratch2/j_rink02/master/lgt/0_data/local_blast/results_Lyzozymes_nA.bed -fo /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates/$i.Lyzozyme.nA.fa; done


# For the Etherases
#prokaryotic database candidates
for i in *; do bedtools getfasta -fi /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy/$i -bed /global/scratch2/j_rink02/master/lgt/0_data/local_blast/results_Etherase.bed -fo /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates/$i.Etherase.fa; done

#noAnt database candidates
for i in *; do bedtools getfasta -fi /global/scratch2/j_rink02/master/lgt/0_data/assemblies_copy/$i -bed /global/scratch2/j_rink02/master/lgt/0_data/local_blast/results_Etherase_nA.bed -fo /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates/$i.Etherase.nA.fa; done

# -fi: <input FASTA>
# -bed BED input file with coordinates
# -fo: Specify an output file name
```
`bedtools getfasta` searches for matching sequences according to the coordinates in the `results.bed` file in every genome. After that, all results can be found in the `dfast_candidates` folder separated by genome.

`cd /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates`

`ls`:
```bash
GAGA-0328.fa			 GAGA_genomes_database.03.nin.fa  GAGA_genomes_database.09.nhr.fa  NCBI-0005.fa
GAGA-0331.fa			 GAGA_genomes_database.03.nsq.fa  GAGA_genomes_database.09.nin.fa  NCBI-0006.fa
GAGA-0378.fa			 GAGA_genomes_database.04.nhr.fa  GAGA_genomes_database.09.nsq.fa  NCBI-0007.fa
GAGA-0382.fa			 GAGA_genomes_database.04.nin.fa  GAGA_genomes_database.10.nhr.fa  NCBI-0008.fa
GAGA-0533.fa			 GAGA_genomes_database.04.nsq.fa  GAGA_genomes_database.10.nin.fa  NCBI-0009.fa
GAGA-0578.fa			 GAGA_genomes_database.05.nhr.fa  GAGA_genomes_database.10.nsq.fa  NCBI-0010.fa
GAGA-0579.fa			 GAGA_genomes_database.05.nin.fa  GAGA_genomes_database.11.nhr.fa  NCBI-0011.fa
GAGA_genomes_database.00.nhr.fa  GAGA_genomes_database.05.nsq.fa  GAGA_genomes_database.11.nin.fa  NCBI-0012.fa
GAGA_genomes_database.00.nin.fa  GAGA_genomes_database.06.nhr.fa  GAGA_genomes_database.11.nsq.fa  NCBI-0013.fa
GAGA_genomes_database.00.nsq.fa  GAGA_genomes_database.06.nin.fa  GAGA_genomes_database.12.nhr.fa  NCBI-0014.fa
GAGA_genomes_database.01.nhr.fa  GAGA_genomes_database.06.nsq.fa  GAGA_genomes_database.12.nin.fa  NCBI-0015.fa
GAGA_genomes_database.01.nin.fa  GAGA_genomes_database.07.nhr.fa  GAGA_genomes_database.12.nsq.fa  NCBI-0016.fa
GAGA_genomes_database.01.nsq.fa  GAGA_genomes_database.07.nin.fa  GAGA_genomes_database.nal.fa	   NCBI-0017.fa
GAGA_genomes_database.02.nhr.fa  GAGA_genomes_database.07.nsq.fa  NCBI-0001.fa			   OUT-0001.fa
GAGA_genomes_database.02.nin.fa  GAGA_genomes_database.08.nhr.fa  NCBI-0002.fa			   OUT-0002.fa
GAGA_genomes_database.02.nsq.fa  GAGA_genomes_database.08.nin.fa  NCBI-0003.fa			   tmp
GAGA_genomes_database.03.nhr.fa  GAGA_genomes_database.08.nsq.fa  NCBI-0004.fa
```

However, the command creates a file for every single genome. To only use files with hits and extracted FASTA sequences for the dfast pipeline, all empty files with no matching results need to be deleted:

```bash
find . -type f -empty -delete
```

Files in the directory `dfast_candidates`after deletion of empty files:
```bash
GAGA-0200.Etherase.nA.fa  GAGA-0331.Lyzozyme.nA.fa  GAGA-0374.Etherase.nA.fa  GAGA-0533.Lyzozyme.fa	GAGA_genomes_database.Etherase.fa     tmp
GAGA-0221.Etherase.fa	  GAGA-0360.Etherase.nA.fa  GAGA-0378.Lyzozyme.fa     GAGA-0533.Lyzozyme.nA.fa	GAGA_genomes_database.Etherase.nA.fa
GAGA-0221.Etherase.nA.fa  GAGA-0361.Etherase.fa     GAGA-0378.Lyzozyme.nA.fa  GAGA-0578.Lyzozyme.fa	GAGA_genomes_database.Lyzozyme.fa
GAGA-0328.Lyzozyme.fa	  GAGA-0361.Etherase.nA.fa  GAGA-0382.Lyzozyme.fa     GAGA-0578.Lyzozyme.nA.fa	GAGA_genomes_database.Lyzozyme.nA.fa
GAGA-0328.Lyzozyme.nA.fa  GAGA-0362.Etherase.fa     GAGA-0382.Lyzozyme.nA.fa  GAGA-0579.Lyzozyme.fa	NCBI-0005.Etherase.nA.fa
GAGA-0331.Lyzozyme.fa	  GAGA-0362.Etherase.nA.fa  GAGA-0396.Etherase.nA.fa  GAGA-0579.Lyzozyme.nA.fa	run_dfast_BLAST_candidates.sh
```

Now, dfast can be run on all remaining files in this directory and all information can be obtained for those additional candidates identified with BLAST.


#### 3.3 Run dfast on additionally identified LGT candidates

##### 3.3.1 Write a dfast job script

`cd /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates`

nano `run_dfast_BLAST_candidates.sh`
```
#$ -S /bin/bash
#$ -N batchdfastjob
#$ -cwd
#$ -w e
#$ -V
#$ -pe smp 20
#$ -l h_vmem=6G
#$ -o /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates/tmp/batch.dfast.job.out
#$ -e global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates/tmp/batch.dfast.job.err
#$ -wd /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates

conda activate dfast

echo "Running on file: $file"

dfast -g $file --force --metagenome --cpu 20 --debug --use_original_name t --minimum_length 100 --database /global/scratch2/databases/dfast/uniprot_bacteria-0.9.ref -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/blast_output_dfast/$file --config /home/j/j_rink02/anaconda3/envs/dfast/bin/custom_config.py
```

Submit the script and parse bash variables to the script using `-v file="file1.fa"``

To submit jobs more easily, use parallel:
```
find . -name "GAGA*" | parallel -I% --max-args 1 qsub -v file="%" run_dfast_BLAST_candidates.sh

find . -name "NCBI*" | parallel -I% --max-args 1 qsub -v file="%" run_dfast_BLAST_candidates.sh

find . -name "OUT*" | parallel -I% --max-args 1 qsub -v file="%" run_dfast_BLAST_candidates.sh
```

##### 3.3.2 Build specific blast databases for dfast

The dfast jobs would take too long for all additionally identified candidates to run against, so we will build specific bacterial databases to significantly shorten the dfast search.

###### Build a database for all Enterobacteria to run dfast for the additional CFA synthases:
```bash
cd /global/scratch2/databases/dfast

# Extract relevant bacteria
 seqkit grep -n -r -p "Sodalis|Serratia|Yersinia|Rahnella|Escherichia" /global/scratch2/databases/dfast/uniprot_bacteria-0.9.fasta > ~/relevant.enterobacteria.CFA.uni90.fa
 ```

 Slightly change the dfast script and use less memory with the smaller database for Enterobacteria:

 `cd /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates`

 nano `run_dfast_CFA_candidates.sh`
 ```
 #$ -S /bin/bash
 #$ -N batchCFAsearch
 #$ -cwd
 #$ -w e
 #$ -V
 #$ -pe smp 15
 #$ -l h_vmem=2G
 #$ -o /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates/tmp/batch.dfast.job.out
 #$ -e global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates/tmp/batch.dfast.job.err
 #$ -wd /global/scratch2/j_rink02/master/lgt/0_data/local_blast/dfast_candidates

 conda activate dfast

 echo "Running on file: $file"

 dfast -g $file --force --metagenome --cpu 15 --debug --use_original_name t --minimum_length 100 --database ~/relevant.enterobacteria.CFA.uni90.fa -o /global/scratch2/j_rink02/master/lgt/2_analysis/gene_annotation/blast_output_dfast/$file --config /home/j/j_rink02/anaconda3/envs/dfast/bin/custom_config.py
 ```

 Submit the CFA jobs:
 ```
 find . -name "GAGA*" | parallel -I% --max-args 1 qsub -hostname=sebb09 -v file="%"  run_dfast_CFA_candidates.sh

 find . -name "NCBI*" | parallel -I% --max-args 1 qsub -hostname=sebb09 -v file="%" run_dfast_CFA_candidates.sh

 find . -name "OUT*" | parallel -I% --max-args 1 qsub -hostname=sebb09 -v file="%" run_dfast_CFA_candidates.sh
 ```
