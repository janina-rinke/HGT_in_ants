## Re-Evaluation of candidates from the first seven genomes
#### Janina Rinke, 14.05.2021

### Re-evaluation of gc-content and evalue
Extraction of `best evalue` and `gc-content`out of the .tsv file from every genome.

```bash
# Get best evalue from the column bestProHit
cd /global/scratch2/j_rink02/master/lgt/0_data
cat filename.tsv | cut -f 5 |cut -f 4 -d ";" |sed 1d

# Explanation:
cut -f 5: only show the fifth column ("BestProHit"), columns are separated by tabs.
cut -f 4 -d ";" : only show fourth column, columns are separated by ";"
-> This only prints out the e-value
sed 1d: remove the first line from the output to remove the header.

# Get gc content for every predicted candidate based on the column name
head -1 filename.tsv | tr '\t' '\n' | cat -n | grep "gc"
cat filename.tsv | cut -f 19 -d |sed 1d

```

Re-evaluated number of good candidates in the first seven genomes


    Genome        nr good candidates    evalue range         gc range    BitDiffSum
    GAGA-0025    2 (piece of one LGT)   8e-16 - 3e-58         ~ 0.22      95 - 235
    GAGA-0063    2 (piece of one LGT)  2e-139 - 3e-160        ~ 0.28     334 - 551
    GAGA-0074            0
    GAGA-0082  3 (2 as piece of one LGT) 0.0e - 1e-95         ~ 0.3      120 - 7467
    GAGA-0083            0
    GAGA-0084    2 large LGTs (Wolbachia)0.0e - 1e-188        ~ 0.3       97 - 45510
    GAGA-0515            1                  5e-93             ~ 0.57        9881
