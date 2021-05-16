## Re-Evaluation of candidates from the first seven genomes
#### Janina Rinke, 14.05.2021

### Re-evaluation of gc-content and evalue
Extraction of `best evalue` and `gc-content`out of the .tsv file from every genome.

```bash
# Get best eval from the column bestProHit
cd /global/scratch2/j_rink02/master/lgt/0_data
egrep -o '[0-9].[0-9].[0-9]E-[0-9]*' filename.tsv

# Get gc content for every predicted candidate based on the column name
head -1 filename.tsv | tr '\t' '\n' | cat -n | grep "gc"

cut -f nr_of_gc_column filename.tsv > newfile_gc.tsv
```

Re-evaluated number of good candidates in the first seven genomes


    Genome        nr good candidates    evalue range         gc range
    GAGA-0025    2 (piece of one LGT)   8e-16 - 3e-58         ~ 0.22
    GAGA-0063    2 (piece of one LGT)  2e-139 - 3e-160        ~ 0.28
    GAGA-0074            0
    GAGA-0082            2              1e-52 - 1e-95         ~ 0.3
    GAGA-0083            0
    GAGA-0084    2 large LGTs (Wolbachia)1e-25 - 1e-188       ~ 0.3
    GAGA-0515            1                  5e-93             ~ 0.57
