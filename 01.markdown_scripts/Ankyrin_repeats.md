# Analysis of Ankyrin Repeat Proteins in GAGA ants

### Plotting the distribution of ankyrin repeats across the GAGA phylogeny
HGT candidates were checked for ankyrin repeat domains at the following website:
```bash
http://smart.embl-heidelberg.de/
```

In case of detection of ANK domains, these proteins were incorporated in the ANK repeat analysis.

After that, all proteins [CDS] were extracted and plotted onto the GAGA phylogeny with the following R-Script:
```bash
./02.r_scripts/11_ANKs_analysis.Rmd
```