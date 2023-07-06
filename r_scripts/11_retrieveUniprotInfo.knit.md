---
title: "Retrieve information for all LGTs from UniProt"
author: "Janina Rinke"
output: 
  pdf_document: default
  html_document:
    toc: true
    toc_float: true
    toc_depth: 2
    collapsed: true
    number_sections: true
    df_print: paged
---
<STYLE TYPE="text/css">
<!--
  td{
    font-family: Arial;
    font-size: 8pt;
    padding:0px;
    cellpadding="0";
    cellspacing="0"
  }
  th {
    font-family: Arial;
    font-size: 8pt;
    height: 20px;
    font-weight: bold;
    text-align: right;
    background-color: #ccccff;
  }
  table {
    border-spacing: 0px;
    border-collapse: collapse;
  }
--->
</STYLE>
---




```r
library(vroom)
#BiocManager::install("GenomicAlignments")
#install.packages("UniprotR")
library(UniprotR)
library(ggplot2)
library(dplyr)
library(stringr)
library(utils)
library(tidyr)
library(viridis)
```

## Load your data

```r
library(readxl)
df<-read_xlsx("/Users/Janina/Documents/MASTER/Masterarbeit/lgt/0_data/LGTs.CDSLevel.Expression.xlsx")
```





























