# *Psidium guajava* evolutionary genomics (assembly → pangenomics)

End-to-end assembly and comparative genomics across *Psidium* and related wild species, including 11 *Psidium guajava*, spanning genome QC, assembly, pangenome graph inference.

Metadata in `samples.info`, with read info and input for [puzzler](https://github.com/merondun/puzzler) assembly in `samples.tsv`.

## Directory map

- [**01_qa_qc_genomescope/**](01_qa_qc_genomescope/) — read/QC and GenomeScope summaries for genome size/heterozygosity context.
- [**02_genome_assembly/**](02_genome_assembly/) — assembly generation and post-processing notes/scripts.
- [**03_purge_dups_sensitivity/**](03_purge_dups_sensitivity/) — assessing purge_dups params on tetraploid samples.
- [**04_whole_chromosome_phylogeny/**](04_whole_chromosome_phylogeny/) — infer chr1 phylogeny for sample ordering
- [**05_whole_genome_alignments/**](05_whole_genome_alignments/) — whole genome alignments using nucmer & dotplots
- [**06_pangenome/**](06_pangenome/) — *psidium guajava* pangenome.
- [**07_annotation/**](07_annotation/) — gene annotation including interproscan. 
- [**08_gene_alignments/**](08_gene_alignments/) — whole genome alignments using genes: jcvi.
- [**09_repeats/**](09_repeats/) — repeat annotation. 
- [**10_orthofinder/**](10_orthofinder/) — orthofinder orthogroup annotation and species tree. 
- [**11_selection_dnds/**](11_selection_dnds/) — dnds along the guava tree. 
- [**12_flesh_candidates/**](12_flesh_candidates/) — flesh color candidates using the pangenome.
- [**13_beast/**](13_beast/) — BEAST time dated tree. 


## Overview

This project generates collapsed haploid assemblies for 20 novel guava samples and wild relatives. 

Included are 11 guava (Psidium guajava) and 9 related species. One guava sample is triploid, and six relatives are tetraploid. 

| Plant or Cultivar Name | Accession | Species                     | ploidy | HiFi (Gb) | HiFi (Mean Length) | HiC (Gb) |
| ---------------------- | --------- | --------------------------- | ------ | --------- | ------------------ | -------- |
| ARAZA FRUIT            | HEUG_001  | Eugenia stipitata           | 2      | 34.3      | 13223.9            | 38.71    |
| PINEAPPLE GUAVA        | HPSI_069  | Feijoa sellowiana           | 2      | 33.13     | 14520.6            | 31.84    |
| AMAZON GUAVA           | HPSI_010  | Psidium acutangulum         | 4      | 35.16     | 17240.2            | 37.16    |
| COSTA RICAN GUAVA      | HPSI_059  | Psidium friedrichsthalianum | 4      | 50.52     | 16823.7            | 35.3     |
| PINK ACID              | HPSI_007  | Psidium guajava             | 2      | 28.21     | 16471.2            | 45.94    |
| PUERTO RICO 2          | HPSI_016  | Psidium guajava             | 2      | 21.9      | 17556.4            | 27.71    |
| KONA 1                 | HPSI_019  | Psidium guajava             | 2      | 22.8      | 18456.7            | 32.43    |
| J.B. WHITE             | HPSI_027  | Psidium guajava             | 2      | 38.46     | 18509.15           | 39.84    |
| KA HUA KULA            | HPSI_035  | Psidium guajava             | 2      | 18.94     | 16196.8            | 36.29    |
| BEAUMONT               | HPSI_037  | Psidium guajava             | 2      | 44.33     | 17738.3            | 12.82    |
| THAI MAROON            | HPSI_041  | Psidium guajava             | 2      | 30.09     | 17442.35           | 22.15    |
| KLOM TOONKLAO          | HPSI_060  | Psidium guajava             | 2      | 32.33     | 17538.7            | 37.89    |
| THAILAND               | HPSI_065  | Psidium guajava             | 2      | 18.03     | 5865.25            | 13.73    |
| GEMA DE DORO           | HPSI_068  | Psidium guajava             | 2      | 18.4      | 6405.2             | 21.84    |
| INDONESIAN SEEDLESS    | HPSI_003  | Psidium guajava             | 3      | 30.67     | 15131.8            | 32.95    |
| PUERTO RICAN GUAVA     | HPSI_072  | Psidium microphyllum        | 4      | 40.83     | 19275.7            | 29.43    |
| ARRAYAN FRUIT          | HPSI_080  | Psidium sartorianum         | 4      | 35.36     | 13464.2            | 34.97    |
| WATER APPLE            | HSYZ_002  | Syzygium aqueum             | 2      | 23.27     | 22708.5            | 15.81    |
| WAX APPLE              | HSYZ_001  | Syzygium samarangense       | 4      | 20.12     | 19098.6            | 1.74     |
| WAX APPLE              | HSYZ_003  | Syzygium samarangense       | 4      | 35.15     | 16641.5            | 6.9      |

Colors:

```
setwd('/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/')
library(tidyverse)

md <- read.table('samples.info',sep='\t',header = TRUE,comment.char = '') %>% as_tibble
cud_10 <- c(
  "#000000", # black
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#9AD0F3", # pale sky blue
  "#BCB0CC"  # pale lavender
)
spcols <- md %>% distinct(Species,Ploidy) %>% mutate(Color = rev(cud_10), Shape = c(24,24,25,25,25,25,23,22,21,21),Group = paste0(Species,' (',Ploidy,'n)'), Group = factor(Group,unique(Group)))
spcols %>% 
  ggplot(aes(y=Group,x=Ploidy,fill=Group,shape=Group))+
  geom_point(size=5)+
  scale_fill_manual(values=spcols$Color,breaks=spcols$Group)+
  scale_shape_manual(values=spcols$Shape,breaks=spcols$Group)+
  scale_x_continuous(breaks = seq(min(spcols$Ploidy), max(spcols$Ploidy), by = 1))+
  theme_bw()
```

Symbology across the paper:

| Group                            | Shape | Color   |
| -------------------------------- | ----- | ------- |
| Psidium guajava (3n)             | 21    | #BCB0CC |
| Psidium guajava (2n)             | 21    | #9AD0F3 |
| Psidium sartorianum (4n)         | 22    | #CC79A7 |
| Psidium friedrichsthalianum (4n) | 22    | #D55E00 |
| Psidium acutangulum (4n)         | 22    | #0072B2 |
| Psidium microphyllum (4n)        | 22    | #F0E442 |
| Feijoa sellowiana (2n)           | 23    | #009E73 |
| Eugenia stipitata (2n)           | 24    | #56B4E9 |
| Syzygium aqueum (2n)             | 25    | #E69F00 |
| Syzygium samarangense (4n)       | 25    | #000000 |

![symbology](/imgs/symbology_map.png)



## Qs & Cs

Questions or comments reach out to Justin Merondun heritabilities [@] gmail.com or make an issue here. 

