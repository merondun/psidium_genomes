# 13_beast/

Estimate a time-dated tree using crown fossil data from [Vasconcelos et al 2017 Myrteae calibration](https://doi.org/10.1016/j.ympev.2017.01.002). Subset genes based on 4-fold degenerate sites and a random 150k codons. 

Output:

![beast](/imgs/20260706_BEAST_Divergence_Dating.png)



Summaries from tracer:

| Run                | Posterior  ESS | Likelihood  ESS | Prior ESS | Tree height  mean | Tree height  95% HPD | Clock rate  mean | Clock rate  95% HPD    |
| ------------------ | -------------- | --------------- | --------- | ----------------- | -------------------- | ---------------- | ---------------------- |
| 150K  Codons       | 4669.3         | 4135.4          | 9001      | 69.6439           | [66.0721, 76.9754]   | 2.46E-04         | [2.2159E-4, 2.6103E-4] |
| 4-Fold  Degenerate | 4864.3         | 4194.6          | 8848.6    | 69.5691           | [66.1613, 76.6403]   | 8.93E-04         | [8.0792E-4, 9.3989E-4] |



___



This will use the 4-fold degenerate site fasta files output from `11_selection_dnds`. These 4-fold gene files have no sequences with > 20% gaps and all 18 samples present with a proper ATG start. They were aligned with codon-aware alignment using MACSE, filtered with codon-aware filtering using clipkit `clipkit --codon -m kpic` with 4-fold sites extracted with [merothon](https://github.com/merondun/merothon#extract-4fold-degenerate-positions-from-msa)  `extract_4fold` .

First we just need to combine the ~4k gene files: 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=20
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# concat and convert
seqkit concat beast_out/*fa > Guava_4Fold.fa 2> seqkit.concat.log
seqret -sequence Guava_4Fold.fa -outseq Guava_4Fold.nex -osformat nexus

# also creat a full genes file
find hyphy_out/ -type f -name "*.fa" -print0 |
while IFS= read -r -d '' f; do
    if [ "$(seqkit stat "$f" | awk 'NR==2{print $4}')" -eq 18 ]; then
        printf '%s\0' "$f"
    fi
done | xargs -0 seqkit concat -o Guava_FullGenes.fa
seqret -sequence Guava_FullGenes.fa -outseq Guava_FullGenes.nex -osformat nexus

# subset 150k codons
seqkit subseq -r 1:450000 Guava_FullGenes.fa > Guava_FullGenes_150k.fa
seqret -sequence Guava_FullGenes_150k.fa -outseq Guava_FullGenes_150k.nex -osformat nexus
```

Submit:

```
realpath raw/*fa > raw_genes.list
sbatch 03_Extract_4Fold.sh raw_genes.list
```

Count variant and invariant sites in all alignments:

```
#!/usr/bin/env python3
import sys
from itertools import zip_longest

def read_fasta(path):
    seqs = []
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    seqs.append("".join(seq))
                    seq = []
            else:
                seq.append(line)
        if seq:
            seqs.append("".join(seq))
    return seqs

def count_sites(seqs):
    invariant = 0
    variant = 0
    aln_len = len(seqs[0])
    for col in zip(*seqs):
        bases = set(col) - set("-")   # ignore gaps
        if len(bases) == 1:
            invariant += 1
        else:
            variant += 1
    return invariant, variant

if __name__ == "__main__":
    for fa in sys.argv[1:]:
        seqs = read_fasta(fa)
        inv, var = count_sites(seqs)
        print(fa)
        print("Invariant sites:", inv)
        print("Variant sites:", var)
        print()
```

Count for fastas:

```bash
Guava_4Fold.fa
Invariant sites: 702731
Variant sites: 106669

Guava_FullGenes_150k.fa
Invariant sites: 427901
Variant sites: 22099
```



Merge the 4-fold degenerate fasta files and then import them into beauti:

* Gamma model, 4 categories, estimated shape, GTR with estimated frequencies
* Strict clock, log normal default prior
* Yule model, tMRCA prior based on [Vasconcelos et al 2017 Myrteae calibration](https://doi.org/10.1016/j.ympev.2017.01.002) paper: crown split between Szygium / rest: 66 Ma. Give it some wiggle room: log normal prior with Offset = 66, Mean = 1.0, Sigma = 1.0, uncheck mean in real space. 
* Under trees, make sure branch lengths are in substitutions/site. 
* 50M chains, log every 5k 

From the paper: 

> Approach A, considered three fossil records: *Myrceugeneloxylon antarticus*, the oldest fossil in Myrteae, was placed on the crown node of Myrteae calibrating it at 66 million years ago (mya). 

ALSO run a BEAST analysis with 150k codons from the full gene alignments for sensitivity. 

### Plot BEAST

```R
#### Plot BEAST annotated trees 
setwd('/project/coffea_pangenome/Guava/selection_dnds/beast')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

# metadata
md <-  read_tsv('~/psidium_genomes/samples.info') %>% mutate(Sample = gsub('_','',Sample))

files = list.files('.',paste0('.*ann'))

counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file) 
  t2 <- drop.tip(iqtree,'HSYZ002')
  gg = ggtree(t2,layout='rectangular') %<+% md
  
  #add label for 95% CIs
  lab = gsub('.ann','',file)
  heights = gg$data$height_0.95_HPD
  df = as.data.frame(do.call(rbind, heights)) #convert the list to a data frame
  df$node_value = 1:nrow(df) # Add node values as a new column
  colnames(df) = c("value1", "value2", "node")
  df = df[, c("node", "value1", "value2")]
  df = df %>% 
    mutate(
      value1 = if (grepl("mu", file)) value1 / 1e6 else value1,
      value2 = if (grepl("mu", file)) value2 / 1e6 else value2,
      lab = paste0(round(value1,1),' - ',round(value2,1))) %>% 
    select(!c(value1,value2))
  
  leg = md %>% select(Group,Color,Shape) %>% unique
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=2) +
    geom_tippoint(aes(fill = Group, shape = Group), size=2)+
    geom_tiplab(aes(label=Cultivar),size=1.5,hjust=-0.5)+
    geom_nodelab(aes(label=lab),size=3,vjust=1) +
    scale_fill_manual(values=md$Color,breaks=md$Group)+
    scale_shape_manual(values=md$Shape,breaks=md$Group)+    
    ggtitle(lab)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=5) +
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 3)+
    xlim(0,max(gg$data$x)*1.5)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')
  ggp
  assign(paste0('p',counter),ggp)
} 

ggarrange(p1,p2,common.legend = TRUE,nrow=1)

pdf('~/symlinks/guava/figures/20260706_BEAST_Divergence_Dating.pdf',height=4,width=7)
ggarrange(p1,p2,common.legend = TRUE,nrow=1)
dev.off()

```

