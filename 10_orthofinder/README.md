# 10_orthofinder/

To investigate shared genes and generate species tree, using orthofinder:

```bash
while read SAMPLE; do \
  infile="outputs/${SAMPLE}/complete.proteins.faa"; \
  outfile="/project/coffea_pangenome/Guava/orthofinder/proteins/${SAMPLE}.pep.fa"; \
  awk 'BEGIN{c=1} /^>/ {printf(">egapxtmp_%06d-P1\n", c++); next} {print}' "$infile" > "$outfile"; \
  echo "Processed $SAMPLE → $outfile"; \
done < Samples.list
```

This is pretty much just to get the initial species tree among our samples. 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=20
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate orthofinder

WD=/project/coffea_pangenome/Guava/orthofinder

t=20
RUN=$1
cd ${WD} 

orthofinder -f ${RUN} -a ${t}

# Afterwards, clean up for long term storage since this generates thousdands of files 
#find proteomes/${REF} -type d \( -name "MultipleSequenceAlignments" -o -name "Gene_Trees" -o -name "WorkingDirectory" -o -name "Orthogroup_Sequences" -o -name "Orthologues" -o -name "Resolved_Gene_Trees" \) -exec rm -rf {} +
```

Plot tree and output nwk:

```R
library(ape)
library(ggtree)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(tidyverse)

# Read in and prune ficus 
t <- read.tree('/project/coffea_pangenome/Guava/orthofinder/proteins/OrthoFinder/Results_Jun29/Species_Tree/SpeciesTree_rooted_node_labels.txt')
plot(t)
t2 <- drop.tip(t,c('HSYZ_001.pep','HSYZ_003.pep'))
plot(t2)
is.rooted(t2)
is.binary(t2)
is.ultrametric(t2)
t2$tip.label <- gsub('.pep','',t2$tip.label)
xt <- 0.1
t2$node.label <- LETTERS[1:17]; ggtree(t2)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label))

# add md 
md <- read.table('~/psidium_genomes/samples.info',sep='\t',header = TRUE,comment.char = '')

tp <- ggtree(t2, layout = "rectangular")  %<+% md
sp_tree <- tp +
  geom_tiplab(aes(label=Cultivar),hjust = -0.1,size=2)+
  geom_tippoint(aes(fill = Group, shape = Group), size=2)+
  scale_fill_manual(values=md$Color,breaks=md$Group)+
  scale_shape_manual(values=md$Shape,breaks=md$Group)+
  xlim(0,max(tp$data$x)*1.3)+
  theme(legend.text = element_text(size = 5),legend.title = element_text(size = 6),
        legend.key.size = unit(0.03, "cm"),    legend.position = 'top')
sp_tree

ggsave('~/symlinks/guava/figures/20260630_species-tree-orthofinder.pdf',sp_tree,height=4,width=5.5)
write.tree(t2,'~/symlinks/guava/selection_dnds/tree.nwk')

```

