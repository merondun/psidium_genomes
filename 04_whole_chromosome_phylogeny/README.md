# 04_whole_chromosome_phylogeny/

Create quick whole-chromosome phylogeny for sample ordering.

Outputs:

![chr_phylogeny](/imgs/20260521_chr1_tree.png)



___



Infer distance matrix, convert to phylip, and infer a tree:

```bash
#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

set -euo pipefail
# source activate pggb
CHR=${1:?Provide chr as first argument}
GENOMES=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm

mkdir -p raw mash_work trees

for SAMPLE in $(cat ../Samples.list); do 
    echo "Extracting ${SAMPLE}"
    samtools faidx ${GENOMES}/${SAMPLE}.fa ${CHR} > raw/${SAMPLE}.fa
done 

# sketch
mash sketch -o mash_work/${CHR} -p 20 raw/*

# distances
mash dist mash_work/${CHR}.msh mash_work/${CHR}.msh > mash_work/${CHR}.mash.dist

# prep phylip
input=mash_work/${CHR}.mash.dist
output=mash_work/${CHR}.mash.phy

# Temporary sorted sample list
samples=$(awk '{print $1; print $2}' "$input" | sort -u)

# Count samples
ntaxa=$(echo "$samples" | wc -l)

# Write number of taxa to output
echo "$ntaxa" > "$output"

# Create associative array for distances
declare -A dist

# Read mash output and store distances
while read -r a b d p shared; do
    dist["$a,$b"]="$d"
    dist["$b,$a"]="$d"
done < "$input"

# Build PHYLIP matrix
for i in $samples; do
    line="$i"
    for j in $samples; do
        if [[ "$i" == "$j" ]]; then
            line="$line 0"
        else
            val="${dist[$i,$j]}"
            if [[ -z "$val" ]]; then
                val="0"   # or NA, but FastME prefers 0
            fi
            line="$line $val"
        fi
    done
    echo "$line" >> "$output"
done

# tree
fastme -T 20 -i mash_work/${CHR}.mash.phy -o trees/${CHR}.nwk

```

Plot:

```
#### Plot whole chr tree
setwd('/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/whole_chromosome_phylogeny')
library(ggtree)
library(phytools)

t <- read.tree('trees/chr1.nwk')
t$tip.label <- gsub('raw/','',t$tip.label)
t$tip.label <- gsub('.fa','',t$tip.label)
md <- read_tsv('~/symlinks/assgv/samples.info') %>% dplyr::rename(label=Accession)
g <- ggtree(t) %<+% md

gp <- g + 
  geom_tiplab(size=3,aes(label=paste0(label,' (',Species,')')))+
  xlim(min(g@data$x),max(g@data$x)*1.4)

ggsave('~/symlinks/assgv/figures/20260521_chr1_tree.pdf',gp,height=4,width=5)

```

