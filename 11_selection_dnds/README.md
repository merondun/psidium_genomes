# 11_selection_dnds

Branch specific dnds.

* Extracts longest‑isoform CDS for all genomes, identifies HSYZ_002 reciprocal‑best‑hit orthologs, builds per‑gene codon alignments, concatenates RBH families, and runs HyPhy (MG94) to estimate branch‑specific dN/dS across the tree.

* Cleans and filters per‑branch estimates, aggregates selection statistics, visualizes branch‑level positive selection on the species tree, and functionally annotates top candidate genes under adaptive evolution.

Output:

![dnds](/imgs/20260706_dNdS_BranchColorsPOSITIVE_zoom.png)



___



Using the same script as above to grab the proteins, we instead grab the gene CDS:

```bash
cd /project/coffea_pangenome/Guava/annotation/only_longest_transcript_per_gene

for i in $(ls *cds.fa | sed 's/.cds.fa//g'); do sed 's/ .*//g' ${i}.cds.fa > ../../selection_dnds/cds/${i}.cds.fa; done 
```

## Identify RBH with HSYZ_002 

Run RBH:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=16
#SBATCH --mem=32Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate isoseq_ann

JOBS=16
WD=/project/coffea_pangenome/Guava/selection_dnds
cd ${WD}

# Submit with sample 
B="${1:?usage: $0 <SAMPLE>}"
first=0
[ "$B" = "HEUG_001" ] && first=1

# dNdS Across tree 
A="HSYZ_002"
Af="cds/${A}.cds.fa"

# output dirs
mkdir -p blast db genes

echo "Blasting ${A} against ${B}"
Bf="cds/${B}.cds.fa"

# Make blast dbs
if [ ! -f "db/${B}.nhr" ]; then
  makeblastdb -in $Bf -dbtype nucl -out db/${B}
fi

if [ ! -f "db/${A}.nhr" ]; then
  makeblastdb -in $Af -dbtype nucl -out db/${A}
fi

# blast each against the other
blastn -query $Af -db db/${B} -out blast/${A}_vs_${B}.tsv -outfmt "6 qseqid sseqid pident length evalue bitscore" -max_target_seqs 1 -evalue 1e-5
blastn -query $Bf -db db/${A} -out blast/${B}_vs_${A}.tsv -outfmt "6 qseqid sseqid pident length evalue bitscore" -max_target_seqs 1 -evalue 1e-5

# sort for best hits 
sort -k1,1 -k6,6nr blast/${A}_vs_${B}.tsv \
  | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' > blast/best_${A}_to_${B}.txt
sort -k1,1 -k6,6nr blast/${B}_vs_${A}.tsv \
  | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' > blast/best_${B}_to_${A}.txt
awk 'NR==FNR {a[$1]=$2; next} {if (a[$2]==$1) print $2"\t"$1}' \
  blast/best_${A}_to_${B}.txt \
  blast/best_${B}_to_${A}.txt \
  > blast/RBH_${A}_${B}.txt

# Export for parallel subshells
export WD Af Bf A B first

# SELF per ida, only if first==1; run once per unique ida in parallel ---
if [ "${first:-0}" -eq 1 ]; then
  cut -f1 "blast/RBH_${A}_${B}.txt" | sort -u | \
  parallel --jobs ${JOBS} '
    ida={}
    mkdir -p ${WD}/genes/${ida}
    self="${WD}/genes/${ida}/${A}.fa"
    # overwrite or create; trimming terminal stop codon
    samtools faidx "'"$Af"'" "$ida" | sed -E "s/(TAA|TAG|TGA)$//" > "$self"
  '
fi

# Extract each idb from Bf in parallel, one file per pair ---
parallel --jobs ${JOBS} --colsep '\t' '
  ida={1}; idb={2}
  mkdir -p ${WD}/genes/${ida}
  out="${WD}/genes/${ida}/${B}.fa"
  samtools faidx "'"$Bf"'" "$idb" | sed -E "s/(TAA|TAG|TGA)$//" > "$out"
' :::: "blast/RBH_${A}_${B}.txt"
```

Run:

```
for i in $(egrep -v 'HSYZ|HEUG' Samples.list); do sbatch -J blast_${i} 01_alignments.sh ${i}; done 
```

Sanity, check number of genes per sample:

```bash
find genes -maxdepth 2 -type f -printf "%f\n" \
    | grep -o -f <(sed 's/$/\\.fa/' Samples.list) \
    | sort | uniq -c
    
  18550 HEUG_001.fa
  19155 HPSI_003.fa
  19111 HPSI_007.fa
  18963 HPSI_010.fa
  19173 HPSI_016.fa
  19158 HPSI_019.fa
  19163 HPSI_027.fa
  19159 HPSI_035.fa
  19153 HPSI_037.fa
  19192 HPSI_041.fa
  19154 HPSI_059.fa
  19139 HPSI_060.fa
  19160 HPSI_065.fa
  19167 HPSI_068.fa
  19471 HPSI_069.fa
  18850 HPSI_072.fa
  19181 HPSI_080.fa
  18550 HSYZ_002.fa
```

## Concatenate RBH Alignments

Concatenate the alignents:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

set -euo pipefail

WD=/project/coffea_pangenome/Guava/selection_dnds
SRC="${WD}/genes"
OUT="${WD}/raw"
JOBS=48
rm -rf raw

mkdir -p "$OUT"

export SRC OUT

find "$SRC" -mindepth 1 -maxdepth 1 -type d | \
parallel -j ${JOBS} '
    d="{1}"
    gene=$(basename "$d")
    out="${OUT}/${gene}.fa"
    tmp=$(mktemp "${OUT}/.${gene}.XXXXXX")

    # require HSYZ_002.fa to exist and be non-empty
    [ -s "${d}/HSYZ_002.fa" ] || { rm -f "$tmp"; exit 0; }

    # concatenate all fasta files with corrected headers
    find "$d" -maxdepth 1 -type f -name "*.fa" -size +0c -printf "%f\n" \
        | LC_ALL=C sort \
        | while read -r f; do
            sample="${f%.fa}"
            # rewrite header and append
            sed "s/^>.*/>${sample}/" "${d}/${f}" >> "$tmp"
        done

    # only keep if something was written
    if [ -s "$tmp" ]; then
        mv "$tmp" "$out"
        echo "Built ${gene}.fa"
    else
        rm -f "$tmp"
    fi
' ::: $(find "$SRC" -mindepth 1 -maxdepth 1 -type d)
```

## Run HyPhy dNdS

Subset those tips from the tree, and run hyphy for dnds estimation:

```bash
find raw -type f -print0 | xargs -0 realpath > genes_all.list
split -n l/8 -d genes_all.list genes_all.list_
ls *list_* | xargs -I {} echo sbatch -J dnds_{} 03_nodespecific_dnds.sh {}
ls *list_* | grep -v '00' | xargs -I {} echo sbatch -J dnds_{} 03_nodespecific_dnds.sh {}
```

Run this:

```bash
#!/bin/bash

#SBATCH --time=3-00:00:00    
#SBATCH --cpus-per-task=48
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann

JOBS=44
LIST="${1:?usage: $0 <genes_list_file>}"
WD=/project/coffea_pangenome/Guava/selection_dnds
TREE="${WD}/tree.nwk"
GENEDIR="${WD}/raw"

mkdir -p ${WD}/hyphy ${WD}/hyphy_out ${WD}/beast_out

cd ${WD}

echo "Working on ${LIST}"

# Ensure required tools exist
for tool in macse pal2nal.pl hyphy jq parallel; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

SCR=${SLURM_TMPDIR:-/tmp}
export WD TREE RUN GENEDIR SCR

process_gene() {
  local fa="$1"           
  local gene
  gene="$(basename "${fa%.fa}")" 

  # Work in a per-gene dir
  (
    outdir="${SCR}/hyphy/${gene}"
    mkdir -p "$outdir"
    cd "$outdir"
    echo "Working on ${gene}"

    # 1) Align with MACSE
    macse -prog alignSequences \
      -seq "${fa}" \
      -out_NT aln_NT.fasta \
      -out_AA aln_AA.fasta > macse.log 2>&1

    # 2) Clean alignment
    macse -prog exportAlignment \
      -align aln_NT.fasta \
      -codonForExternalFS --- \
      -codonForFinalStop --- \
      -codonForInternalFS --- \
      -codonForInternalStop --- \
      -charForRemainingFS - \
      -out_NT aln_NT.clean.fasta 2>&1

    # 3) Normalization, ensure underscores stripped
    sed 's/_//g' aln_NT.clean.fasta > aln_NT.clean.fa

    # 4) Trim gappy codon positions
    clipkit aln_NT.clean.fa -o aln_NT.clipkit.fa --codon -m kpic > clipkit.log 2>&1

    # 4b) Create a "strict" codon alignment that drops sequences with >20% gaps
    seqkit fx2tab aln_NT.clipkit.fa \
    | awk -F'\t' '
        {
          seq = $2
          gsub(/[^-]/, "", gapped)
        }
        {
          total = length($2)
          gaps = gsub(/-/, "", $2)
          if (total > 0 && gaps/total <= 0.20)
              print $1
        }
    ' | seqkit grep -f - aln_NT.clipkit.fa > aln_NT.clipkit.strict.fa
    
    # 5) Extract present taxa from fasta headers
    grep '^>' aln_NT.clipkit.fa | sed 's/^>//' > tips.txt

    ntips=$(wc -l < tips.txt)
    if [ "$ntips" -lt 15 ]; then
      echo "Skipping ${gene}: fewer than 15 taxa after filtering" > skip.txt
      cp skip.txt "${WD}/hyphy_out/${gene}.skip.txt"
      exit 0
    fi

    # For BEAST, Require 18 samples with ATG aligned and extract 
    python /home/justin.merondun/merothon/merothon/scripts/Extract_4Fold.py -i aln_NT.clipkit.strict.fa -m 18 -o extract
    cp extract.4fold.fa "${WD}/beast_out/${gene}.fa"

    # 6) Prune master tree to present taxa
    cp "${TREE}" tree.full.nwk
    nw_prune -v tree.full.nwk $(tr '\n' ' ' < tips.txt) > tree.pruned.nwk 2> prune.log

    # make sure pruning succeeded and tree is non-empty
    if [ ! -s tree.pruned.nwk ]; then
      echo "Skipping ${gene}: pruned tree is empty" > skip.txt
      cp skip.txt "${WD}/hyphy_out/${gene}.skip.txt"
      exit 0
    fi

    # 7) Run HyPhy on pruned tree
    hyphy "${WD}/FitMG94.bf" \
      --alignment extract.trimmed.fa \
      --tree tree.pruned.nwk \
      --output fitmg94.json \
      --type local \
      --lrt Yes > hyphy.log 2>&1

    # 8) Extract branch results
    jq -r '
      .input["number of sites"] as $L
      | ["Branch","dN","dS","N_exp","S_exp","LB","MLE","UB","LRT","FDR"],
        (.["branch attributes"]["0"]
          | to_entries
          | map([
              .key,
              (.value["dN"] // "NA"),
              (.value["dS"] // "NA"),
              ((.value["dN"] // 0) * $L),
              ((.value["dS"] // 0) * $L),
              (.value["Confidence Intervals"]["LB"] // "NA"),
              (.value["Confidence Intervals"]["MLE"] // "NA"),
              (.value["Confidence Intervals"]["UB"] // "NA"),
              (try .value["LRT"]["LRT"] catch "NA"),
              (try .value["LRT"]["FDR"] catch "NA")
            ])
        )[] | @tsv
    ' fitmg94.json > fitmg94.tsv

    # 9) Append gene 
    awk -F'\t' -v OFS='\t' -v Gene="$gene" 'NR==1 { print $0, "Gene"; next }
            { print $0, Gene }
            ' fitmg94.tsv > merged.tsv
            
    cp merged.tsv "${WD}/hyphy_out/${gene}.tsv"
    cp tree.pruned.nwk "${WD}/hyphy_out/${gene}.tree.nwk"
    cp extract.trimmed.fa "${WD}/hyphy_out/${gene}.fa"
    
  )
}

export -f process_gene

# Run in parallel
parallel --will-cite -j ${JOBS} process_gene :::: "${LIST}"

```

Check after:

```
find hyphy_out/ -type f -name '*tsv' | wc -l
find hyphy_out/ -type f -name '*skip*' | wc -l
```

Compile:

```bash
# Compile results From the directory containing hyphy_out/
awk -v OFS='\t' '
  FNR==1 {
    # Extract gene from filename (strip directory and .tsv)
    gene = FILENAME
    sub(/.*\//, "", gene)
    sub(/\.tsv$/, "", gene)
  }

  # Print header once, from the first file, and add Gene column
  FNR==1 && NR==1 { print $0, "Gene"; next }

  # Skip headers in subsequent files
  FNR==1 { next }

  # Print data line + gene ID
  NF { print $0, gene }
' hyphy_out/*.tsv > Node_dNdS_20260701.tsv
```

## Summarize Selection & Plot Tree

Visualization of dnds along the guava tree. 

Outputs

- Filtered results table: `Node_dNdS_20260701-RInput.tsv.gz`
- Phylogenetic tree with selection strength

```R
# dNdS along the tree 
setwd('/project/coffea_pangenome/Guava/selection_dnds')
library(tidyverse)
library(scales)
library(stringr)
library(data.table)
library(ggrepel)
library(ape)
library(ggtree)
library(ggpubr)
library(viridis)
library(meRo)

# read in dnds
dnds <- read_tsv('Node_dNdS_20260701.tsv') %>% select(!Gene...11) %>% dplyr::rename(Gene = Gene...12)

# Read in trees
tr <- read.tree('tree.nwk')
ggtree(tr) + geom_nodelab() +geom_tiplab()

# Apply some sanity thresholds like # syn mutations 
omega_cap <- 5    # ceiling for dnds
min_syn <- 1; max_syn <- 5000     # require at least 1 syn mutation for stability 
max_non <- 5000
max_ub <- 100 # maximum omega upper bound, above indicates model instability - these will be REMOVED
wstats_clean <- dnds %>%
  # Remove extrmee 
  filter(S_exp >= min_syn & S_exp < max_syn & N_exp < max_non & UB < max_ub) %>% 
  mutate(omega = pmin(MLE, omega_cap),
         pos = FDR < 0.10 & MLE > 1,
         puri = FDR < 0.10 & MLE < 1) 
hist(wstats_clean$omega)
fwrite(wstats_clean,'Node_dNdS_20260701-RInput.tsv.gz')
wstats_clean <- fread('Node_dNdS_20260701-RInput.tsv.gz') %>% as_tibble

# First, aggregate: compare # purified genes to total for each branch 
# Build per-ID counts
counts <- wstats_clean %>%
  mutate(
    pos  = as.integer(FDR < 0.10 & MLE > 1),
    puri = as.integer(FDR < 0.10 & MLE < 1)
  ) %>%
  group_by(Branch) %>%
  summarise(pos = sum(pos), puri = sum(puri), total = n(), .groups = "drop")
counts

##### Plot per-branch variation #####
# per-branch summaries
summarize_branch_metrics <- function(df) {
  df %>%
    mutate(
      ci_w = pmax(0, UB - LB),
      pos = FDR < 0.10 & MLE > 1,
      puri = FDR < 0.10 & MLE < 1,
    ) %>%
    group_by(Branch) %>%
    summarise(
      med_omega = median(omega, na.rm = TRUE),
      pos_frac  = mean(pos,  na.rm = TRUE),
      total_pos = sum(pos, na.rm = TRUE),
      puri_frac = mean(puri, na.rm = TRUE),
      total_puri = sum(puri, na.rm = TRUE),
      ci_w_med  = median(ci_w, na.rm = TRUE),
      n_genes   = n(),
      .groups   = "drop"
    ) 
}

# apply 
branch_summ <- summarize_branch_metrics(wstats_clean)
branch_summ 

##### Plot dNdS on the tree: lwd and color PURIFYING #####
# import md
md <-  read_tsv('~/psidium_genomes/samples.info') %>% mutate(label = gsub('_','',Sample))

# Limits for the plot, in case we want to do discrete bins  
limits <- range(branch_summ$pos_frac, na.rm = TRUE)
dummy <- data.frame(
  pos_frac = limits,          # two rows: min and max
  x = c(0,0),
  y = c(0,0)
)
# Build tree data and join branch summaries
td <- fortify(tr) %>% 
  dplyr::left_join(
    branch_summ, by = c("label" = "Branch"))

# Base tree
trm <- left_join(td,md)
p0 <- ggtree(tr) %<+% trm
range(p0@data$pos_frac,na.rm=TRUE)
p <- ggtree(tr) %<+% trm +
  geom_tree(aes(color = pos_frac, linewidth = pos_frac)) +
  geom_tiplab(size = 3.0) +
  geom_nodelab(aes(label=total_pos),hjust=3)+
  theme_tree() +
  scale_color_viridis(
    option = 'plasma',
    name = "Fraction Positive"
  ) +
  scale_linewidth_continuous(
    name  = "Frac FDR<0.1 & dnds>1",
    range = c(0.5, 2.5),limits=limits)

# Expand x-axis to create room for tip labels
xmax <- max(p$data$x, na.rm = TRUE)
p <- p + ggplot2::xlim(0, xmax * 1.3)
p

# zoom
zoomnode <- p$data$node[p$data$label == "F"]
subtree <- tree_subset(tr, zoomnode, levels_back = 0)
g_sub <- ggtree(subtree)
subdf <- g_sub$data
subdf2 <- left_join(subdf, trm %>% select(label, pos_frac, total_pos), by="label")

zoom <- ggtree(subtree) %<+% subdf2 +
  geom_tree(aes(color = pos_frac, linewidth = pos_frac)) +
  geom_tiplab(size = 3.0) +
  geom_nodelab(aes(label=total_pos),hjust=3)+
  theme_tree() +
  scale_color_viridis(
    option = 'plasma',
    name = "Fraction Positive"
  ) +
  scale_linewidth_continuous(
    name  = "Frac FDR<0.1 & dnds>1",
    range = c(0.5, 2.5),limits=limits)+
  xlim(c(0,max(g_sub@data$x)*1.3))
zoom

ggarrange(p,zoom,common.legend = TRUE)
ggsave('~/symlinks/guava/figures/20260706_dNdS_BranchColorsPOSITIVE_zoom.pdf',
       ggarrange(p,zoom),height=4,width=14)


# confirm
library(meRo)
branch_summ %>% arrange(desc(pos_frac))
# # A tibble: 33 × 8
# Branch  med_omega pos_frac total_pos puri_frac total_puri ci_w_med n_genes
# <chr>       <dbl>    <dbl>     <int>     <dbl>      <int>    <dbl>   <int>
#   1 HPSI065     0.218  0.00985         2     0.246         50    0.781     203
# 2 N           0.233  0.00694         3     0.238        103    0.758     432
# 3 HPSI035     0.220  0.00676         2     0.264         78    0.758     296
# 4 I           0.202  0.00543         1     0.234         43    0.757     184
# 5 HPSI027     0.193  0.00521         2     0.310        119    0.701     384
# 6 G           0.114  0.00432         6     0.207        288    0.747    1390
# 7 HPSI003     0.156  0.00358         3     0.228        191    0.730     839
# 8 K           0.162  0.00356         4     0.212        238    0.739    1124
# 9 H           0.153  0.00318         3     0.181        171    0.745     944
#10 L           0.127  0.00317         2     0.238        150    0.759     631
##### What are the positive genes? ######
write.table(wstats_clean %>% filter(pos == TRUE),file='TopGenes_Pos_20260706.tsv',quote=F,sep='\t',row.names=F)

```

## Functional Annotation

Using the interproscan results, merge our top candidate positive selection genes to see if there are any interesting candidates. 

```
setwd("/project/coffea_pangenome/Guava/selection_dnds")
#Positive selection: functional annotation of genes 
library(tidyverse)
library(scales)
library(stringr)
library(GO.db)
library(AnnotationDbi)
library(data.table)	
library(ggrepel)

# Read in genes 
genes <- read_tsv("TopGenes_Pos_20260706.tsv")
wstats_clean <- fread("Node_dNdS_20260701-RInput.tsv.gz") %>% as_tibble
id <- fread("/project/coffea_pangenome/Guava/annotation/interproscan/HSYZ_002.interpro.tsv") %>% as_tibble

func <- left_join(genes,id %>% dplyr::select(-Sample))

interest <- c('HPSI065','N','HPSI035','I','HPSI027','G','HPSI003','K','H','L')
targs <- func  %>% filter(Branch %in% interest & ONTOLOGY == "MF")
targ_func <- targs %>% 
  group_by(Branch) %>% 
  count(TERM) %>% 
  ggplot(aes(y=TERM,x=n,fill=Branch))+
  geom_col()+
  scale_fill_manual(values=viridis(7,option='turbo'))+
  theme_bw()

targs %>% distinct(Branch,Gene) %>% count(Gene)
# # A tibble: 16 × 2
# Gene                   n
# <chr>              <int>
#   1 egapxtmp_000403-R1     1
# 2 egapxtmp_000853-R1     1
# 3 egapxtmp_001841-R1     1
# 4 egapxtmp_003126-R1     1
# 5 egapxtmp_003151-R1     1
# 6 egapxtmp_003607-R1     1
# 7 egapxtmp_006508-R1     1
# 8 egapxtmp_008937-R1     2
# 9 egapxtmp_009693-R1     2
# 10 egapxtmp_010295-R1     1
# 11 egapxtmp_015631-R1     1
# 12 egapxtmp_019199-R1     1
# 13 egapxtmp_019718-R1     1
# 14 egapxtmp_021820-R1     1
# 15 egapxtmp_022971-R1     1
# 16 egapxtmp_027114-R1     1

targs  %>% filter(Gene %in% c('egapxtmp_008937-R1','egapxtmp_009693-R1'))

ggsave('~/symlinks/guava/figures/20260706_function_annotation_positive_top10_clades.pdf',targ_func,height=6,width=9)

```

