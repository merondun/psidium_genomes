# 09_repeats/

* Run earlgrey v7.2.2 to annotate TEs and other repeats, plot high level summaries. 

![repeats](/imgs/20260630_RepeatsHighLevelSummary.png)



____



Run earlgrey to annotate repeats:

```bash
#!/bin/bash

#SBATCH --time=14-00:00:00   
#SBATCH --cpus-per-task=16
#SBATCH --mem=512Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=16

#module load miniconda
#source activate earlgrey722

SAMPLE=$1
WD=/90daydata/coffea_pangenome/scratch/repeats
ASM=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm

cd ${WD}

FASTA="${ASM}/${SAMPLE}.fa" 

echo -e "\e[43m~~~~ Starting repeat annotation for ${SAMPLE} ~~~~\e[0m"
# Run earlgrey with eudicotyledons repeatmasker search time, generating soft-masked genome and run helitrons. 
earlGrey -r eudicotyledons -d yes -t ${t} -g ${FASTA} -q yes -s ${SAMPLE} -o /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/repeats
```

And copy:

```bash
DIR=/project/coffea_pangenome/Guava/assemblies/20250101_JustinAssemblies/repeats
REP=/project/coffea_pangenome/Guava/repeats

for SAMPLE in $(cat Samples.list) ; do 
cp ${DIR}/${SAMPLE}_EarlGrey/${SAMPLE}_summaryFiles/* ${REP}/
done 
```

Add accession to each output:

```bash
for SAMPLE in $(ls *.familyLevelCount.txt | sed 's/.familyLevelCount.txt//g'); do 
echo "${SAMPLE}"
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}.familyLevelCount.txt > ${SAMPLE}.families.out
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}.highLevelCount.txt > ${SAMPLE}.summary.out
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}_divergence_summary_table.tsv > ${SAMPLE}.divergence.out
done 

mergem *families.out > Repeat_Families.txt
mergem *summary.out > Repeat_Summaries.txt
mergem *divergence.out > Divergence_Summaries.txt
```

## Plot Repeat Variation

```
setwd('/project/coffea_pangenome/Guava/repeats/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
library(vegan)
library(broom)
library(ggrepel)
library(caper)

# Add metadata information
md <- read_tsv('~/psidium_genomes/samples.info')

##### High Level #####
high_level <- read_tsv('Repeat_Summaries.txt') 
names(high_level) <- c('Classification','Copies','Count','Proportion','Gen','Distinct_Classifications','Sample')
fl <- left_join(high_level,md) %>% filter(!grepl('HSYZ_001|HSYZ_003',Sample))
fl$Sample <- factor(fl$Sample,levels=md$Sample)
fl <- fl %>%
  arrange(Sample) %>%
  mutate(Group = factor(Group, levels = unique(Group)))

# aggregate nested
fl <- fl %>%
  mutate(ClassSimple = gsub("-nested", "", Classification))
fl_sum <- fl %>%
  group_by(Sample, Group, ClassSimple) %>%
  summarise(Proportion = sum(Proportion), .groups="drop") %>% 
  filter(ClassSimple != "Total Interspersed Repeat")

fl_sum$ClassSimple <- factor(fl_sum$ClassSimple,levels=c('Non-Repeat','Unclassified','Other (Simple Repeat, Microsatellite, RNA)','DNA','Penelope','Rolling Circle','LTR','LINE','SINE'))
cols <- fl_sum %>% dplyr::select(ClassSimple) %>% distinct %>% mutate(Color = brewer.pal(9,'Paired'))

# ltr labs
fl_labels <- fl_sum %>%
  filter(ClassSimple == "LTR") %>%
  mutate(
    label = paste0(round(Proportion, 1), "%"),
    text_color = ifelse(Proportion > 8, "black", "black")
  )


# Plot landscape 
all <- fl_sum %>% 
  #mutate(Coverage = Coverage / 1e6) %>% 
  pivot_longer(c(Proportion)) %>%
  #filter( !(name == 'Distinct_Classifications' & (Classification == 'Unclassified' | Classification == "Other (Simple Repeat, Microsatellite, RNA)")) ) %>% 
  ggplot(aes(y = Sample, x = value, fill = ClassSimple)) +
  geom_bar(stat = 'identity', position = position_stack()) +
  #add LTR percent labels
  geom_text(
    data = fl_labels,
    aes(y = Sample, x = Proportion, label = label),
    position = position_stack(vjust = 0.5),
    color = fl_labels$text_color,
    size = 2.5
  ) +
  theme_bw() +
  facet_grid(Group ~ name, scales = 'free', space = 'free_y') +
  scale_fill_manual(values = cols$Color, breaks = cols$ClassSimple) +
  theme(strip.text.y = element_text(angle = 0),
        legend.position='top',
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5)) +
  ylab('') + xlab('Distinct Classifications') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

all
ggsave('~/symlinks/guava/figures/20260630_RepeatsHighLevelSummary.pdf',
       all,dpi=300,height=5.5,width=5)

fl_sum %>% dplyr::select(Sample,Group,ClassSimple,Proportion) %>% pivot_wider(names_from = ClassSimple,values_from=Proportion)
# Sample   Group                        DNA  LINE   LTR `Non-Repeat` Other (Simple Repeat…¹ Penelope `Rolling Circle`    SINE Unclassified
# <fct>    <fct>                      <dbl> <dbl> <dbl>        <dbl>                  <dbl>    <dbl>            <dbl>   <dbl>        <dbl>
#   1 HPSI_003 Psidium guajava (3n)        1.91 1.91   19.0         44.0                   3.19        0            0.297 0.00530         30.4
# 2 HPSI_027 Psidium guajava (2n)        1.86 1.88   19.9         44.1                   3.21        0            0.387 0.00317         29.3
# 3 HPSI_060 Psidium guajava (2n)        1.79 1.71   18.1         43.6                   2.95        0            0.394 0.00455         32.1
# 4 HPSI_041 Psidium guajava (2n)        1.75 1.76   20.2         43.9                   2.68        0            0.399 0.00373         29.9
# 5 HPSI_019 Psidium guajava (2n)        1.77 1.69   19.3         43.5                   3.24        0            0.357 0.00878         30.8
# 6 HPSI_007 Psidium guajava (2n)        1.86 1.36   27.0         43.4                   3.46        0            0.443 0.00903         23.2
# 7 HPSI_068 Psidium guajava (2n)        2.00 2.01   21.2         45.9                   1.94        0            0.408 0.00703         27.2
# 8 HPSI_065 Psidium guajava (2n)        1.93 1.94   19.5         45.6                   1.62        0            0.323 0.00983         29.6
# 9 HPSI_037 Psidium guajava (2n)        1.93 1.86   22.3         44.5                   3.53        0            0.366 0.00237         26.2
# 10 HPSI_035 Psidium guajava (2n)        1.93 1.91   20.9         45.2                   2.96        0            0.369 0.0131          27.4
# 11 HPSI_016 Psidium guajava (2n)        1.93 2.03   23.7         44.4                   4.20        0            0.421 0.00661         24.1
# 12 HPSI_080 Psidium sartorianum (4n)    2.22 1.71   24.4         46.0                   2.90        0            0.606 0.00687         23.0
# 13 HPSI_059 Psidium friedrichsthalian…  2.16 1.60   27.6         44.0                   1.90        0            0.569 0.00864         22.8
# 14 HPSI_010 Psidium acutangulum (4n)    2.62 1.75   27.0         43.7                   2.05        0            1.06  0.00263         22.7
# 15 HPSI_072 Psidium microphyllum (4n)   2.44 1.84   22.1         45.9                   2.90        0            0.908 0.0108          24.5
# 16 HPSI_069 Feijoa sellowiana (2n)      1.31 1.19   17.9         57.9                   5.70        0            0.456 0.00287         16.1
# 17 HEUG_001 Eugenia stipitata (2n)      3.77 0.816  23.3         52.0                   6.28        0            0.343 0.00241         14.1
# 18 HSYZ_002 Syzygium aqueum (2n)        2.33 2.65   24.6         45.6                   3.48        0            1.53  0.00233         20.6

write.csv(fl_sum %>% dplyr::select(Sample,Group,ClassSimple,Proportion),'~/psidium_genomes/09_repeats/Repeat_proportions_coverage_summarized.csv',quote = F,row.names = F)

```

