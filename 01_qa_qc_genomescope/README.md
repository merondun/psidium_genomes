# 01_qa_qc_genomescope/

* Runs GenomeScope2 to estimate genome kmers.
* Extracts k‑mer pairs and generates a smudgeplot to infer ploidy.

Outputs include smudgeplots with ploidy estimates and genomescope histos:



![Smudges_histos](/imgs/Combined_Histos_Smudges.png)



___



Leveraging HiFi and HiC data. First, admin: assign chromosomes based on chromosome-level guava. Name chrs from that assembly (instead of CMXXXXX.X) with:

```bash
awk '/^>/ {print ">chr" ++i; next} {print}' GCA_016432845.1_guava_v11.23_genomic.fna > renamed_sequences.fna
for i in $(seq 1 11); do echo chr$i >> chrs.list; done
seqtk subseq renamed_sequences.fna chrs.list > GCA_016432845.1_guava_v11.23_genomic.chrs.fa
chr1 40370300
chr2 38465871
chr3 50577630
chr4 48287879
chr5 44769864
chr6 42823316
chr7 35363649
chr8 33379055
chr9 32331238
chr10 37879573
chr11 37021556
```

And also the Syzygium genome:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/216/185/GCA_031216185.1_ASM3121618v1/GCA_031216185.1_ASM3121618v1_genomic.fna.gz

awk '/^>/ {print ">chr" ++i; next} {print}' GCA_031216185.1_ASM3121618v1_genomic.fna > renamed_sequences.fna
for i in $(seq 1 11); do echo chr$i >> chrs.list; done
seqtk subseq renamed_sequences.fna chrs.list > GCA_031216185.1_ASM3121618v1_genomic.chrs.fa

```



# Read Stats

Number of reads and bases, read length for HiFi counted with `seqkit stats`.

```bash
#!/bin/bash

#SBATCH --time=08:00:00 
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1 
#SBATCH --mem=12Gb
#SBATCH --partition=short 

ID=$1

RAWDATA=/project/coffea_pangenome/Guava/RawData
OUTDIR=/project/coffea_pangenome/Guava/Reads_Concatenated

cat ${RAWDATA}/HiFi/${ID}*gz > ${OUTDIR}/${ID}.HiFi.fastq.gz
cat ${RAWDATA}/HiC/${ID}*HiC_R1*gz > ${OUTDIR}/${ID}.HiC.R1.fastq.gz
cat ${RAWDATA}/HiC/${ID}*HiC_R2*gz > ${OUTDIR}/${ID}.HiC.R2.fastq.gz

# Get stats on bases / reads, will output to slurm 
bbduk.sh in=${OUTDIR}/${ID}.HiFi.fastq.gz 2> ${ID}.HiFi.log
bbduk.sh in=${OUTDIR}/${ID}.HiC.R1.fastq.gz 2> ${ID}.HiC.log

# Extract 
HIFI_BASES=$(grep 'Input:' ${ID}.HiFi.log | awk '{print $4}')
HIFI_READS=$(grep 'Input:' ${ID}.HiFi.log | awk '{print $2}')
HIC_BASES=$(grep 'Input:' ${ID}.HiC.log | awk '{print $4}')
HIC_READS=$(grep 'Input:' ${ID}.HiC.log | awk '{print $2}')

echo -e "${ID}\t${HIFI_BASES}\t${HIFI_READS}\t${HIC_BASES}\t${HIC_READS}" > ${ID}.info
```

## Smudgeplot & Genomescope

Run smudge & genome scope from hifi data:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --cpus-per-task=24
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome
 
set -euo pipefail

module load miniconda
source activate puzzler201

SAMPLE=${1:?Provide SAMPLE id as first argument}

# Read tsv
IFS=$'\t,' read -r _ PLOIDY HIFI < <(
    awk -F'[\t,]' -v sample="$SAMPLE" '$1 == sample {print $0}' smudge_info.tsv
)

# Where our analyses are run 
WD=/project/coffea_pangenome/Guava/GenomeSizeEstimation/202501_Justin_Kmer
echo -e "=======================================================================\nParameters for sample: ${SAMPLE} \nPLOIDY: ${PLOIDY}\nHIFI: ${HIFI}\n=======================================================================\n"

# Create directory to save all results 
mkdir -p ${WD}/results_${SAMPLE}
cd ${WD}/results_${SAMPLE}
echo "Working on ${SAMPLE} for file ${HIFI}" 

# This is a kmer value, in general 31 works well for a range of eukaryotes
K=31

# First, generate your kmer table, this is memory and time intensive 
FastK -v -t1 -k${K} -M70 -T24 -NFastK_Table_k${K}_${SAMPLE} ${HIFI}

# Now, generate a coverate histogram from that table 
Histex -G FastK_Table_k${K}_${SAMPLE} > ${SAMPLE}_k${K}.hist

# Run genomescope2 on that histogram table 
genomescope2 --input ${SAMPLE}_k${K}.hist --output . --ploidy ${PLOIDY} --kmer_length ${K} --name_prefix ${SAMPLE}_k${K}

conda deactivate
source activate smudge 

# Create directory to save all results
cd ${WD}/results_${SAMPLE}

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot hetmers -L 12 -t 4 -tmp . -o ${SAMPLE}_kmerpairs --verbose FastK_Table_k${K}_${SAMPLE}

# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot all --format pdf -o ${SAMPLE}_smudge ${SAMPLE}_kmerpairs.smu

# Afterwards, you will want to clean up that massive file 
#rm FastK*
```

## Plot 

Plot:

```R
setwd('~/psidium_genomes/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(ggtext)

md <- read_tsv('samples.info') %>% mutate(HapCoverage = HiFi_Gb / 0.4)
md

info_dat <- NULL
hist_dat <- NULL
files <- list.files('01_qa_qc_genomescope/genomescope/',pattern = '*_k31.info')
for (file in files) {
  id = gsub('_k31.info','',file)
  cat('Processing: ',id,'\n')
  info <- read_tsv(paste0('01_qa_qc_genomescope/genomescope/',file))
  info_dat <- rbind(info_dat,info)
  hist <- read_tsv(paste0('01_qa_qc_genomescope/genomescope/',id,'_k31.hist'),col_names=F)
  hist <- hist %>% mutate(Accession = id)
  hist_dat <- rbind(hist_dat,hist)
}

names(hist_dat) <- c('Bin','Coverage','Accession')
hm <- left_join(hist_dat,md)
hm$Accession <- factor(hm$Accession,levels=unique(md$Accession))

# limits
xlims <- hm %>% group_by(Accession,Species) %>% summarize(xmin=5,xmax=HapCoverage*Ploidy) %>% distinct

hm2 <- hm %>%
  left_join(xlims, by = c("Accession", "Species")) %>%
  mutate(ID = paste0(Accession, " \n(<i>", Species, "</i> ",Ploidy,"N)"))
spord <- hm2 %>% select(ID,Order) %>% arrange(Order) %>% distinct(ID)
hm3 <- hm2 %>% 
  mutate(ID = factor(ID, levels = spord$ID)) %>%
  filter(Bin >= xmin, Bin <= xmax)

hist_plot <- hm3 %>%
  ggplot(aes(x = Bin, y = Coverage, fill = Group)) +
  geom_bar(stat='identity') +
  facet_wrap(~ ID, scales = "free", nrow = 5, ncol = 4) +
  scale_fill_manual(values = md$Color, breaks = md$Group) +
  theme_bw(base_size=5) +
  labs(x='Coverage Peak',y='')+
  theme(
    strip.text = element_markdown(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
hist_plot

ggsave('~/psidium_genomes/01_qa_qc_genomescope/20260603_Genomescope_Histograms.pdf',hist_plot,height=4.5,width=7,dpi=300)

```



| Lab ID  | Accession | Species                     | Ploidy | HiFi Kmer Peak | HiFi Gb Output | Genome_Size_Collapsed (1N) | Genome_Size (1N * Ploidy) | Comments                         |
| ------- | --------- | --------------------------- | ------ | -------------- | -------------- | -------------------------- | ------------------------- | -------------------------------- |
| HPSI_69 | HPSI_069  | Acca sellowiana             | 2      | 49             | 33.13          | 338                        | 676                       |                                  |
| HEUG_1  | HEUG_001  | Eugenia stipitata           | 2      | 49             | 34.3           | 350                        | 700                       |                                  |
| HPSI_10 | HPSI_010  | Psidium acutangulum         | 2      | 21             | 35.16          | 837                        | 1674                      |                                  |
| HPSI_59 | HPSI_059  | Psidium friedrichsthalianum | 4      | 31             | 50.52          | 407                        | 1628                      | 4N?                              |
| HPSI_3  | HPSI_003  | Psidium guajava             | 2      | 34             | 30.67          | 451                        | 902                       | 2N very likely, very unlikely 3N |
| HPSI_7  | HPSI_007  | Psidium guajava             | 2      | 31             | 28.21          | 455                        | 910                       |                                  |
| HPSI_16 | HPSI_016  | Psidium guajava             | 2      | 25             | 21.9           | 438                        | 876                       |                                  |
| HPSI_19 | HPSI_019  | Psidium guajava             | 2      | 24             | 22.8           | 475                        | 950                       |                                  |
| HPSI_27 | HPSI_027  | Psidium guajava             | 2      | 8              | 7.58           | 474                        | 948                       |                                  |
| HPSI_35 | HPSI_035  | Psidium guajava             | 2      | 21             | 18.94          | 451                        | 902                       |                                  |
| HPSI_37 | HPSI_037  | Psidium guajava             | 2      | 50             | 44.33          | 443                        | 886                       |                                  |
| HPSI_41 | HPSI_041  | Psidium guajava             | 2      | 33             | 30.09          | 456                        | 912                       |                                  |
| HPSI_60 | HPSI_060  | Psidium guajava             | 2      | 37.5           | 32.33          | 431                        | 862                       |                                  |
| HPSI_65 | HPSI_065  | Psidium guajava             | 2      | 21             | 18.03          | 429                        | 858                       |                                  |
| HPSI_68 | HPSI_068  | Psidium guajava             | 2      | 21             | 18.4           | 438                        | 876                       |                                  |
| HPSI_72 | HPSI_072  | Psidium microphyllum        | 4      | 26             | 40.83          | 393                        | 1572                      | 4N?                              |
| HPSI_80 | HPSI_080  | Psidium sartorianum         | 4      | 22             | 35.36          | 402                        | 1608                      | 4N?                              |
| HSYZ_2  | HSYZ_002  | Syzygium aqueum             | 2      | 29             | 23.27          | 401                        | 802                       |                                  |
| HSYZ_1  | HSYZ_001  | Syzygium samarangense       | 4      | 11             | 20.12          | 457                        | 1828                      | 4N?                              |
| HSYZ_3  | HSYZ_003  | Syzygium samarangense       | 4      | 22             | 35.15          | 399                        | 1596                      | 4N?                              |

