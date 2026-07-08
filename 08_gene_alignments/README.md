# 08_gene_alignments/

* Extracts longest‑isoform CDS/pep/gff3 per gene for all samples, converts them to BED, and prepares formatted FASTA files for JCVI synteny.

* Runs pairwise genomic alignments with JCVI (orthologs, anchors, dotplots), filters anchors, labels inter‑chromosomal links, and generates chromosome layout + karyotype plots for multi‑sample comparative synteny.

Output:

![jcvi alignment](/imgs/20260610_gene_alignments_min10.png)



___



## Prep

Extract cds and essential files for single longest transcript per gene:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=6
#SBATCH --mem=32Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

WD=/project/coffea_pangenome/Guava/annotation/
mkdir -p ${WD}/only_longest_transcript_per_gene
cd ${WD}

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sample>"
    exit 1
fi
# module load miniconda
# source activate isoseq_ann

SAMPLE=$1

TARGET=/project/coffea_pangenome/Guava/assemblies/unmasked/${SAMPLE}.fa

agat_sp_keep_longest_isoform.pl --gff gtfs/${SAMPLE}.gtf -o only_longest_transcript_per_gene/${SAMPLE}.longest_transcript_per_gene.gtf
gffread only_longest_transcript_per_gene/${SAMPLE}.longest_transcript_per_gene.gtf \
    -o only_longest_transcript_per_gene/${SAMPLE}.gff3 --keep-genes -O \
    -g ${TARGET} \
    -w only_longest_transcript_per_gene/${SAMPLE}.cds.fa \
    -y only_longest_transcript_per_gene/${SAMPLE}.pep.fa
```

Submit serial via sample:

```bash
cat Samples.list | xargs -I {} sbatch -J ann_{} 01_Extra_Single_Transcript.sh {} 
```

Convert to bed with jcvi:

```bash
WD=/project/coffea_pangenome/Guava/annotation/jcvi
cd ${WD}
# source activate jcvi
for i in $(cat  ../Samples.list); do 
	echo "Working on ${i}"
	python -m jcvi.formats.gff bed --type=transcript --key=ID ${i}.gff3 -o ${i}.bed
	python -m jcvi.formats.fasta format ${i}.cds.fa ${i}.cds
done 
```

## Alignments

Synteny with JCVI:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate 
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_pair> <second_pair>"
    exit 1
fi

reference=$1  #reference
target=$2 #target
echo "Working on $target aligned to $reference"

WD=/project/coffea_pangenome/Guava/annotation/jcvi
cd ${WD}

rm ${reference}*last* ${reference}*lifted* ${reference}*pdf

# Align and extract. Note the formatting must be EXACT with file names, cscore = RBH according to jcvi documentation
python -m jcvi.compara.catalog ortholog ${reference} ${target} --cscore=.99 
#--no_strip_names 

# for dotplots
magick -density 300 ${reference}.${target}.pdf ../output_jcvi/${reference}_${target}.png

# for filtering anchors 
python -m jcvi.compara.synteny screen --minsize=10 --simple ${reference}.${target}.anchors ${reference}.${target}.anchors.minsize10

```

Submit:

```bash
cat Pairs.list | xargs -L 1 sbatch 03_JCVI_Alignment.sh
```

For inter-chromosomal colors:

```bash
#!/usr/bin/env bash

# Usage:
#   ./label_interchrom.sh Pairs.list
#
# inputs:
#   ${ref}.bed
#   ${tgt}.bed
#   ${ref}.${tgt}.anchors.simple
#
# Produces:
#   ${ref}.${tgt}.anchors.cols.simple   (annotated with g* for inter-chr)

pairs_file="$1"

while read -r ref tgt; do
    echo "Processing pair: $ref  $tgt"

    ref_bed="${ref}.bed"
    tgt_bed="${tgt}.bed"
    anchor_file="${ref}.${tgt}.anchors.simple"
    out_file="${ref}.${tgt}.anchors.cols.simple"

    if [[ ! -f "$ref_bed" || ! -f "$tgt_bed" || ! -f "$anchor_file" ]]; then
        echo "  Missing files for pair ($ref, $tgt). Skipping."
        continue
    fi

    # Build chromosome lookup tables
    awk '{chr[$4]=$1} END {for (k in chr) print k, chr[k]}' "$ref_bed" > ref.chr.map
    awk '{chr[$4]=$1} END {for (k in chr) print k, chr[k]}' "$tgt_bed" > tgt.chr.map

    # Annotate anchors
    awk '
        BEGIN {
            while (getline <"ref.chr.map" > 0) { refchr[$1]=$2 }
            while (getline <"tgt.chr.map" > 0) { tgtchr[$1]=$2 }
        }
        {
            r1=$1; r2=$2; t1=$3; t2=$4;
            inter = (refchr[r1] != tgtchr[t1] || refchr[r2] != tgtchr[t2])

            if (inter)
                print "g* "$0
            else
                print $0
        }
    ' "$anchor_file" > "$out_file"

    echo "  -> wrote $out_file"

done < "$pairs_file"

# Cleanup temporary maps
```



## Plot 

Create a layout file:

```
#	y,	xstart,	xend,	rotation,	color,	label,	va,	bed
	0.95,	0.2,	0.95,	0,	,	HPSI_003,	top,	HPSI_003.bed
	0.89,	0.2,	0.95,	0,	,	HPSI_027,	top,	HPSI_027.bed
	0.84,	0.2,	0.95,	0,	,	HPSI_060,	top,	HPSI_060.bed
	0.79,	0.2,	0.95,	0,	,	HPSI_041,	top,	HPSI_041.bed
	0.74,	0.2,	0.95,	0,	,	HPSI_019,	top,	HPSI_019.bed
	0.68,	0.2,	0.95,	0,	,	HPSI_007,	top,	HPSI_007.bed
	0.63,	0.2,	0.95,	0,	,	HPSI_068,	top,	HPSI_068.bed
	0.58,	0.2,	0.95,	0,	,	HPSI_065,	top,	HPSI_065.bed
	0.53,	0.2,	0.95,	0,	,	HPSI_037,	top,	HPSI_037.bed
	0.47,	0.2,	0.95,	0,	,	HPSI_035,	top,	HPSI_035.bed
	0.42,	0.2,	0.95,	0,	,	HPSI_016,	top,	HPSI_016.bed
	0.37,	0.2,	0.95,	0,	,	HPSI_080,	top,	HPSI_080.bed
	0.32,	0.2,	0.95,	0,	,	HPSI_059,	top,	HPSI_059.bed
	0.26,	0.2,	0.95,	0,	,	HPSI_010,	top,	HPSI_010.bed
	0.21,	0.2,	0.95,	0,	,	HPSI_072,	top,	HPSI_072.bed
	0.16,	0.2,	0.95,	0,	,	HPSI_069,	top,	HPSI_069.bed
	0.11,	0.2,	0.95,	0,	,	HEUG_001,	top,	HEUG_001.bed
	0.05,	0.2,	0.95,	0,	,	HSYZ_002,	bottom,	HSYZ_002.bed
#       edges
e,      0,      1,      HPSI_003.HPSI_027.anchors.cols.simple
e,      1,      2,      HPSI_027.HPSI_060.anchors.cols.simple
e,      2,      3,      HPSI_060.HPSI_041.anchors.cols.simple
e,      3,      4,      HPSI_041.HPSI_019.anchors.cols.simple
e,      4,      5,      HPSI_019.HPSI_007.anchors.cols.simple
e,      5,      6,      HPSI_007.HPSI_068.anchors.cols.simple
e,      6,      7,      HPSI_068.HPSI_065.anchors.cols.simple
e,      7,      8,      HPSI_065.HPSI_037.anchors.cols.simple
e,      8,      9,      HPSI_037.HPSI_035.anchors.cols.simple
e,      9,      10,     HPSI_035.HPSI_016.anchors.cols.simple
e,      10,     11,     HPSI_016.HPSI_080.anchors.cols.simple
e,      11,     12,     HPSI_080.HPSI_059.anchors.cols.simple
e,      12,     13,     HPSI_059.HPSI_010.anchors.cols.simple
e,      13,     14,     HPSI_010.HPSI_072.anchors.cols.simple
e,      14,     15,     HPSI_072.HPSI_069.anchors.cols.simple
e,      15,     16,     HPSI_069.HEUG_001.anchors.cols.simple
e,      16,     17,     HEUG_001.HSYZ_002.anchors.cols.simple
```

And the chrs.txt file. I took this ordering from the chromsyn ordering. 

```
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11
```

And create karyoplot:

```bash
python -m jcvi.graphics.karyotype \
	--format png --font Arial --seed 1 \
	-o ../output_jcvi/20260630_gene_alignments_interchr_min10.pdf \
	chrs.txt chr_layout.txt --basepair
```



## Interproscan Annotation

Annotate the proteins for select samples using interproscan:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=8
#SBATCH --mem=52Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load interproscan/5.78-109.0

SAMPLE="${1:?usage: $0 <SAMPLE>}"
DATA_DIR=/project/coffea_pangenome/Guava/annotation/only_longest_transcript_per_gene
echo "Running on ${SAMPLE}"

interproscan.sh \
  -i ${DATA_DIR}/${SAMPLE}.pep.fa \
  -f tsv \
  --goterms \
  --pathways \
  --cpu 8
  
# Clean annotations
module unload interproscan/5.78-109.0
module load miniconda
source activate r
Rscript 02_clean_interpro.R ${SAMPLE}
```

Ensure that this script is available in the `cwd` which will simplify the annotations:

```R
#!/usr/bin/env Rscript

# Simplify interproscan annotations
# mamba create -n r r-base=4.4.1 r-tidyverse r-ggpubr bioconda::bioconductor-go.db bioconda::bioconductor-annotationdbi r-data.table
library(tidyverse)
library(GO.db)
library(AnnotationDbi)
library(data.table)	

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript 02_clean_interpro.R SAMPLE_ID\n")
}

sample <- args[1]

cat("Processing sample:", sample, "\n")

base_dir <- "/project/coffea_pangenome/Guava/annotation/interproscan"

input_file  <- file.path(base_dir, paste0(sample, ".pep.fa.tsv"))
output_file <- file.path(base_dir, paste0(sample, ".interpro.tsv"))

cat("Input :", input_file, "\n")
cat("Output:", output_file, "\n")


id <- fread(input_file) %>% as_tibble %>%
  transmute(Gene = V1, source_db = V4, interpro_id = V12, go_raw = V14) %>%
  distinct()

gene2go_long <- id %>%
  filter(!is.na(go_raw), go_raw != "-", str_detect(go_raw, "GO:")) %>%
  # split on |
  separate_rows(go_raw, sep = "\\|") %>%
  # extract the GO:######## part only 
  mutate(GO = str_extract(go_raw, "GO:\\d{7}")) %>%
  filter(!is.na(GO)) %>%
  distinct(Gene, GO)

gene2go_long

# Map ontology
go_ontology <- AnnotationDbi::select(
  GO.db,
  keys = unique(gene2go_long$GO),
  columns = c("ONTOLOGY", "TERM"),
  keytype = "GOID"
) %>%
  as_tibble() %>%
  dplyr::rename(GO = GOID) %>%
  distinct(GO, ONTOLOGY, TERM)

gene2go_annot <- gene2go_long %>%
  inner_join(go_ontology, by = "GO") %>%
  filter(!is.na(ONTOLOGY)) %>% 
  mutate(Sample = sample) %>% 
  relocate(Sample)

gene2go_annot
write.table(
  gene2go_annot,
  output_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

cat("Done.\n")
```

