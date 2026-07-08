# 05_whole_genome_alignments/

## BUSCO-level Synteny

This workflow runs chromsyn to place BUSCO genes onto chromosomes and summarize synteny using BUSCO anchors. It generates plotting inputs (BUSCO tables, telomere tracks, and repeat/telomere-window scores), merges them into a chromsyn report (PDF/XLSX), and summarizes BUSCO counts and total syntenic block lengths/heatmaps in R.

Output:

![chromsyn](/imgs/chromsyn.png)



___



```bash
WD=/project/coffea_pangenome/Guava/Synteny_BUSCO
GENOMES=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm/chrs
```

Run the chromsyn workflow:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=20

# Check if the correct number of arguments is provided
set -euo pipefail

module load miniconda
source activate chromsyn

FASTA="${1:?usage: $0 <FASTA>}"
TARGET=$(basename ${FASTA} .fa)
FILE=$(realpath ${FASTA})
WD=/project/coffea_pangenome/Guava/Synteny_BUSCO
GENOMES=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm/chrs

echo "Working on ${TARGET}, file ${FASTA}"
export PYTHONWARNINGS="ignore::SyntaxWarning"

# Prep busco db 
BUSCO_DB=/project/coffea_pangenome/Software/Merondun/busco_downloads
LINEAGE=embryophyta_odb12
if [ -d "${BUSCO_DB}/lineages/${LINEAGE}" ]; then
        echo "BUSCO db ${LINEAGE} already present at ${BUSCO_DB}/lineages/${LINEAGE} – skipping download"
else
        busco --download ${LINEAGE} --download_path ${BUSCO_DB}
fi

mkdir -p ${WD}/work ${WD}/plotting_inputs
cd ${WD}/work

# Generate inputs
TELO_DIR=/project/coffea_pangenome/Software/Merondun/telociraptor/code
if [ -f ${TARGET}.telomeres.tdt]; then
        echo "Telociraptor output exists for ${TARGET} – skipping"
else
        python ${TELO_DIR}/telociraptor.py seqin=${FILE} basefile=${FILE} i=-1 tweak=F telonull=T
fi

# busco 
if [ -f ${TARGET}.busco5.tsv]; then
        echo "BUSCO already ran on ${TARGET} - skipping"
else
        busco -f -o run_${TARGET} -i ${FILE} -l ${BUSCO_DB}/lineages/${LINEAGE} --cpu ${t} -m genome
        cp -v run_${TARGET}/run_${LINEAGE}/full_table.tsv ${TARGET}.busco5.tsv
        rm -rf run_${TARGET}*
fi

# repeat scores
if [ -f ${TARGET}.tidk.tsv]; then
        echo "TIDK already ran on ${TARGET} - skipping"
else
        tidk search --dir search --output ${TARGET} -s AACCCT ${FILE}
        cp -v search/${TARGET}_telomeric_repeat_windows.tsv ${TARGET}.tidk.tsv
fi

# Copy outputs
cp ${TARGET}.tidk.tsv ${TARGET}.gaps.tsv ${TARGET}.busco5.tsv ${TARGET}.telomeres.tdt ${TARGET}.contigs.tdt ${WD}/plotting_inputs/
```

Merge the outputs and plot:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

WD=/project/coffea_pangenome/Guava/Synteny_BUSCO
cd ${WD}

> busco.fofn > gaps.fofn > sequences.fofn > tidk.fofn

for i in $(cat Samples.list); do 
    echo -e "${i} ${WD}/plotting_inputs/${i}.chr.busco5.tsv" >> busco.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.chr.gaps.tdt" >> gaps.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.chr.telomeres.tdt" >> sequences.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.chr.tidk.tsv" >> tidk.fofn
done 

Rscript ~/symlinks/software/chromsyn/chromsyn.R labelsize=1.5 opacity=0.4 pdfheight=8
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -sOutputFile=output.pdf chromsyn.pdf
```



## Reference-based Dotplots

Assess collinearity just using reference-based dotplots, using the published guava genome `GCA_016432845.1_guava_v11.23` as the reference - except for the wax apples, use our water apple as a reference since they are very diverged: `HSYZ_002.fa` 

```bash
#!/bin/bash

#SBATCH --time=0-02:00:00   
#SBATCH --cpus-per-task=4
#SBATCH --mem=16Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate puzzler201

t=16
set -euo pipefail
SAMPLE=${1:?Provide SAMPLE as first argument}
GENOMES=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm/chrs
WD=/project/coffea_pangenome/Guava/dotplots

# Default REF
REFERENCE=/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa

if [[ "${SAMPLE}" == "HSYZ_001" || "${SAMPLE}" == "HSYZ_002" || "${SAMPLE}" == "HSYZ_003" ]]; then
    REFERENCE="/project/coffea_pangenome/Guava/Assemblies/GCA_031216185.1_ASM3121618v1_genomic.chrs.fa"
fi

mkdir -p ${WD}/chrs_pafs

echo "WGA for ${SAMPLE}"
mashmap -r ${REFERENCE} -q ${GENOMES}/${SAMPLE}.chr.fa -t ${t} -s 10000 --perc_identity 85 -o ${WD}/chrs_pafs/${SAMPLE}.paf
Rscript ~/apptainer/paf2dotplot.R ${WD}/chrs_pafs/${SAMPLE}.paf -r 1e6 -m 1e4 -p 4 -c 1 -i chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11 --identity-lower-color 100


```

## Minimap & NGenomeSyn

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --cpus-per-task=5
#SBATCH --mem=24Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <first_pair> <second_pair> <asmX>"
    exit 1
fi 

#Submit as cat MM_Pairs.list  | xargs -L 1 sbatch 02_Guava_Minimap.sh
#module load miniconda
#source activate wga

# Use a file with $SAMPLE_1 \t $SAMPLE_2 \t $ASM5/10/20 within Pairs.list, submit
ref=$1  #reference
qry=$2 #query
sens=$3 #minimap -asm 
GENOME_DIR=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm/chrs/
echo "Working on $ref aligned to $qry with -${sens}"

MIN_LEN=5000

mkdir -p minimap
cd minimap

#align
#minimap2 -x ${sens} -t 10 ${GENOME_DIR}/${ref}.chr.fa  ${GENOME_DIR}/${qry}.chr.fa > ${ref}-${qry}.paf

#filter len
awk -v l=${MIN_LEN} '
  # Only keep long alignments
  ($4 - $3 >= l) {

    # Extract divergence (dv:f:xxx)
    dv = -1
    for(i=12; i<=NF; i++){
      if($i ~ /^dv:f:/){
        split($i,a,":")
        dv = a[3]
      }
    }

    # If dv missing, skip
    if(dv < 0) next

    # Identity check
    ident = 1 - dv
    if(ident < 0.95) next

    # Only primary mappings
    if($0 !~ /tp:A:P/) next

    # Output 6-column link format for NGenomeSyn:
    # ref_chr, ref_start, ref_end, qry_chr, qry_start, qry_end
    printf "%s\t%d\t%d\t%s\t%d\t%d\n", $1, $3, $4, $6, $8, $9
  }
' ${ref}-${qry}.paf > ${ref}-${qry}.links
```

NGenomeSyn:

Create the chrs files:

```bash
GENOME_DIR=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm/
mkdir -p chrsizes
for s in $(cat Samples.list); do 
    awk '{OFS="\t"}{print $1, "1", $2}' ${GENOME_DIR}/chrs/${s}.chr.fa.fai > chrsizes/${s}.len;
done 
```

Config file:

```
###  # Detailed parameters please see  at the  NGenomeSyn_manual_English.pdf 

SetParaFor = global
########## global and required parameters ######
GenomeInfoFile1=./chrsizes/HPSI_003.len
GenomeInfoFile2=./chrsizes/HPSI_027.len
GenomeInfoFile3=./chrsizes/HPSI_060.len
GenomeInfoFile4=./chrsizes/HPSI_041.len
GenomeInfoFile5=./chrsizes/HPSI_019.len
GenomeInfoFile6=./chrsizes/HPSI_007.len
GenomeInfoFile7=./chrsizes/HPSI_068.len
GenomeInfoFile8=./chrsizes/HPSI_065.len
GenomeInfoFile9=./chrsizes/HPSI_037.len
GenomeInfoFile10=./chrsizes/HPSI_035.len
GenomeInfoFile11=./chrsizes/HPSI_016.len
GenomeInfoFile12=./chrsizes/HPSI_080.len
GenomeInfoFile13=./chrsizes/HPSI_059.len
GenomeInfoFile14=./chrsizes/HPSI_010.len
GenomeInfoFile15=./chrsizes/HPSI_072.len
GenomeInfoFile16=./chrsizes/HPSI_069.len
GenomeInfoFile17=./chrsizes/HEUG_001.len
GenomeInfoFile18=./chrsizes/HSYZ_002.len
GenomeInfoFile19=./chrsizes/HSYZ_001.len
GenomeInfoFile20=./chrsizes/HSYZ_003.len

LinkFileRef1VsRef2=./links/HPSI_003-HPSI_027.intra.links
LinkFileRef2VsRef3=./links/HPSI_027-HPSI_060.intra.links
LinkFileRef3VsRef4=./links/HPSI_060-HPSI_041.intra.links
LinkFileRef4VsRef5=./links/HPSI_041-HPSI_019.intra.links
LinkFileRef5VsRef6=./links/HPSI_019-HPSI_007.intra.links
LinkFileRef6VsRef7=./links/HPSI_007-HPSI_068.intra.links
LinkFileRef7VsRef8=./links/HPSI_068-HPSI_065.intra.links
LinkFileRef8VsRef9=./links/HPSI_065-HPSI_037.intra.links
LinkFileRef9VsRef10=./links/HPSI_037-HPSI_035.intra.links
LinkFileRef10VsRef11=./links/HPSI_035-HPSI_016.intra.links
LinkFileRef11VsRef12=./links/HPSI_016-HPSI_080.intra.links
LinkFileRef12VsRef13=./links/HPSI_080-HPSI_059.intra.links
LinkFileRef13VsRef14=./links/HPSI_059-HPSI_010.intra.links
LinkFileRef14VsRef15=./links/HPSI_010-HPSI_072.intra.links
LinkFileRef15VsRef16=./links/HPSI_072-HPSI_069.intra.links
LinkFileRef16VsRef17=./links/HPSI_069-HEUG_001.intra.links
LinkFileRef17VsRef18=./links/HEUG_001-HSYZ_002.intra.links
LinkFileRef18VsRef19=./links/HSYZ_002-HSYZ_001.intra.links
LinkFileRef19VsRef20=./links/HSYZ_001-HSYZ_003.intra.links

LinkFileRef1VsRef2=./links/HPSI_003-HPSI_027.inter.links
LinkFileRef2VsRef3=./links/HPSI_027-HPSI_060.inter.links
LinkFileRef3VsRef4=./links/HPSI_060-HPSI_041.inter.links
LinkFileRef4VsRef5=./links/HPSI_041-HPSI_019.inter.links
LinkFileRef5VsRef6=./links/HPSI_019-HPSI_007.inter.links
LinkFileRef6VsRef7=./links/HPSI_007-HPSI_068.inter.links
LinkFileRef7VsRef8=./links/HPSI_068-HPSI_065.inter.links
LinkFileRef8VsRef9=./links/HPSI_065-HPSI_037.inter.links
LinkFileRef9VsRef10=./links/HPSI_037-HPSI_035.inter.links
LinkFileRef10VsRef11=./links/HPSI_035-HPSI_016.inter.links
LinkFileRef11VsRef12=./links/HPSI_016-HPSI_080.inter.links
LinkFileRef12VsRef13=./links/HPSI_080-HPSI_059.inter.links
LinkFileRef13VsRef14=./links/HPSI_059-HPSI_010.inter.links
LinkFileRef14VsRef15=./links/HPSI_010-HPSI_072.inter.links
LinkFileRef15VsRef16=./links/HPSI_072-HPSI_069.inter.links
LinkFileRef16VsRef17=./links/HPSI_069-HEUG_001.inter.links
LinkFileRef17VsRef18=./links/HEUG_001-HSYZ_002.inter.links
LinkFileRef18VsRef19=./links/HSYZ_002-HSYZ_001.inter.links
LinkFileRef19VsRef20=./links/HSYZ_001-HSYZ_003.inter.links

#### LINK COLORS
SetParaFor=Link1
fill=grey
stroke=grey
SetParaFor=Link2
fill=grey
stroke=grey
SetParaFor=Link3
fill=grey
stroke=grey
SetParaFor=Link4
fill=grey
stroke=grey
SetParaFor=Link5
fill=grey
stroke=grey
SetParaFor=Link6
fill=grey
stroke=grey
SetParaFor=Link7
fill=grey
stroke=grey
SetParaFor=Link8
fill=grey
stroke=grey
SetParaFor=Link9
fill=grey
stroke=grey
SetParaFor=Link10
fill=grey
stroke=grey
SetParaFor=Link11
fill=grey
stroke=grey
SetParaFor=Link12
fill=grey
stroke=grey
SetParaFor=Link13
fill=grey
stroke=grey
SetParaFor=Link14
fill=grey
stroke=grey
SetParaFor=Link15
fill=grey
stroke=grey
SetParaFor=Link16
fill=grey
stroke=grey
SetParaFor=Link17
fill=grey
stroke=grey
SetParaFor=Link18
fill=grey
stroke=grey
SetParaFor=Link19
fill=grey
stroke=grey
SetParaFor=Link20
fill=grey
stroke=green
SetParaFor=Link21
fill=grey
stroke=green
SetParaFor=Link22
fill=grey
stroke=green
SetParaFor=Link23
fill=grey
stroke=green
SetParaFor=Link24
fill=grey
stroke=green
SetParaFor=Link25
fill=grey
stroke=green
SetParaFor=Link26
fill=grey
stroke=green
SetParaFor=Link27
fill=grey
stroke=green
SetParaFor=Link28
fill=grey
stroke=green
SetParaFor=Link29
fill=grey
stroke=green
SetParaFor=Link30
fill=grey
stroke=green
SetParaFor=Link31
fill=grey
stroke=green
SetParaFor=Link32
fill=grey
stroke=green
SetParaFor=Link33
fill=grey
stroke=green
SetParaFor=Link34
fill=grey
stroke=green
SetParaFor=Link35
fill=grey
stroke=green
SetParaFor=Link36
fill=grey
stroke=green
SetParaFor=Link37
fill=grey
stroke=green
SetParaFor=Link38
fill=grey
stroke=green


############ canvas and figure[optional] ######
#body=1200                 ### size of canvas with width and height. plot region:  (up/down/left/right)=(55,25,100,120)
up=200
down=200
left=100
right=120
#CanvasHeightRitao=0.25     ## adjust height of the plot
#CanvasWidthRitao=1.0      ## adjust width of the plot
NoPng=1                   ## No OutPut the png File


############# adjust genome setting [optional]##########
#SetParaFor = Genome1       ## GenomeALL/GenomeX
#ZoomChr=1.0               ## adjust chr length, 1 for equal; >1 for enlarge; <1  for
#RotateChr=30              ## rotate the chr with 30 degrees
#ShiftX=0
#ShiftY=0                  ## move the start of chr to (X,Y)
#MoveToX                   ## MoveToY   ## similar to ShiftX and ShiftY


#ChrWidth=20               ## chr width
#LinkWidth=180             ## link height between this genome and next genome
#ChrSpacing=10             ## spacing width of chr/scaffolds
#NormalizedScale=0         ## custom scale for the geome relative the  default.
#SpeRegionFile=Spe.bed     ## input file for highlighted regions[chr start end key1=value1] in the genome.
#ZoomRegion                ## Zoom the specific Region,format (ZoomRegion=chr2:1000:5000)


#GenomeNameRatio
#GenomeName
## GenomeName  GenomeNameSizeRatio  GenomeNameColor  GenomeNameShiftX GenomeNameShiftY
## ChrNameShow ChrNameShiftX ChrNameShiftY ChrNameSizeRatio ChrNameColor ChrNameRotate
## ShowCoordinates=1     ## Show Coordinates with other para [ScaleNum=10 ScaleUpDown ScaleUnit LabelUnit  LablefontsizeRatio  RotateAxisText NoShowLabel ]


#SetParaFor = Genome2       ##

############# adjust link setting [optional]##########
#SetParaFor=Link1           ### LinkALL/LinkX  setting for link X
#StyleUpDown=UpDown        ## UpDown DownUp UpUp DownDown line
#Reverse=1                 ## reverse links
#HeightRatio=1.0           ## ratio of link height relative to the default
## other attributes:  fill|stroke|stroke-opacity|fill-opacity|stroke-width


```

```bash
~/symlinks/software/NGenomeSyn/bin/NGenomeSyn -InConf ngenomesyn.config -OutPut test.svg
```

