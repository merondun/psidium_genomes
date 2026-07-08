# 03_purge_dups_sensitivity/

Assess purging thresholds for a primary assembly, using 4 different coverage thresholds. Evalulate with BUSCO and juicer outputs:

* Create HiFi-only assembly
* Compute coverage and cutoffs for purge dups
* Run purge_dups trials
  * Extract purged and haplotig FASTAs
  * Map Hi‑C reads
  * Run HapHiC scaffolding
  * Generate Juicer/Hi-C maps
  * Run BUSCO completeness check

Outputs, for the 4 tetraploid samples. I selected **agg2**, which consistently showed the lowest duplicated BUSCO counts. agg3 was also very similar, but all performed better than default, so I kept agg2 for consistency. 

Compared to default, this modifies the cutoffs:

**c (lower haploid coverage boundary)** is increased by **20%**

**d (upper haploid coverage boundary)** is increased by **20%**

**e (lower diploid coverage boundary)** is decreased by **15%**



| Sample   | Label   | SizeBP    | Sequences | Contigs | Gaps | ContigN50 | ScafN50  | GC      | Merqury_Complete | Merqury_QV | BUSCO_Complete | BUSCO_singlecopy | BUSCO_duplicated |
| -------- | ------- | --------- | --------- | ------- | ---- | --------- | -------- | ------- | ---------------- | ---------- | -------------- | ---------------- | ---------------- |
| HPSI_010 | agg1    | 420729839 | 62        | 90      | 28   | 23152470  | 38055516 | 0.40533 | 51.432           | 71.005     | 97.8           | 94.8             | 3                |
| HPSI_010 | agg2    | 416559177 | 60        | 90      | 30   | 23152470  | 38055516 | 0.4054  | 51.1519          | 70.9618    | 97.8           | 95.6             | 2.2              |
| HPSI_010 | agg3    | 416525142 | 60        | 91      | 31   | 23152470  | 38055516 | 0.4054  | 51.1515          | 70.9614    | 97.8           | 95.6             | 2.2              |
| HPSI_010 | default | 434884293 | 68        | 93      | 25   | 18026470  | 37774841 | 0.40518 | 52.3239          | 71.0326    | 97.8           | 89.2             | 8.6              |
| HPSI_059 | agg1    | 456255667 | 94        | 158     | 64   | 10590830  | 22622641 | 0.40333 | 60.6143          | 70.7455    | 98.8           | 88.3             | 10.5             |
| HPSI_059 | agg2    | 436839903 | 85        | 156     | 71   | 12064372  | 23605113 | 0.4038  | 59.4367          | 70.6656    | 98.7           | 92.1             | 6.6              |
| HPSI_059 | agg3    | 434939542 | 85        | 164     | 79   | 10590830  | 23479345 | 0.40381 | 59.3421          | 70.6467    | 98.7           | 92.8             | 5.9              |
| HPSI_059 | default | 473313350 | 102       | 164     | 62   | 13058537  | 23588982 | 0.40322 | 61.5924          | 70.7308    | 98.8           | 84.3             | 14.5             |
| HPSI_072 | agg1    | 400064285 | 38        | 84      | 46   | 17255591  | 35843552 | 0.3993  | 43.7769          | 69.6914    | 97.2           | 94.4             | 2.8              |
| HPSI_072 | agg2    | 399870954 | 36        | 85      | 49   | 16156023  | 35831952 | 0.3993  | 43.776           | 69.6893    | 97.2           | 94.4             | 2.8              |
| HPSI_072 | agg3    | 399734829 | 36        | 87      | 51   | 16156023  | 35831952 | 0.39929 | 43.7695          | 69.6878    | 97.2           | 94.5             | 2.7              |
| HPSI_072 | default | 400381121 | 38        | 78      | 40   | 17045874  | 35852492 | 0.39934 | 43.7896          | 69.6948    | 97.2           | 94.3             | 2.9              |
| HPSI_080 | agg1    | 530717058 | 284       | 311     | 27   | 14773535  | 34407057 | 0.41065 | 58.8413          | 69.5647    | 98.9           | 62.4             | 36.4             |
| HPSI_080 | agg2    | 428814935 | 280       | 305     | 25   | 14773535  | 32194574 | 0.41221 | 52.1419          | 69.251     | 98.7           | 90.9             | 7.8              |
| HPSI_080 | agg3    | 422060271 | 280       | 306     | 26   | 14773535  | 32194574 | 0.41243 | 51.6786          | 69.2433    | 98.7           | 93               | 5.6              |
| HPSI_080 | default | 618548367 | 291       | 321     | 30   | 15599390  | 31108673 | 0.40997 | 64.5209          | 69.7973    | 98.9           | 42.3             | 56.6             |
| HSYZ_001 | agg1    | 412490973 | 43        | 56      | 13   | 20179189  | 23709579 | 0.40788 | 60.776           | 71.1286    | 98.5           | 91.4             | 7.1              |
| HSYZ_001 | agg2    | 406825562 | 42        | 56      | 14   | 20179189  | 23709579 | 0.40795 | 60.3202          | 71.0685    | 98.5           | 92.6             | 5.9              |
| HSYZ_001 | agg3    | 400281073 | 42        | 56      | 14   | 20179189  | 23709579 | 0.40802 | 59.8039          | 71.0636    | 98.5           | 93.8             | 4.7              |
| HSYZ_001 | default | 415619550 | 45        | 56      | 11   | 20179189  | 23709579 | 0.40788 | 60.9895          | 71.1614    | 98.5           | 90.5             | 8                |
| HSYZ_003 | agg1    | 440387253 | 142       | 156     | 14   | 20460975  | 29016061 | 0.40675 | 59.514           | 68.0678    | 96.9           | 85.9             | 11               |
| HSYZ_003 | agg2    | 433010379 | 131       | 148     | 17   | 19379737  | 29002929 | 0.40676 | 59.0551          | 67.9945    | 96.9           | 87.1             | 9.9              |
| HSYZ_003 | agg3    | 432674658 | 131       | 152     | 21   | 18574881  | 29002929 | 0.40673 | 59.0323          | 67.9911    | 96.9           | 87.1             | 9.9              |
| HSYZ_003 | default | 444164814 | 150       | 159     | 9    | 20460975  | 29016061 | 0.40677 | 59.78            | 68.1049    | 97             | 84.3             | 12.7             |



___



Run this, submitting with sample as positional:



```bash
#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=36
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

set -euo pipefail
# module load miniconda
# source activate puzzler200
JARFILE="${CONDA_PREFIX}/lib/python3.12/site-packages/HapHiC/utils/juicer_tools.1.9.9_jcuda.0.8.jar"
RUN_JUICERTOOLS="java -Xmx16G -jar ${JARFILE}"
FILTER_BAM="${CONDA_PREFIX}/lib/python3.12/site-packages/HapHiC/utils/filter_bam"

SAMPLE=${1:?Provide SAMPLE id as first argument}
HIFI=/project/coffea_pangenome/Guava/Reads_Concatenated/${SAMPLE}.HiFi.fastq.gz
HIC_R1=/project/coffea_pangenome/Guava/Reads_Concatenated/${SAMPLE}.HiC.R1.fastq.gz
HIC_R2=/project/coffea_pangenome/Guava/Reads_Concatenated/${SAMPLE}.HiC.R2.fastq.gz
WD=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies
REFERENCE=/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa
BUSCO_DB=/project/coffea_pangenome/Software/Merondun/busco_downloads/lineages/embryophyta_odb12
t=36

# ---- Working dirs ----
ASM_DIR=${WD}/${SAMPLE}/01b_hifiasm_hifi
TRIALS_DIR=${WD}/${SAMPLE}/02b_purge_trials

mkdir -p "${TRIALS_DIR}" "${ASM_DIR}"

# ---- First, assembly without HiC and ignore hom_cov ----
if [[ ! -f "${ASM_DIR}/asm.p_ctg.gfa" ]]; then
  echo "[INFO] Creating HiFi only assembly..."
  cd ${WD}/${SAMPLE}/01b_hifiasm_hifi
  hifiasm --primary -t ${t} -o asm ${HIFI}
else
  echo "[INFO] Reusing existing asm.p_ctg.gfa."
fi

DRAFT="${WD}/${SAMPLE}/02b_purge_trials/p_ctg.purged.fa"

cd ${TRIALS_DIR}

echo -e "\033[43m~~~ Starting Purge_Dups parameter sweep for ${SAMPLE} ~~~\033[0m"

# ---- Coverage stats and cutoffs ----
if [[ ! -f PB.base.cov || ! -f cutoffs ]]; then
  echo "[INFO] Computing coverage & cutoffs..."
  awk '/^S/{print ">"$2;print $3}' ${ASM_DIR}/asm.p_ctg.gfa > ${DRAFT}
  minimap2 -x map-hifi -t "${t}" ${DRAFT} ${HIFI} | gzip -c > hifi.paf.gz
  pbcstat hifi.paf.gz
  calcuts PB.stat > cutoffs
else
  echo "[INFO] Reusing existing PB.cov and cutoffs."
fi

# paths to reuse in the loop
COV="${TRIALS_DIR}/PB.base.cov"
CUT="${TRIALS_DIR}/cutoffs"
SELF_PAF_GZ=${TRIALS_DIR}/self.paf.gz

# ---- Ensure self-alignment (split + minimap2) exists ----
if [[ ! -f "${SELF_PAF_GZ}" ]]; then
  echo "[INFO] Self PAF not found; generating split and self-alignment..."
  split_fa "${DRAFT}" > pri.split.fa
  minimap2 -x asm5 -DP -t "${t}" pri.split.fa pri.split.fa | gzip -c > ${TRIALS_DIR}/self.paf.gz
else
  echo "[INFO] Reusing existing ${SELF_PAF_GZ}"
fi

# ---- derive 4 cutoff sets from the base cutoffs ----
CUTDIR="${TRIALS_DIR}/cutoff_sets"
mkdir -p "${CUTDIR}"

make_cutoff_file_aggressive() {
  local infile=$1
  local outfile=$2
  local c_scale=$3
  local d_scale=$4
  local e_scale=$5

  read -r a b c d e f < "${infile}"

  local c2 d2 e2
  c2=$(awk -v x="${c}" -v s="${c_scale}" 'BEGIN{printf "%d", x*s}')
  d2=$(awk -v x="${d}" -v s="${d_scale}" 'BEGIN{printf "%d", x*s}')
  e2=$(awk -v x="${e}" -v s="${e_scale}" 'BEGIN{printf "%d", x*s}')

  (( c2 < b )) && c2=$b
  (( d2 <= c2 )) && d2=$((c2+1))
  (( e2 <= d2 )) && e2=$((d2+1))
  (( e2 >= f )) && e2=$((f-1))

  echo "${a} ${b} ${c2} ${d2} ${e2} ${f}" > "${outfile}"
}

[[ -f "${CUTDIR}/default.cutoffs" ]] || cp "${CUT}" "${CUTDIR}/default.cutoffs"
[[ -f "${CUTDIR}/agg1.cutoffs" ]] || make_cutoff_file_aggressive "${CUT}" "${CUTDIR}/agg1.cutoffs" 1.10 1.10 0.90
[[ -f "${CUTDIR}/agg2.cutoffs" ]] || make_cutoff_file_aggressive "${CUT}" "${CUTDIR}/agg2.cutoffs" 1.20 1.20 0.85
[[ -f "${CUTDIR}/agg3.cutoffs" ]] || make_cutoff_file_aggressive "${CUT}" "${CUTDIR}/agg3.cutoffs" 1.25 1.25 0.80

run_purge_trial() {
  local label=$1
  local cutfile=$2
  local flags=$3

  local OUTDIR="${TRIALS_DIR}/${label}"
  mkdir -p "${OUTDIR}"
  cd "${OUTDIR}"

  echo "[INFO] purge_dups ${label}"

  [[ -f used.cutoffs ]] || cp "${cutfile}" used.cutoffs

  if [[ ! -f dups.bed ]]; then
    purge_dups ${flags} -c "${COV}" -T used.cutoffs "${SELF_PAF_GZ}" > dups.bed
  fi

  if [[ ! -f purged.fa || ! -f hap.fa ]]; then
    get_seqs dups.bed "${DRAFT}" > "get_seqs.${label}.out" 
    [[ -f purged.fa ]] || mv *.purged.fa purged.fa 2>/dev/null || true
    [[ -f hap.fa    ]] || mv *.hap.fa hap.fa 2>/dev/null || true
  fi
}

# Run all 4 
run_purge_trial default "${CUTDIR}/default.cutoffs" ""
run_purge_trial agg1   "${CUTDIR}/agg1.cutoffs"   ""
run_purge_trial agg2   "${CUTDIR}/agg2.cutoffs"   ""
run_purge_trial agg3   "${CUTDIR}/agg3.cutoffs"   ""

run_postpurge_eval() {
  local label=$1
  local PURGED_FA="${TRIALS_DIR}/${label}/purged.fa"

  local EVAL_DIR="${WD}/${SAMPLE}/03_eval_${label}"
  local HIC_DIR="${EVAL_DIR}/hic"
  local HAPHIC_DIR="${EVAL_DIR}/haphic"
  local JUICER_DIR="${EVAL_DIR}/juicer"
  local BUSCO_OUT="${EVAL_DIR}/busco"

  mkdir -p "${HIC_DIR}" "${HAPHIC_DIR}" "${JUICER_DIR}" "${BUSCO_OUT}"

  # ---- Hi-C alignment ----
  if [[ ! -f "${HIC_DIR}/filtered.bam" ]]; then
    echo "[INFO] Mapping HiC reads: ${label}"
    cd "${HIC_DIR}"
    bwa-mem2 index "${PURGED_FA}"

    bwa-mem2 mem -5SP -t ${t} "${PURGED_FA}" "${HIC_R1}" "${HIC_R2}" | \
      samblaster | \
      samtools view - -@ ${t} -S -h -b -F 3340 | \
      ${FILTER_BAM} - 1 --nm 3 --threads ${t} | \
      samtools view - -b -@ ${t} -o filtered.bam
  fi

  # ---- HapHiC ----
  if [[ ! -f "${HAPHIC_DIR}/haphic.complete" ]]; then
    echo "[INFO] Running HapHiC: ${label}"
    rm -rf "${HAPHIC_DIR}/run"
    mkdir -p "${HAPHIC_DIR}/run"
    cd "${HAPHIC_DIR}/run"

    set +euo pipefail
    haphic pipeline "${PURGED_FA}" "${HIC_DIR}/filtered.bam" "11" \
      --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t}
    set -euo pipefail

    touch "${HAPHIC_DIR}/haphic.complete"
    rm -rf 01.cluster 02.reassign 03.sort 
  fi

  # ---- Juicer / JBAT ----
  if [[ -f "${HAPHIC_DIR}/haphic.complete" && ! -f "${JUICER_DIR}/haphic_JBAT.hic" ]]; then
    echo "[INFO] Running Juicer: ${label}"
    cd "${JUICER_DIR}"

    mashmap -r "${REFERENCE}" -q "${PURGED_FA}" -t ${t} -s 10000 --perc_identity 85 \
      -o asm_to_ref.paf

    haphic refsort "${HAPHIC_DIR}/run/04.build/scaffolds.raw.agp" asm_to_ref.paf \
      > alignment.agp

    samtools faidx "${HAPHIC_DIR}/run/04.build/scaffolds.fa"
    samtools faidx "${PURGED_FA}"


    juicer pre \
              -a -q 1 \
              -o haphic_JBAT \
              "${HIC_DIR}/filtered.bam" \
              alignment.agp \
              "${PURGED_FA}.fai" > haphic_JBAT.log 2>&1

    grep PRE_C_SIZE haphic_JBAT.log | awk '{print $2" "$3}' > chrom.sizes

    ${RUN_JUICERTOOLS} pre \
      -r 5000000,4000000,3000000,2000000,1500000,1000000,750000,500000,250000,100000,50000 \
      haphic_JBAT.txt haphic_JBAT.hic chrom.sizes
  fi

  # ---- BUSCO ----
  if [[ ! -f "${BUSCO_OUT}/busco.complete" ]]; then
    echo "[INFO] Running BUSCO: ${label}"
    cd "${BUSCO_OUT}"

    busco -i "${PURGED_FA}" \
      -l "${BUSCO_DB}" \
      -m genome \
      -c ${t} \
      -o "${SAMPLE}.${label}" \
      -f \
      --offline

    cp "${BUSCO_OUT}/${SAMPLE}.${label}/short_summary.specific"*.txt \
       "${BUSCO_OUT}/${SAMPLE}.${label}.busco.txt"

    rm -rf \
      "${BUSCO_OUT}/${SAMPLE}.${label}/run_embryophyta_odb12/hmmer_output" \
      "${BUSCO_OUT}/${SAMPLE}.${label}/run_embryophyta_odb12/miniprot_output" \
      "${BUSCO_OUT}/${SAMPLE}.${label}/run_embryophyta_odb12/busco_sequences"

    touch "${BUSCO_OUT}/busco.complete"
  fi
}

for label in default agg1 agg2 agg3; do
  [[ -f "${TRIALS_DIR}/${label}/purged.fa" ]] || { echo "[ERROR] Missing purged.fa for ${label}"; exit 1; }
  run_postpurge_eval "${label}"
done

echo -e "label\tbusco_file" > "${TRIALS_DIR}/busco_summary_files.tsv"
for label in default agg1 agg2 agg3; do
  echo -e "${label}\t${WD}/${SAMPLE}/03_eval_${label}/busco/${SAMPLE}.${label}.busco.txt"
```

Search outputs for BUSCO: 

```bash
find . -name "*busco.txt" -exec grep -H "C:" {} \;
./03_eval_agg1/busco/HSYZ_003.agg1.busco.txt:   C:96.9%[S:85.9%,D:11.0%],F:0.8%,M:2.3%,n:2026,E:2.8%       
./03_eval_agg2/busco/HSYZ_003.agg2.busco.txt:   C:96.9%[S:87.1%,D:9.9%],F:0.8%,M:2.3%,n:2026,E:2.8%        
./03_eval_agg3/busco/HSYZ_003.agg3.busco.txt:   C:96.9%[S:87.1%,D:9.9%],F:0.8%,M:2.3%,n:2026,E:2.8%        
./03_eval_default/busco/HSYZ_003.default.busco.txt:     C:97.0%[S:84.3%,D:12.7%],F:0.8%,M:2.2%,n:2026,E:2.8%       
```

Size:

```bash
find . -name "*scaffolds.fa" -exec seqkit stats {} \;
file                                             format  type  num_seqs      sum_len  min_len      avg_len     max_len
./03_eval_agg1/haphic/run/04.build/scaffolds.fa  FASTA   DNA        101  440,391,353   17,894  4,360,310.4  56,435,745
file                                             format  type  num_seqs      sum_len  min_len      avg_len     max_len
./03_eval_agg2/haphic/run/04.build/scaffolds.fa  FASTA   DNA        101  433,013,379   17,894  4,287,261.2  63,489,989
file                                             format  type  num_seqs      sum_len  min_len      avg_len     max_len
./03_eval_agg3/haphic/run/04.build/scaffolds.fa  FASTA   DNA        101  432,677,658   17,894  4,283,937.2  63,463,210
file                                                format  type  num_seqs      sum_len  min_len      avg_len     max_len
./03_eval_default/haphic/run/04.build/scaffolds.fa  FASTA   DNA        101  444,169,714   17,894  4,397,719.9  56,435,745
```

Probably agg2 or agg3. Check juicer:

```bash
mkdir -p juicerpurge

for d in 03_*; do
    if [ -d "$d/juicer" ]; then
        id=$(echo "$d" | sed 's/03_eval_//')
        
        # Copy .hic
        if [ -f "$d/juicer/haphic_JBAT.hic" ]; then
            cp "$d/juicer/haphic_JBAT.hic" "juicerpurge/${id}.haphic_JBAT.hic"
        fi

        # Copy .assembly
        if [ -f "$d/juicer/haphic_JBAT.assembly" ]; then
            cp "$d/juicer/haphic_JBAT.assembly" "juicerpurge/${id}.haphic_JBAT.assembly"
        fi
    fi
done
```

Assess coverage / completeness

```R
setwd('/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/polypurge')
library(tidyverse)
library(readr)
library(stringr)
library(patchwork)
library(viridis)
library(ggtext)

t <- read_tsv('20260514_PurgingStats.tsv')
t$Label <- factor(t$Label,levels=c('default','agg1','agg2','agg3'))
t %>% 
  pivot_longer(c(Merqury_Complete,BUSCO_duplicated,SizeBP)) %>% 
  ggplot(aes(x=Sample,y=value,col=Label))+
  scale_color_manual(values=viridis(4))+
  geom_point()+
  facet_grid(name~.,scales='free')+
  theme_bw()

# coverage
c <- read_tsv('20260514_CovStats.tsv')
tot <- c %>% filter(contig == 'total') %>% select(sample,label,total = mean_cov)
c %>% 
  left_join(tot) %>% 
  filter(length > 5e6) %>% 
  ggplot(aes(x=label,y=mean_cov,fill=label))+
  geom_point(data=tot,aes(x=label,y=total),pch=8,size=3)+
  geom_boxplot()+
  facet_grid(sample~.,scales='free')+
  theme_bw()
```



