# 06_pangenome/

For the pangenome, we will first start at the chromosome-level and build a collapsed pangenome for **chrN** using the 12 collapsed guava assemblies. 

First, rename the chrs so that they include sample and chromosome IDs, following the [panSN spec](https://github.com/pangenome/PanSN-spec) for non-haplotype phased:

for `Chr01`, looks like:

```
zcat chr01.fa.gz | grep '>' | head
>H6#chr01
>HART001#chr01
>HART030#chr01
>HART032#Chr01
>HART033#Chr01
```

Submit:

```bash
#!/bin/bash
#SBATCH --time=0-05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

set -euo pipefail

indir=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm
outdir=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/pangenome/01_assemblies

mkdir -p "$outdir"

for chrnum in $(seq 1 11); do
    chr="chr${chrnum}"
    outfile="$outdir/${chr}.fa"
    : > "$outfile"

    for fa in "$indir"/*.fa; do
        sample=$(basename "$fa" .fa)

        awk -v chr="$chr" -v sample="$sample" '
            /^>/ {
                header=$0
                sub(/^>/, "", header)

                split(header, fields, /[[:space:]]+/)

                if (fields[1] == chr) {
                    print ">" sample "#" chr
                    keep=1
                } else {
                    keep=0
                }
                next
            }
            keep { print }
        ' "$fa" >> "$outfile"
    done

    bgzip -f "$outfile"
    samtools faidx "$outfile" 
done
```

## pggb Optimization

* Builds PGGB pangenomes across parameter grids (p, s). 

* Generates per‑run summary TSVs by extracting VCF and ODGI statistics.
* Plots variation across s and p. 

**Graph QC metrics** for memory in case I forget later: 

**length**
 • total bp in graph
 • too big = over‑alignment
 • too small = over‑collapse

**nodes**
 • fragmentation indicator
 • fewer = cleaner graph
 • many = graph shattered into tiny pieces

**edges**
 • graph connectivity
 • too many = tangled
 • too few = under‑aligned

------

Best graph is....

• **lowest fragmentation:** fewest nodes, largest avg node length (length/nodes)
 • **reasonable complexity:** edges stable, length not inflated
 • **stable variation:** SNP/indel/other counts plateau across runs

Total bp input for each chr using the index files:

```bash
# summarize total inputs
echo -e "chr\tinput_size" > chr_cumulative.tsv
for chr in $(ls *fai | sed 's/.fa.gz.fai//g'); do 
	size=$(awk '{sum += $2} END {print sum}' ${chr}.fa.gz.fai)
    echo -e "${chr}\t${size}" >> chr_cumulative.tsv
done 

```

Run:

```bash
#!/bin/bash

#SBATCH --time=5-00:00:00   
#SBATCH --cpus-per-task=24
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load miniconda
source activate pggb

CHR=${1:?Provide chr as first argument}
GDIR=/project/coffea_pangenome/Guava/pangenome/01_assemblies
[ -f ${GDIR}/${CHR}.fa.gz.fai ] || samtools faidx ${GDIR}/${CHR}.fa.gz

mkdir -p trials 

# -p, --map-pct-id PCT        percent identity for mapping/alignment [default: 90]
# -s, --segment-length N      segment length for mapping [default: 5000]
# -k, --min-match-len N       filter exact matches below this length [default: 23]

for P in 90 95 98; do 
	for S in 5000 25000 50000; do 
		for K in 23; do 

		OUTDIR="output_${CHR}_p${P}_s${S}_k${K}"
		id="${CHR}_p${P}_s${S}_k${K}"
		GFA=$(find ${OUTDIR} -name *gfa)
		if [[ -f ${GFA} ]]; then
			echo "Skipping PGGB build for ${id}, final GFA ${GFA} already exists."
		else
		echo "Building pangenome for -p ${P} and -s ${S} and -k ${K} for ${CHR}" 
		pggb build -i ${GDIR}/${CHR}.fa.gz \
			-o ${OUTDIR} \
			-t 48 \
			-p ${P} \
			-s ${S} \
			-k ${K} \
			--poa-params asm10 \
			-V 'HPSI_007:1000'
		fi

		TSV=trials/${CHR}_p${P}_s${S}_k${K}.tsv

		third_col=$(awk 'NR==2 {print $3}' "$TSV" 2>/dev/null)

		if [[ -n "$third_col" ]]; then
			echo "Skipping PGGB summary for ${id}, summary tsv ${TSV} already exists."
		else
			# header
			echo "Creating summary tsv for ${id}"
			echo -e "id\trecords\tSNPs\tindels\tothers\tts\ttrvs\tts_tv\tsize_og\tsize_gfa\tsize_vcf\tlength\tnodes\tedges\tpaths\tsteps" > trials/${CHR}_p${P}_s${S}_k${K}.tsv

			# VCF stats
			stats=$(find ${OUTDIR} -name ${CHR}*HPSI_007.vcf.stats)
			records=$(grep -P "\tnumber of records" "$stats" | awk '{print $6}')
			snps=$(grep -P "\tnumber of SNPs" "$stats" | awk '{print $6}')
			indels=$(grep -P "\tnumber of indels" "$stats" | awk '{print $6}')
			others=$(grep -P "\tnumber of others" "$stats" | awk '{print $6}')

			ts=$(grep "^TSTV" "$stats" | awk '{print $3}')
			tv=$(grep "^TSTV" "$stats" | awk '{print $4}')
			ts_tv=$(grep "^TSTV" "$stats" | awk '{print $5}')

			og=$(find ${OUTDIR} -name ${CHR}*.og)
			gfa=$(find ${OUTDIR} -name ${CHR}*.gfa)
			vcf=$(find ${OUTDIR} -name ${CHR}*.HPSI_007.vcf)
			size_og=$(stat -c%s "$og")
			size_gfa=$(stat -c%s "$gfa")
			size_vcf=$(stat -c%s "$vcf")

			# odgi stats
			odgi_stats=$(odgi stats -i "$og" -S | tail -n1)

			length=$(echo "$odgi_stats" | awk '{print $1}')
			nodes=$(echo "$odgi_stats" | awk '{print $2}')
			edges=$(echo "$odgi_stats" | awk '{print $3}')
			paths=$(echo "$odgi_stats" | awk '{print $4}')
			steps=$(echo "$odgi_stats" | awk '{print $5}')

			echo -e "${id}\t${records}\t${snps}\t${indels}\t${others}\t${ts}\t${tv}\t${ts_tv}\t${size_og}\t${size_gfa}\t${size_vcf}\t${length}\t${nodes}\t${edges}\t${paths}\t${steps}" >> trials/${CHR}_p${P}_s${S}_k${K}.tsv
			plot=$(find "${OUTDIR}" -name "${CHR}*smooth.final.og.viz_multiqc.png")
			if [[ -n "$plot" ]]; then
				cp "$plot" "trials/${CHR}_p${P}_s${S}_k${K}.smooth.final.og.viz_multiqc.png"
			else
				echo "Warning: No plot found in ${OUTDIR}"
			fi
		fi
		done 
	done 
done 

```

Which failed

```
for CHR in {1..11}; do
  for P in 90 95 98; do
    for S in 5000 25000 50000; do
      for K in 23; do
        echo "chr${CHR}_p${P}_s${S}_k${K}.smooth.final.og.viz_multiqc.png"
      done
    done
  done
done > expected.txt
ls *.smooth.final.og.viz_multiqc.png > actual.txt
grep -v -F -f actual.txt expected.txt
grep -v -F -f actual.txt expected.txt | sed 's/_.*//g' | sort | uniq  > ../RERUN.list
```



Summarize: `mergem *tsv > 20260706_pggb_optimization.tsv` 

## Plot Parameterization

```
setwd('/project/coffea_pangenome/Guava/pangenome')
library(tidyverse)
library(ggpubr)

l <- read_tsv('01_assemblies/chr_cumulative.tsv')
p <- read_tsv('02_pggb/trials/20260612_pggb_optimization.tsv') %>% 
  separate(id, into=c('chr','p','s','k'), remove =F ) %>% 
  mutate(p = as.numeric(gsub('p','',p)),
         s = as.numeric(gsub('s','',s)),
         k = as.numeric(gsub('k','',k)))
p
lp <- left_join(p,l)
lp <- lp %>% 
  mutate(compression = 1 - length / input_size,
         chrn = as.numeric(gsub('chr','',chr))) %>% 
  arrange(chrn) %>% 
  mutate(chr = factor(chr, levels = unique(chr)))

comp <- lp %>% 
  pivot_longer(c(compression,nodes)) %>% 
  ggplot(aes(x=value, y=chr, shape=as.factor(p), fill=s, size=as.factor(k)))+
  geom_point(position=position_jitter(height=0.2))+
  scale_shape_manual(values=c(21,22,24))+
  facet_grid(.~name,scales='free')+
  scale_fill_continuous(low='yellow',high='red')+
  theme_bw()
comp

# pareto
lp %>% 
  ggplot(aes(x=compression, y=nodes, shape=as.factor(p), fill=s, size=k))+
  geom_point(position=position_jitter(height=0.2))+
  scale_shape_manual(values=c(21,22,24))+
  facet_grid(.~chr,scales='free')+
  scale_fill_continuous(low='yellow',high='red')+
  theme_bw()

# heatmap
comp_heat <- lp %>% 
  ggplot(aes(x=s, y=p, fill=compression))+
  geom_tile()+
  facet_grid(.~chr,scales='free')+
  scale_fill_continuous(low='yellow',high='red')+
  theme_bw()
node_heat <- lp %>% 
  ggplot(aes(x=s, y=p, fill=nodes))+
  geom_tile()+
  facet_grid(.~chr,scales='free')+
  scale_fill_continuous(low='yellow',high='red')+
  theme_bw()
ggarrange(comp_heat, node_heat)

norm <- lp %>% 
  group_by(chr) %>% 
  mutate(min_node = min(nodes),
         norm_nodes = nodes / min_node,
         min_comp = min(compression),
         norm_comp = compression / min_comp)
norm %>% 
  pivot_longer(c(norm_comp, norm_nodes)) %>% 
  ggplot(aes(x=value, y=chr, shape=as.factor(p), fill=s, size=k))+
  geom_point(position=position_jitter(height=0.2))+
  scale_shape_manual(values=c(21,22,24))+
  facet_grid(.~name,scales='free')+
  scale_fill_continuous(low='yellow',high='red')+
  theme_bw()

```





