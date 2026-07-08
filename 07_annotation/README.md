# 07_annotation/

Annotate using egapx with existing short read data, using the same species when available, or if not, a member of the same genus. 

* Assigns each genome a species‑matched short‑read dataset (via `lookup.tsv`) and uses eGAPx to generate transcriptome‑guided structural annotations per sample.

* Builds YAML configs on‑the‑fly, runs eGAPx under SLURM, then extracts final GTFs and gene counts for downstream comparative analyses.

These are the accessions I will use for annotation:

| Sample                | Run                                                          | Link                                                   | Tissue        |
| --------------------- | ------------------------------------------------------------ | ------------------------------------------------------ | ------------- |
| Psidium guajava       | [SRR34021006](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR34021006) | https://www.ncbi.nlm.nih.gov/sra/SRX29216563[accn]     | bark          |
| Psidium guajava       | [SRR26158591](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR26158591) | https://www.ncbi.nlm.nih.gov/sra/SRX21870914[accn]     | root          |
| Psidium guajava       | [SRR26158593](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR26158593) | https://www.ncbi.nlm.nih.gov/sra/SRX21870912[accn]     | mature leaf   |
| Psidium guajava       | [SRR26158595](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR26158595) | https://www.ncbi.nlm.nih.gov/sra/SRX21870910[accn]     | young leaf    |
| Psidium guajava       | [SRR26158598](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR26158598) | https://www.ncbi.nlm.nih.gov/sra/SRX21870907[accn]     | immature leaf |
| Psidium guajava       | [SRR26156677](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR26156677) | https://www.ncbi.nlm.nih.gov/sra/SRX21869069[accn]     | flower bud    |
| Psidium guajava       | [SRR11746930](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR11746930) | https://www.ncbi.nlm.nih.gov/sra/SRX8300260[accn]      | fruit         |
| Feijoa sellowiana     | [SRR26156476](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR26156476) | https://www.ncbi.nlm.nih.gov/sra/?term=Acca+sellowiana | leaf          |
| Eugenia uniflora      | [SRR500513](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR500513) | https://www.ncbi.nlm.nih.gov/sra/SRX149432[accn]       | unknown       |
| Eugenia uniflora      | [SRR29331828](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR29331828) | https://www.ncbi.nlm.nih.gov/sra/SRX24847979[accn]     | leaf          |
| Eugenia uniflora      | [SRR29331821](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR29331821) | https://www.ncbi.nlm.nih.gov/sra/SRX24847986[accn]     | unknown       |
| Syzygium samarangense | [SRR24058007](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR24058007) | https://www.ncbi.nlm.nih.gov/sra/SRX19859363[accn]     | leaf          |
| Syzygium samarangense | [SRR2029397](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR2029397) | https://www.ncbi.nlm.nih.gov/sra/SRX1029964[accn]      | fruit         |
| Syzygium aqueum       | [SRR24331231](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR24331231) | https://www.ncbi.nlm.nih.gov/sra/SRX20124125[accn]     | leaf          |
| Syzygium aqueum       | [SRR24331240](https://trace.ncbi.nlm.nih.gov/Traces?run=SRR24331240) | https://www.ncbi.nlm.nih.gov/sra/SRX20124116[accn]     | stem          |

And I will create a `lookup.tsv` file which assigns each sample one of these `$SAMPLE` and an ncbi taxon id:

| Sample   | Read File                   | TaxonID |
| -------- | --------------------------- | ------- |
| HPSI_003 | psidium.reads               | 120290  |
| HPSI_027 | psidium.reads               | 120290  |
| HPSI_060 | psidium.reads               | 120290  |
| HPSI_041 | psidium.reads               | 120290  |
| HPSI_019 | psidium.reads               | 120290  |
| HPSI_007 | psidium.reads               | 120290  |
| HPSI_068 | psidium.reads               | 120290  |
| HPSI_065 | psidium.reads               | 120290  |
| HPSI_037 | psidium.reads               | 120290  |
| HPSI_035 | psidium.reads               | 120290  |
| HPSI_016 | psidium.reads               | 120290  |
| HPSI_080 | psidium.reads               | 1924220 |
| HPSI_059 | psidium.reads               | 681482  |
| HPSI_010 | psidium.reads               | 1835410 |
| HPSI_072 | psidium.reads               | 120290  |
| HPSI_069 | feijoa.reads                | 260130  |
| HEUG_001 | eugenia.reads               | 1453397 |
| HSYZ_002 | syzygium_aqueum.reads       | 219865  |
| HSYZ_001 | syzygium_samarangense.reads | 260143  |
| HSYZ_003 | syzygium_samarangense.reads | 260143  |

Then submit this with e.g.: `sbatch 01_Annotation_eGAPx.sh HPSI_003` 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load egapx/0.5.2
module unload java
module load java/21

# Variables for genome to annotate, and which isoseq reads to use 
ID=${1:?Provide ID as first argument}

META=lookup.tsv
WD=/project/coffea_pangenome/Guava/annotation
GENOMES=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies/primary_asm
RUN="${ID}_srr"
echo "WORKING ON ${RUN}"

# Extract readfile + taxid
readfile=$(awk -v s="$ID" '$1==s{print $2}' "$META")
taxid=$(awk -v s="$ID" '$1==s{print $3}' "$META")

if [[ -z "$readfile" || -z "$taxid" ]]; then
    echo "ERROR: Could not find $ID in $META"
    exit 1
fi

echo "Using readfile: $readfile"
echo "Using taxid: $taxid"


# Extract SRRs from the readfile
if [[ ! -f "$readfile" ]]; then
    echo "ERROR: readfile $readfile not found"
    exit 1
fi

# Format SRR list as YAML array
SRR_BLOCK=$(sed 's/^/  - /' "$readfile")


# Create YAML
YAML=${WD}/${RUN}.yaml

cat > "$YAML" <<EOF
genome: ${GENOMES}/${ID}.fa
taxid: ${taxid}
short_reads:
${SRR_BLOCK}
cmsearch:
  enabled: false
trnascan:
  enabled: false
EOF

echo "YAML created at: $YAML"

# Run 
egapx.py ${YAML} -e slurm -w ${WD}/work/${RUN} -o ${WD}/out/${RUN}
```

Afterwards, go to the output folder and extract the gtfs, just count genes for sanity first:

```bash
mkdir -p ../gtfs
for SAMPLE in */; do ID=$(echo $SAMPLE | sed 's@_srr/@@g'); cp ${SAMPLE}/complete.genomic.gtf ../gtfs/${ID}.gtf; ngenes=$(awk '$3 == "gene"' ../gtfs/${ID}.gtf | wc -l); echo -e "${ID}\t${ngenes}"; done
```

Genes:

```
HEUG_001        22493
HPSI_003        25714
HPSI_007        25584
HPSI_010        26278
HPSI_016        25675
HPSI_019        25691
HPSI_027        25628
HPSI_035        25695
HPSI_037        25686
HPSI_041        25709
HPSI_059        27365
HPSI_060        25677
HPSI_065        25607
HPSI_068        25694
HPSI_069        25962
HPSI_072        26053
HPSI_080        26570
HSYZ_001        27796
HSYZ_002        27293
HSYZ_003        29084
```

