# 02_genome_assembly/

Primary collapsed assembly proceeds with [puzzler v2.0.1](https://github.com/merondun/puzzler) from HiFi and HiC. 

Primary outputs being a genome, of course, which we can summarize with HiC contact maps and merqury kmer histograms:



![contacts_merq](/imgs/20260603_contacts_merq.png)



Output assembly statistics:



| Sample   | Order | Cultivar            | Species                     | Ploidy | SizeBP   | Sequences | Contigs | Gaps | ContigN50 | ScafN50  | GC      | Merqury_Complete | Merqury_QV | WithinChrsBP | PropWithinChrs | GapsInChrs | BUSCO_Complete | BUSCO_singlecopy | FCS_Contaminants | Contaminant_BP |
| -------- | ----- | ------------------- | --------------------------- | ------ | -------- | --------- | ------- | ---- | --------- | -------- | ------- | ---------------- | ---------- | ------------ | -------------- | ---------- | -------------- | ---------------- | ---------------- | -------------- |
| HPSI_003 | 1     | INDONESIAN SEEDLESS | Psidium guajava             | 3      | 4.19E+08 | 141       | 216     | 75   | 10653233  | 39209275 | 0.39663 | 90.2215          | 65.3423    | 4.13E+08     | 0.9862         | 75         | 98.6           | 96.7             | 1                | 22629          |
| HPSI_027 | 2     | J.B. WHITE          | Psidium guajava             | 2      | 4.12E+08 | 31        | 97      | 66   | 11308119  | 37557431 | 0.39483 | 99.2097          | 71.3273    | 4.1E+08      | 0.9967         | 66         | 98.7           | 96.8             | 297              | 20678787       |
| HPSI_060 | 3     | KLOM TOONKLAO       | Psidium guajava             | 2      | 4.22E+08 | 124       | 212     | 88   | 10090123  | 38175560 | 0.39675 | 96.3477          | 64.4828    | 4.15E+08     | 0.983          | 86         | 98.7           | 96.7             | 179              | 11568321       |
| HPSI_041 | 4     | THAI MAROON         | Psidium guajava             | 2      | 4.17E+08 | 19        | 90      | 71   | 11489848  | 38552022 | 0.39453 | 97.5031          | 69.6728    | 4.16E+08     | 0.9976         | 71         | 98.7           | 96.9             | 172              | 25850705       |
| HPSI_019 | 5     | KONA 1              | Psidium guajava             | 2      | 4.22E+08 | 122       | 181     | 59   | 12687117  | 38950305 | 0.39659 | 91.2954          | 66.4767    | 4.16E+08     | 0.9857         | 59         | 98.7           | 96.9             | 4                | 125527         |
| HPSI_007 | 6     | PINK ACID           | Psidium guajava             | 2      | 4.21E+08 | 96        | 172     | 76   | 14773278  | 38146420 | 0.39622 | 87.1893          | 67.9077    | 4.18E+08     | 0.9913         | 76         | 98.1           | 96.3             | 8                | 323892         |
| HPSI_068 | 7     | GEMA DE DORO        | Psidium guajava             | 2      | 4E+08    | 51        | 129     | 78   | 7626385   | 37322223 | 0.39492 | 99.4289          | 69.0656    | 3.98E+08     | 0.9934         | 77         | 98.6           | 96.7             | 358              | 9720875        |
| HPSI_065 | 8     | THAILAND            | Psidium guajava             | 2      | 4.06E+08 | 64        | 142     | 78   | 10404437  | 37400934 | 0.3951  | 97.3133          | 69.2993    | 4.04E+08     | 0.995          | 78         | 98.5           | 96.7             | 286              | 5558207        |
| HPSI_037 | 9     | BEAUMONT            | Psidium guajava             | 2      | 4.18E+08 | 136       | 195     | 59   | 10921675  | 37132908 | 0.39797 | 93.2802          | 63.5379    | 4.1E+08      | 0.9817         | 59         | 98.7           | 96.9             | 346              | 28858006       |
| HPSI_035 | 10    | KA HUA KULA         | Psidium guajava             | 2      | 4.19E+08 | 342       | 392     | 50   | 9364672   | 37251746 | 0.39826 | 96.0124          | 66.2485    | 4.08E+08     | 0.9758         | 50         | 98.7           | 96.8             | 1809             | 49557435       |
| HPSI_016 | 11    | PUERTO RICO 2       | Psidium guajava             | 2      | 4.16E+08 | 141       | 193     | 52   | 12671186  | 37528474 | 0.39786 | 96.9959          | 67.3806    | 4.1E+08      | 0.9846         | 50         | 98.8           | 96.9             | 68               | 3669622        |
| HPSI_080 | 12    | ARRAYAN FRUIT       | Psidium sartorianum         | 4      | 3.98E+08 | 29        | 68      | 39   | 15045748  | 36498344 | 0.40323 | 50.8278          | 72.5602    | 3.97E+08     | 0.9953         | 38         | 98.6           | 95.6             | 250              | 18101487       |
| HPSI_059 | 13    | COSTA RICAN GUAVA   | Psidium friedrichsthalianum | 4      | 4.33E+08 | 22        | 127     | 105  | 10030556  | 38360436 | 0.40324 | 59.3255          | 70.6442    | 4.31E+08     | 0.9959         | 104        | 98.6           | 92.8             | 42               | 1735417        |
| HPSI_010 | 14    | AMAZON GUAVA        | Psidium acutangulum         | 4      | 4.15E+08 | 20        | 64      | 44   | 23152470  | 38132351 | 0.405   | 51.1515          | 70.9614    | 4.15E+08     | 0.9983         | 44         | 97.7           | 95.5             | 33               | 1327083        |
| HPSI_072 | 15    | PUERTO RICAN GUAVA  | Psidium microphyllum        | 4      | 3.99E+08 | 14        | 72      | 58   | 16156023  | 36228238 | 0.3992  | 43.7695          | 69.6878    | 3.98E+08     | 0.999          | 58         | 97.1           | 94.4             | 15               | 1189782        |
| HPSI_069 | 16    | PINEAPPLE GUAVA     | Feijoa sellowiana           | 2      | 3.22E+08 | 257       | 281     | 24   | 20123824  | 28504751 | 0.40464 | 91.4559          | 63.1786    | 3.11E+08     | 0.9662         | 24         | 98.8           | 96.4             | 239              | 27806303       |
| HEUG_001 | 17    | ARAZA FRUIT         | Eugenia stipitata           | 2      | 3.28E+08 | 333       | 391     | 58   | 12310540  | 29282433 | 0.39906 | 88.2646          | 62.2609    | 3.16E+08     | 0.9636         | 58         | 98.7           | 97.3             | 62               | 6590391        |
| HSYZ_002 | 18    | WATER APPLE         | Syzygium aqueum             | 2      | 3.86E+08 | 41        | 85      | 44   | 18257371  | 35742225 | 0.40971 | 99.2297          | 67.8222    | 3.83E+08     | 0.9925         | 43         | 99.3           | 96.8             | 0                | 0              |
| HSYZ_001 | 19    | WAX APPLE           | Syzygium samarangense       | 4      | 3.93E+08 | 15        | 55      | 40   | 20179189  | 40382743 | 0.40796 | 59.5331          | 70.982     | 3.92E+08     | 0.9973         | 40         | 98.5           | 95.1             | 0                | 0              |
| HSYZ_003 | 20    | WAX APPLE           | Syzygium samarangense       | 4      | 4.19E+08 | 33        | 80      | 47   | 18574881  | 35938324 | 0.40773 | 58.1989          | 70.447     | 4.17E+08     | 0.996          | 43         | 96.8           | 89               | 77               | 4116427        |

Puzzler inputs `samples.tsv` 

| sample   | runtime | container | wd                                                           | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs                                                     | reference | hom_cov                                         | fcs_database | fcs_taxid                                                    | busco_lineage | busco_database |
| -------- | ------- | --------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | --------- | ----------------------------------------------- | ------------ | ------------------------------------------------------------ | ------------- | -------------- |
| HPSI_069 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_069.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_069.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_069.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 98        | /project/coffea_pangenome/Software/Merondun/fcs | 260130       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HEUG_001 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HEUG_001.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HEUG_001.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HEUG_001.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 98        | /project/coffea_pangenome/Software/Merondun/fcs | 1453397      | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_007 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_007.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_007.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_007.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 62        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_016 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_016.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_016.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_016.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 50        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_019 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_019.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_019.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_019.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 48        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_027 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_027.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_027.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_027.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 88        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_035 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_035.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_035.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_035.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 42        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_037 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_037.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_037.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_037.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 100       | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_041 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_041.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_041.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_041.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 66        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_060 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_060.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_060.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_060.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 75        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_065 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_065.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_065.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_065.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 42        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_068 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_068.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_068.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_068.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 42        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HSYZ_002 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_002.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_002.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_002.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 58        | /project/coffea_pangenome/Software/Merondun/fcs | 219865       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_003 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_003.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_003.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_003.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 69        | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_010 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_010.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_010.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_010.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 84        | /project/coffea_pangenome/Software/Merondun/fcs | 1835410      | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_059 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_059.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_059.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_059.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 124       | /project/coffea_pangenome/Software/Merondun/fcs | 681482       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_072 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_072.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_072.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_072.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 104       | /project/coffea_pangenome/Software/Merondun/fcs | 120290       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HPSI_080 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_080.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_080.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HPSI_080.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 88        | /project/coffea_pangenome/Software/Merondun/fcs | 1924220      | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HSYZ_001 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_001.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_001.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_001.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 44        | /project/coffea_pangenome/Software/Merondun/fcs | 260143       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |
| HSYZ_003 | conda   | NA        | /project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_003.HiFi.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_003.HiC.R1.fastq.gz | /project/coffea_pangenome/Guava/Reads_Concatenated/HSYZ_003.HiC.R2.fastq.gz | 11/project/coffea_pangenome/Guava/Assemblies/GCA_016432845.1_guava_v11.23_genomic.chrs.fa | 88        | /project/coffea_pangenome/Software/Merondun/fcs | 260143       | embryophyta_odb12/project/coffea_pangenome/Software/Merondun/busco_downloads |               |                |

Submit:

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --cpus-per-task=16
#SBATCH --mem=512Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=16

#module load miniconda
#source activate puzzler201

SAMPLE=$1

puzzler -s ${SAMPLE} -m samples.tsv --threads 16 --mem 512
```



## Remove haplotype chromosomes

```
cd primary_asm
mkdir -p init
mv HPSI_010* init/
cd init
egrep -v 'chr10_unloc|chr11_unloc' HPSI_010.fa.fai | cut -f1 > HPSI_010.keep.list
samtools faidx HPSI_010.fa $(cat HPSI_010.keep.list) > HPSI_010.nohaps.fa

# also compare the chrs... 
for CHR in chr10 chr11; do 
    echo "Estimating distance for ${CHR}"
    samtools faidx HPSI_010.fa ${CHR} > ${CHR}.hap.fa
    samtools faidx HPSI_010.fa ${CHR}_unloc1 > ${CHR}_unloc1.hap.fa
done 


# Sketch
files=$(ls *.hap.fa | tr '\n' ' ')
echo "Calculating distances for ${SAMPLE} using files: ${files}"
mash sketch -o HPSI_010 ${files}

# Then calculate all pairwise distances
mash dist HPSI_010.msh HPSI_010.msh > HPSI_010_pairwise_distances.txt
sed 's/.hap.fa//g' HPSI_010_pairwise_distances.txt \
    | awk '$1 != $2 {
        pair = ($1 < $2 ? $1 FS $2 : $2 FS $1)
        if (!seen[pair]++) print
    }'
chr10_unloc1    chr10   0.0877506       0       86/1000
chr11   chr10   0.140258        3.47813e-94     27/1000
chr11_unloc1    chr10   0.191731        2.13941e-31     9/1000
chr11   chr10_unloc1    0.197292        1.68092e-29     8/1000
chr11_unloc1    chr10_unloc1    0.197292        1.64542e-30     8/1000
chr11_unloc1    chr11   0.0662649       0       142/1000
```

### HSYZ_001 Haplotypes

Split haplotype contigs to maximize unique sequence and contacts, while minimizing haplotypic redundancies. 

After re-working the scaffolds, do this:

```
module load miniconda
source activate puzzler201

SAMPLE=HSYZ_001
WD=/project/coffea_pangenome/Guava/Assemblies/20250101_JustinAssemblies
MAP_FILE=tetrasamps.tsv
t=4
MEM=32
IFS=$'\t,' read -r _ RUNTIME CONTAINER WD HIFI HIC_R1 HIC_R2 NUM_CHRS REFERENCE HOM_COV FCS_DB FCS_TAXID BUSCO_LINEAGE BUSCO_DB < <(
    awk -F'[\t,]' -v sample="$SAMPLE" '$1 == sample {print $0}' "${MAP_FILE}"
)


##########
# If juicebox review file exists and reference provided, get ready for chromosome naming 
mkdir -p ${WD}/${SAMPLE}/05_postjuicebox
echo -e "\033[43m~~~ Extracting post-curation assembly and mapping to reference for ${SAMPLE} ~~~\033[0m"
cd ${WD}/${SAMPLE}/05_postjuicebox
${PUZZLER} juicer post \
    -o haphic-post_JBAT \
    ${WD}/juicer_files/${SAMPLE}_JBAT.review.assembly \
    ${WD}/${SAMPLE}/04_juicer/haphic_JBAT.liftover.agp \
    ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa 2> ${SAMPLE}.juicer.post.log
mv haphic-post_JBAT.FINAL.fa ${WD}/${SAMPLE}/05_postjuicebox/post_juicer_asm.fa
# Renaming: to reference chromosomes 
${PUZZLER} mashmap -r ${REFERENCE} -q post_juicer_asm.fa -t ${t} -s 10000 --perc_identity 85 -o asmpost_to_paf.paf 2> mashmap.postjuicer.log
${PUZZLER} samtools faidx post_juicer_asm.fa
${PUZZLER} samtools faidx ${REFERENCE}
${PUZZLER} map_chromosomes --paf asmpost_to_paf.paf --fai post_juicer_asm.fa.fai --out map.txt --min_size 0.1 &> mapping_renaming.log

#########
echo -e "\033[43m~~~ Renaming chromosomes for ${SAMPLE} ~~~\033[0m"

cd ${WD}/${SAMPLE}/05_postjuicebox
# This file has $haphic_scafID \t $ref_chr_ID \t $strand: ONLY WORKS IF THE SCAFFOLDS IN REF ARE NAMED E.G. >chr1 or >Chr1 
awk '{OFS="\t"}{print $1, $2, $5, $3, $4}' map.txt | grep -E 'chr|Chr' > chromosome_naming_map.txt 

# Exclude those chromsomes from map, and then add e.g. scaf101405 for the remaining non-chromosome scaffolds 
awk '{print $1}' chromosome_naming_map.txt > exclude_chr_scafIDs.txt
grep -vwf exclude_chr_scafIDs.txt post_juicer_asm.fa.fai | awk '{OFS="\t"}{print $1, "scaffold" NR, "+"}' >> chromosome_naming_map.txt

# Must be the same number of scaffolds
if [ $(cat chromosome_naming_map.txt | wc -l) -eq $(cat post_juicer_asm.fa.fai | wc -l) ]; then
    echo -e "\033[43m~~~ Scaffold sanity check passed for renaming, proceeding! ~~~\033[0m"
else
    echo -e "\033[41m~~~ Not same number of scaffolds, stop and inspect ${WD}/${SAMPLE}/05_postjuicebox/chromosome_naming_map.txt ~~~\033[0m"
    exit 1
fi # exit scaffold check 

# Orient chromosomes in the same strand direction 
awk '{if ($3 == "-") print $1}' chromosome_naming_map.txt > reverse_list.txt
awk '{if ($3 == "+") print $1}' chromosome_naming_map.txt > positive_list.txt
${PUZZLER} seqtk subseq post_juicer_asm.fa reverse_list.txt > to_revercomp.fa
${PUZZLER} seqtk subseq post_juicer_asm.fa positive_list.txt > positive.fa
${PUZZLER} seqkit seq --line-width 0 -t DNA -v -r -p to_revercomp.fa > revcomp.fa
cat positive.fa revcomp.fa > orient.fa
rm reverse_list.txt positive_list.txt to_revercomp.fa positive.fa revcomp.fa

# Detect if there are duplicate scaffolds
awk '{print $2}' chromosome_naming_map.txt | sort | uniq -d > duplicates.txt

##### THIS WILL DEAL WITH DUPLICATES, IN CASE THERE ARE 2 SCAFFOLDS CORRESPONDING TO A CHR! eg. Chr1_unloc1
if [ -s duplicates.txt ]; then

    echo -e "\033[41m~~~ Multiple scaffolds corresponding to single Chr for ${SAMPLE}, Renaming them e.g. Chr1_unloc1, Chr1_unloc2.. ~~~\033[0m"
    awk '
    {
        chr=$2
        len=$5
        data[chr][len] = $0
        lengths[chr][len] = len
        count[chr]++
    }
    END {
        PROCINFO["sorted_in"] = "@ind_num_asc"
        for (chr in data) {
            n = 0
            asorti(lengths[chr], sorted_lengths, "@val_num_desc")
            for (i in sorted_lengths) {
                # First one keeps chr as-is, subsequent ones get chr_unlocN
                label = (n == 0) ? chr : chr "_unloc" n
                split(data[chr][sorted_lengths[i]], line, "\t")
                line[2] = label
                print line[1], line[2], line[3], line[4], line[5]
                n++
            }
        }
    }
    ' OFS='\t'  chromosome_naming_map.txt > chromosome_renamed_map.txt
    cp chromosome_naming_map.txt chromosome_naming_map.ORIGINAL.txt
    mv chromosome_renamed_map.txt chromosome_naming_map.txt

else 

    echo -e "\033[43m~~~ Single scaffolds corresponding to a single Chr for ${SAMPLE} ~~~\033[0m"
fi # Exit duplicate scaffold check 

# Finally, rename the chromosomes
awk '$2 ~ /^[Cc]hr/ {print $0}' chromosome_naming_map.txt | sort -k2,2V | awk '{print $2}' > sorted_chr.txt
awk '$2 !~ /[Cc]hr/ {print $0}' chromosome_naming_map.txt | sort -k2,2V | awk '{print $2}' >> sorted_chr.txt

${PUZZLER} seqkit replace --line-width 0 -p "(.*)" -r "{kv}" -k chromosome_naming_map.txt orient.fa > haphic_renamed_unord.fa 2> seqkit_renaming.log
${PUZZLER} samtools faidx haphic_renamed_unord.fa
xargs ${PUZZLER} samtools faidx haphic_renamed_unord.fa < sorted_chr.txt > final_asm.fa
cp ${WD}/${SAMPLE}/05_postjuicebox/final_asm.fa ${WD}/primary_asm/${SAMPLE}.fa
${PUZZLER} samtools faidx ${WD}/primary_asm/${SAMPLE}.fa

count1=$(grep '>' haphic_renamed_unord.fa | wc -l)
count2=$(grep '>' final_asm.fa | wc -l)
count3=$(grep '>' post_juicer_asm.fa | wc -l)

if [ "$count1" -ne "$count2" ] || [ "$count1" -ne "$count3" ]; then
    echo -e "\033[41m~~~ Scaffold number mismatch between post_juicer_asm.fa and final_asm.fa in ${SAMPLE}/05_postjuicebox/ INSPECT ~~~\033[0m"
fi

######## exclude
mv final_asm.fa init_asm.fa
samtools faidx init_asm.fa
cut -f1 init_asm.fa.fai | grep -Ev 'chr3_unloc1$|chr10_unloc1$' > keep.list
samtools faidx init_asm.fa $(cat keep.list) > final_asm.fa
samtools faidx final_asm.fa
cp ${WD}/${SAMPLE}/05_postjuicebox/final_asm.fa ${WD}/primary_asm/${SAMPLE}.fa
```

### HSYZ_003 Haplotypes

```
SAMPLE=HSYZ_003
cd ${WD}/${SAMPLE}/05_postjuicebox

######## exclude
mv final_asm.fa init_asm.fa
samtools faidx init_asm.fa
cut -f1 init_asm.fa.fai | grep -Ev 'chr3_unloc1$' > keep.list
samtools faidx init_asm.fa $(cat keep.list) > final_asm.fa
samtools faidx final_asm.fa
cp ${WD}/${SAMPLE}/05_postjuicebox/final_asm.fa ${WD}/primary_asm/${SAMPLE}.fa
```

### HPSI_080 Haplotypes

```
SAMPLE=HPSI_080
cd ${WD}/${SAMPLE}/05_postjuicebox

######## exclude
mv final_asm.fa init_asm.fa
samtools faidx init_asm.fa
cut -f1 init_asm.fa.fai | grep -Ev 'chr8_unloc1$|chr9_unloc1$' > keep.list
samtools faidx init_asm.fa $(cat keep.list) > final_asm.fa
samtools faidx final_asm.fa
cp ${WD}/${SAMPLE}/05_postjuicebox/final_asm.fa ${WD}/primary_asm/${SAMPLE}.fa
```

## Merqury Histograms

- collect Merqury k‑mer copy‑number histograms for all accessions and combine into a unified dataset 

Just grab the merqury histos:

```bash
for i in $(cat Samples.list); do cp ${i}/08_merqury/merq.${i}.spectra-cn.hist merqury_spectra/;done
```

And plot:

```R
setwd('~/psidium_genomes/02_genome_assembly/merqury_spectra/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(viridis)
library(ggtext)

md <- read.table('../../samples.info',sep='\t',header = TRUE,comment.char = '') %>% as_tibble
grpcol <- md %>% distinct(Group, Color) %>% deframe()

hist_dat <- NULL
files <- list.files('.',pattern = '*hist')
for (file in files) {
  id = gsub('merq.','',gsub('.spectra-cn.hist','',file))
  cat('Processing: ',id,'\n')
  hist <- read_tsv(file) %>% mutate(Sample = id)
  hist_dat <- rbind(hist_dat,hist)
}

hm <- left_join(hist_dat,md) %>% arrange(Order)
hm$Sample <- factor(hm$Sample,levels=unique(md$Sample))

# limits
xlims <- hm %>% group_by(Sample,Group) %>% summarize(xmin=5,xmax=(HiFi_Gb/(SizeBP/1e9))+(HiFi_Gb/(SizeBP/1e9))*0.5) %>% distinct

hm2 <- hm %>%
  left_join(xlims, by = c("Sample", "Group")) %>%
  mutate(ID = paste0(Sample, "(", Group, "; ",round(Merqury_Complete,1),'%)')) %>% 
  filter(kmer_multiplicity >= xmin, kmer_multiplicity <= xmax ) %>% 
  mutate(Copies = factor(Copies,levels = c("read-only", "1"     ,    "2"        , "3" ,        "4"    ,     ">4")))
ord <- hm2 %>% arrange(Sample) %>%  distinct(Sample,ID)
hm2$ID <- factor(hm2$ID,levels=ord$ID)

hist_plot <- hm2 %>%
  ggplot(aes(x = kmer_multiplicity, y = Count, fill = Copies)) +
  geom_bar(stat='identity',alpha=0.6) +
  facet_wrap(~ ID, scales = "free", nrow = 5, ncol = 4) +
  scale_fill_manual(values = c('grey40',viridis(5))) +
  theme_bw(base_size=5) +
  labs(x='Coverage Peak',y='')+
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
hist_plot

ggsave('../20260603_Merqury_Histograms.pdf',hist_plot,height=4,width=6,dpi=300)

```



## Clean up for NCBI

NCBI identified some contam scaffolds, and mtDNA embedded within scaffold.

This will remove contaminated scaffolds entirely, listed in $SAMPLE.list:

```bash
for i in $(cat exclude_whole_scaffold.list); do 
    echo "exluding for ${i}"
    cut -f1 ../${i}.fa.fai | grep -vwf ${i}.list > ${i}.keep
    BEFORE=$(cat ../${i}.fa.fai | wc -l)
    TOTAL=$(cat ${i}.list | wc -l)
    samtools faidx ../${i}.fa $(cat ${i}.keep) > ${i}.fa
    samtools faidx ${i}.fa 
    AFTER=$(cat ${i}.fa.fai | wc -l)
    echo "for ${i} there are ${BEFORE} before and ${AFTER} after, removed ${TOTAL}"
done
```

And for masking within a scaffold:

From excel:

```bash
for i in $(cat mask_scaffold.list); do grep -w ${i} bad_regions.bed | awk '{OFS="\t"}{print $2, $3, $4}' > ${i}.bed; done 
```



```bash
for i in $(cat mask_scaffold.list); do 
    echo "masking for ${i}"
    bedtools maskfasta -fi ../${i}.fa -fo ${i}.fa -bed ${i}.bed

    awk '{OFS="\t"}{print $1, $2}' ../${i}.fa.fai > ${i}.pre.genome
    bedtools slop -i ${i}.bed -g ${i}.pre.genome -b 15 > ${i}.slop15.bed
    bedtools getfasta -fi ../${i}.fa -bed ${i}.slop15.bed -fo ${i}.before.fa
    bedtools getfasta -fi ${i}.fa -bed ${i}.slop15.bed -fo ${i}.after.fa
done
```

And rename the scaffolds, after preparing a lookup table:

```
cat sp_lookup.tsv 
HPSI_003        Psidium guajava
HPSI_027        Psidium guajava
HPSI_060        Psidium guajava
HPSI_041        Psidium guajava
HPSI_019        Psidium guajava
HPSI_007        Psidium guajava
HPSI_068        Psidium guajava
HPSI_065        Psidium guajava
HPSI_037        Psidium guajava
HPSI_035        Psidium guajava
HPSI_016        Psidium guajava
HPSI_080        Psidium sartorianum
HPSI_059        Psidium friedrichsthalianum
HPSI_010        Psidium acutangulum
HPSI_072        Psidium microphyllum
HPSI_069        Feijoa sellowiana
HEUG_001        Eugenia stipitata
HSYZ_002        Syzygium aqueum
HSYZ_001        Syzygium samarangense
HSYZ_003        Syzygium samarangense
```

Submit: 

```python
#!/usr/bin/env python3

import argparse
import os
import re
import sys


def read_lookup(path):
    lookup = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            sample = parts[0]
            organism = " ".join(parts[1:])
            lookup[sample] = organism
    return lookup


def infer_sample_from_filename(fasta_path):
    base = os.path.basename(fasta_path)
    sample = re.sub(r'\.(fa|fasta|fna)(\.gz)?$', '', base, flags=re.IGNORECASE)
    return sample


def sanitize_seqid(text):
    # Allowed: letters, digits, underscore, dot, hyphen, colon, asterisk, #.
    text = re.sub(r'[^A-Za-z0-9._\-:#*]', '_', text)
    return text


def parse_chr_name(header_core):
    """
    Convert Chr01 → 1, Chr1 → 1, chr28 → 28.
    Returns None if not chromosome primary scaffold.
    """
    m = re.fullmatch(r'Chr0*([1-9][0-9]*)', header_core, flags=re.IGNORECASE)
    if m:
        return str(int(m.group(1)))
    return None


def parse_unlocalized(header_core):
    """
    Detect unlocalized chromosome scaffolds.
    Examples:
       chr3_unloc1
       Chr05_unloc10
    Returns chromosome number as string, or None.
    """
    m = re.fullmatch(r'Chr0*([0-9]+)_unloc[0-9]+', header_core, flags=re.IGNORECASE)
    if m:
        return str(int(m.group(1)))
    return None


def rename_header(orig_header, sample, organism):
    header_core = orig_header.strip().split()[0]

    # 1) Primary chromosome sequence (Chr01, Chr1, etc.)
    chr_name = parse_chr_name(header_core)
    if chr_name is not None:
        seqid = sanitize_seqid(f"{sample}_{chr_name}")
        return (
            f">{seqid} "
            f"[organism={organism}] "
            f"[isolate={sample}] "
            f"[location=chromosome] "
            f"[chromosome={chr_name}]"
        )

    # 2) Unlocalized scaffold belonging to a chromosome: chr3_unloc1
    unloc_chr = parse_unlocalized(header_core)
    if unloc_chr is not None:
        seqid = sanitize_seqid(f"{sample}_{header_core}")
        return (
            f">{seqid} "
            f"[organism={organism}] "
            f"[isolate={sample}] "
            f"[location=chromosome] "
            f"[chromosome={unloc_chr}]"
        )

    # 3) Everything else = unplaced scaffold
    safe_core = sanitize_seqid(header_core)
    seqid = sanitize_seqid(f"{sample}_{safe_core}")
    return (
        f">{seqid} "
        f"[organism={organism}] "
        f"[isolate={sample}]"
    )


def process_fasta(infile, outfile, lookup, sample=None):
    if sample is None:
        sample = infer_sample_from_filename(infile)

    if sample not in lookup:
        sys.exit(f"ERROR: sample '{sample}' not found in lookup file")

    organism = lookup[sample]

    with open(infile) as fin, open(outfile, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                orig_header = line[1:].rstrip("\n")
                fout.write(rename_header(orig_header, sample, organism) + "\n")
            else:
                fout.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="Rename FASTA headers for NCBI genome submission."
    )
    parser.add_argument("-l", "--lookup", required=True,
                        help="Lookup table: sample_id <tab/space> organism name")
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA")
    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA")
    parser.add_argument("-s", "--sample", default=None,
                        help="Sample ID override; default inferred from FASTA filename")
    args = parser.parse_args()

    lookup = read_lookup(args.lookup)
    process_fasta(args.input, args.output, lookup, sample=args.sample)


if __name__ == "__main__":
    main()
```

It will rename your scaffolds like:

```bash
grep '>' HPSI_003.ncbi.fa | head 
>HPSI_003_1 [organism=Psidium guajava] [isolate=HPSI_003] [location=chromosome] [chromosome=1]
>HPSI_003_2 [organism=Psidium guajava] [isolate=HPSI_003] [location=chromosome] [chromosome=2]
>HPSI_003_3 [organism=Psidium guajava] [isolate=HPSI_003] [location=chromosome] [chromosome=3]
>HPSI_003_4 [organism=Psidium guajava] [isolate=HPSI_003] [location=chromosome] [chromosome=4]
>HPSI_003_5 [organism=Psidium guajava] [isolate=HPSI_003] [location=chromosome] [chromosome=5]
>HPSI_003_chr5_unloc1 [organism=Psidium guajava] [isolate=HPSI_003] [location=chromosome] [chromosome=5]
>HPSI_003_chr5_unloc2 [organism=Psidium guajava] [isolate=HPSI_003] [location=chromosome] [chromosome=5]
```



