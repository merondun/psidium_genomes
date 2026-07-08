# 12_flesh_candidates/

Identify potential flesh color candidates using the  guava pangenome. Under construction.



___



Using the pggb outputs:

```bash
VCF=chr2.fa.gz.8b44a75.ec291a3.44557c4.smooth.final.HPSI_007.vcf

grep -vx '^HPSI_007$' ../../pinkflesh.tsv | awk '{print $1"\tpink"}' > flesh.groups.tsv
awk '{print $1"\twhite"}' ../../whiteflesh.tsv >> flesh.groups.tsv

bcftools norm -m -any "$VCF" -Ou \
  | bcftools +fill-tags -Ou -- -S flesh.groups.tsv -t AC,AN,AF,NS \
  | bcftools view \
      -i 'AF_pink=0 && AF_white=1 && NS_pink=6 && NS_white=3' \
      -Oz -o fixed_white_vs_pink.pggb.vcf.gz
bcftools index fixed_white_vs_pink.pggb.vcf.gz

bcftools view fixed_white_vs_pink.pggb.vcf.gz \
  -i 'abs(strlen(ALT)-strlen(REF))>=50' \
  -Oz -o fixed_white_vs_pink.SVlike.pggb.vcf.gz
bcftools index fixed_white_vs_pink.SVlike.pggb.vcf.gz

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF_pink\t%AF_white\n' fixed_white_vs_pink.SVlike.pggb.vcf.gz \
| awk -v OFS='\t' '
{
  chrom=$1
  pos=$2
  id=$3
  ref=$4
  alt=$5
  afp=$6
  afw=$7

  ref_len=length(ref)
  alt_len=length(alt)
  start=pos-1

  if (ref_len >= alt_len) {
    end=start+ref_len
    svtype="DEL_in_white"
    svlen=alt_len-ref_len
  } else {
    end=start+1
    svtype="INS_in_white"
    svlen=alt_len-ref_len
  }

  name=chrom ":" pos ":" svtype ":" svlen
  print chrom,start,end,name,svtype,svlen,ref_len,alt_len,afp,afw,id
}' \
| sort -k1,1 -k2,2n > fixed_white_vs_pink.SVlike.bed
```

Afterwards, take the HPSI_007 annotation:

```bash
awk -v OFS='\t' '$3=="gene" {
  id="."
  if (match($9,/ID=([^;]+)/,m)) id=m[1]
  name=id
  if (match($9,/Name=([^;]+)/,n)) name=n[1]
  print $1,$4-1,$5,name,$7,$9
}' HPSI_007.gff3 \
| sort -k1,1 -k2,2n > HPSI_007.genes.bed
```

