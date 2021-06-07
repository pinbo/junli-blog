---
title: Use SnpEff to annotate wheat variations
author: Junli Zhang
date: '2021-06-07'
slug: use-snpeff-to-annotate-wheat-variations
categories:
  - Method
tags:
  - tutorial
  - vcf
  - snpEff
---

1. Install SnpEff based on the directions here: https://pcingola.github.io/SnpEff/download/
```sh
# Go to home dir
cd
# Download latest version
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# Unzip file
unzip snpEff_latest_core.zip
```
2. Download the wheat annotation database (hexaploid wheat, they also have tetraploid annotations)
```sh
cd  snpEff
java -jar snpEff.jar download Triticum_aestivum
```

3. Run snpEff
```sh
java -Xmx20g -jar /path/to/path/to/snpEff/snpEff.jar -v -stats summary.html Triticum_aestivum your.vcf.gz > your.Ta.ann.vcf
```

4. Convert vcf to an variant table with annotation
```sh
cat your.vcf.gz.Ta.ann.vcf | java -jar /path/to/snpEff/SnpSift.jar filter "( QUAL > 100 )" | /path/to/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /path/to/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT QUAL AC AN DP MQ "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].IMPACT" "ANN[*].HGVS_P" "GEN[*]" > your-ann-table.txt
```

5. Convert vcf alleles from numbers to ATGC
Get the python script [here](https://github.com/pinbo/myscripts/blob/master/Python/convert_vcf_calls_to_SNP_and_add_Blosum62_score.py).
```sh
python ./convert_vcf_calls_to_SNP_and_add_Blosum62_score.py your-ann-table.txt converted-your-ann-table.txt
```