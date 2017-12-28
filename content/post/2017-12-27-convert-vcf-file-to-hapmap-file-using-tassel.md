---
title: Convert vcf file to hapmap file using Tassel
author: Junli Zhang
date: '2017-12-27'
slug: convert-vcf-file-to-hapmap-file-using-tassel
categories:
  - Research
tags:
  - Tips
---

I was working on exon capture data recently. I found Tassel could load vcf file and output hapmap file easily. If the vcf file has SnpEff prediction, it is also easy to extract fields of interest. I am not sure why some variations could have multiple effects. Maybe there are more than two variations among lines. Anyway, I just want to know what is the SNP effect. So I used [this command](http://snpeff.sourceforge.net/SnpSift.html#Extract) to extract all the variations and their SnpEff:

```sh
java -jar SnpSift.jar extractFields example_snpeff.vcf CHROM POS REF ALT  ANN[0].EFFECT ANN[0].BIOTYPE ANN[0].GENE ANN[0].GENEID ANN[0].HGVS_C ANN[0].HGVS_P > all_SNPeff.txt
```

Then I can combine the two files (the one from Tassel and the one from SnpSift) and get a complete genotying file with SnpEff.