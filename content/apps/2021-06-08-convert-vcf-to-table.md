---
title: Convert a vcf file to SNP table
author: Junli Zhang
date: '2021-06-08'
slug: vcf-to-snp
categories:
  - Method
tags:
  - vcf
  - SNP
---

This is a simple tool to convert a vcf file to SNP table. 

<label for="infile">Select a .vcf file (NOT compressed gz file)</label>:<br>
<input type="file" id="infile"><br><br>
<input type="button" id="start" value="Start Converting">

<button id="download-btn" onclick="download()">Download converted SNPs</button><br>
<output id="output" style="display:none"></output>

<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/vcf2snp.js"></script>

### Help
Output includes:

> **CHROM**: chromosome  
> **POS**: position  
> **REF**: reference allele  
> **ALT**: alternative alleles (mutations compared to the reference)  
> **QUAL**: variant quality (the more reads and more counts of ALT alleles, the higher quality)  
> **AC**: "Allele count in genotypes, for each ALT allele, in the same order as listed"  
> **AN**: "Total number of alleles in called genotypes"; AN/2 = number of lines with data  
> **DP**: total reads; "Approximate read depth; some reads may have been filtered"  
> **MQ**: Mapping Quality

Then genotying data of each line: **N** is missing and **H** is heterozygous.