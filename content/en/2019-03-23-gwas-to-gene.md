---
title: GWAS to gene
author: Junli Zhang
date: '2019-03-23'
slug: gwas-to-gene
categories: []
tags: []
---

I am a little disappointed in **Scientific Reports** after reading a newly published Genome-wide association study (GWAS) in wheat for grain size: [Genome-wide association study revealed that the *TaGW8* gene was associated with kernel size in Chinese bread wheat](https://www.nature.com/articles/s41598-019-38570-2). This study mapped a significant quantitative trait locus (QTL) for kernel length (KL) and thousand kernel weight (TKW) in a haplotype block of 28.67 Mb on chromosome 7BS. This block contains 247 annotated genes and one of them is a homolog of *OsGW8* (*TaGW8*). Since *OsGW8* affects grain weight in rice, the authors concluded that *TaGW8* is the gene responsible for the TKW QTL but without any validation. Later they characterized TaGW8 and showed that the two alleles of TaGW had significant differences for multiple traits in the association mapping panel. However, any genes in this block that are polymorphic would give you the same conclusion. If they just say *TaGW8* is the candidate gene, I would not have so many complaints.

### GWAS results need validation

GWAS gives us a fast and easy way to discover QTLs from collective populations. Although the QTL confidence region is smaller than a biparental population with similar size, on average the QTL confidence region is still around 10 Mb in wheat, where tens to hundreds of genes could be present. So it is almost impossible to get the gene directly from GWAS without fine mapping and validation.

One big concern of GWAS is false positives, so validation in other populations, especially biparental populations, is necessary. Biparental populations in which the QTL is validated can also be used directly for fine mapping and cloning of the gene.

### Current protocol we are using for gene cloning from GWAS results

1. Validate GWAS in biparental populations
2. Create near-isogenic lines (NILs): use BCxFn or heterogeneous inbred families (HIF) from heterozygous F5 or higher generations
3. Screen and evaluate recombinants from NILs or HIFs to reduce the confidence region
4. Predict candidate genes from a fined QTL confidence region according to gene annotations, polymorphisms, gene expression styles etc
5. Validate (or disprove) with mutants
6. Complementation test of mutants

Without mutant validation and complementation validation, you can only say that you get a candidate gene, but you cannot say you get the gene. Even this gene is involved in the pathway controlling the trait of interest, we cannot rule out that other genes in your fined confidence region are not involved at all. Actually, the validation of one candidate is only to prove it is controlling the trait, but it did not prove that other genes in the final confidence region (usually too small to break by recombinations) are not controlling the traits. So we also need evidence that the other genes in the final confidence region are unlikely controlling the traits. Always keep in mind that there might be more than one candidate gene underneath your QTL.

