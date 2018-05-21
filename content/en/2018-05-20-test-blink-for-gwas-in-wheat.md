---
title: Test BLINK for GWAS in wheat
author: Junli Zhang
date: '2018-05-20'
slug: test-blink-for-gwas-in-wheat
categories:
  - Research
tags:
  - GWAS
  - BLINK
---

I have been using [Tassel](http://www.maizegenetics.net/tassel) for GWAS in wheat for several years. I tried [GAPIT](http://www.zzlab.net/GAPIT/index.html) too before, but for some reason, GAPIT is quite slow when I have over 20 K markers and 260 lines. So most of the time, I still use Tassel. In recent years, new GWAS methods have been reported, but I have not had a chance to try them yet. [Miao et al (2008)](https://www.biorxiv.org/content/early/2018/04/29/310391) compared different GWAS software/methods on the trait-variant association in crops under varying genetic architectures and found that FarmCPU provides the greatest statistical power for moderately complex traits. A few people also recommended FarmCPU to me. So I went to the [FarmCPU webiste](http://www.zzlab.net/FarmCPU/index.html), but I found Zhiwu's group has developed a new software called ["BLINK"](http://www.zzlab.net/blink/index.html), which supposed to outperform FarmCPU in both speed and accuracy. I tested it with my wheat 90K SNP data and kernel shape data, but I did not get any significant SNPs with the message "The signal from LM is too weak!". I read the [preprint paper of BLINK](https://www.biorxiv.org/content/early/2017/11/30/227249) and found that BLINK only allows examination of SNPs with P values of 0.01 after a multiple test correction. I think that is why BLINK is not working in wheat. I used the R version of BLINK, in which I could change the primary criteria, but I am not confident about this modification. I will contact the authors to confirm.

Wheat is a highly inbred species with very long LD blocks, but most of these software developed for GWAS are for maize, an outcrossing crop with a very short LD decay distance. And a lot of GWAS are also reported in outcrossing species. When I first touched GWAS in wheat, using the Bonferroni criteria or the FDR criteria, I could not find any significant markers. Later, after discussing with my wheat colleagues, I think the problem is that wheat has a long distance of LD decay, which is good for such a huge genome species because I can use a small number of markers to find MTAs, but the bad thing is that the associations are usually not that significant. What I do now is to see consistency rather than just highly significant P values. A few steps for GWAS in my research:

1. Estimate the overall LD decay distance and critical `$r^2$`. In wheat, the average LD is about 10 cM depending on the populations. In [our latest GWAS paper](https://doi.org/10.1007/s00122-018-3111-9) with a collection of 262 elite spring wheat materials from North America, the LDs ranged from 0.55 cM to 12.8 cM. In wheat, 1 cM is about 1 Mb on average.

2. Estimate the number of highly informative markers based on the critical `$r^2$` value using Haploview. In our latest GWAS paper, we used 22,226 SNPs, but there are only 1090 highly informative markers. So if we use Bonferroni correction at `$\alpha = 0.1$`, the significant raw P value should be `$0.1/1090$`, about `$1e-4$`, which is the P value we used for "highly significant SNPs".

3. Significant SNPs were identified using raw P values based on consistency across environments, since we have at least 3 trials for field experiments. In our latest GWAS paper, we defined "significant SNPs" as those that have raw P values < 0.05 in at least 3 environments and raw P values < 0.01 in at least 1 environment (we have 4 to 10 trials for different traits).

4. To better reduce the false positives, we also validated the significant SNPs using 8 nested association mapping populations (NAMs), which also allows us to have donor materials for gene clonings and breeding.