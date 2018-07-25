---
title: Write another script for designing sgRNAs in wheat
author: Junli Zhang
date: '2018-06-03'
slug: write-another-script-for-designing-sgrnas-in-wheat
categories:
  - Research
tags:
  - record
---

Yanpeng talked to me last week and mentioned that he hopes I could write a script to help design gRNA selection in wheat. There is one website, [CRISPR-P 2.0](http://crispr.hzau.edu.cn/CRISPR2/), for designing gRNAs in plants, but it does not include wheat. Beyond that, wheat is a hexaploid wheat and usually, we want to design gRNAs that can target all 3 homeologs in wheat's A, B, D genomes. So even the wheat genome is publicly available, the CRISPR-P 2.0 may still need some big update to handle polyploid species.

The script is now finished, and it is called [CRISPR-wheat](https://github.com/pinbo/CRISPR-wheat), which is available in my [github](https://github.com/pinbo). I used the same software, [BatMis](https://code.google.com/archive/p/batmis/), to find potential mismatches in the wheat genome. After working on CRISPR sgRNAs for a while, I think you can just use  [CRISPR-P 2.0](http://crispr.hzau.edu.cn/CRISPR2/) if they included wheat genome. The website does not need to update others, because they can use the gene on one of the A, B, D genome of wheat to find sgRNAs, and then you can check the hits on the genome to know whether they hit the other two sub-genomes.

More details of [CRISPR-wheat](https://github.com/pinbo/CRISPR-wheat) is in the README. I also listed detailed steps on how to prepare all the needed files. Hope it helps.