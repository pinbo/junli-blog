---
title: Maybe a trick of designing CAPS for backcrossing
author: Junli Zhang
date: '2018-03-23'
slug: maybe-a-trick-of-designing-caps-for-backcrossing
categories:
  - Research
tags:
  - Primer design
---

Today I was screening BC2F1 plants for further backcrossing. I tested about 10 plants and the majority of them are heterozygous, so I wonder whether it is because of incomplete digestion, but then I found the control (the gene donor) was digested completely, so I know the donor allele is the one that gets digested. If I can see the lower band, then it means it has the gene, even under incomplete digestion! Therefore I get one trick:

**When using CAPS/dCAPS marker in marker-assisted backcrossing, make sure the donor allele is the one that gets digested!**

![CAPS](/images/20180323_CAPS.svg)

In the figure above, L1 to L3 are some BCnF1 plants, P1 is the recurrent parent, and P2 is the donor plant. Design your CAPS marker to make the donor allele get digested.

### Another trick from Jorge
Jorge mentioned another trick on designing CAPS/dCAPS primer: if possible, include an extra common cut site in the PCR product to use as an internal control for the well, so you know the digestion in that well works. For example in the figure below: the 3rd band is a control band, all the wells should have this band!
![CAPS trick2](/images/20180323_CAPS_trick2.svg)

### Another trick on dCAPS primers
I prefer agarose gel to polyacrylamide gel electrophoresis (PAGE). But dCAPS primer usually only has 20 bp long, so is the fragment length difference. To make the difference longer, we need to make our primer longer. Two ways:

1. Make the primer itself longer (30 bp for example), but this only works when the GC content is low near the SNP.

1. Manually add a tail (for example the FAM tail) and this will increase 20 more bp difference, which will be easily separated on a 2% agarose gel.