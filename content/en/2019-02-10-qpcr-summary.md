---
title: qPCR summary
author: Junli Zhang
date: '2019-02-10'
slug: qpcr-summary
categories:
  - Research
tags:
  - qPCR
---

I just finished a few qPCR experiments. Here I just want to summarize some experiences and lessons before I forgot them.

## A few basic rules in designing qPCR primers
You do not need fancy software to design qPCR primers. I just use [Primer3](http://primer3.ut.ee/) to design qPCR primers like regular primers. Here are my rules.

- Optimum Tm at 60 °C
- Tm difference < 2 degrees between forward and reverse primers
- Amplicon size 60 bp to 200 bp
- If possible, one primer should span an exon-exon junction to avoid amplifying genomic DNA; if the intron is big, two primers on each side of the intron will not amplifying the gDNA either.
- Target specific, at least 2-bp difference in the first 4 bps from the 3' end. If the gene homeologs/paralogs are too conserved, we have to design gene-specific primers on the 5' and 3' UTRs. Just make sure you know the UTR size.

## Test primer efficiency
Before using the primers, make sure to test the primer efficiency. Only primer pairs with efficiency 90% to 110% are good for using.

1. make a series of dilution: depending on your gene's expression level, you can dilute 2x, 4x or 10x on each dilution. When the Ct get to over 34, I found qPCR is not reliable, so sometimes I need to do a regular PCR with my cDNA to increase the concentration of my gene, 10 to 15 cycles of regular PCR is enough. Then I dilute the PCR product 10x as 1x template, then dilute 4x for each dilution to make 5 dilutions.

1. Check the melting point after the qPCR run. If the melting peak is not unique, your primers are not specific.

1. If the melting point is good, you can go ahead to calculate the efficiency. Here is [a good article](https://biosistemika.com/blog/qpcr-efficiency-over-100/) for you to learn how to do that.

## Other tricks
- Always have negative controls (water)
- Always use filter tips to avoid contamination
- Have more than 3 biological replicates
- I usually do not do technical replicates when I have enough biological replicates.
- Choose a good internal control, ACTIN is not always good, depending on your tissue and stage.
- Sometimes you need to set the threshold manually, for example, when to compare two different runs for the same set of samples.

## Update 2019-03-22

I repeated my qPCR using the same set of primers for several times, and later I found I got contaminations, because my negative control (water) also got strong amplification (Ct < 34). So I have to order new primers, replace all reagents, clean my bench and pipettes with 10% bleach, and even changed new clothes and prepared the qPCR in a hood. So here is something more to consider:

- Do not open your PCR tubes in the preparation area to avoid contamination around your bench.
- Always use filter tips and change gloves often too.
- You may also need to consider setting another negative control: reverse transcription control (with or without reverse transcriptase), which will tell you the genomic DNA contamination and whether your reverse transcription is good.

