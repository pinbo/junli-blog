---
title: Analyze multiplexed long-read amplicon sequencing
author: Junli Zhang
date: '2025-04-03'
slug: analyze-multiplexed-long-read-amplicon-sequencing
categories:
  - NGS
  - method
tags:
  - long-read
  - ampliconSeq
  - minimap2
---

[Azanta/Genwiz](https://www.genewiz.com/) has Oxford Nanopore Technology (ONT) long-read sequencing service for plasmid ([plasmid-EZ](https://www.genewiz.com/public/services/next-generation-sequencing/whole-plasmid-sequencing-plasmid-ez)) and long PCR amplicons (>500 bp, [PCR-EZ](https://web.genewiz.com/pcr-ez)). If you want to check the variations for a long gene (>500 bp) for multiple genotypes, you can do two rounds of PCR with the same barcodes and adaptor used for the CRISPR editing check. Then mix them and sequence all the long amplicons with Genwiz [PCR-EZ](https://web.genewiz.com/pcr-ez) (or [plasmid-EZ](https://www.genewiz.com/public/services/next-generation-sequencing/whole-plasmid-sequencing-plasmid-ez)).

The results will be a fastq file with ONT long reads. You can do the analysis using the following steps:

1. Demultiplex the fastq files [here](/apps/demultiplex-reads-with-barcodes-on-both-ends/).
2. Prepare a fasta file with your gene references and map the demultiplexed fastq files [here](/apps/map-long-reads-with-minimap2/) with `minimap2`.
3. Download the bam files and view them in [IGV](https://igv.org/). 

**Hint**: ONT reads are very noisy, so only consistent variations across reads are real.