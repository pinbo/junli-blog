---
title: Check CRISPR editing events with PCR amplicon NGS
author: Junli Zhang
date: '2021-06-21'
slug: check-crispr-editing-with-pcr-amplicon-ngs
categories:
  - Method
tags:
  - CRISPR
  - NGS
  - Amplicon sequencing
---

I wrote a [method](/en/process-ngs-of-crispr-editing) to check CRISPR editing in transgenic plants with NGS, but that method needs a whole sequening lane, which is usually an overkill for CRISPR checking. We found [MGH](https://dnacore.mgh.harvard.edu/new-cgi-bin/site/pages/crispr_sequencing_pages/crispr_sequencing_submission.jsp) and [Genewiz](https://www.genewiz.com/en/Public/Services/Next-Generation-Sequencing/Amplicon-Sequencing-Services) provide PCR amplicon sequencing service. Each sample is about $70, but we can mix our PCR amplicons with barcodes using the same two-round PCR method in my last [post](/en/process-ngs-of-crispr-editing). Because the service provider will add sequencing adapters and their barcodes to the submitted samples,  we need differnt adapters and barcodes this time. I have designed new PCR adapters and barcodes for PCR amplicon sequencing. You can download the protocol [here](/files/Check-CRISPR-editing-events-in-transgenic-wheat-with-NGS.pdf) and the barcodes [here](/files/PCR-amplicon-NGS-barcodes.xlsx).

The results are fastq files and you can analyze by the following steps:

1. Check quality and filtering with fastp [here](/apps/filter-fastq-files-with-fastp/).
2. Demultiplex the filtered fastq files [here](/apps/demultiplex-a-fastq-file/).
3. Map reads to your PCR template and get SNP and small indel calls [here](/apps/make-bam-files-with-bwa-and-samtools/) with [BWA-MEM](http://bio-bwa.sourceforge.net/) and [SAMTOOLS](http://www.htslib.org/).
4. Check large indels (>15bp) and structure variations [here](/apps/make-bams-and-indel-calls-with-subread/) with [Subread](http://subread.sourceforge.net/).
5. View the downloaded bam files with [**IGV**](https://software.broadinstitute.org/software/igv/download). You need to install **IGV** first. The web version does not work well.