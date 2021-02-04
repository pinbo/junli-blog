---
title: Make bam files with BWA and SAMTOOLS
author: Junli Zhang
date: '2021-02-03'
slug: make-bam-files-with-bwa-and-samtools
categories:
  - Research
tags:
  - CRISPR
  - NGS
  - bwa
  - samtools
---

<h4>I. Choose reference file (a fasta file)</h4>
<input id="reference" type="file">

<h4>II. Choose demultiplexed fastq(.gz) files</h4>
<input id="fastq" type="file" multiple>

<p id="demo1"></p>
<p id="demoRef"></p>
<p id="demoFq" style="display:none;"></p>

<h4>III. Map reads and create bam files</h4>

After loading the template fasta file and all the fastq files, now we will use tool `bwa` and `samtools` to create indexed bam files for viewing in the software [IGV](https://software.broadinstitute.org/software/igv/download).

<button onclick="makeSam()">Step 1: Map reads with BWA</button>
<p id="bwa"></p>
<button onclick="makeBam()">Step 2: Make sorted bam files</button>
<p id="sort"></p>
<button onclick="downloadBam()">Step 3: Download indexed bam files</button>
<p id="download"></p>
<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/bwa-samtools.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>