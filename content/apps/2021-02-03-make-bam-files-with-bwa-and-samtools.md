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

## Before start
BWA needs [SIMD](https://v8.dev/features/simd) for vector calculation. Please enable it in your web brower first (just need to do once).

- chromium based browsers (Google Chrome and new Microsoft Edge): go to URL [chrome://flags/](chrome://flags/), search `WebAssembly SIMD support`, and select "Enabled";

- Firefox: go to URL [about:config](about:config), search `javascript.options.wasm_simd`, then choose `true`;

## Get started

This tool is a WebAssembly implementation of [BWA](http://bio-bwa.sourceforge.net/) and [SAMTOOLS](http://www.htslib.org/). It runs commands like this:
```sh
## map reads to the templates with bwa
bwa index <your-references.fa>
bwa mem <your-references.fa> <R1.fastq.gz> <R2.fastq.gz> > out.sam
## make sorted bam files with samtools
samtools sort out.sam > out.bam
samtools index out.bam
```
This tool is for **paired end** fastq files. Please run the 3 steps below to get the indexed bams from fastq files. 

<label for="suffix">The suffix of your read 1 (R1) fastq files, default is "_R1_001.fastq.gz" for files names like xxx_R1_001.fastq.gz and xxx_R2_001.fastq.gz.</label><br>
<input id="suffix" name="LeftAdapter" value="_R1_001.fastq.gz" size="40"><br>

<h4>I. Choose reference file (a fasta file)</h4>
<input id="reference" type="file">

<h4>II. Choose demultiplexed fastq(.gz) files</h4>
<input id="fastq" type="file" multiple>

<p id="demo1"></p>
<p id="demoRef" style="display:none;"></p>
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