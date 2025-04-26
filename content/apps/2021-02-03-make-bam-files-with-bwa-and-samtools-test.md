---
title: CRISPR editing check with BWA and SAMTOOLS v2
author: Junli Zhang
date: '2022-05-30'
slug: make-bam-files-with-bwa-and-samtools-test
categories:
  - Research
tags:
  - CRISPR
  - NGS
  - bwa
  - samtools
---

This tool supports both paired end fastq files (for example, xxx_R1_001.fastq.gz and xxx_R2_001.fastq.gz) and single end fastq files (please make sure they have suffix1 below).  
Please run the 3 steps below to get the indexed bams from a list of fastq files.
<p id=recommend" style="color:darkviolet;">Recommend using private browser windows to avoid troubles caused by cookies and caches (open from the menu at the topright corner)</p>
<p id=recommend2" style="color:red;">No spaces are allowed in input file names!</p>

**Provide the suffix of your fastq files**:  
<label for="suffix1">Read 1:</label>
<input id="suffix1" value="_R1_001.fastq.gz" size="40"><br>
<label for="suffix2">Read 2:</label>
<input id="suffix2" value="_R2_001.fastq.gz" size="40"><br>

<h4>I. Choose reference file (a fasta file)</h4>
<input id="reference" type="file">

<h4>II. Choose demultiplexed fastq(.gz) files</h4>
<input id="fastq" type="file" multiple>

<p id="indexErr" style="color:red;"></p>
<p id="demoRef" style="display:none;"></p>
<p id="demoFq" style="display:none;"></p>

<h4>III. Map reads and create bam files</h4>

After loading the template fasta file and all the fastq files, now we will use `bwa` and `samtools` to create indexed bam files for viewing in the software [IGV](https://software.broadinstitute.org/software/igv/download).

<div id="options" style="font-size:90%;color:blue;">

Options:<br>
<input size="2" id="mq" value="0" type="text"> min mapping quality when calling variants  

</div>

<button onclick="analyzeBam()">Map reads and Make bam files</button>
<p id="bwa"  style="color:tomato;font-style: italic;"></p>
<p id="sort" style="color:tomato;font-style: italic;"></p>
<button id="download-btn" onclick="downloadBam()" style="visibility:hidden">Download indexed bam files</button>
<p id="download" style="color:tomato;font-style: italic;"></p>
<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/bwa-samtools-v5.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>

## Help

This tool is a WebAssembly implementation of [BWA](http://bio-bwa.sourceforge.net/) and [SAMTOOLS](http://www.htslib.org/). It runs commands like this:
```sh
## map reads to the templates with bwa
bwa index <your-references.fa>
bwa mem <your-references.fa> <R1.fastq.gz> <R2.fastq.gz> > out.sam
## make sorted bam files with samtools
samtools sort out.sam > out.bam
samtools index out.bam
## I added editcall.c to bwa tool sets to call variants
bwa editcall -f your-reference.fa -o calledSNPs.txt out.sam
```

Visit the GitHub page for more details: [https://github.com/pinbo/bwa-samtools-web](https://github.com/pinbo/bwa-samtools-web).

### Enable SIMD for your browser

BWA needs [SIMD](https://v8.dev/features/simd) for vector calculation. Please enable it in your web brower first (just need to do once).

- chromium based browsers (Google Chrome and new Microsoft Edge): go to URL [chrome://flags/](chrome://flags/), search `WebAssembly SIMD support`, and select "Enabled"; **seems enabled by default since 2021.**

- Firefox: go to URL [about:config](about:config), search `javascript.options.wasm_simd`, then choose `true`; **seems enabled by default since 2022.**

- Other browers were not checked.

## updates

- 2022-05-08: add `editcall` [nim version](https://github.com/pinbo/editcall) to call indels and inversions from `bwa mem` sam file.
- 2022-05-30: add `editcall` [c version](https://github.com/pinbo/practice_c) to call SNPs, indels and inversions from `bwa mem` sam file.
- 2022-05-30: removed exactSNP results due to its wrong call of SNPs.
- 2022-12-27: support single end reads now.
- 2023-11-02: bwa-mem now put the reads with mapq=0 (multi-map) in the same chromosome as its mate if its mate has mapq > 0.

## Acknowledgement

- Thank the authors of **bwa**, **samtools** and **Subread**.
- Thank [Robert Aboukhalil](https://github.com/robertaboukhalil) for his development of [aioli and biowasm](https://github.com/biowasm). I successfully compiled **bwa**, **samtools** and **Subread** by learning the patch and compiling scripts from [biowasm](https://github.com/biowasm/biowasm).