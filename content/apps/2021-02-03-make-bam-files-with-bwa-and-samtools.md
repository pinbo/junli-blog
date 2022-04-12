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

**Before start**: Please enable [SIMD](https://v8.dev/features/simd) in your web brower [see [Help](#enable-simd-for-your-browser) below].

This tool is for **paired end** fastq files (for example, xxx_R1_001.fastq.gz and xxx_R2_001.fastq.gz).  
Please run the 3 steps below to get the indexed bams from a list of fastq files.
<p id=recommend" style="color:darkviolet;">Recommend using private browser windows to avoid troubles caused by cookies and caches (open from the menu at the topright corner)</p>

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

<button onclick="analyzeBam()">Map reads and Make bam files</button>
<p id="bwa"  style="color:tomato;font-style: italic;"></p>
<p id="sort" style="color:tomato;font-style: italic;"></p>
<button id="download-btn" onclick="downloadBam()" style="visibility:hidden">Download indexed bam files</button>
<p id="download" style="color:tomato;font-style: italic;"></p>
<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/bwa-samtools.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>

## Help

This tool is a WebAssembly implementation of [BWA](http://bio-bwa.sourceforge.net/), [SAMTOOLS](http://www.htslib.org/) and [subread exactSNP](http://subread.sourceforge.net/). It runs commands like this:
```sh
## map reads to the templates with bwa
bwa index <your-references.fa>
bwa mem <your-references.fa> <R1.fastq.gz> <R2.fastq.gz> > out.sam
## make sorted bam files with samtools
samtools sort out.sam > out.bam
samtools index out.bam
## call variants with subread exactSNP (it has bug calling from bam, so we call SNPs from sam)
exactSNP -i out.sam -g your-reference.fa -o calledSNPs.vcf
```

`bwa mem` can only detect small indels. For indels > 15 bp, please try [subread](http://subread.sourceforge.net/) or my web app [Make bams and indel calls with Subread](/apps/make-bams-and-indel-calls-with-subread).

Visit the GitHub page for more details: [https://github.com/pinbo/bwa-samtools-web](https://github.com/pinbo/bwa-samtools-web).

### Enable SIMD for your browser

BWA needs [SIMD](https://v8.dev/features/simd) for vector calculation. Please enable it in your web brower first (just need to do once).

- chromium based browsers (Google Chrome and new Microsoft Edge): go to URL [chrome://flags/](chrome://flags/), search `WebAssembly SIMD support`, and select "Enabled"; **seems enabled by default since 2021.**

- Firefox: go to URL [about:config](about:config), search `javascript.options.wasm_simd`, then choose `true`;

- Other browers do not seem to support this yet.

## Acknowledgement

Thank [Robert Aboukhalil](https://github.com/robertaboukhalil) for his development of [aioli and biowasm](https://github.com/biowasm). I successfully compiled **bwa**, **samtools** and **Subread** by learning the patch and compiling scripts from [biowasm](https://github.com/biowasm/biowasm).