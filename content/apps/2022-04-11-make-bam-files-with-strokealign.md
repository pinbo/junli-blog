---
title: Make bam files with Strobealign
author: Junli Zhang
date: '2022-04-11'
slug: make-bam-files-with-strobealign
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

After loading the template fasta file and all the fastq files, now we will use `strobeAlign` and `samtools` to create indexed bam files for viewing in the software [IGV](https://software.broadinstitute.org/software/igv/download).
<div id="options" style="font-size:90%;color:blue;">
Alignment parameters (change unless you know what you are doing):  
<input size="2" id="match" value="4" type="text"> matching score  
<input size="2" id="mismatch" value="8" type="text"> mismatch penalty  
<input size="2" id="gapopen" value="12" type="text"> gap open penalty  
<input size="2" id="gapext" value="1" type="text"> gap extension penalty

<button onclick="analyzeBam()">Map reads and Make bam files</button>
<p id="bwa"  style="color:tomato;font-style: italic;"></p>
<p id="sort" style="color:tomato;font-style: italic;"></p>
<button id="download-btn" onclick="downloadBam()" style="visibility:hidden">Download indexed bam files</button>
<p id="download" style="color:tomato;font-style: italic;"></p>
<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/strobealignweb.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>
</div>
## Help

This tool is a WebAssembly implementation of [strobeAlign](https://github.com/ksahlin/StrobeAlign/), [SAMTOOLS](http://www.htslib.org/) and [subread exactSNP](http://subread.sourceforge.net/). It runs commands like this:
```sh
## map reads to the templates with strobealign
strobealign <your-references.fa> <R1.fastq.gz> <R2.fastq.gz> > out.sam
## make sorted bam files with samtools
samtools sort out.sam > out.bam
samtools index out.bam
## call variants with subread exactSNP (exactSNP has a bug calling sorted bams, so use sams here)
exactSNP -i out.sam -g your-reference.fa -o calledSNPs.vcf
```

**strobeAlign** can detect SNPs and relative large indels from CRIPSR editing. For better SNP and small indels (<10bp), 
[bwa mem](/apps/make-bam-files-with-bwa-and-samtools/) seems more accurate. For structure variations like inversions, please try [subread](http://subread.sourceforge.net/) or my web app [Make bams and indel calls with Subread](/apps/make-bams-and-indel-calls-with-subread).

Visit the GitHub page for more details: [https://github.com/pinbo/bwa-samtools-web](https://github.com/pinbo/bwa-samtools-web).

### Enable SIMD for your browser

strobeAlign needs [SIMD](https://v8.dev/features/simd) for vector calculation. Please enable it in your web brower first (just need to do once).

- chromium based browsers (Google Chrome and new Microsoft Edge): go to URL [chrome://flags/](chrome://flags/), search `WebAssembly SIMD support`, and select "Enabled"; **seems enabled by default since 2021.**

- Firefox: go to URL [about:config](about:config), search `javascript.options.wasm_simd`, then choose `true`;

- Other browers do not seem to support this yet.

## Acknowledgement

Thank [Robert Aboukhalil](https://github.com/robertaboukhalil) for his development of [aioli and biowasm](https://github.com/biowasm). I successfully compiled **strobeAlign**, **samtools** and **Subread** by learning the patch and compiling scripts from [biowasm](https://github.com/biowasm/biowasm).