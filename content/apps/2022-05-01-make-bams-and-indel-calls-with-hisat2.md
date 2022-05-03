---
title: CRISPR NGS editing check with HISAT2
author: Junli Zhang
date: '2022-05-01'
slug: make-bams-and-indel-calls-with-hisat2
categories:
  - Method
tags:
  - bam
  - indel
  - hisat2
  - NGS
  - CRISPR
---

**Before start**: Please enable [SIMD](https://v8.dev/features/simd) in your web brower [see [Help](#enable-simd-for-your-browser) below].

This tool is for **paired end** fastq files (for example, xxx_R1_001.fastq.gz and xxx_R2_001.fastq.gz). It uses [HISAT2](http://daehwankimlab.github.io/hisat2/) to find big deletions (treated as introns).
Please run the 3 steps below to get the indexed bams from a list of fastq files.
In the end, you will get a summary file of called indels and all the bam files to be viewed in [IGV](https://software.broadinstitute.org/software/igv/download).  

<p id=recommend2" style="color:red;">No spaces are allowed in input file names!</p>

**Provide the suffix of your fastq files**:  
<label for="suffix1">Read 1:</label>
<input id="suffix1" value="_R1_001.fastq.gz" size="40"><br>
<label for="suffix2">Read 2:</label>
<input id="suffix2" value="_R2_001.fastq.gz" size="40"><br>

<h4>I. Load reference file (a fasta file)</h4>

<input id="reference" type="file">

<h4>II. Load demultiplexed fastq(.gz) files (adapters need to be trimmed)</h4>
<input id="fastq" type="file" multiple>

<p id="indexErr" style="color:red;"></p>
<p id="demoRef" style="display:none;"></p>
<p id="demoFq" style="display:none;"></p>

<h4>III. Map reads and create bam files</h4>

After loading the template fasta file and all the fastq files, now we will use [hisat2](http://daehwankimlab.github.io/hisat2/) to map reads to your templates and use [exactSNP](http://subread.sourceforge.net/) to call variations.
<div id="options" style="font-size:90%;color:blue;">
Addition parameters for HISAT2 &nbsp;&nbsp;&nbsp;&nbsp;<input size="20" id="hisat2" value="" type="text">  
Addition parameters for exactSNP <input size="20" id="exactSNP" value="" type="text">
</div>
<button onclick="makeAll()">Map reads and Make bam files</button>
<p id="bam" style="color:Tomato;font-style: italic;"></p>
<p id="bamErr" style="color:red;font-style: italic;"></p>
<button id="download-btn" onclick="downloadBam()" style="visibility:hidden">Download indexed bam files</button>
<p id="download" style="color:Tomato;font-style: italic;"></p>

**Running Summary** [you can also check the debug information with the brower developer tool (Ctrl+Shift+I for Chrome and Firefox)]:
<textarea id="stderr" name="stderr" rows="14" cols="85" style="font-family: monospace;font-size: 12px;" placeholder="Software running informaiton will be shown here"></textarea><br>

<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/hisat2web.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>

## Help

This tool is a WebAssembly implementation of [hisat2](http://daehwankimlab.github.io/hisat2/). It runs commands like this:

```sh
## index the references
hisat2-build-s <your-references.fa> my_index
## make sam files and big indels will be treated as splice
hisat2-align-s -x my_index -1 xxx_R1_001.fastq -2 xxx_R2_001.fastq -S out.sam --pen-noncansplice 0
```

Although hisat2 can soft clip unmapped fragments, it seems it can map more reads when using adapter-trimmed reads. So please trim adapters when demulitiplexing your fastq files and when filtering the fastq files with fastp.

Visit the GitHub page for more details: [https://github.com/pinbo/bwa-samtools-web](https://github.com/pinbo/bwa-samtools-web).

### Enable SIMD for your browser

hisat2 needs [SIMD](https://v8.dev/features/simd) for vector calculation. Please enable it in your web brower first (just need to do once).

- chromium based browsers (Google Chrome and Microsoft Edge): **seems enabled by default since 2021.**; 

- Firefox: go to URL [about:config](about:config), search `javascript.options.wasm_simd`, then choose `true`;

- Other browers do not seem to support this yet.

## Acknowledgement

- Thank authors of [hisat2](http://daehwankimlab.github.io/hisat2/).
- Thank [Robert Aboukhalil](https://github.com/robertaboukhalil) for his development of [aioli and biowasm](https://github.com/biowasm).