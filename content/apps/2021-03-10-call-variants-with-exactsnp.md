---
title: Call variants with exactSNP
author: Junli Zhang
date: '2021-03-10'
slug: call-variants-with-exactsnp
categories:
  - Method
tags:
  - bam
  - indel
  - subread
  - NGS
  - CRISPR
  - SNP
---

This tool is for calling SNPs and small indels from bam files with subread tool `exactSNP`. All the operations are within your memory.
<p id=recommend" style="color:darkviolet;">Recommend using private browser windows to avoid troubles caused by cookies and caches (open from the menu at the topright corner)</p>
<p id=recommend2" style="color:red;">No spaces are allowed in input file names!</p>

<h4>I. Load the reference (a fasta file)</h4>
<input id="reference" type="file"><br>
<p id="demoRef" style="display:none;"></p>

<h4>II. Load bam files</h4>
<input id="fastq" type="file" multiple>
<p id="demoFq" style="display:none;"></p>

<h4>III. Call variants</h4>

One vcf file will be created for each bam. A summary tab-delimited text file merging all vcf files will be created too.

<button onclick="makeAll()">Call Variants</button>
<p id="bam" style="color:tomato;font-style: italic;"></p>
<p id="indexErr" style="color:red;"></p>
<button id="download-btn" onclick="downloadVar()" style="visibility:hidden">Download variant summary</button>
<p id="download" style="color:tomato;font-style: italic;"></p>

<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/exactSNP.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>

## Help

This tool is a WebAssembly implementation of [Subread exactSNP](http://subread.sourceforge.net/). It runs commands like this:

```sh
exactSNP -b -i mapping_results.bam -g mm10.fa -o calledSNPs.vcf
```

Visit the GitHub page for more details: [https://github.com/pinbo/bwa-samtools-web](https://github.com/pinbo/bwa-samtools-web).

## Acknowledgement

- [Subread](http://subread.sourceforge.net/).
- Thank [Robert Aboukhalil](https://github.com/robertaboukhalil) for his development of [aioli and biowasm](https://github.com/biowasm).