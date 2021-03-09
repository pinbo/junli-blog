---
title: Index bams with SAMTOOLS
author: Junli Zhang
date: '2021-03-09'
slug: index-bams-with-samtools
categories:
  - Method
tags:
  - bam
  - indel
  - subread
  - NGS
  - CRISPR
---

This tool is for indexing small bam files with `samtools index`. All the operations are within your memory.
<p id=recommend" style="color:darkviolet;">Recommend using private browser windows to avoid troubles caused by cookies and caches (open from the menu at the topright corner)</p>

<h4>I. Load bam files</h4>
<input id="fastq" type="file" multiple>
<p id="demoFq" style="display:none;"></p>

<h4>II. Make index</h4>

`.bai` indexes will be created with `samtools index`.

<button onclick="makeAll()">Index bams</button>
<p id="bam" style="color:pink;font-style: italic;"></p>
<p id="indexErr" style="color:red;"></p>
<button id="download-btn" onclick="downloadBam()" style="visibility:hidden">Download indexed bam files</button>
<p id="download" style="color:pink;font-style: italic;"></p>

<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/samtools-index.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>

## Help

This tool is a WebAssembly implementation of [samtools](http://www.htslib.org/). It runs commands like this:

```sh
## index the bam
samtools index <bam file>
```

Visit the GitHub page for more details: [https://github.com/pinbo/bwa-samtools-web](https://github.com/pinbo/bwa-samtools-web).

## Acknowledgement

- [samtools](http://www.htslib.org/).
- Thank [Robert Aboukhalil](https://github.com/robertaboukhalil) for his development of [aioli and biowasm](https://github.com/biowasm).