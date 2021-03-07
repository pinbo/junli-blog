---
title: Make bams and indel calls with Subread
author: Junli Zhang
date: '2021-03-06'
slug: make-bams-and-indel-calls-with-subread
categories:
  - Method
tags:
  - bam
  - indel
  - subread
  - NGS
  - CRISPR
---

This tool is for **paired end** fastq files (for example, xxx_R1_001.fastq.gz and xxx_R2_001.fastq.gz).  
Please run the 3 steps below to get the indexed bams from a list of fastq files.  
**Please make sure you have at least 1.5 Gb of free memory.**
<p id=recommend" style="color:darkviolet;">Recommend using private browser windows to avoid troubles caused by cookies and caches (open from the menu at the topright corner)</p>

**Provide the suffix of your fastq files**:  
<label for="suffix1">Read 1:</label>
<input id="suffix1" value="_R1_001.fastq.gz" size="40"><br>
<label for="suffix2">Read 2:</label>
<input id="suffix2" value="_R2_001.fastq.gz" size="40"><br>

<h4>I. Load reference file (a fasta file)</h4>
<input id="reference" type="file">

<h4>II. Load demultiplexed fastq(.gz) files</h4>
<input id="fastq" type="file" multiple>

<p id="indexErr" style="color:red;"></p>
<p id="demoRef" style="display:none;"></p>
<p id="demoFq" style="display:none;"></p>

<h4>III. Map reads and create bam files</h4>

After loading the template fasta file and all the fastq files, now we will use [subread](http://subread.sourceforge.net/) to create indexed bam files for viewing in the software [IGV](https://software.broadinstitute.org/software/igv/download) and call indels and structure variations.

<button onclick="makeAll()">Map reads and Make bam files</button>
<p id="bam" style="color:pink;font-style: italic;"></p>
<button id="download-btn" onclick="downloadAll()" style="visibility:hidden">Download indexed bam files</button>
<p id="download" style="color:pink;font-style: italic;"></p>

**Running Summary** [you can also check the debug information with the brower developer tool (Ctrl+Shift+I for Chrome and Firefox)]:
<textarea id="stderr" name="stderr" rows="27" cols="100" style="font-family: monospace;font-size: 10px;" placeholder="Software running informaiton will be shown here"></textarea><br>

<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/subread.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>

## Help

This tool is a WebAssembly implementation of [subread](http://subread.sourceforge.net/). It runs commands like this:

```sh
## index the references
subread-buildindex -M 1000 -o my_index <your-references.fa>
## make sorted bam files and vcf files for indels (<=16 bp) 
#  and structure variaitons (indels > 16 bp and inversions etc)
subread-align -i my_index -r xxxx_R1_001.fastq.gz -R xxxx_R2_001.fastq.gz
              -o xxxx.bam -I 16  -sv --sortReadsByCoordinates
```

Visit the GitHub page for more details: [https://github.com/pinbo/bwa-samtools-web](https://github.com/pinbo/bwa-samtools-web).

## Acknowledgement

- [subread](http://subread.sourceforge.net/) is a great software for calling big indels.
- Thank [Robert Aboukhalil](https://github.com/robertaboukhalil) for his development of [aioli and biowasm](https://github.com/biowasm).