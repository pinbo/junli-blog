---
title: Filter multiple fastq files with FASTP
author: Junli Zhang
date: '2021-12-22'
slug: filter-multiple-fastq-files-with-fastp
categories:
  - Research
tags:
  - fastq
  - fastp
  - NGS
---

FASTP is a tool designed to provide fast all-in-one preprocessing for FastQ files.  
Please check out the fastp github page for some detailed information:  
https://github.com/OpenGene/fastp
<p id=recommend" style="color:darkviolet;">Recommend using private browser windows to avoid troubles caused by cookies and caches (open from the menu at the topright corner)</p>

**Provide the suffix of your paired fastq files**:  
<label for="suffix1">Read 1:</label>
<input id="suffix1" value="_R1_001.fastq.gz" size="40"><br>
<label for="suffix2">Read 2:</label>
<input id="suffix2" value="_R2_001.fastq.gz" size="40"><br>

<h4>Step 1: load fastq files and set filtering options</h4>
<div id="options" style="font-size:90%;color:blue;">
<label for="fastq">Fastq files to filter (choose 1 file for single end (SE) sequencing or 2 files for paired end (PE) sequencing or 1 file for interleaved PE)</label><br>
<input id="fastq" type="file" multiple><br>
<p id="demoFq" style="display:none;"></p>

<input id="basequality" name="basequality" value="15" size="2">
<label for="basequality">The minimum base quality</label>

<input type="checkbox" id="interleaved" name="interleaved" value="--interleaved_in">
<label for="interleaved">Interleaved PE? (Read1 and read2 are in one fastq file)</label>

<input type="checkbox" id="merge" name="merge" value="-m">
<label for="merge">Merge read1 and read2 if there are overlaps? (3 output: merged, non-merged-R1, non-merged-R2)</label>

<label for="addopt">Additional filtering options. Click the help button below to see a list of options</label><br>
<textarea id="addopt" name="addopt" rows="3" cols="60" placeholder="Additional Options for the tool: -g -x etc"></textarea><br>

<button onclick="printHelp()">Print Help of fastp</button>
<button onclick="clearHelp()">Clear Help</button><br>
<p id="help"></p>
</div>
<h4>Step 2: Start filtering</h4>
<button onclick="makeAll()">Start Filtering</button>
<input type="checkbox" id="interleaved_out" name="interleaved_out" value="--stdout">
<label for="interleaved_out" style="font-size:90%;color:blue;">Interleaved output?</label><br>

<div id="download-btn" style="display:none">
    <h4>Step 3: Download filtered files and summary html</h4>
    <button onclick="download()">Download the filtered fastq file(s)</button><br><br>
</div>
<p id="error" style="color:red;"></p>
<pre><code id="stdout"></code></pre>

**Running Log**

<textarea id="stderr" name="stderr" rows="28" cols="85" style="font-family: monospace;font-size: 12px;" placeholder="Software running informaiton will be shown here"></textarea><br>


## Acknowledgement
1. **fastp**: https://github.com/OpenGene/fastp
2. **WebAssembly of fastp** was compiled from source based on the patch file of https://github.com/biowasm/biowasm/tree/main/tools/fastp with additional option `--interleaved_out` added.

<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/jszip.min.js"></script>
<script src="/libs/pako_deflate.min.js"></script>
<script src="/libs/fastp-multiplex.js"></script>