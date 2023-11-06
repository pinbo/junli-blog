---
title: Simple NGS assembly with velvet
author: Junli Zhang
date: '2023-11-05'
slug: simple-NGS-assembly-with-velvet
categories:
  - Method
tags:
  - assembly
  - NGS
  - velvet
---

<!--<script src="https://biowasm.com/cdn/v3/aioli.js"></script>
<script type="module" src="/libs/velvet.js"></script>
 <script type="module">
  import { run, printHelp, clearHelp } from "/libs/velvet.js";
  window.run = run;
  window.printHelp = printHelp;
  window.clearHelp = clearHelp;
</script> -->

[Velvet](https://github.com/dzerbino/velvet) is software to perform dna De novo assembly from short reads by manipulating de Bruijn graphs. Some introduction can be found here:
[Introduction to de novo assembly with Velvet](https://www.melbournebioinformatics.org.au/tutorials/tutorials/assembly/assembly-background/). It is an old software and recently [SPAdes](https://github.com/ablab/spades) is the most popular de novo assembler. I used `Velvet` here simply because it is easy to be exported to WebAssembly and it is good enough for simple projects like gene/plasmid assembly.



<p id=recommend" style="color:darkviolet;">Recommend using private browser windows to avoid troubles caused by cookies and caches (open from the menu at the topright corner)</p>
<p id=recommend2" style="color:red;">No spaces are allowed in input file names!</p>

<h4>Step 1: load fastq files and set up parameters</h4>
<div id="options" style="font-size:90%;color:blue;">
<label for="fastq">Load fastq files (choose 1 file for single end (SE) sequencing or 2 files for paired end (PE) sequencing or 1 file for interleaved PE)</label><br>
<input id="fastq" type="file" multiple><br>
<p id="demoFq" style="display:none;"></p>

<input type="checkbox" id="interleaved" name="interleaved" value="--interleaved_in">
<label for="interleaved">Interleaved PE? (Read1 and read2 are in one fastq file)</label>

<select name="fileFormat" id="fileFormat">
  <option value="-fastq.gz">fastq.gz</option>
  <option value="-fastq">fastq</option>
  <!-- <option value="-fasta.gz">fasta.gz</option>
  <option value="-fasta">fasta</option>
  <option value="-sam">sam</option>
  <option value="-bam">bam</option> -->
</select>
<label for="fileFormat">File format</label>

<input id="minCov" name="minCov" value="3" size="2">
<label for="minCov">Minimum coverage to output</label>

<input id="insertSize" name="insertSize" value="300" size="2">
<label for="insertSize">Average insert size (just a guess based on the gel picture or other quality test)</label>

<input id="expCov" name="expCov" value="30" size="2">
<label for="expCov">Expected region coverage (just a guess)</label>

<input id="hashLength" name="hashLength" value="31" size="2">
<label for="hashLength">Hash length (an odd integer <= 31 OR m,M,s where m and M are odd integers <= 31 and s is a step (even number))</label>

<textarea id="addopt1" name="addopt2" rows="1" cols="60" placeholder="Additional Options for program velveth"></textarea><br>

<textarea id="addopt2" name="addopt2" rows="1" cols="60" placeholder="Additional Options for program velvetg"></textarea><br>

<!-- <button id="printHelp" onclick="printHelp()">Print Help of fastp</button>
<button id="clearHelp" onclick="clearHelp()">Clear Help</button><br> -->
<button id="printHelp">Print Help of velvet</button>
<button id="clearHelp">Clear Help</button><br>

<p id="help"></p>
</div>
<h4>Step 2: Start Assembling</h4>
<!-- <button id="run" onclick="run()">Start</button> -->
<button id="run">Start</button>

<div id="download-btn" style="display:none">
    <h4>Step 3: Download the assembled contig</h4>
    <button id="download" >Download the assembled contigs</button><br><br>
    <!-- <a id="downloadLink" download>Download Assembly</a> -->
</a>
</div>
<p id="error" style="color:red;"></p>
<pre><code id="stdout"></code></pre>

## Acknowledgements
1. **velvet**: https://github.com/dzerbino/velvet
2. [**biowasm project**](https://github.com/biowasm/biowasm) for easy exporting C/C++ tools to WebAssembly.

<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/velvet-v2.js"></script>