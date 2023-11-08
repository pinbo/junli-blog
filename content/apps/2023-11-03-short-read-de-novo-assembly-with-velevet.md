---
title: Short read de novo assembly with Velevet
author: Junli Zhang
date: '2023-11-03'
slug: short-read-de-novo-assembly-with-velevet
categories:
  - Method
tags:
  - biowasm
  - velvet
  - assembly
  - Aioli_v3
---

[Velvet](https://github.com/dzerbino/velvet) is software to perform dna De novo assembly from short reads by manipulating de Bruijn graphs. Some introduction can be found here:
[Introduction to de novo assembly with Velvet](https://www.melbournebioinformatics.org.au/tutorials/tutorials/assembly/assembly-background/). It is an old software and recently [SPAdes](https://github.com/ablab/spades) is the most popular de novo assembler. I used `Velvet` here simply because it is easy to be exported to WebAssembly and it is good enough for simple projects like gene/plasmid assembly.

You can read the [Velvet manual](https://github.com/dzerbino/velvet/blob/master/Manual.pdf) or the [SPAdes instructions](https://github.com/ablab/spades#assembling-long-illumina-paired-reads-2x150-and-2x250) for some parameter setting suggestion.

<p id=recommend" style="color:darkviolet;">Recommend using private browser windows to avoid troubles caused by cookies and caches (open from the menu at the topright corner)</p>
<p id=recommend2" style="color:red;">No spaces are allowed in input file names!</p>

<h4>Step 1: load fastq files and set assembling parameters</h4>
<div id="options" style="font-size:90%;color:blue;">
<label for="fastq">Fastq files to load (1 file for single end (SE) sequencing OR 2 files for paired end (PE) sequencing OR 1 file for interleaved PE)</label><br>
<input id="fastq" type="file" multiple><br>
<p id="demoFq" style="display:none;"></p>

<input type="checkbox" id="interleaved" name="interleaved" value="--interleaved_in">
<label for="interleaved">Interleaved PE? (Read1 and read2 are in one fastq file)</label>

<input id="minCov" name="minCov" value="3" size="3">
<label for="minCov">Minimum coverage to output ("auto" or an integer)</label>

<input id="insertSize" name="insertSize" value="300" size="3">
<label for="insertSize">Average insert size (a guess based on the gel picture or other quality test)</label>

<input id="expCov" name="expCov" value="30" size="3">
<label for="expCov">Expected region coverage ("auto" or a guess based on data size)</label>

<input id="hashLength" name="hashLength" value="31" size="3">
<label for="hashLength">Hash/k-mer length (an odd integer from 21 to about half of the read length)</label>

<label for="addopt1">Additional parameters. Click the help button below to see a list of options</label><br>
<textarea id="addopt1" name="addopt2" rows="1" cols="60" placeholder="Additional Options for program velveth"></textarea><br>

<label for="addopt2">Additional parameters. Click the help button below to see a list of options</label><br>
<textarea id="addopt2" name="addopt2" rows="1" cols="60" placeholder="Additional Options for program velvetg"></textarea><br>

<!-- <button id="printHelp" onclick="printHelp()">Print Help of fastp</button>
<button id="clearHelp" onclick="clearHelp()">Clear Help</button><br> -->
<button id="printHelp1">Print Help for velveth</button>
<button id="printHelp2">Print Help for velvetg</button>
<button id="clearHelp">Clear Help</button><br>

<p id="help"></p>
</div>
<h4>Step 2: Start Assembling</h4>
<!-- <button id="run" onclick="run()">Start</button> -->
<button id="run">Start</button>

<div id="download-btn" style="display:none">
    <h4>Step 3: Download assembly</h4>
    <!-- <button id="download" >Download the assembled contigs</button><br><br> -->
    <a id="downloadLink" download="contigs.fa">**Download Assembly**</a><br>
    <a id="downloadlog" download="log.txt">Download Log File</a><br>
    <a id="downloadstats" download="stats.txt">Download Stats File</a>
</a>
</div>
<p id="error" style="color:red;"></p>
<pre><code id="stdout"></code></pre>



## Acknowledgements
1. **velvet**: https://github.com/dzerbino/velvet
2. [**biowasm project**](https://github.com/biowasm/biowasm) for easy exporting C/C++ tools to WebAssembly.

<script src="https://biowasm.com/cdn/v3/aioli.js"></script>
<script type="module" src="/libs/velvet.js"></script>
<!-- <script type="module">
  import { run, printHelp, clearHelp } from "/libs/velvet.js";
  window.run = run;
  window.printHelp = printHelp;
  window.clearHelp = clearHelp;
</script> -->
