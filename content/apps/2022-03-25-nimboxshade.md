---
title: Shade sequence alignments with nimBoxshade
author: Junli Zhang
date: '2022-03-25'
slug: nimboxshade
categories:
  - Research
tags:
  - Boxshade
  - multiple sequence alignment
  - nim
draft: false
---

This is a duplicate to the app [boxshade page](/apps/boxshade/). It is just my testing of [Nim language](https://nim-lang.org/). It uses [nimBoxshade](https://github.com/pinbo/nimBoxshade) to create good looking printouts from multiple-aligned protein and DNA sequences.

<h4>Step 1: load alignment file (fasta format)</h4>
<div id="options" style="font-size:90%;color:blue;">
<input id="snpfile" type="file"><br>
<!-- 
<textarea id="paste" name="paste" rows="6" cols="85" placeholder="OR paste your sequences here in fasta format"></textarea><br>
<p id="demoFq" style="display:none;"></p> -->



<!-- Output file name (without extension) <input size="20" id="output" value="" type="text">   -->

<input type="checkbox" id="ruler" value="1" checked> add ruler to the alignment?  
<input type="checkbox" id="seqnum" value="1" checked> add sequence position number?  
<input type="checkbox" id="consensus" value="1"> create consensus line?  
<input type="checkbox" id="dna" value="1"> input is aligned DNA sequences?  
<input size="2" id="fraction" value="0.5" type="text"> the fraction of sequences that must agree for a consensus (0-1)  
<input size="2" id="outlen" value="60" type="text"> output width

<p id="help"></p>
</div>
<h4>Step 2: Start shading</h4>
<button onclick="process()">Start shading</button>


<div id="download-btn" style="display:none">
    <h4>Step 3: Download output</h4>
    <button id="download" onclick="download()">Download the formatted alignment</button><br><br>
</div>
<p id="error" style="color:red;"></p>
<pre><code id="stdout"></code></pre>


## Help

nimBoxshade now only support aligned fasta file as the input file and output a RTF file that can be opened with Word.
Please see the github page for more information: https://github.com/pinbo/nimboxshade

The amino acid similarity matrix is from [this paper](https://doi.org/10.1186/1471-2105-10-394).

<!-- <script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script> -->
<script src="/tools/aioli/latest/aioli.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/nimboxshadeweb.js"></script>