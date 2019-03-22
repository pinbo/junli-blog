---
title: Call SNPs from multiple sequence alignment
author: Junli Zhang
date: '2019-03-21'
slug: call-snps-from-multiple-sequence-alignment
tags:
  - MSA
  - javascript
---

Please paste your multiple sequence alignment (MSA) in a fasta format in the textbox below. Make sure all sequences have the same length. The program will call SNPs and indels from your input. The script is in my github repository [msa2snp](https://github.com/pinbo/msa2snp).


**Please paste your sequences below** (or click "Example input" button to paste an example)

<script src='/libs/msa2snp.js'></script>

<textarea rows="10" cols="75" id="input"></textarea>
<br />

<button onclick="callsnps()">Submit</button>
<button onclick="paste_example()">Example input</button>
<button onclick="clearseq()">Clear</button>

<p>Output below</p>
<textarea rows="10" cols="75" id="output" ></textarea>
<br />
