---
title: 'Extract cDNA sequences'
author: Junli Zhang
date: '2024-07-21'
slug: extract-sequences
categories:
  - tools
tags:
  - DNA
  - samtools
---

Extract cDNAs from CS and Kronos.

Please paste gene IDs (e.g. TraesCS5A02G391700) or specific cDNA IDs (e.g. TraesCS5A02G391700.2) below.

*Each line is a gene*.

<textarea rows="10" cols="75" id="input"></textarea>
<br />

<button id="run">Submit</button>
<button id="clearseq">Clear</button>
<button id="example">Example</button>

<select id="box1">
  <option value="CS_cDNA_HC_v1.1">CS IWGSC cDNA v1.1 HC</option>
  <option value="Kronos_cDNA_v1.0">Kronos cDNA v1.0</option>
  <option value="CS_cDNA_LC_v1.1">CS IWGSC cDNA v1.1 LC</option>
</select>

<p>Output below</p>
<textarea rows="10" cols="75" id="output" ></textarea>
<br />


<script src="https://biowasm.com/cdn/v3/aioli.js"></script>
<script type="module" src="/libs/extract-sequences.js"></script>