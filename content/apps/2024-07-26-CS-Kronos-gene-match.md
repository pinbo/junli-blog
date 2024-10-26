---
title: 'Gene ID conversion between Kronos and CS'
author: Junli Zhang
date: '2024-07-26'
slug: CS-Kronos-gene-match
categories:
  - tools
tags:
  - wheat
  - Kronos
  - CS
  - sqlite3
---

Match the genes in CS RefSeq v1.1 annotation HC with Kronos v1.0 annotation (best hit with >97% identity and alignment length > 100 bp).

I did not do more filtering than that, so some genes may have multiple hits or not on the same chromosome. 

<p style="color:red";>Select with your own eyes based on chromosome matches and %identity.</p>

Please paste gene IDs (e.g. TraesCS5A02G391700) below.
*Each line is a gene*.

**Gene ID conversion**:
<select id="box1">
  <option value="K2C">Kronos to CS</option>
  <option value="C2K">CS to Kronos</option>
</select>

<textarea rows="10" cols="75" id="input"></textarea>
<br />

<button id="run">Submit</button>
<button id="clearseq">Clear</button>
<button id="example">Example</button>

**Output below**

<!-- <textarea rows="10" cols="75" id="output" ></textarea> -->
<!-- <br /> -->
<p id="alert" style="color:blue";></p>
<button id="copytable">Copy to clipboard</button>
<a download="geneID-conversion-Kronos-CS.csv" href="#" onclick="return ExcellentExport.csv(this, 'datatable');" style="color:Tomato;">or Export to CSV</a>
<!-- The button used to copy the text -->

<table id="datatable" style="font-size: 11px;" align="left">
<thead>
    <tr>
        <th>Kronos transcript ID</th>
        <th>CS transcript ID</th>
        <th>%identity</th>
        <th>alignment length</th>
        <th>Kronos CDS length</th>
        <th>CS CDS length</th>
    </tr>
</thead>
    <tbody id="tbody"></tbody>
</table>

<script src="/tools/sqljs/v1.10.3/sql-wasm.js"></script>
<script type="module" src="/libs/geneID-conversion-Kronos-CS.js"></script>
<script src="/libs/excellentexport.min.js"></script>
<script src="/libs/pako_inflate.min.js"></script>
