---
title: 'Find wheat potential homeologs'
author: Junli Zhang
date: '2024-07-23'
slug: wheat-homeolog
categories:
  - tools
tags:
  - wheat
  - homeolog
  - sqlite3
---

Find potential wheat homeologs (best hit with >90% identity and alignment >60% of the CDS length) and their functions based on Arabidopsis (At) and rice (Os) blast results (top 1 hit).

Please paste gene IDs (e.g. TraesCS5A02G391700) below. *Each line is a gene*.

To start a new job, click "**Clear**" button below, and resubmit (faster than refresh the page).

**Database** to search:
<select id="box1">
  <option value="Kronos_cDNA_v1.0">Kronos cDNA v1.0</option>
  <option value="CS_cDNA_HC_v1.1">CS IWGSC cDNA v1.1 HC</option>
</select>  
<input type="checkbox" id="check1" name="check1" />
<label for="check1">Output At/Os best hits only</label>

<textarea rows="10" cols="75" id="input" placeholder="paste gene names here: one gene per line"></textarea>
<br />

<button id="run">Submit</button>
<button id="clearseq">Clear</button>
<button id="example">Example</button>

**Output below**

<!-- <textarea rows="10" cols="75" id="output" ></textarea> -->
<!-- <br /> -->
<p id="alert" style="color:blue";></p>
<button id="copytable">Copy to clipboard</button>
<a download="wheat_homeolog_and_function.csv" href="#" onclick="return ExcellentExport.csv(this, 'datatable');" style="color:Tomato;">or Export to CSV</a>
<!-- The button used to copy the text -->
<table id="datatable" style="font-size: 10px;" align="left">
<thead id="thead">
    <tr>
        <th>WheatGeneID</th>
        <th>Best Wheat matches</th>
        <th>Wheat %identity</th>
        <th>Best At matches</th>
        <th>At %identity</th>
        <th>At description</th>
        <th>Best Os matches</th>
        <th>Os %identity</th>
        <th>Os description</th>
    </tr>
</thead>
    <tbody id="tbody"></tbody>
</table>

<script src="/tools/sqljs/v1.10.3/sql-wasm.js"></script>
<!-- <script src="/tools/sqlite3/3.46.0/sqlite3.js"></script> -->
<script type="module" src="/libs/get-wheat-homeologs-v6.js"></script>
<script src="/libs/excellentexport.min.js"></script>
<script src="/libs/pako_inflate.min.js"></script>
