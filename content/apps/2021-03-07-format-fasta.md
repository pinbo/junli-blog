---
title: Format FASTA to fixed length
author: Junli Zhang
date: '2021-03-07'
slug: format-fasta
categories:
  - Method
tags:
  - fasta
---

This tool format a fasta file to fixed line length, which is required by some software, such as `subread`.

<label for="fasta">Select a FASTA file</label>:
<input type="file" id="fasta"/><br>
<label for="line_width">Line length to format (bp)</label>
<input id="line_width" value="60" size="4"><br>
</form>
<button id="download-btn" onclick="download()">Download formatted fasta file</button>
<pre id="formatted" style="font-size:11px;"></pre>

<script src="/libs/format_fasta.js"></script>
<script src="/libs/FileSaver.min.js"></script>


## Acknowledgement

- The originial script is from: https://www.biostars.org/p/185103/