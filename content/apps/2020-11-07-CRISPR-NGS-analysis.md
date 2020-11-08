---
title: 'CRISPR Editing Analysis'
author: Junli Zhang
date: '2020-11-07'
slug: CRISPR-editing-check
categories:
  - Research
tags:
  - javascript
  - CRISPR
  - NGS
---

**See Help below on your first use**

<script src="/libs/pako_inflate.min.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/crisjs-functions.js"></script>
<dev class="crisjs">
    <label for="left">Your left flanking sequence:</label>
    <input id="left" name="LeftSequence" placeholder="a short sequence left of the gRNA" size="40"><br>
    <label for="right">Your right flanking sequence:</label>
    <input id="right" name="RightSequence" placeholder="a short sequence right of the gRNA" size="40"><br>
    <label for="grna">Your guide RNA sequence:</label>
    <input id="grna" name="guideRNA" placeholder="the gRNA sequence" size="40"><br>
    <label for="wt">Wild type sequence</label><br>
    <textarea id="wt" name="Wildtype" rows="5" cols="80" placeholder="The target gene sequences (long enough to include both the left and right flanking sequences)"></textarea>
    <button onclick="clearseq()">Clear</button>
    <button onclick="putExample()">Example Input</button><br>
    <label for="files">Choose fastq or fastq.gz files for analysis</label><br>
    <input type="file" id="files" name="files[]" multiple /><br>
    <button onclick="analyze()">Start Analyze</button>
    <select id="box1">
      <option value="F">Forward Strand (R1)</option>
      <option value="R">Reverse Strand (R2)</option>
    </select>
    <label for="box1">(Read1 or Read2? Forward strand is to check read1; reverse strand is for read2)</label><br><br>
    <button style="visibility:hidden" id="download-btn" onclick="download()"> Download output csv file</button><br>
    <p id="demo"></p>
    <output id="output" style="display:none;"></output>
</dev>

## Help
This program summarize CRISPR editing results based on the idea of [CRIS.py](https://github.com/patrickc01/CRIS.py). You can read the paper [here](https://www.nature.com/articles/s41598-019-40896-w). The basic idea is summarized in its Fig. 1. 
<img src="/images/CRISjs.png" alt="CRISjs dialog" width=100%>

All three sequences (the left and right flanking sequences and the gRNA sequences) should be on the **same strand as the template**. If you are working on polyploid species, your left or right flanking sequences should be unique to your template (subgenome) if the fastq files include homolog/homeologs.

### Step

1. Input all the required sequences. Make sure both the left and the right flanking sequences are in the same read or merged read in the fastq file.
2. Select fastq or fastq.gz files by clicking "**Choose files**".
3. Start analyzing by clicking the "**Start Analyze**" button. 
4. Download the summary file by clicking the "**Download output csv file**" button.

You can click the "**Example input**" button to get the example inputs (then select the example fastq.gz files https://github.com/pinbo/CRISjs/tree/main/example-input).

If you have a lot of fastq files, you may use the standalone command line program [CRISgo](https://github.com/pinbo/CRISgo).

