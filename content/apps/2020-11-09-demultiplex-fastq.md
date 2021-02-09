---
title: 'Demultiplex a Fastq File'
author: Junli Zhang
date: '2020-11-09'
slug: Demultiplex-a-Fastq-File
categories:
  - Research
tags:
  - javascript
  - CRISPR
  - NGS
---

**Input is an interleaved fastq file**  
You can use tool ["Filter fastq files with FASTP"](/apps/filter-fastq-files-with-fastp) to make an interleaved fastq if you have two fastq files (R1 and R2).


<label for="left">Your left adapter sequence:</label>
<input id="left" name="LeftAdapter" placeholder="adapter added to your forward PCR primer" size="40"><br>
<label for="right">Your right adapter sequence:</label>
<input id="right" name="RightAdapter" placeholder="adapter added to your reverse PCR primer" size="40"><br>
<label for="leftBarcodeLen">Left barcode length:</label>
<input id="leftBarcodeLen" name="leftBarcodeLen" placeholder="left barcode length" size="20" value="8"><br>
<label for="rightBarcodeLen">Right barcode length:</label>
<input id="rightBarcodeLen" name="rightBarcodeLen" placeholder="right barcode length" size="20" value="8"><br>
<button onclick="clearseq()">Clear</button>
<button onclick="putExample()">Example Adapters</button><br>
<label for="barcode">Select the barcode file (a file with 3 columns separated by tab: sample ID, left barcode sequence, right barcode sequence)</label><br>
<input type="file" id="barcode" name="barcode" /><br>
<label for="fastq">Select the fastq or fastq.gz file to demultiplex</label><br>
<input type="file" id="fastq" name="fastq" /><br>

<button style="display:none" onclick="readBarcode()">Read Barcodes</button>
<button onclick="startAnalyze()">Start Demultiplex</button>
<p id="demo1"></p>
<p id="demo2"></p>
<button style="visibility:hidden" id="download-btn" onclick="download()"> Download Demultiplexed Files</button><br>
<p id="demo3"></p>
<output id="output" style="display:none"></output>

<script src="/libs/pako.min.js"></script>
<script src="/libs/FileSaver.min.js"></script>
<script src="/libs/demultiplex.js"></script>
<script src="/libs/jszip.min.js"></script>

## Help

We are using the CRISPR sequencing service at [MGH DNA Core](https://dnacore.mgh.harvard.edu/new-cgi-bin/site/pages/crispr_sequencing_main.jsp). We usually add our own barcodes to pool multiple samples as one sample for submission. Over there, they will add the sequencing adapters and barcodes for Illumina NGS. The returned data is an interleaved fastq file (R1 and R2 are in the same file). Therefore, we need to demultiplex it before checking mutations.

We do two round of PCRs to add barcodes to the PCR amplicons. The final PCR amplicons are like this:

![NGS-PCR](/images/NGS-PCR.png)

The program firstly look for left adapter and right adapter to get the orientation of the read: left adapter = R1 and right adapter = R2. Then check the barcode combination to sort them into different samples.

The final output is a zipped file containing all fastq.gz files. You can check the editing events in the samples with the [**CRISPR Editing Analysis**](/apps/crispr-editing-check) program.