---
title: 'CRISPR Editing Analysis v2'
author: Junli Zhang
date: '2020-11-15'
slug: CRISPR-editing-check-v2
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
<script src="/libs/crisjs-functions-v2.js"></script>
<dev class="crisjs">
<label for="reference">Select the reference (a fasta file with all the gene sequences)</label><br>
	<input type="file" id="reference" name="reference" /><br>
	<label for="sequences">Please input the gene name, Left flanking, Right flanking, and gRNA sequences, separated by SPACE:</label><br>
	<textarea id="sequences" name="sequences" rows="5" style="width: 100%; max-width: 100%;" placeholder="template ID, Left flanking, Right flanking, and gRNA sequences"></textarea>
	<br>
	<button onclick="clearseq()">Clear</button>
	<button onclick="putExample()">Example Input</button><br>
	<label for="reference">Select the fastq or fastq.gz files for analyzing</label><br>
	<input type="file" id="files" name="files[]" multiple /><br>
	<button onclick="analyze()">Start Analyze</button>
	<select id="box1">
		<option value="F">Forward Strand</option>
		<option value="R">Reverse Strand</option>
	</select>
	<label for="box1">(Read1 or Read2? Forward strand is to check read1; reverse strand is for read2)</label><br>
	<progress style="visibility:hidden" id="progress" value="0" max="100"></progress><br>
	<dev id="download" style="visibility:hidden">
		<button id="download-btn" onclick="download()"> Download output csv file</button><br>
		<label for="download-name">Download file name (optional)</label>
		<input id="download-name" name="download-name" placeholder="CRISPR-eidting-summary.csv" size="60"/><br>
	</dev>
	<p id="demo"></p><br>
	<output id="output" style="display:none"></output>
	<output id="template" style="display:none"></output>
</dev>

## Help
This program summarize CRISPR editing results based on the idea of [CRIS.py](https://github.com/patrickc01/CRIS.py). You can read the paper [here](https://www.nature.com/articles/s41598-019-40896-w). The basic idea is summarized in its Fig. 1. 
<img src="/images/CRISjs.png" alt="CRISjs dialog" width=100%>

All three sequences (the left and right flanking sequences and the gRNA sequences) should be on the **same strand as the template**. If you are working on polyploid species, your left or right flanking sequences should be unique to your template (subgenome) if the fastq files include homolog/homeologs.

**GitHub page**: https://github.com/pinbo/CRISjs

### Differences with version 1
Version 1 can be found [here](/apps/crispr-editing-check).
1. Version 2 now get the template from a fasta file containing all sequences.
2. All 3 sequences needed now are in one input form, so it is easy to copy and paste.

### Steps

1. Input all the required sequences. Make sure both the left and the right flanking sequences are in the same read or merged read in the fastq file.
2. Select fastq or fastq.gz files by clicking "**Choose files**".
3. Start analyzing by clicking the "**Start Analyze**" button. 
4. Download the summary file by clicking the "**Download output csv file**" button.

You can click the "**Example input**" button to get the example inputs (then select the example fastq.gz files https://github.com/pinbo/CRISjs/tree/main/example-input).

If you have a lot of fastq files, you may use the standalone command line program [CRISgo](https://github.com/pinbo/CRISgo).

### Output

The program will output a csv files reporting the percentage of intact gRNAs (no editing) and percentage of indels. It also gives the sequences and locations of the top 2 mutations and their percentages.

The first line is "Intact reference": the unedited sequences from the left flanking to the right flanking sequences. The 3rd line is the header of the output table. Explanations here:

> **fastq_file**: the fastq file name  
> **number_of_matched_reads**: number of reads that have both flanking sequences  
> **number_of_reads_with_intact_gRNA**: matched reads with intact gRNA sequence  
> **%intact_gRNA**: % of reads with intact gRNA sequence  
> **total_indel**: number of reads with an indel  
> **%indel**: %number of reads with an indel  
> **number_of_reads_with_leftSeq**: number of reads with the left flanking sequence  
> **number_of_reads_with_rightSeq**: number of reads with the right flanking sequence  
> **nleftSeq/nrightSeq**: number_of_reads_with_leftSeq / number_of_reads_with_rightSeq  
> **#1_mutation**: the most frequenct mutations (SNP or indel)  
> **#1_count**: number of reads with #1_mutation  
> **#1_\%**: #1_count / number_of_matched_reads * 100  
> **#1_seq**: the #1 mutation sequence from the left flanking to the right flanking  
> **#1_ref**: the reference allele  
> **#1_alt**: the mutation allele  
> **#1_bp_left_of_PAM**: distance from the PAM sequence  
> **#2_mutation**: the 2nd most frequenct mutations (SNP or indel)  
> **#2_count**: number of reads with #2_mutation  
> **#2_\%**: #2_count / number_of_matched_reads * 100  
> **#2_seq**: the #2 mutation sequence from the left flanking to the right flanking  
> **#2_ref**: the reference allele  
> **#2_alt**: the mutation allele  
> **#2_bp_left_of_PAM**: distance from the PAM sequence



