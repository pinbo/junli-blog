---
title: Process NGS of CRISPR editing
author: Junli Zhang
date: '2020-07-26'
slug: process-ngs-of-crispr-editing
categories:
  - Method
tags:
  - CRISPR
  - NGS
---

CRISPR editing is hard to predict: editing frequency varies, editing events could be multiple. So the best way to check the editing efficiency is to use next-generation sequencing (NGS). Sequencing libraries can be prepared with PCR since you know where are your targets. You just need two rounds of PCR to add barcodes and adapters to the ends of your PCR amplicons.


## First round PCR
Just regular PCR (annealing temperature is the same as without the adapters) with sequencing primers. Make sure your PCR amplicons are less than 300 bp for PE150 or 500 bp for PE250, so there are overlaps between read 1 and read 2 of the NGS results.

```
Forward (Read1): 5’-ACACTCTTTCCCTACACGACGCTCTTCCGATCT -your forward primer -3’  
Reverse (Read2): 5’-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-Your-reverse-primer -3’

94 C 5 min
35 cycles of 
    94 C 30 s
    58 C 30 s (depending on your primers’ Tm)
    72 C 40 s
72 C 5 min
12 C forever
```
## Second round PCR
Use 1 ul of 1st round PCR (you can also dilute 10x if the PCR product is strong). Then use the primers below:

```
Index 1 (XXXXXXXX: plate barcoces: i5 barcodes for example)
5'-AATGATACGGCGACCACCGAGATCTACAC-XXXXXXXX-ACACTCTTTCCCTACACGACGCTCTTCCGATCT -3'
Index 2 (NNNNNNNN: plate well barcode: i7 barcodes for example)
5'-CAAGCAGAAGACGGCATACGAGAT-NNNNNNNN-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3'

94 C 5 min
11 cycles of 
    94 C 30 s
    65 C 30 s
    72 C 40 s
72 C 5 min
12 C forever
```

After the 2nd round PCR, your amplicons will look like this (136 bp longer than your PCR target):
```
5' - AATGATACGGCGACCACCGAGATCTACAC-XXXXXXXX-ACACTCTTTCCCTACACGACGCTCTTCCGATCT
- PCR target - 
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC-NNNNNNNN-ATCTCGTATGCCGTCTTCTGCTTG - 3'
```

Use 5 ul to check on an agarose gel to make sure the amplifications are good and check whether primer dimers are strong. Then pools samples together for AMPure XP beads cleaning or gel recovery.

Then send for sequencing. We are using MiSeq NANO 500 (PE250: paired-end 250 bp, $555 per run, 1 million reads).

Sequencing results will be like (all adapters will be trimmed, only PCR amplicons):

**Read 1 (R1)**: Left primer until 150 bp (you use SE150 or PE150) or 250bp (if you use SE250 or PE250)

**Read 2 (R2)**: Right primer until 150 bp (you use SE150 or PE150) or 250bp (if you use SE250 or PE250)

## NGS result analysis

### Method 1: use `bwa mem` to map reads on the references

You need to have `bwa` and `samtools` installed on your system. Also `bcftools` if you want to call SNPs automatically.

1. Make a fasta file with all your gene sequences (whole genomic sequences or just amplicons)
2. Index the fasta file with command `bwa index your-reference.fasta`
3. Download the `bwa mem` script [here](/files/run_aln_v2.sh) and put it to somewhere.
4. Put all your fastq.gz files in a folder and open a terminal there (or `cd` there), then run  
`for i in *_R1*.fastq.gz; do base=${i%_R1*}; path/to/run_aln_v2.sh path/to/your-reference.fasta ${base}_R1*.fastq.gz ${base}_R2*.fastq.gz out_$base; done`
5. You can check which gene has more than 10 reads in the output file `log.txt`
6. Call SNPs/indels with bcftools for each bam file: `for i in *.bam; do bcftools mpileup -Ou -f path/to/your-reference.fasta $i | bcftools call -mv -Oz -o indels-calls-${i}.vcf.gz; done`  
Or call SNPs/indels for all bams at once: `bcftools mpileup -Ou -f path/to/your-reference.fasta *.bam | bcftools call -mv -Oz -o indels-calls-all.vcf.gz`
7. Then use [TASSEL](https://www.maizegenetics.net/tassel) or my script  [convert_vcf_calls_to_SNP_v3](https://github.com/pinbo/myscripts/blob/master/Go/convert_vcf_calls_to_SNP_v3.go) or other software to convert the vcf file to SNP tables.

I am not sure how big of indels `bcftools` can call. So you'd better check the bam files with [**IGV**](http://software.broadinstitute.org/software/igv/).

### Method 2: use `CRISgo` to do e-PCR directly on the fastq.gz files

This method is very simple and the idea is from the [CRISpy paper](https://www.nature.com/articles/s41598-019-40896-w). `CRISpy` uses two flanking sequences around the guide RNA to select reads that match both flanking sequences, basically an e-PCR. I rewrote [CRISgo](https://github.com/pinbo/CRISgo) with that idea. Just like `CRISpy`, `CRISgo` searches two flanking sequences of gRNA in fastq files (plain text files), and check whether the gRNA was edited by comparing the matched sequences with the wild type sequence, but `CRISgo` is much faster and gives more information about the indels, such as position and the actual base changes. Therefore, you do not need to do alignments manually to see what the indels are. Here I will use **PE250** results as an example to show you how to use `CRISgo`.

![CRISgo](/images/20200726-CRISgo.png)

1. Choose flanking sequences around the guide RNA. The two flanking sequences should meet the following criteria:
  - Inside PCR primers, do not use only the PCR primers, because your primers might amplify off-target regions. Overlapping with the PCR primers will allow checking a big region to test potential big deletions.
  - The combinations of the two flanking sequences should be unique to your target. If you working on polyploid crops like wheat, you need to make sure they are sub-genome specific. So if your primers can amplify all homeologs, do an alignment and select unique sequences to your target.
  - Length usually is 10-20bp, but has no limitations. Too long may exclude some reads due to PCR or sequencing errors.
2. If your PCR amplicon is less than 500 bp (2x of the length of reads), then your R1 and R2 will have overlap and you can use software such as [fastp](https://github.com/OpenGene/fastp) to remove adapters and connect two reads to get the whole PCR amplicons.
    - Download the [fastp](https://github.com/OpenGene/fastp#get-fastp) to your system.
    - Put all your fastq.gz files in a folder and open a terminal there (or `cd` there), then run  
    `for i in *_R1*.fastq.gz; do base=${i%_R1*}; path/to/fastp -i ${base}_R1*.fastq.gz -I ${base}_R2*.fastq.gz -o out_${base}-fail2merge_R1_001.fastq.gz -O out_${base}-fail2merge_R2_001.fastq.gz -m --merged_out out_${base}-merged_R1_001.fastq.gz; done`
    - Then run command `mkdir merged_reads; mv *merged_R1_001.fastq.gz merged_reads/`
    - Download the latest `CRISgo` release [here](https://github.com/pinbo/CRISgo/releases/latest). We will use `CRISgov5.1` as an example. 
    - Run the following command `cd merged_reads; path/to/CRISgov5.1 reference-file output_file_prefix gene_name left_flanking_sequence right_flanking_sequence gRNA_sequence`. **CRISgov5.1** will process all the fataq.gz file in the current folder and output a csv file. Read more in the [github page](https://github.com/pinbo/CRISgo).

If your guide RNA is complete in R1 or R2 or both, you can choose the two flanking sequences so both of them are within 250 bp from at least one end, then you do not need to merge R1 and R2. But merging R1 and R2 will allow you to check a bigger region if your PCR amplicon is more than 250 bp.