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
<input type="checkbox" id="check2" name="check2" />
<label for="check2">Find wheat genes that match given At/Os genes (e.g. At/Os genes -> wheat genes)</label>

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
<table id="datatable" style="font-size: 11px;" align="left">
<thead id="thead">
    <tr>
        <th>WheatGeneID</th>
        <th>Best Wheat matches</th>
        <th>Wheat %identity</th>
        <th>Best At matches</th>
        <th>At %identity</th>
        <th>At align length</th>
        <th>At description</th>
        <th>Best Os matches</th>
        <th>Os %identity</th>
        <th>At align length</th>
        <th>Os description</th>
    </tr>
</thead>
    <tbody id="tbody"></tbody>
</table>

<div id="gap1" style="margin-top: 300px;"></div>

**Update**

- 2024-09-18: modify the `blastp` method (`-seg yes`) to match Ensembl blast output (only affect some top hits of Arabidopsis).
- 2024-09-18: add some low confidence genes that are hits of high confidence genes. For example, the B homeolog of PLATZ-A1 (TraesCS6A02G156600) is a low confidence gene.
- 2024-11-01: add alignment length from BLAST for At and Os hits. Without the alignment length, we cannot tell which wheat gene is best At/Os homolog.


**Methods**

Here are the commands I used for preparing homeologs and the best hits in Arabidopsis and rice. Arabidopsis and rice seequnces were downloaded from [Ensembl Plants](https://plants.ensembl.org/index.html). Kronos cDNAs were downloaded from [Zenodo](https://zenodo.org/records/11106422). CS IWGSC annotation v1.1 HC cDNAs were downloaded from [Wheat URGI](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/).

``` sh
## homeolog search by self blast
### blast self
blastn -task blastn -db ../blastdb/Kronos.v1.0.all.cds.fa -query ../blastdb/Kronos.v1.0.all.cds.fa -outfmt "6 std qlen slen" -perc_identity 90 -word_size 20 -num_threads 40 -out out_Kronos_v1.0_cdna_self_wordsize20.txt &
blastn -task blastn -query /Users/galaxy/blastdb/IWGSC_v1.1_HC_20170706_cds.fasta -db /Users/galaxy/blastdb/IWGSC_v1.1_HC_20170706_cds.fasta -outfmt "6 std qlen slen" -perc_identity 90 -word_size 20 -num_threads 40 -out out_CS_v1.1_HC_self_wordsize20.txt &

### organize results: self3, use 0.6 length as cut point, due to splice variation
gawk '$4>$13*0.6 {split($1,aa,"."); split($2,bb,"."); qq=aa[1]; ss=bb[1]; if(!(qq"\t"ss in cc)) {cc[qq"\t"ss]++; printf("%s\t%s\t%.f\t%s\n",qq,ss,$3,$4)} }' out_CS_v1.1_HC_self_wordsize20.txt > filtered_CS_v1.1_HC_self3.txt
gawk '$4>$13*0.6 {split($1,aa,"."); split($2,bb,"."); qq=aa[1]; ss=bb[1]; if(!(qq"\t"ss in cc)) {cc[qq"\t"ss]++; printf("%s\t%s\t%.f\t%s\n",qq,ss,$3,$4)} }' out_Kronos_v1.0_cdna_self_wordsize20.txt > filtered_Kronos_self3.txt

## blast Os and At
# update 2024-09-18: add '-seg yes'
### Kronos
blastp -db ../blastdb/Arabidopsis_thaliana.TAIR10.pep.all.fa -query ../blastdb/Kronos.v1.0.all.pep.fa -outfmt "6 std qlen slen stitle" -max_target_seqs 6 -word_size 3 -num_threads 40 -out out_Kronos_v1.0_against_Arabidopsis_TAIR10_pep_wordsize3.txt -seg yes &
blastn -task blastn -db /Users/galaxy/blastdb/Oryza_sativa.IRGSP-1.0.cds.all.fa -query ../blastdb/Kronos.v1.0.all.cds.fa -outfmt "6 std qlen slen stitle" -max_target_seqs 6 -word_size 15 -num_threads 40 -out out_Kronos_v1.0_against_rice_IRGSP-1.0_cdna_wordsize15.txt &

gawk 'bb[$1]<1{bb[$1]=1; print}' out_Kronos_v1.0_against_Arabidopsis_TAIR10_pep_wordsize3.txt > top1hit_out_Kronos_v1.0_against_Arabidopsis_TAIR10_pep_wordsize3.txt
sed -i 's/ gene:/\t/g;s/ gene_symbol:/\t/g;s/ description:/\t/g;s/ \[Source/\t/g' top1hit_out_Kronos_v1.0_against_Arabidopsis_TAIR10_pep_wordsize3.txt

gawk 'bb[$1]<1{bb[$1]=1; print}' out_Kronos_v1.0_against_rice_IRGSP-1.0_cdna_wordsize15.txt > top1hit_out_Kronos_v1.0_against_rice_IRGSP-1.0_cdna_wordsize15.txt
sed -i 's/ gene:/\t/g;s/ gene_biotype:/\t/g; s/ gene_symbol:/\t/g;s/ description:/\t/g' top1hit_out_Kronos_v1.0_against_rice_IRGSP-1.0_cdna_wordsize15.txt

### CS
blastp -db ../blastdb/Arabidopsis_thaliana.TAIR10.pep.all.fa -query ../blastdb/Triticum_aestivum.IWGSC.pep.all.fa -outfmt "6 std qlen slen stitle" -max_target_seqs 6 -word_size 3 -num_threads 40 -out out_CS_v1.1_against_Arabidopsis_TAIR10_pep_wordsize3.txt -seg yes &
blastn -task blastn -db /Users/galaxy/blastdb/Oryza_sativa.IRGSP-1.0.cds.all.fa -query /Users/galaxy/blastdb/IWGSC_v1.1_HC_20170706_cds.fasta -outfmt "6 std qlen slen stitle" -max_target_seqs 6 -word_size 11 -num_threads 40 -out out_CS_v1.1_against_rice_IRGSP-1.0_cdna_wordsize11.txt &

gawk 'bb[$1]<1{bb[$1]=1; print}' out_CS_v1.1_against_Arabidopsis_TAIR10_pep_wordsize3.txt > top1hit_out_CS_v1.1_against_Arabidopsis_TAIR10_pep_wordsize3.txt
sed -i 's/ gene:/\t/g;s/ gene_symbol:/\t/g;s/ description:/\t/g;s/ \[Source/\t/g' top1hit_out_CS_v1.1_against_Arabidopsis_TAIR10_pep_wordsize3.txt

gawk 'bb[$1]<1{bb[$1]=1; print}' out_CS_v1.1_against_rice_IRGSP-1.0_cdna_wordsize11.txt > top1hit_CS_v1.1_against_rice_IRGSP-1.0_cdna_wordsize11.txt
sed -i 's/ gene:/\t/g; s/ gene_biotype:/\t/g; s/ gene_symbol:/\t/g; s/ description:/\t/g'  top1hit_CS_v1.1_against_rice_IRGSP-1.0_cdna_wordsize11.txt

## then I prepared a sqlite3 database for the webtool

```

**Acknowledgment**

- [IWGSC](https://www.wheatgenome.org/) for CS Refseq v1 assembly and annotation
- [Krasileva lab](https://krasilevalab.org/) for Kronos v1 assembly and annotation
- [Ensembl Plants](https://plants.ensembl.org/index.html) for hosting many plant genomes
- [sqlite3](https://www.sqlite.org/) for database preparation
- [sql.js](https://sql.js.org/) for using sqlite3 in the browser

<script src="/tools/sqljs/v1.10.3/sql-wasm.js"></script>
<!-- <script src="/tools/sqlite3/3.46.0/sqlite3.js"></script> -->
<script type="module" src="/libs/get-wheat-homeologs-v6.js"></script>
<script src="/libs/excellentexport.min.js"></script>
<script src="/libs/pako_inflate.min.js"></script>
