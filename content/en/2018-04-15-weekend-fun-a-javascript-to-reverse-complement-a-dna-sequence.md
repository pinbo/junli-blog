---
title: 'Weekend fun: a javascript to Reverse Complement a DNA sequence'
author: Junli Zhang
date: '2018-04-15'
slug: weekend-fun-a-javascript-to-reverse-complement-a-dna-sequence
categories:
  - Learning
tags:
  - javascript
---

This weekend I was having fun to write a javascript script to convert DNA sequences to the reverse, complement, or reverse-complement sequences. Good websites on doing this are http://reverse-complement.com and http://www.bioinformatics.org/sms2/rev_comp.html. Below is just my simple practice.

<script src='/libs/simple-bioinformatics.js'></script>

<p id="demo">Please paste your sequences below.</p>

<textarea rows="10" cols="75" id="input"></textarea>
<br />

<button onclick="revcomp()">Submit</button>
<button onclick="clearseq()">Clear</button>
<select id="box1">
    <option value="RC">Reverse-complement</option>
    <option value="R">Reverse</option>
    <option value="C">Complement</option>
</select>
<select id="box2">
    <option value="fasta">FASTA mode</option>
    <option value="single_line">Single line mode</option>
</select>

<p>Output below</p>
<textarea rows="10" cols="75" id="output" ></textarea>
<br />
