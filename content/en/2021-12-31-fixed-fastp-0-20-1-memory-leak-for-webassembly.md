---
title: Fixed fastp-0.20.1 memory leak for WebAssembly
author: Junli Zhang
date: '2021-12-31'
slug: fixed-fastp-0-20-1-memory-leak-for-webassembly
categories:
  - programming
tags:
  - cpp
  - webassembly
  - Valgrind
  - Memcheck
---

I compiled [fastp v0.20.1](https://github.com/OpenGene/fastp/releases/tag/v0.20.1) for WebAssembly application based on the patch file from [Biowasm](https://github.com/biowasm/biowasm/blob/main/tools/fastp/patches/0.20.1.patch). It worked well if I only process one pair of fastq files (single run). However, after converting my [web app](/apps/filter-multiple-fastq-files-with-fastp/) for processing multiple files with a loop, the memory increased about 200 Mb after each cycle, so it can only process about 10 paired-end fastq files before exceeding the 2 Gb memory limit.

Firstly, I thought it was a WebAssembly problem because this did not happen when running on my Linux Desktop. I searched and got this [article](https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html#memory-management) but did not understand it. So I had to use a dummy way: terminate the web worker after processing each file using this JavaScript file [fastp-multiplex-v3.js](/libs/fastp-multiplex-v3.js). It works but is a little slower. I knew this was not the perfect way to solve the problem and I really wanted to find out why the memory increased after each cycle. Later it occurred to me that this could be a memory leak problem of `fastp` itself, although I thought this should not be, because all the authors should have checked the memory leak of their software. I checked it anyway just in case with the program [`Valgrind`](https://valgrind.org/docs/manual/quick-start.html). And there were! The major leak is from the function `process` of files "peprocessor.cpp" and "seprocessor.cpp": there is only `initPackRepository()` at the beginning but not `destroyPackRepository()` at the end of the function, which causes 240 Mb of memory leak! Now I will NOT trust any developers and I will always check by myself. For all the fixes of fastp v0.20.1, please see my GitHub page [fastp-0.20.1-JZ](https://github.com/pinbo/fastp-0.20.1-JZ).

Here are the commands I used for the memory check: 

- For single-end fastq files: `valgrind --tool=memcheck --leak-check=yes /path/to/fastp -i 1_R1_001.fastq.gz -o test1.fq.gz`
- For paired-end fastq files: `valgrind --tool=memcheck --leak-check=yes /path/to/fastp -i 1_R1_001.fastq.gz -I 1_R2_001.fastq.gz -o test1.fq.gz -O test2.fq.gz`

I found [`Valgrind`](https://valgrind.org/docs/manual/quick-start.html) is really useful because it can point out which line in which file the bug is from. The quick start page of  [`Valgrind`](https://valgrind.org/docs/manual/quick-start.html) is a must-see. Several important notes are listed below:

- Modify the Makefile to make sure it has `-g` for including exact line numbers in the report and `-O0` (or `-O1`) for more accurate error messages;
- It's worth fixing errors in the order they are reported, as later errors can be caused by earlier errors;
- The first "by" in each error block usually is the causes and later "by"s are just due to the first one;
- Message "Mismatched free() / delete / delete []": possibly need to use `delete []` other than `delete`;
- Possibly you can first go to the end of the report "HEAP SUMMARY" and find out the most serious memory leak such as"80,000,000 bytes in 1 blocks are possibly lost in loss record 71 of 72" (the 2nd "by" is the cause in this case).