#!/bin/bash

ref=$1 # indexed reference fasta file
input_r1=$2 # Read 1
input_r2=$3 # Read 2
output=$4 # output name

## steps

# 1.  make sure the reference is indexed
#bwa index $ref

# 2. align reads to reference
bwa mem $ref $input_r1  $input_r2 |samtools view -Sb - >${output}.bam

# 3. sort bam file
samtools sort ${output}.bam -o ${output}.sorted.bam
rm ${output}.bam

# 4. index sorted bam file
samtools index ${output}.sorted.bam


# 5. check coverage
echo $output >> log.txt
samtools depth ${output}.sorted.bam | awk '$3>10' | cut -f1 | sort | uniq >> log.txt
