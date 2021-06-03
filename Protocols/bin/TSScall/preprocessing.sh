#!/usr/bin/env bash
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 2020-05-27
#
# Part of the PINTS project
# Copyright (C) 2020 Li Yao at the Yu Lab
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# This file is intended to generate the following inputs for TSSCalls:
# 1. FORWARD_BEDGRAPH describing the coverage of 5' ends from Start-seq reads at single-nucleotide resolution (fwd)
# 2. REVERSE_BEDGRAPH describing the coverage of 5' ends from Start-seq reads at single-nucleotide resolution (rev)
# 3. chrom.size is a text file describing the size of each chromosome/contig in the reference genome. (chromosome\tlength)
set -o errexit;
input_file=$1;
exp_design=$2;
output_prefix=$3;
filter_strs="U13369\|chrM\|Mycoplasma\|EBV";
machine_info=$(hostname -f)
if [ "$machine_info" == "cbsuhy01" ]; then
  proot="/local/storage/ly349"
else
  proot="/fs/cbsuhy01/storage/ly349"
fi;

# generate chrom size table from the bam file for following procedures
samtools view -H "${input_file}" | grep ^@SQ | sed 's/@SQ\tSN://g' | sed 's/LN://g' | sort -k 1,1 -k2,2n > chrom.size

# generate TSSs beds and whole bedGraph
if [ "$exp_design" == "R_5" ] || [ "$exp_design" == 'GROcap' ] || [ "$exp_design" == 'PROcap' ] || [ "$exp_design" == 'GROseq' ]
then
    # designs like GROcap
    bedtools bamtobed -i "${input_file}" 2> job.err | awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
        awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | \
        gzip > "${output_prefix}".dirty.bed.gz
elif [ "$exp_design" == "R_5_r" ] || [ "$exp_design" == 'PROseq' ]
then
    # designs like PROseq
    bedtools bamtobed -i "${input_file}" 2> job.err | awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
        awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,"-"}; ($6 == "-") {print $1,$3-1,$3,$4,$5,"+"}' | \
        gzip > "${output_prefix}".dirty.bed.gz
elif [ "$exp_design" == "R1_5" ]
then
    # designs like RAMPAGE
    # extract read1
    # samtools view -hbf 64 "${input_file}" > "${output_prefix}"_SR.bam; # SR for Single Reads
    python "${proot}"/projects/peakcalling/remote_debug/scripts/bam_extractor.py -i "${input_file}" \
     -o "${output_prefix}"_SR.bam -m PE2SE_R1

    bedtools bamtobed -i "${output_prefix}"_SR.bam 2> job.err | awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
        awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | \
        gzip > "${output_prefix}".dirty.bed.gz
elif [ "$exp_design" == "R2_5" ] || [ "$exp_design" == 'CoPRO' ]
then
    # designs like CoPRO
    # extract read2
    # samtools view -hbf 128 "${input_file}" > "${output_prefix}"_SR.bam;
    python ${proot}/projects/peakcalling/remote_debug/scripts/bam_extractor.py -i "${input_file}" \
     -o "${output_prefix}"_SR.bam -m PE2SE_R2

    bedtools bamtobed -i "${output_prefix}"_SR.bam 2>> job.err | awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
        awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | \
        gzip > "${output_prefix}".dirty.bed.gz
elif [ "$exp_design" == 'mNETseq' ]
then
    # designs like mNETseq
    # extract read2 and reverse it
    # samtools view -hbf 128 "${input_file}" > "${output_prefix}"_SR.bam;
    python ${proot}/projects/peakcalling/remote_debug/scripts/bam_extractor.py -i "${input_file}" \
     -o "${output_prefix}"_SR.bam -m PE2SE_R2 -r

    bedtools bamtobed -i "${output_prefix}"_SR.bam 2>> job.err | awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
        awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | \
        gzip > "${output_prefix}".dirty.bed.gz
else
    exit 1;
fi

# outputs for MACS2
zcat "${output_prefix}".dirty.bed.gz | grep '^chr' | grep ${filter_strs} -v | \
    awk '{if (index($1, "_")==0) print $0}' - | sort-bed - | gzip > "${output_prefix}".bed.gz

# generate TSSs bedGraphs
bedtools genomecov -bg -i "${output_prefix}".bed.gz -g chrom.size -strand + | sort -k1,1 -k2,2n | awk 'BEGIN{FS="\t";OFS="\t"}{printf "%s\t%d\t%d\t%d\n",$1,$2,$3,$4}' > "${output_prefix}"_TSS_plus.bedGraph
bedtools genomecov -bg -i "${output_prefix}".bed.gz -g chrom.size -strand - | sort -k1,1 -k2,2n | awk 'BEGIN{FS="\t";OFS="\t"}{printf "%s\t%d\t%d\t%d\n",$1,$2,$3,$4}' > "${output_prefix}"_TSS_minus.bedGraph
