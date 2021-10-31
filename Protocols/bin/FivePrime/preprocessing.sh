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

# This file is intended to generate bam file for FivePrime
# basically, it converts pe to se, and makes sure the precise end locates at 5' of reads
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

# generate TSSs beds and whole bedGraph
if [ "$exp_design" == "R_5" ] || [ "$exp_design" == 'GROcap' ] || [ "$exp_design" == 'PROcap' ] || [ "$exp_design" == 'GROseq' ]
then
    # designs like GROcap
    # outputs for Tfit
    cp "${input_file}" "${output_prefix}"_SR.bam
elif [ "$exp_design" == "R_5_r" ] || [ "$exp_design" == 'PROseq' ]
then
    python "${proot}"/projects/peakcalling/remote_debug/scripts/bam_extractor.py -i "${input_file}" \
     -o "${output_prefix}"_SR.bam -m SE -r
elif [ "$exp_design" == "R1_5" ]
then
    # designs like RAMPAGE
    # extract read1
    # samtools view -hbf 64 "${input_file}" > "${output_prefix}"_SR.bam; # SR for Single Reads
    python "${proot}"/projects/peakcalling/remote_debug/scripts/bam_extractor.py -i "${input_file}" \
     -o "${output_prefix}"_SR.bam -m PE2SE_R1
elif [ "$exp_design" == "R2_5" ] || [ "$exp_design" == 'CoPRO' ]
then
    # designs like CoPRO
    # extract read2
    # samtools view -hbf 128 "${input_file}" > "${output_prefix}"_SR.bam;
    python "${proot}"/projects/peakcalling/remote_debug/scripts/bam_extractor.py -i "${input_file}" \
     -o "${output_prefix}"_SR.bam -m PE2SE_R2

elif [ "$exp_design" == 'mNETseq' ]
then
    # designs like mNETseq
    # extract read2 and reverse it
    # samtools view -hbf 128 "${input_file}" > "${output_prefix}"_SR.bam;
    python "${proot}"/projects/peakcalling/remote_debug/scripts/bam_extractor.py -i "${input_file}" \
     -o "${output_prefix}"_SR.bam -m PE2SE_R2 -r

else
    exit 1;
fi
