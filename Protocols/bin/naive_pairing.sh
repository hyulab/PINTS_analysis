#!/usr/bin/env bash
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 2019-06-11
#
# PINTS: Peak Identifier for Nascent Transcripts Sequencing
# Copyright (C) 2019 Li Yao at the Yu Lab
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
#
# This program takes three parameters:
# 1. input bed
# 2. output prefix
# 3. pairing distance (optional, by default, 300bp)
set -o errexit

if [ -z "$3" ]; then
  max_pairing_distance=300;
else
  max_pairing_distance=$3;
fi
if [ -f "$1" ]; then
  awk 'BEGIN{FS="\t";OFS="\t"}{if ($6=="+") print $1,$2,$3,$4,$5,$6}' "$1" > "$2".plus.bed
  awk 'BEGIN{FS="\t";OFS="\t"}{if ($6=="-") print $1,$2,$3,$4,$5,$6}' "$1" > "$2".minus.bed
  bedtools intersect -v -a "$2".plus.bed -b "$2".minus.bed > "$2".naive.bed
  bedtools intersect -v -b "$2".plus.bed -a "$2".minus.bed >> "$2".naive.bed
  bedtools sort -i "$2".naive.bed > "$2".temp.bed
  mv "$2".temp.bed "$2".naive.bed
  # -w: Base pairs added upstream and downstream of each entry in A when searching for overlaps in B.
  # -Sm: Only report hits in B that overlap A on the opposite strand.
  bedtools window -w "$max_pairing_distance" -Sm -a "$2".naive.bed -b "$2".naive.bed | \
    awk 'BEGIN{FS="\t";OFS="\t"}{min=$2;max=$3; if($8<$2)min=$8; if($9>$3)max=$9; print $1,min,max}' | \
     sort -k1,1 -k2,2n | bedtools merge -i stdin > "$2".paired.bed
  rm "$2".plus.bed
  rm "$2".minus.bed
  rm "$2".naive.bed
fi
