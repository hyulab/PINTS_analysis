#!/bin/bash
set -o errexit;
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 2020-05-27
#

# This file is modified from the original script from FivePrime
# it skips the downsampling process
# the original script takes 7 parameters and the new script takes 6,
# here are the map
# original_arg  new_arg description
# 1             1       The input bam file
# 2             2       The name of the output directory (will create this directory if it does not exist).
# 3             -       The number of reads to downsample to (most of the analysis in our paper used 20000000)
# 4             3       Min Value
# 5             4       Min density rise
# 6             5       Min pos with data
# 7             6       Min Sum
INPUT=$1
OUTPUT_DIR=$2
MinVal=$3
min_density_rise=$4
min_pos_with_data=$5
min_sum=$6

machine_info=$(hostname -f)
if [ "$machine_info" == "cbsuhy01" ]; then
  proot="/local/storage/ly349"
else
  proot="/fs/cbsuhy01/storage/ly349"
fi;

CTSS=${OUTPUT_DIR}/paraclu.ctss

echo "Make directory"
if [ ! -d "$OUTPUT_DIR" ]; 
then
    mkdir $OUTPUT_DIR
fi
echo $OUTPUT_DIR

echo "Make CTSS"
# ./make_ctss3.sh $DS $CTSS
${proot}/projects/peakcalling/tools/FivePrime-master/PeakCallingPipeline/make_ctss3.sh $1 $CTSS

echo "Run Paraclu"
PARACLU_OUT=${OUTPUT_DIR}/paraclu
${proot}/projects/peakcalling/tools/FivePrime-master/PeakCallingPipeline/RunParaclu.sh $CTSS $MinVal $PARACLU_OUT $min_density_rise $min_pos_with_data $min_sum



