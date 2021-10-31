#!/usr/bin/env bash
# This script downloads both pl and mn bigwig files from encode server,
# it takes three required parameters
# 1. url(s) to the pl bw file(s)
# 2. url(s) to the mn bw file(s)
# 3. output name
# If multiple files are provided, they will be merged as one sample
# The bigwig (wig) files from ENCODE CAGE and RAMPAGE experiments were
# generated by STAR, which performs RPM normlization on the counts
# This script will try to convert the normalized counts back to approximate raw counts
# by first finding the smallest value (denote as `s`) in the wig track, 
# and assume it stands for one raw read (s=1_000_000*1/X, where 
# X is the total uniquely mapped reads in the lib), so for all other normalized 
# values (denote as `c`), they can be converted back by c / s
# Note: the min count in the raw library is not necessary to be 1 (assume it's k), in that
# case all converted "raw" counts will be `1/k` to the real raw counts
set -o errexit;
export PATH=/programs/kentUtils/bin:$PATH

if (($# == 0)); then
  echo "Please pass argumensts -p <pl1><pl2>... -m <mn1><mn2>.. -n name"
  exit 2
fi

while getopts ":p:m:n:" opt; do
  case $opt in
    p)
      PLS=($OPTARG)
      ;;
    m)
      MNS=($OPTARG)
      ;;
    n)
      NAME=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

echo "Script started"
date
source /home/ly349/miniconda3/etc/profile.d/conda.sh
conda activate wiggletools

if [ ! -f promoters_1kb_tss_centered.bed ];
then
    cp /fs/cbsuhy01/storage/ly349/refs/annotations/human/hg38/segmentation/promoters_1kb_tss_centered.bed .
fi

if [ ! -f hg38.chrom.sizes ];
then
    # fetch chromosome size
    fetchChromSizes hg38 > hg38.chrom.sizes
fi

ITER=0
pl_intermediate_arr=()
mn_intermediate_arr=()
echo "Started downloading files"
date
for pls in ${!PLS[@]}; do 
    plf=${PLS[$pls]}
    mnf=${MNS[${ITER}]}
    
    # download files to local
    if [ ! -f $(basename $plf) ];
    then
        wget $plf
    fi;
    bigWigToBedGraph $(basename $plf) ${NAME}_${ITER}_pl.bg
    if [ ! -f $(basename $mnf) ];
    then
        wget $mnf
    fi;
    bigWigToBedGraph $(basename $mnf) ${NAME}_${ITER}_mn.bg
    
    min=$(cat ${NAME}_${ITER}_pl.bg ${NAME}_${ITER}_mn.bg | awk 'NR==1{min = $4 + 0; next} {if ($4 < min) min = $4;} END {print min}')
    echo $min

    awk -v scale="$min" 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,int($4/scale)}' ${NAME}_${ITER}_pl.bg | grep -v "_" | grep -v "chrEBV" > ${NAME}_${ITER}_pl.converted.bg
    awk -v scale="$min" 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,-1*int($4/scale)}' ${NAME}_${ITER}_mn.bg | grep -v "_" | grep -v "chrEBV" > ${NAME}_${ITER}_mn.converted.bg
    
    # convert file back to bigwig
    bedGraphToBigWig ${NAME}_${ITER}_pl.converted.bg hg38.chrom.sizes ${NAME}_${ITER}_pl.bw
    bedGraphToBigWig ${NAME}_${ITER}_mn.converted.bg hg38.chrom.sizes ${NAME}_${ITER}_mn.bw
    pl_intermediate_arr+=" ${NAME}_${ITER}_pl.bw"
    mn_intermediate_arr+=" ${NAME}_${ITER}_mn.bw"
    ITER=$(expr $ITER + 1)
done
echo "Files downloads"
date

echo "Started merging files"
date
wiggletools sum $pl_intermediate_arr > ${NAME}_pl.wig
wigToBigWig ${NAME}_pl.wig hg38.chrom.sizes ${NAME}_pl.bw
wiggletools sum $mn_intermediate_arr > ${NAME}_mn.wig
wigToBigWig ${NAME}_mn.wig hg38.chrom.sizes ${NAME}_mn.bw
echo "Files merged"
date

echo "Clean up"
rm ${NAME}_*.bg
rm ${NAME}_*.wig
rm *.bigWig  # downloaded bigwigs
