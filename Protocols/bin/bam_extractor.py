#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 4/11/20
import pysam
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-bam", required=True, help="Input file (bam format)")
    parser.add_argument("-o", "--output-bam", required=True, help="Save the output to")
    parser.add_argument("-m", "--modes", nargs="+", required=True, help="Modes: PE2SE_R1, PE2SE_R2, SE")
    parser.add_argument("-f", "--chr-filters", nargs="*", help="Chromosome with names in this list will be ignored")
    parser.add_argument('-r', '--reverse-complement', action='store_true', dest='seq_reverse_complement',
                        required=False, default=False,
                        help='Set this switch if reads in this library represent the reverse complement of nascent '
                             'RNAs, like PROseq')

    args = parser.parse_args()
    filters = set(args.chr_filters) if args.chr_filters is not None else set()
    supported_modes = ("PE2SE_R1", "PE2SE_R2", "SE")

    PE2SE_R1 = False
    PE2SE_R2 = False
    SE = False
    for mode in args.modes:
        if mode == "PE2SE_R1":
            PE2SE_R1 = True
        elif mode == "PE2SE_R2":
            PE2SE_R2 = True
        elif mode == "SE":
            SE = True

    assert PE2SE_R1 + PE2SE_R2 != 2, "Cannot set PE2SE_R1 and PE2SE_R2 at the same time"
    if SE:
        assert args.seq_reverse_complement, "If you want to use this script with se libraries, rc must be set."
    source_bam = pysam.AlignmentFile(args.input_bam)
    with pysam.AlignmentFile(args.output_bam, "wb", template=source_bam) as dest_bam:
        for read in source_bam:
            if read.reference_name in filters:
                continue
            if PE2SE_R1:
                if read.is_read1:
                    read.is_read1 = False
                    read.is_read2 = False
                    if args.seq_reverse_complement:
                        read.is_reverse = not read.is_reverse
                else:
                    continue
            elif PE2SE_R2:
                if read.is_read2:
                    read.is_read1 = False
                    read.is_read2 = False
                    if args.seq_reverse_complement:
                        read.is_reverse = not read.is_reverse
                else:
                    continue
            elif SE:
                if args.seq_reverse_complement:
                    read.is_reverse = not read.is_reverse
            dest_bam.write(read)
