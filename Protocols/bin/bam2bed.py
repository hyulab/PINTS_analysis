#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Li Yao
# @Date: 2/16/21
import multiprocessing
import os
import argparse

_BRICKS = {
    "bamtobed": "bedtools bamtobed -split -i",
    "flip_strand": "awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{if ($6 == \"+\") print $1,$2,$3,$4,$5,\"-\"; else print $1,$2,$3,$4,$5,\"+\"}'",
    "extract_read1": "samtools view -bf 64",
    "extract_read2": "samtools view -bf 128"
}
_COMMANDS_POOL = {
    "SE": _BRICKS["bamtobed"] + " %s | sed 's/U13369.1/chrrRNA/g' - > %s",
    "SEr": _BRICKS["bamtobed"] + " %s | " + _BRICKS["flip_strand"] + " | sed 's/U13369.1/chrrRNA/g' - > %s",
    "PE1": _BRICKS["extract_read1"] + " %s | " + _BRICKS["bamtobed"] + " stdin | sed 's/U13369.1/chrrRNA/g' - > %s",
    "PE1r": _BRICKS["extract_read1"] + " %s | " + _BRICKS["bamtobed"] + " stdin | " + _BRICKS[
        "flip_strand"] + " | sed 's/U13369.1/chrrRNA/g' - > %s",
    "PE2": _BRICKS["extract_read2"] + " %s | " + _BRICKS["bamtobed"] + " stdin | sed 's/U13369.1/chrrRNA/g' - > %s",
    "PE2r": _BRICKS["extract_read2"] + " %s | " + _BRICKS["bamtobed"] + " stdin | " + _BRICKS[
        "flip_strand"] + " | sed 's/U13369.1/chrrRNA/g' - > %s",
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-bams", nargs="+", required=True, help="Input file (bam format)")
    parser.add_argument("-o", "--save-to", required=True, help="Save the output to")
    parser.add_argument("-l", "--layouts", nargs="+", required=True, help="Layouts: SE, SEr, PE1, PE1r, PE2, PE2r")
    parser.add_argument("-p", "--processes", required=False, default=16, type=int,
                        help="Subprocesses allowed to create")

    args = parser.parse_args()

    assert len(args.input_bams) == len(args.layouts), "# of input bams must match # of layouts % (%d, %d)" % (
        len(args.input_bams), len(args.layouts))
    assert os.path.exists(args.save_to), "The output path doesn't exist"
    jobs = []
    for i, bam_file in enumerate(args.input_bams):
        base_name = os.path.splitext(os.path.split(bam_file)[1])[0]
        jobs.append(_COMMANDS_POOL[args.layouts[i]] % (bam_file,
                                                       os.path.join(args.save_to, base_name + ".bed")))

    rcs = []
    with multiprocessing.Pool(args.processes) as pool:
        rcs = pool.map(os.system, jobs)

    if all([rc == 0 for rc in rcs]):
        print("All files converted successfully.")
