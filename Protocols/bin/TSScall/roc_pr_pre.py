#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 6/1/20
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
import sys
import os
import re
import argparse
import socket
import numpy as np
import pandas as pd

home_dir = "/fs/cbsuhy01/storage/ly349" if socket.gethostname().find(
    "cbsuhy01") == -1 else "/local/storage/ly349"
sys.path.append(os.path.join(home_dir, "projects/peakcalling/remote_debug/helpers"))
from run import run_command
from uuid import uuid4
from multiprocessing import Pool
from decimal import Decimal, getcontext


def _1dgrid(par, par_range, log_to, threads=16):
    commands = []
    par_output_mapping = []
    for threshold in par_range:
        output_prefix = str(uuid4())
        commands.append(par.format(threshold=threshold, output=output_prefix))
        par_output_mapping.append((threshold, output_prefix))

    with Pool(threads) as pool:
        results = pool.map(run_command, commands)

    for result in results:
        assert result[2] == 0, result[1]
    log_df = pd.DataFrame(par_output_mapping, columns=("Threshold", "Output_prefix"))
    log_df.to_csv(log_to, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--input-pl", required=True, help="Input file (forward strand)")
    parser.add_argument("-m", "--input-mn", required=True, help="Input file (reverse strand)")
    parser.add_argument("-c", "--chrom-size", required=True, help="Chromosome sizes")
    parser.add_argument("-t", "--trails", default=101, type=int)
    parser.add_argument("-s", "--TSScall-script", required=False,
                        default=os.path.join(home_dir, "projects/peakcalling/tools/TSScall/TSScall.py"),
                        help="Path to TSScall")
    parser.add_argument("-b", "--bidirectional-script", required=False,
                        default=os.path.join(home_dir,
                                             "projects/peakcalling/remote_debug/helpers/TSScall/TSScall_bidirectional_peaks.py"),
                        help="Path to Bidirectional pair generator")
    parser.add_argument("-a", "--annotation", required=False,
                        default=os.path.join(home_dir, "refs/annotations/human/hg38/gencode.v24.annotation.gtf"),
                        help="Path to gene annotations (gtf)")

    args = parser.parse_args()
    reg = r"hg(\d+)"
    matches = re.finditer(reg, args.input_pl)
    genome_release = None
    for i, match in enumerate(matches, start=1):
        genome_release = match.group()

    assert genome_release is not None, "Cannot extract genome release from file name"

    cmd = "python %s --fdr {threshold} -a %s --detail_file {output}.txt " \
          "--cluster_bed {output}.cluster.bed %s %s %s {output}.bed" % (
              args.TSScall_script, args.annotation, args.input_pl, args.input_mn, args.chrom_size)

    n_small_pvals = int(args.trails * 0.25)
    n_exp_pvals = int(args.trails * 0.25)
    small_pvals = np.linspace(0, 0.1, n_small_pvals).tolist()
    exp_pvals = [10 ** i for i in range(-1 * n_exp_pvals, -1)]
    # make sure enough precision
    getcontext().prec = n_exp_pvals + 1
    actual_values = len(set(small_pvals + exp_pvals))
    n_general_pvals = args.trails - actual_values + 1
    general_pvals = np.linspace(0.1, 1, n_general_pvals).tolist()
    search_range = set(small_pvals + exp_pvals + general_pvals)
    _1dgrid(cmd, search_range, log_to="TSScall.csv")
    cmd = "for txt in *.txt; do python %s " \
          "--include-convergent --input ${txt} --save-to ${txt/txt/bid.bed}; done" % args.bidirectional_script
    run_command(cmd)
