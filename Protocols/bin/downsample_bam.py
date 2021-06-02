#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@gmail.com)
# Created on: 3/1/20
# @todo: samtools-based downsampling engine

import pysam
import argparse
import logging
import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from multiprocessing import Pool
from run import run_command

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
logging.basicConfig(format='%(name)s - %(asctime)s - %(levelname)s: %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO,
                    handlers=[
                        logging.FileHandler(os.path.join(os.getcwd(), 'Downsampling_%s.log' % timestamp)),
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("PINTS - Downsampler")


def _de_picard(bam, out, rate, pe_tag=0, picard_executable="/programs/picard-tools-2.19.2/picard.jar"):
    seed = int(datetime.timestamp(datetime.now()))
    std, err, rc = run_command("java -jar %s DownsampleSam I=%s O=%s RANDOM_SEED=%d PROBABILITY=%f STRATEGY=Chained" % (
        picard_executable, bam, out, seed, rate))
    assert rc == 0, "Picard error."
    t = None
    if rc == 0:
        pysam.index(out)
        tmp = pysam.AlignmentFile(out)
        count = tmp.count()
        count = count / 2 if pe_tag else count
        t = (bam, out, seed, rate, count)
        logger.info("%s\t%s\t%d\t%f\t%d" % t)
        tmp.close()
    else:
        logger.error(err)
    return t


def _de_samtools(bam, out, rate, pe_tag=0, samtools_executable="~/bin/samtools"):
    seed = int(datetime.timestamp(datetime.now()))
    std, err, rc = run_command("samtools view -bs %d.%s %s > %s" % (
        samtools_executable, seed, rate, bam, out))
    t = None
    if rc == 0:
        pysam.index(out)
        tmp = pysam.AlignmentFile(out)
        count = tmp.count()
        count = count / 2 if pe_tag else count
        t = (bam, out, seed, rate, count)
        logger.info("%s\t%s\t%d\t%f\t%d" % t)
        tmp.close()
    else:
        logger.error(err)
    return t


def _generate_clean_bam(bf, sample_label, layout, save_to, mapq, trace, filters=[]):
    logger.info("Generating clean BAM file for sample %s" % sample_label)
    total_counts = 0
    counts_pass_filter = 0
    pe_tag = 0
    if not os.path.exists(bf + ".bai"):
        pysam.index(bf)
    bam = pysam.AlignmentFile(bf)

    if len(filters) == 0 and mapq == 0:
        total_counts = bam.count()
        counts_pass_filter = total_counts
        cb_fn = bf
    else:
        cb_fn = os.path.join(save_to, sample_label + "_clean.bam")
        if os.path.exists(cb_fn):
            if not os.path.exists(cb_fn + ".bai"):
                pysam.index(cb_fn)
            try:
                clean_bam = pysam.AlignmentFile(cb_fn)
            except OSError as e:
                logger.error(bf)
                logger.error(e)
                sys.exit(1)
            assert clean_bam.check_index(), "%s is missing index" % cb_fn
            total_counts = bam.count()
            counts_pass_filter = clean_bam.count()
        else:
            clean_bam = pysam.AlignmentFile(cb_fn, "wb", template=bam)

            for read in bam:
                total_counts += 1
                watch = 0
                for keyword in filters:
                    if read.reference_name.find(keyword) != -1:
                        watch = 1
                if watch:
                    continue
                if read.mapq < mapq:
                    continue
                if read.is_secondary:
                    read.is_secondary = False
                counts_pass_filter += 1
                clean_bam.write(read)
        clean_bam.close()
        bam.close()
        if not os.path.exists(cb_fn + ".bai"):
            pysam.index(cb_fn)
    if layout == "pe":
        pe_tag = 1
        total_counts /= 2
        counts_pass_filter /= 2
    # counts[i] = counts_pass_filter
    # raw_counts[i] = total_counts
    return trace, cb_fn, counts_pass_filter, total_counts, pe_tag


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-bams", nargs="+", required=True)
    parser.add_argument("-l", "--layout-labels", nargs="+", required=True)
    parser.add_argument("-c", "--downsample-to", type=int, required=False,
                        help="If this one is set, the input bam will be downsampled to this depth. "
                             "Only applicable when there's one bam provided.")
    parser.add_argument("-k", "--keep-all", action="store_true", default=False,
                        help="Keep all reads before downsampling, filters and mapq will be applied after downsampling.")
    parser.add_argument("-f", "--filters", nargs="*")
    parser.add_argument("-r", "--repeats", type=int, default=3)
    parser.add_argument("-m", "--mapq", type=int, default=0, help="Min mapq, if you want to keep all reads. Default: 0")
    parser.add_argument("-t", "--thread", type=int, default=1, help="Threads. Default: 1")
    parser.add_argument("-o", "--save-to", required=False)

    args = parser.parse_args()
    logger.info("Parameters: %s" % " ".join(sys.argv))
    raw_counts = np.ones(len(args.input_bams), dtype=int)
    counts = np.ones(len(args.input_bams), dtype=int)
    sample_labels = tuple(map(lambda x: os.path.splitext(os.path.split(x)[1])[0], args.input_bams))
    pe_tags = np.zeros(len(args.input_bams), dtype=int)

    filters = []
    if args.filters is not None and len(args.filters) > 0:
        filters = set(args.filters)

    clean_bams = dict()

    for bf in args.input_bams:
        if not os.path.exists(bf):
            logger.error("Cannot read file %s" % bf)
            import sys

            sys.exit(-1)

    if not args.keep_all:
        para_args = []
        for i, bf in enumerate(args.input_bams):
            # bf, sample_label, layout, save_to, mapq, filters=[]
            para_args.append((bf, sample_labels[i], args.layout_labels[i], args.save_to, args.mapq, i, filters))

        depths = []
        with Pool(args.thread) as pool:
            depths = pool.starmap(_generate_clean_bam, para_args)

        for res in depths:
            index, cb_fn, counts_pass_filters, total_counts, pe_tag = res
            counts[index] = counts_pass_filters
            raw_counts[index] = total_counts
            pe_tags[index] = pe_tag
            clean_bams[index] = cb_fn
    else:
        for i, bf in enumerate(args.input_bams):
            counts_pass_filter = 0
            pe_tag = 0
            if not os.path.exists(bf + ".bai"):
                pysam.index(bf)
            with pysam.AlignmentFile(bf) as bam:
                total_counts = bam.count()
            if args.layout_labels[i] == "pe":
                pe_tags[i] = 1
                total_counts /= 2
            counts[i] = total_counts
            raw_counts[i] = total_counts
            clean_bams[i] = bf

    if len(args.input_bams) > 1:
        bottom_line = counts.min()
    elif len(args.input_bams) == 1 and args.downsample_to > 0:
        bottom_line = args.downsample_to
    else:
        sys.exit("You need to provide either multiple bam files or one bam file with a downsampling target")

    sampling_rates = bottom_line / counts
    logger.info("Coverages:")
    for i, sample in enumerate(sample_labels):
        norm_tag = "" if counts[i] > bottom_line else "*"
        logger.info("%s: %d (total mapped: %d, %f%s)" % (sample, counts[i], raw_counts[i], sampling_rates[i], norm_tag))

    logs = []
    para_args = []
    for r in range(args.repeats):
        for i, bf in clean_bams.items():
            save_to = os.path.join(args.save_to, sample_labels[i] + "_ds%d.bam" % (r + 1))
            if not os.path.exists(save_to):
                para_args.append((bf, save_to, sampling_rates[i], pe_tags[i]))
    with Pool(args.thread) as pool:
        logs = pool.starmap(_de_picard, para_args)
        # logs.append(_de_picard(bf, save_to, sampling_rates[i], pe_tags[i]))
    log_filename = "Downsampling_%s.csv" % timestamp
    pd.DataFrame(logs, columns=("Input", "Output", "Seed", "Sampling rate", "# Reads")).to_csv(log_filename,
                                                                                               index=False)

    if args.keep_all:
        para_args = []
        for r in range(args.repeats):
            for i, bf in clean_bams.items():
                saved_to = os.path.join(args.save_to, sample_labels[i] + "_ds%d.bam" % (r + 1))
                save_to = os.path.join(args.save_to, sample_labels[i] + "_ds%d_clean.bam" % (r + 1))
                para_args.append((saved_to, sample_labels[i] + "_ds%d" % (r + 1), args.layout_labels[i], args.save_to,
                                  args.mapq, i, filters))
        with Pool(args.thread) as pool:
            logs = pool.starmap(_generate_clean_bam, para_args)
    logger.info("Logs wrote to %s" % log_filename)
