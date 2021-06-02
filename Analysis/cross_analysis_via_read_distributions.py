#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Li Yao
# @Date: 2/28/21
import logging
import argparse
import os
import json
import pybedtools
import numpy as np
import pandas as pd
from .assay_specificity import get_reads_from_abundant_transcripts
from .utils import run_command
from configparser import ConfigParser

logging.basicConfig(format='%(name)s - %(asctime)s - %(levelname)s: %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO,
                    handlers=[
                        logging.FileHandler(os.path.join(os.getcwd(), 'AssayCrossAnalysisReads.log')),
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("PINTS - AssayCrossAnalysisViaReads")


def _ref_region_mcnemar(raw_df, condition_a, condition_b, threshold=np.log(6)):
    """
    Packed function for McNemar's test for comparing coverage among reference regions

    Parameters
    ----------
    raw_df : pd.DataFrame
        DataFrame where cells store coverages in a region (row) by an assay (col)
    condition_a : str
        Name of a column in the DataFrame
    condition_b : str
        Name of a column in the DataFrame
    threshold : float
        Read coverages higher than this value will be considered as covered

    Returns
    -------
    test_statistic : float

    pvalue : float
        p-value of the null hypothesis of equal marginal distributions
    df_matched_pairs : pd.DataFrame
        DataFrame records the matched pairs
    """
    from statsmodels.stats.contingency_tables import mcnemar
    new_df = raw_df.loc[:, (condition_a, condition_b)].copy()
    df_matched_pairs = pd.DataFrame(0, columns=("T", "F"),
                                    index=("T", "F"))
    a = new_df[condition_a].apply(lambda x: "T" if x > threshold else "F")
    b = new_df[condition_b].apply(lambda x: "T" if x > threshold else "F")
    x = pd.DataFrame({"A": a, "B": b})
    for (a, b), sdf in x.groupby(["A", "B"]):
        df_matched_pairs.loc[a, b] = sdf.shape[0]

    nd = df_matched_pairs.loc["T", "F"] + df_matched_pairs.loc["F", "T"]
    min_count = df_matched_pairs.min().min()
    correction = True if min_count < 5 else False
    exact = True if nd >= 20 else False
    bunch = mcnemar(df_matched_pairs, exact=exact, correction=correction)
    logger.info(
        f"McNemar's test for {condition_a} vs {condition_b} (Continuous correction: {correction}, Exact test: {exact}), p-vale: {bunch.pvalue}")
    return bunch.statistic, bunch.pvalue, df_matched_pairs


def rRNA_spec(beds_for_per_rrna, beds_for_coverage, pos_ref_bed, abundant_rna, save_to, global_registry, local_registries, tmp_dir="."):
    """
    Zoom-in study to show the importance of depleting rRNA

    Parameters
    ----------
    beds_for_per_rrna : dict
        dict of bed files, each bed file represents all mapped reads from a corresponding bam file
    beds_for_coverage : list or tuple
        List of bed files, each bed file contains a wildcard
    pos_ref_bed : str
        Path to a bed file, which defines true enhancer elements
    abundant_rna : str
        Path to a bed file, which defines genomic regions for abundant RNAs
    save_to : str
        Path to write outputs
    global_registry : int
        Global registration for outputs
    local_registries :
        Local registrations (order) for outputs
    tmp_dir : str
        Temporary files will be written into this path

    Returns
    -------
    r1 : str
        Path to a csv file which stores # of reads from abundant transcripts
    r2_1, r2_2 : list of str
        - Path to a csv file which stores % of true enhancer coverage
        - Path to a csv file which stores pairwise enhancer coverage
    """
    from .assay_sensitivity import bidirectional_coverage_among_ref_regions

    r1 = get_reads_from_abundant_transcripts(bed_files=beds_for_per_rrna,
                                             abundant_rna_file=abundant_rna,
                                             save_to=save_to,
                                             global_registry=global_registry,
                                             local_registry=local_registries[0],
                                             tmp_dir=tmp_dir)
    outputs = map(run_command, [f"wc -l {v}" for k, v in beds_for_per_rrna.items()])
    read_counts = []
    for out in outputs:
        stdout, stderr, rc = out
        assert rc == 0
        rc, _ = stdout.split(" ")
        rc = int(rc)
        read_counts.append(rc)
    df = pd.read_csv(r1, index_col=0)
    df = df.append(pd.DataFrame([read_counts, ], columns=df.columns))
    df.to_csv(r1)

    r2 = bidirectional_coverage_among_ref_regions(beds_for_coverage, pos_ref_bed, save_to, global_registry,
                                                  local_registries[1:],
                                                  detectable_threshold=(5,), n_samples=3, n_threads=16)
    return r1, r2[0], r2[1]


def main(beds_for_per_rrna, beds_for_coverage, data_save_to, true_enhancers,
         abundant_rna, data_prefix):
    analysis_summaries = {
        "rrna_zoomin": [],
    }
    analysis_summaries["rrna_zoomin"].extend(rRNA_spec(beds_for_per_rrna=beds_for_per_rrna,
                                                        beds_for_coverage=beds_for_coverage,
                                                        pos_ref_bed=true_enhancers, abundant_rna=abundant_rna,
                                                        save_to=data_save_to, global_registry=global_registry,
                                                        local_registries=("BruUVPerrRNA",
                                                                            "BruUVTECov",
                                                                            "BruUVTECovPer")))
    with open(os.path.join(data_save_to, f"{data_prefix}_summary.json"), "w") as fh:
        json.dump(analysis_summaries, fh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-save-to", required=True, help="Save output data to")
    parser.add_argument("--te-bed", required=False,
                        help="True enhancer loci in a bed file")
    parser.add_argument("--abundant-rna", required=False,
                        help="Abundant RNA loci defined in a bed file")
    parser.add_argument("--data-prefix", required=False, default="cross_analysis", help="Prefix for all outputs")
    parser.add_argument("--global-registry", type=int, help="Global registry", default=7)
    parser.add_argument("--config-file", default="config.conf", help="Configuration file")
    parser.add_argument("--tmp-dir", default="./tmp", help="All temporary files will be written into this folder")

    args = parser.parse_args()

    assert os.path.exists(args.config_file)

    if not os.path.exists(args.tmp_dir):
        os.mkdir(args.tmp_dir)

    cfg = ConfigParser()
    cfg.optionxform = str
    cfg.read(args.config_file)
    collapsed_assays = cfg.get("assays", "plot_order_simplified").split("|")

    global tmp_dir, n_threads, n_samples, global_registry
    tmp_dir = args.tmp_dir
    pybedtools.set_tempdir(tmp_dir)
    n_threads = int(cfg.get("global", "n_threads"))
    n_samples = int(cfg.get("assays", "n_downsamples"))
    global_registry = args.global_registry

    all_reads_beds = dict()
    ds_beds = dict()

    import socket

    server1_home_dir = "/fs/cbsuhy01/storage/ly349/" if socket.gethostname().find(
        "cbsuhy01") == -1 else "/local/storage/ly349/"
    server2_home_dir = "/fs/cbsuhy02/storage/ly349/" if socket.gethostname().find(
        "cbsuhy02") == -1 else "/local/storage/ly349/"
    bioq_dir = os.path.join(server2_home_dir, "BioQueue/workspace")

    from .utils import load_bioq_datasets

    all_reads_beds = load_bioq_datasets("rRNA_zoom_in_bed", bioq_dir, cfg_file=args.config_file)
    ds_beds = load_bioq_datasets("rRNA_zoom_in_ds_bed", bioq_dir, cfg_file=args.config_file)
    main(beds_for_per_rrna=all_reads_beds, beds_for_coverage=ds_beds,
         data_save_to=args.data_save_to, fig_save_to=args.fig_save_to, abundant_rna=args.abundant_rna,
         true_enhancers=args.te_bed, data_prefix=args.data_prefix)
