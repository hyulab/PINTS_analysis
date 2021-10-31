#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 2/28/21
import logging
import argparse
import os
import json
import pybedtools
import numpy as np
import pandas as pd
from .assay_specificity import get_reads_from_abundant_transcripts
from .utils import run_command, get_file_hash, msg_wrapper
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


def rRNA_spec(beds_for_per_rrna, beds_for_coverage, pos_ref_bed, abundant_rna, save_to, global_registry, 
              local_registries, tmp_dir="."):
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
    local_registries : list/tuple of str
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


@msg_wrapper(logger)
def capping_bias(ds_beds, enhancer_file, save_to, global_registry, local_registry, n_reps=3):
    """
    Evaluating potential bias introduced by capturing only capped RNAs

    Parameters
    ----------
    ds_beds : dict
        key : str, cellLine_cappingStatus
        value : str, path to downsampled libraries with `%d` placed for downsample trails
    enhancer_file : str
        Path to a bed file, which defines enhancer loci
    save_to : str
        Save output to the directory
    global_registry : str or numeric
        Global registry
    local_registry : str
        Local registry
    n_reps : int
        Number of downsamplings performed for each library

    Returns
    -------
    final_result_file : str

    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        enhancer_bed = pybedtools.BedTool(enhancer_file)
        sub_dfs = []
        for k, v in ds_beds.items():
            for i in range(1, n_reps + 1):
                bed_file = pybedtools.BedTool(v % i)
                tmp_df = enhancer_bed.intersect(bed_file, c=True).to_dataframe()
                tmp_df.index = tmp_df["chrom"] + ":" + tmp_df["start"].map(str) + "-" + tmp_df["end"].map(str)
                tmp_df.drop(columns=["chrom", "start", "end", "name"], inplace=True)
                tmp_df.rename(columns={"score": f"{k}_{i}", }, inplace=True)
                sub_dfs.append(tmp_df)
        merged_counts_df = pd.concat(sub_dfs, axis=1)
        merged_counts_df["Capped"] = merged_counts_df.loc[:, ("K562_Capped_1", "K562_Capped_2", "K562_Capped_3")].mean(
            axis=1)
        merged_counts_df["Capped + Uncapped"] = merged_counts_df.loc[:,
                                                ("K562_RppH_1", "K562_RppH_2", "K562_RppH_3")].mean(axis=1)
        merged_counts_df.to_csv(final_result_file)
    return final_result_file

    
@msg_wrapper(logger)
def paused_polymerase_efficiency(library_beds, pause_indexes, transcript_segmentation,
                                 save_to, global_registry, local_registry, tmp_dir=".",
                                 target_transcripts=set()):
    """
    Analysis of whether run-on assays can efficiently unleash and capture paused polymerase

    Parameters
    ----------
    library_beds : dict
        key : name of the library
        value : path to the alignment in bed format
    pause_indexes : dict
        key : name of the library used for calculating pause index
        value : path to the pause index file (tab-separated)
    transcript_segmentation : str
        Path to a bed-like file which contains the segmentation info for transcripts
            * chromosome
            * start
            * end
            * transcript id (needs to be matched with the IDs used in pause_indexes files)
            * segmentation type: Promoter, Gbody, PAS
            * strand
    save_to : str
        Path to write outputs
    global_registry : int
        Global registration for outputs
    local_registry : str
        Local registration (order) for output
    tmp_dir : str
        Path to write tmp files
    target_transcripts : set
        If this set is not empty, then only transcripts in this set will be considered.

    Returns
    -------
    final_result_file : str
        Path to a csv file which stores results
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        # load pausing indexes
        pause_index_dfs = dict()
        for k, v in pause_indexes.items():
            tdf = pd.read_csv(v, sep="\t")
            tdf = tdf.loc[:, (tdf.columns[0], "treatment_tags/ Pausing Ratio")]
            tdf.columns = ("Transcript ID", k)
            tdf.set_index("Transcript ID", inplace=True)
            pause_index_dfs[k] = tdf
        pause_index_mdf = pd.concat(list(pause_index_dfs.values()), axis=1)
        # count reads among genes, especially promoter regions
        gene_segmentation_bed = pybedtools.BedTool(transcript_segmentation)
        gk = get_file_hash(transcript_segmentation)[:7]
        promoter_rc_mapping = dict()
        for k in library_beds:
            library_bed = library_beds[k]
            lk = get_file_hash(library_bed)[:7]
            expected_cache = os.path.join(tmp_dir, f"paused_polymerase_efficiency_{gk}_{lk}.csv")
            if not os.path.exists(expected_cache):
                lib_count = gene_segmentation_bed.intersect(library_bed, c=True).to_dataframe(
                    names=("chr", "start", "end", "transcript_id", "transcript_segment", "strand", "count"))
                lib_count = lib_count.loc[lib_count["transcript_segment"] == "Promoter", ("transcript_id", "count")]
                lib_count.set_index("transcript_id", inplace=True)
                lib_count.to_csv(expected_cache)
            else:
                logger.info(f"Promoter read count for {lk} exists, loading...")
                lib_count = pd.read_csv(expected_cache, index_col="transcript_id")
            promoter_rc_mapping[k] = lib_count.to_dict()["count"]

        assay_promoter_signal_df = pause_index_mdf.copy()
        for col in library_beds.keys():
            assay_promoter_signal_df[col] = assay_promoter_signal_df.index.map(promoter_rc_mapping[col])

        super_pre_dfs = []
        for pi_base in pause_indexes.keys():
            if len(target_transcripts) > 0:
                xdf = assay_promoter_signal_df.loc[assay_promoter_signal_df.index.isin(target_transcripts), :].copy()
            else:
                xdf = assay_promoter_signal_df.copy()
            xdf["PI base"] = pi_base
            xdf["PI"] = pause_index_mdf.loc[xdf.index, pi_base]
            super_pre_dfs.append(xdf)
        super_merged_df = pd.concat(super_pre_dfs)
        super_merged_df.to_csv(final_result_file)
    return final_result_file


def main(beds_for_per_rrna, beds_for_coverage, cap_analysis_beds, ro_data_dict, data_save_to, true_enhancers,
         abundant_rna, data_prefix):
    analysis_summaries = {
        "capping_analysis": [],
        "ro_paused_polymerase": [],
        "rrna_zoomin": [],
    }
    analysis_summaries["capping_analysis"].append(
        capping_bias(ds_beds=cap_analysis_beds,
                     enhancer_file=true_enhancers,
                     save_to=data_save_to,
                     global_registry=global_registry,
                     local_registry="EnhancerCappingCoPRO")
    )

    bed_for_ro_pausing = ro_data_dict["library_beds"]
    pause_indexes = ro_pause_files["pause_indexes"]
    transcript_segmentation = ro_pause_files["transcript_segmentation"]
    per_sample_expressed = []
    for expression_profile in ro_pause_files["transcript_expressions"]:
        df = pd.read_csv(expression_profile, sep="\t")
        per_sample_expressed.append(set(df.loc[df.TPM > 5, "transcript_id"].values))

    consistently_expressed = set([t for t in per_sample_expressed[0].intersection(*per_sample_expressed)])
    logger.info(f"Consistently expressed transcripts for capturing efficiency analysis: {len(consistently_expressed)}")
    
    analysis_summaries["ro_paused_polymerase"].append(
        paused_polymerase_efficiency(library_beds=bed_for_ro_pausing,
                                     pause_indexes=pause_indexes,
                                     transcript_segmentation=transcript_segmentation,
                                     save_to=data_save_to,
                                     global_registry=global_registry,
                                     local_registry="ROPauseCountIndex",
                                     tmp_dir=tmp_dir,
                                     target_transcripts=consistently_expressed)
    )

    analysis_summaries["rrna_zoomin"].extend(
        rRNA_spec(beds_for_per_rrna=beds_for_per_rrna,
                  beds_for_coverage=beds_for_coverage,
                  pos_ref_bed=true_enhancers, 
                  abundant_rna=abundant_rna, 
                  save_to=data_save_to, 
                  global_registry=global_registry,
                  local_registries=("BruUVPerrRNA", "BruUVTECov", "BruUVTECovPer")
                  )
    )

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
    cap_analysis_beds = dict()

    ro_pause_files = {
        "library_beds": dict(),
        "pause_indexes": dict(),
        "transcript_segmentation": "",
        "transcript_expressions": [],
    }

    import socket

    server1_home_dir = "/fs/cbsuhy01/storage/ly349/" if socket.gethostname().find(
        "cbsuhy01") == -1 else "/local/storage/ly349/"
    server2_home_dir = "/fs/cbsuhy02/storage/ly349/" if socket.gethostname().find(
        "cbsuhy02") == -1 else "/local/storage/ly349/"
    bioq_dir = os.path.join(server2_home_dir, "BioQueue/workspace")

    from .utils import load_bioq_datasets

    all_reads_beds = load_bioq_datasets("rRNA_zoom_in_bed", bioq_dir, cfg_file=args.config_file)
    ds_beds = load_bioq_datasets("rRNA_zoom_in_ds_bed", bioq_dir, cfg_file=args.config_file)
    cap_analysis_beds = load_bioq_datasets("copro_lib_ec", bioq_dir, cfg_file=args.config_file)
    ro_pause_files["library_beds"] = load_bioq_datasets("RO_pause_libraries", bioq_dir, cfg_file=args.config_file)
    ro_pause_files["pause_indexes"] = load_bioq_datasets("pause_indexes", bioq_dir, cfg_file=args.config_file)
    ro_pause_files["transcript_segmentation"] = cfg.get("references", "gencode_transcript_segmentation")
    ro_pause_files["transcript_expressions"] = (cfg.get("references", "K562_transcript_expression1"),
                                                    cfg.get("references", "K562_transcript_expression2"),)
    main(beds_for_per_rrna=all_reads_beds, beds_for_coverage=ds_beds,
         cap_analysis_beds=cap_analysis_beds, ro_data_dict=ro_pause_files,
         data_save_to=args.data_save_to, fig_save_to=args.fig_save_to, abundant_rna=args.abundant_rna,
         true_enhancers=args.te_bed, data_prefix=args.data_prefix)
