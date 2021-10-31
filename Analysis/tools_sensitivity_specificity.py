#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 8/24/20
import argparse
import logging
import os
import sys
import json
import pickle
import pybedtools
import pyBigWig
import os
from cornerstone.file_parser_mem.gtf import parse_gtf
from copy import deepcopy
import numpy as np
import pandas as pd
from datetime import datetime
from pybedtools import BedTool
from configparser import ConfigParser
from collections import defaultdict, namedtuple

bidirectional_coverage_tracks = namedtuple("BidirectionalCoverageTracks", field_names=("pl", "mn", "digest"))

logging.basicConfig(format='%(name)s - %(asctime)s - %(levelname)s: %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO,
                    handlers=[
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("PINTS - ToolSensitivityResources")

from .utils_roc import *
from .utils import run_command, bin_scores, midpoint_generator, get_file_hash, nfs_mapping


def _filter_candidate_regions_by_exp(input_bed, pl_covs_npy, mn_covs_npy, save_to, window_coverage_threshold=0,
                                     direction=1, cache=True):
    """
    Filter candidate region by requiring signals from input_bam in the same region

    Parameters
    ----------
    input_bed : str
            full path of the input bed file
    pl_covs_npy : dict
        key: chromosome
        value: Path to a npy file for pl coverage
    mn_covs_npy : dict
        key: chromosome
        value: Path to a npy file for mn coverage
    save_to : str
            full path of the output file
    window_coverage_threshold : float
            require regions have signals more than window_coverage_threshold, optional
    direction : int
            direction, 0 for restricted for one strand, 1 for only requiring signal on one strand, 2 for two strands.

    Returns
    -------

    """
    import numpy as np
    assert direction == 0 or direction == 1 or direction == 2, "direction must be 0, 1 or 2"
    pl_covs = dict()
    mn_covs = dict()
    for k, v in pl_covs_npy.items():
        pl_covs[k] = np.load(v)
        mn_covs[k] = np.load(mn_covs_npy[k])

    pls = 0
    mns = 0
    for k in pl_covs:
        pls += np.nansum(pl_covs[k])
        mns += np.nansum(mn_covs[k])

    fh = open(input_bed, "r")
    nfh = open(save_to, "w")
    for line in fh:
        items = line.strip().split("\t")
        if items[0] not in pl_covs_npy:
            logger.warning("Element record %s not in coverage tracks" % line.replace("\n", ""))
            continue
        pl_value = np.sum(pl_covs[items[0]][int(items[1]):int(items[2])])
        mn_value = np.sum(mn_covs[items[0]][int(items[1]):int(items[2])])
        mn_value = mn_value if mn_value >= 0 else mn_value * -1
        ttd = int(pl_value > window_coverage_threshold) + int(mn_value > window_coverage_threshold)
        # restricted single, single or divergent
        if (direction == 0 and ttd == 1) or (direction == 1 and ttd > 0) or (direction == 2 and ttd > 1):
            nfh.write(line)
    fh.close()
    nfh.close()


def _enhancer_set(raw_enhancers, raw_nonenhancers, bam_plc, bam_mnc, output_dir, output_prefix, condition,
                  gencode_annotation, genome_size, size_threshold=250):
    """
    Build the ref set for enhancers

    Parameters
    ----------
    raw_enhancers : str
        Path to the true enhancer set in bed format
    raw_nonenhancers : str
        Path to the non-enhancer set in bed format
    bam_plc : dict
        Dictionary of per base coverage per chromosome (positive strand)
    bam_mnc : dict
        Dictionary of per base coverage per chromosome (negative strand)
    output_dir : str
        Save outputs to this directory
    output_prefix : str
        Prefix for outputs
    condition : str
        Choices: single or bidirectional. Default, bidirectional.
    gencode_annotation : str
        Path to GENCODE annotation file (gtf)
    genome_size : str
        Path to a tab-separated file, which store the size of each chromosome
    size_threshold : int
        Min length a "true" TRE needs to have. Default, 250bp

    Returns
    -------
    tp_enhancers : str
        Path to the compiled ref enhancer file for the assay
    tn_enhancers : str
        Path to the compiled non-enhancer file for the assay
    """
    direction = 1 if condition == "single" else 2
    enhancer_df = pd.read_csv(raw_enhancers, sep="\t", header=None)
    nonenhancer_df = pd.read_csv(raw_nonenhancers, sep="\t", header=None)
    logger.info(f"Size of raw input (enhancer): {enhancer_df.shape[0]}")
    logger.info(f"Size of raw input (non-enhancer): {nonenhancer_df.shape[0]}")
    # filter by size
    size_selected_NE = nonenhancer_df.loc[nonenhancer_df[2] - nonenhancer_df[1] >= size_threshold, :]
    logger.info(f"{size_selected_NE.shape[0]} of {nonenhancer_df.shape[0]} raw non-enhancers passed size filter")

    bed_size_selected_NE = os.path.join(output_dir, output_prefix + "_NE_size.bed")
    size_selected_NE.to_csv(bed_size_selected_NE, sep="\t", index=False, header=False)
    bed_signal_selected_E = os.path.join(output_dir, output_prefix + "_E_size_signal.bed")
    bed_signal_selected_NE = os.path.join(output_dir, output_prefix + "_NE_size_signal.bed")

    # filter by signal
    _filter_candidate_regions_by_exp(bed_size_selected_NE, bam_plc, bam_mnc,
                                     bed_signal_selected_NE, direction=direction)
    _filter_candidate_regions_by_exp(raw_enhancers, bam_plc, bam_mnc,
                                     bed_signal_selected_E, direction=direction)
    ne_signal_selected = pd.read_csv(bed_signal_selected_NE, sep="\t", header=None)
    e_signal_selected = pd.read_csv(bed_signal_selected_E, sep="\t", header=None)
    logger.info("%d of %s raw enhancers passed signal filter" % (e_signal_selected.shape[0], enhancer_df.shape[0]))
    logger.info("%d of %s size-filtered nonenhancers passed signal filter" % (
        ne_signal_selected.shape[0], size_selected_NE.shape[0]))
    ne_signal_selected.sort_values(by=[0, 1], inplace=True)
    e_signal_selected.sort_values(by=[0, 1], inplace=True)
    # merge
    ne_ss_bed = BedTool(ne_signal_selected.to_csv(sep="\t", index=False, header=False), from_string=True)
    e_ss_bed = BedTool(e_signal_selected.to_csv(sep="\t", index=False, header=False), from_string=True)
    ne_ss_merged = ne_ss_bed.merge()
    e_ss_merged = e_ss_bed.merge()
    logger.info(f"After merging overlapped elements, there are {len(e_ss_merged)} enhancers")
    logger.info(f"After merging overlapped elements, there are {len(ne_ss_merged)} nonenhancers")

    tp_enhancers = os.path.join(output_dir, output_prefix + f"_{get_file_hash(raw_enhancers)[:7]}_TP_enhancers.bed")
    e_ss_bed.moveto(tp_enhancers)

    ga = parse_gtf(gencode_annotation)
    ga["start"] -= 1
    ga_transcripts = ga.loc[
        ga["feature"] == "transcript", ("seqname", "start", "end", "transcript_id", "gene_name", "strand")]
    bed_transcripts = BedTool(ga_transcripts.to_csv(sep="\t", index=False, header=False), from_string=True)

    # promoter here is defined as TSS +/- 500bp [-500, TSS, +500]
    bed_promoters = bed_transcripts.flank(g=genome_size, l=500, r=0, s=True).slop(g=genome_size, l=0, r=500, s=True)
    tn_enhancers = os.path.join(output_dir, output_prefix + f"_{get_file_hash(raw_nonenhancers)[:7]}_TN_enhancers.bed")
    bed_tn_enhancers = ne_ss_merged.intersect(b=bed_promoters, v=True)
    logger.info(f"{len(bed_tn_enhancers)} of {len(ne_ss_merged)} nonenhancers passed promoter filer")
    bed_tn_enhancers.moveto(tn_enhancers)
    return tp_enhancers, tn_enhancers


def _filter_by_CAGE(df, pl_cage, mn_cage, inclusion_zone=1000, directioned=0):
    to_ditch = []
    CAGE_pl = pyBigWig.open(pl_cage)
    CAGE_mn = pyBigWig.open(mn_cage)
    for nr, row in df.iterrows():
        if row[5] == "+":
            upstream_offset = row[1] - inclusion_zone
            downstream_offset = row[1] + inclusion_zone
        else:
            upstream_offset = row[1] + inclusion_zone
            downstream_offset = row[1] - inclusion_zone
        pl_stats = 0
        mn_stats = 0
        try:
            pl_stats = CAGE_pl.stats(row[0], int(np.min((int(row[1]), int(upstream_offset)))),
                                     int(np.max((int(row[1]), int(upstream_offset)))))[0]
            mn_stats = CAGE_mn.stats(row[0], int(np.min((int(row[1]), int(downstream_offset)))),
                                     int(np.max((int(row[1]), int(downstream_offset)))))[0]
        except Exception as e:
            logger.error(row, e)
        if pl_stats is None:
            pl_stats = 0
        if mn_stats is None:
            mn_stats = 0
        no_signal = (pl_stats == 0 and mn_stats == 0)
        divergent = (pl_stats > 0 and mn_stats > 0)
        single = (pl_stats > 0 and mn_stats == 0) or (mn_stats > 0 and pl_stats == 0)
        if directioned == 1:  # single + divergent, sa
            if no_signal:
                to_ditch.append(nr)
        elif directioned == 2:  # divergent
            if no_signal or single:
                to_ditch.append(nr)
        else:  # single, ea
            if divergent or no_signal:
                to_ditch.append(nr)
    CAGE_pl.close()
    CAGE_mn.close()
    return to_ditch


def _promoter_set(bam_plc, bam_mnc, output_dir, output_prefix, sampling_sizes, condition, sampling_seed,
                  gencode_annotation, total_RNA, polyA_RNA, CAGE_pl, CAGE_mn, overlapping_gene, fpkm_threshold=1.):
    """
    Build the ref set for promoters

    Parameters
    ----------
    bam_plc : dict
        Dictionary of per base coverage per chromosome (positive strand)
    bam_mnc : dict
        Dictionary of per base coverage per chromosome (negative strand)
    output_dir : str
        Save outputs to this directory
    output_prefix : str
        Prefix for outputs
    sampling_sizes : int
        How many promoters should be drawn from the promoter pool.
    condition : str
        Choices: single or bidirectional. Default, bidirectional.
    sampling_seed : int
        Seed for generating random number and draw samples
    gencode_annotation : str
        Path to GENCODE annotation file (gtf)
    total_RNA : str
        Path to a csv file, which has two columns:
            transcript_id
            fpkm: FPKM values estimated from total RNA-seq
    polyA_RNA : str
        Path to a csv file, which has two columns:
            transcript_id
            fpkm: FPKM values estimated from poly-A+ RNA-seq
    CAGE_pl : str
        Path to a bigwig file, which contains CAGE signals on the forward strand
    CAGE_mn : str
        Path to a bigwig file, which contains CAGE signals on the reverse strand
    overlapping_gene : str
        Path to a txt file, which defines overlapping genes for Homo sapiens. By default, the list from OverGeneDB.
    fpkm_threshold : float
        Only promoters with expression levels (FPKM) higher than this threshold will be considered. Default, 1.

    Returns
    -------
    final_fn : str

    """
    direction = 1 if condition == "single" else 2
    gencode = parse_gtf(gencode_annotation)
    gencode["start"] -= 1
    gencode_transcripts = gencode.loc[
        np.logical_and(gencode["feature"] == "transcript", gencode["gene_type"] == "protein_coding"),
        ("gene_id", "transcript_id", "gene_name", "seqname", "start", "end", "strand")]

    gencode_promoter = deepcopy(gencode_transcripts)
    logger.info(
        "Extracted %d promoters from GENCODE annotations for protein coding transcripts" % gencode_promoter.shape[0])
    tmp = gencode_promoter.apply(
        lambda x: (x["start"] - 500, x["start"] + 500) if x["strand"] == "+" else (x["end"] - 500, x["end"] + 500),
        axis=1, result_type="expand")
    gencode_promoter["ES"] = tmp[0]
    gencode_promoter["EE"] = tmp[1]
    gencode_promoter = gencode_promoter.loc[:,
                       ("seqname", "ES", "EE", "transcript_id", "gene_id", "strand", "gene_name")]
    gencode_promoter.columns = ("seqname", "start", "end", "transcript_id", "gene_id", "strand", "gene_name")
    gencode_promoter.to_csv(os.path.join(output_dir, output_prefix + "_gencode_promoter.bed"), sep="\t", header=False,
                            index=False)

    # Select promoters having specific signals
    _filter_candidate_regions_by_exp(os.path.join(output_dir, output_prefix + "_gencode_promoter.bed"),
                                     bam_plc, bam_mnc, os.path.join(output_dir, output_prefix + "_signal_ok.bed"),
                                     direction=direction)
    signal_ok = pd.read_csv(os.path.join(output_dir, output_prefix + "_signal_ok.bed"), sep="\t", header=None)
    logger.info(f"{signal_ok.shape[0]}/{gencode_promoter.shape[0]} transcripts passed signal filter")

    # Select stable product
    to_ditch_div = _filter_by_CAGE(signal_ok, CAGE_pl, CAGE_mn, directioned=direction)
    cage_filtered_div_df = signal_ok.drop(to_ditch_div)
    logger.info(f"{cage_filtered_div_df.shape[0]}/{signal_ok.shape[0]} transcripts passed CAGE filter")

    total_df = pd.read_csv(total_RNA)
    polyA_df = pd.read_csv(polyA_RNA)
    logger.info(f"Loaded expression for {total_df.shape[0]} transcripts from total RNA-seq")
    logger.info(f"Loaded expression for {polyA_df.shape[0]} transcripts loaded from polyA RNA-seq")
    total_promoters = total_df.set_index("transcript_id").join(cage_filtered_div_df.set_index(3), how="inner")
    polyA_promoters = polyA_df.set_index("transcript_id").join(cage_filtered_div_df.set_index(3), how="inner")
    total_promoters.reset_index(inplace=True)
    polyA_promoters.reset_index(inplace=True)

    total_promoters = total_promoters.loc[total_promoters.fpkm >= fpkm_threshold, (0, 1, 2, "fpkm", "index", 5, 4, 6)]
    polyA_promoters = polyA_promoters.loc[polyA_promoters.fpkm >= fpkm_threshold, (0, 1, 2, "fpkm", "index", 5, 4, 6)]
    logger.info(f"{total_promoters.shape[0]} transcripts passed expression filter (total RNA-seq)")
    logger.info(f"{polyA_promoters.shape[0]} transcripts passed expression filter (polyA RNA-seq)")
    total_fn = os.path.join(output_dir, output_prefix + "_total_expressed.bed")
    polyA_fn = os.path.join(output_dir, output_prefix + "_polyA_expressed.bed")
    total_promoters.to_csv(total_fn, sep="\t", index=False, header=False)
    polyA_promoters.to_csv(polyA_fn, sep="\t", index=False, header=False)

    promoter_merged_div_df = total_promoters.set_index("index").join(polyA_promoters.set_index("index"), how="inner",
                                                                     lsuffix="_total", rsuffix="_polyA")
    logger.info(f"{promoter_merged_div_df.shape[0]} transcripts passed consistency filter (RNA-seq)")

    promoter_merged_div_df.reset_index(inplace=True)
    promoter_merged_div_df = promoter_merged_div_df.loc[:, promoter_merged_div_df.columns[:8]]
    promoter_merged_div_df.columns = ("TranscriptID", "Chromosome", "Start", "End", "FPKM", "Strand", "GeneID", "Gene")
    ov_genes_df = pd.read_csv(overlapping_gene, sep="\t", index_col=False)
    ov_genes = ov_genes_df.Gene1.tolist()
    ov_genes.extend(ov_genes_df.Gene2.tolist())

    to_ditch_ov_genes = []
    current_gene_set = promoter_merged_div_df.Gene.tolist()
    for og in ov_genes:
        if og in current_gene_set:
            to_ditch_ov_genes.append(og)

    promoter_merged_div_df = promoter_merged_div_df.set_index("Gene").drop(to_ditch_ov_genes).reset_index()
    logger.info(f"{promoter_merged_div_df.shape[0]} transcripts passed overlapping gene filter")
    np.random.seed(sampling_seed)  # set seed for reproduction
    # divergent
    assert promoter_merged_div_df.shape[0] >= sampling_sizes, "sampling size is larger than # promoters passed filters"
    gene_choices = np.random.choice(promoter_merged_div_df.index, sampling_sizes, replace=False)
    promoter_merged_div_df = promoter_merged_div_df.loc[gene_choices, ("Chromosome", "Start", "End")]
    logger.info(f"FINAL: {promoter_merged_div_df.shape[0]} promoters were selected")
    final_fn = os.path.join(output_dir, output_prefix + "_promoters.bed")
    promoter_merged_div_df.to_csv(final_fn, sep="\t", header=None, index=False)
    return final_fn


def _generate_contrast_curves_for_ref_sets(assay, corroborative_bws, enhancer_def, promoter_def, negative_def,
                                           chromosome_size, ref_region_extension=1000):
    tp_enhancers = BedTool(midpoint_generator(BedTool(enhancer_def))).slop(l=ref_region_extension,
                                                                           r=ref_region_extension,
                                                                           g=chromosome_size)
    # tp_enhancers = BedTool(enhancer_def).filter(lambda x: len(x) > 100).saveas()
    tp_promoters = BedTool(midpoint_generator(BedTool(promoter_def))).slop(l=ref_region_extension,
                                                                           r=ref_region_extension,
                                                                           g=chromosome_size)
    tn = BedTool(midpoint_generator(BedTool(negative_def))).slop(l=ref_region_extension,
                                                                 r=ref_region_extension,
                                                                 g=chromosome_size)
    result_df = None
    results = []
    for k, v in corroborative_bws.items():
        with pyBigWig.open(nfs_mapping(v)) as bw_obj:
            enhancer_mat, enhancer_m, (enhancer_u, enhancer_l) = bin_scores(regions=tp_enhancers,
                                                                            score_bw=bw_obj)
            promoter_mat, promoter_m, (promoter_u, promoter_l) = bin_scores(regions=tp_promoters,
                                                                            score_bw=bw_obj)
            tn_mat, tn_m, (tn_u, tn_l) = bin_scores(regions=tn, score_bw=bw_obj)
            marker_label = k.split("_")[1]
            results.append(pd.concat([pd.DataFrame({"mean": enhancer_m,
                                                    "upper": enhancer_u,
                                                    "lower": enhancer_l,
                                                    "marker": marker_label,
                                                    "set": "Enhancer",
                                                    "assay": assay}),
                                      pd.DataFrame({"mean": promoter_m,
                                                    "upper": promoter_u,
                                                    "lower": promoter_l,
                                                    "marker": marker_label,
                                                    "set": "Promoter",
                                                    "assay": assay}),
                                      pd.DataFrame({"mean": tn_m,
                                                    "upper": tn_u,
                                                    "lower": tn_l,
                                                    "marker": marker_label,
                                                    "set": "Non-enhancer",
                                                    "assay": assay}),
                                      ]))

    if len(results) > 0:
        result_df = pd.concat(results)
    return result_df


def _build_assay_specific_true_false_sets(assay, random_seed, pl_cov_bw, mn_cov_bw, output_dir, true_tre_set,
                                          false_tre_set, gencode_annotation, genome_size, total_RNA, polyA_RNA,
                                          CAGE_pl, CAGE_mn, overlapping_gene, fpkm_threshold=1.,
                                          condition="bidirectional", track_sig=""):
    """
    Build assay specific true/non-enhancer set

    Parameters
    ----------
    assay : str
        Name of the assay
    random_seed : int
        Seed for generating random number
    pl_cov_bw : str
        Full path to input bigwig file, plus strand
    mn_cov_bw : str
        Full path to input bigwig file, minus strand
    output_dir : str
        Save outputs to this directory
    true_tre_set : str
        Path to the true enhancer set in bed format
    false_tre_set : str
        Path to the non-enhancer set in bed format
    gencode_annotation : str
        Path to GENCODE annotation file (gtf)
    genome_size : str
        Path to a tab-separated file, which store the size of each chromosome
    total_RNA : str
        Path to a csv file, which has two columns:
            transcript_id
            fpkm: FPKM values estimated from total RNA-seq
    polyA_RNA : str
        Path to a csv file, which has two columns:
            transcript_id
            fpkm: FPKM values estimated from poly-A+ RNA-seq
    CAGE_pl : str
        Path to a bigwig file, which contains CAGE signals on the forward strand
    CAGE_mn : str
        Path to a bigwig file, which contains CAGE signals on the reverse strand
    overlapping_gene : str
        Path to a txt file, which defines overlapping genes for Homo sapiens. By default, the list from OverGeneDB.
    fpkm_threshold : float
        Only promoters with expression levels (FPKM) higher than this threshold will be considered. Default, 1.
    condition : str
        Choices: single or bidirectional. Default, bidirectional.

    Returns
    -------

    """
    from pints.io_engine import get_coverage_bw
    pl_covs, mn_covs, rc = get_coverage_bw(bw_pl=pl_cov_bw, bw_mn=mn_cov_bw, chromosome_startswith="chr",
                                           output_dir=output_dir, output_prefix=assay)

    tpe, tne = _enhancer_set(raw_enhancers=true_tre_set, raw_nonenhancers=false_tre_set, bam_plc=pl_covs,
                             bam_mnc=mn_covs, output_dir=output_dir, 
                             output_prefix=f"{assay}_{condition}_{track_sig}",
                             condition=condition, gencode_annotation=gencode_annotation, genome_size=genome_size)

    n_tp = 0
    with open(tpe, "r") as fh:
        for _ in fh:
            n_tp += 1

    tpp = _promoter_set(bam_plc=pl_covs, bam_mnc=mn_covs, output_dir=output_dir, 
                        output_prefix=f"{assay}_{condition}_{track_sig}_{get_file_hash(true_tre_set)[:7]}",
                        sampling_sizes=n_tp, condition=condition, sampling_seed=random_seed,
                        gencode_annotation=gencode_annotation, total_RNA=total_RNA, polyA_RNA=polyA_RNA,
                        CAGE_pl=CAGE_pl, CAGE_mn=CAGE_mn, overlapping_gene=overlapping_gene,
                        fpkm_threshold=fpkm_threshold)
    # merging
    tek = get_file_hash(true_tre_set)[:7]
    nek = get_file_hash(false_tre_set)[:7]
    cmd = f"cat {tpe} {tpp} | sort -k1,1 -k2,2n | "
    cmd += "awk 'BEGIN{OFS=\"\t\";FS=\"\t\"}{print $1,$2,$3 }' > "
    cmd += f"{output_dir}/{assay}_{condition}_{track_sig}_{tek}_positives.bed"
    stdout, stderr, rc = run_command(cmd)
    contrast_df = _generate_contrast_curves_for_ref_sets(assay=assay, corroborative_bws=cfg["corroborative_bws"],
                                                         enhancer_def=tpe,
                                                         promoter_def=tpp,
                                                         negative_def=tne,
                                                         chromosome_size=genome_size)
    contrast_df.to_csv(f"{output_dir}/{assay}_{condition}_{track_sig}_{tek}_{nek}_contrast.csv")
    assert rc == 0, stderr


def generate_rocs(roc_dedicated_jobs, save_to, positive_set, negative_set, gencode_annotation, genome_size,
                  tool_name_mapping, total_RNA, polyA_RNA, CAGE_pl, CAGE_mn, overlapping_gene, local_registry,
                  coverage_tracks, fpkm_threshold=1, element_type="bidirectional", comment=""):
    """
    Generate ROCs

    Parameters
    ----------
    roc_dedicated_jobs : pd.DataFrame
        A dataframe stores info about peak calling jobs for ROC profiling.
    save_to : str
        Path to write outputs
    positive_set : str
        Path to the true enhancer set in bed format
    negative_set : str
        Path to the non-enhancer set in bed format
    gencode_annotation : str
        Path to GENCODE annotation file (gtf)
    genome_size : str
        Path to a tab-separated file, which store the size of each chromosome
    tool_name_mapping : dict
        key: conventional name
        value: official name
    total_RNA : str
        Path to a csv file, which has two columns:
            transcript_id
            fpkm: FPKM values estimated from total RNA-seq
    polyA_RNA : str
        Path to a csv file, which has two columns:
            transcript_id
            fpkm: FPKM values estimated from poly-A+ RNA-seq
    CAGE_pl : str
        Path to a bigwig file, which contains CAGE signals on the forward strand
    CAGE_mn : str
        Path to a bigwig file, which contains CAGE signals on the reverse strand
    overlapping_gene : str
        Path to a txt file, which defines overlapping genes for Homo sapiens. By default, the list from OverGeneDB.
    local_registry : tuple or list
        Local registrations (order) for outputs
    coverage_tracks : dict
        key: name of the assay
        value: instance of `bidirectional_coverage_tracks`
            `pl` : Path to a bigwig file storing signals on the forward strand
            `mn` : Path to a bigwig file storing signals on the reverse strand
            `digest` : hash digest
    fpkm_threshold : float
        Only promoters with expression levels (FPKM) higher than this threshold will be considered. Default, 1.
    element_type : str
        Choices: single or bidirectional. Default, bidirectional.
    comment : str
        abc

    Returns
    -------

    """
    final_result_file1 = os.path.join(save_to,
                                      "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                      local_registry=local_registry[0]))
    final_result_file2 = os.path.join(save_to,
                                      "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                      local_registry=local_registry[1]))
    if all([os.path.exists(f) for f in (final_result_file1, final_result_file2)]):
        logger.info(f"Final output files exist, skip...")
    else:
        logger.info(f"Expected files {final_result_file1} and {final_result_file2} don't exist, running pipeline...")
        pre_merge_sdfs = []
        tf_tre_folder = os.path.join(tmp_dir, "assay_specific_tf_tres")
        if not os.path.exists(tf_tre_folder):
            os.mkdir(tf_tre_folder)
        ref_contrasts = []
        for assay, sdf in roc_dedicated_jobs.groupby(by=["Assay", ]):
            k0 = coverage_tracks[assay].digest  # fingerprint for signal tracks
            k1 = get_file_hash(positive_set)[:7]  # fingerprint for positive loci
            expected_true_tre_file = os.path.join(tf_tre_folder, f"{assay}_{element_type}_{k0}_{k1}_positives.bed")
            k2 = get_file_hash(negative_set)[:7]  # fingerprint for negative loci
            expected_false_tre_file = os.path.join(tf_tre_folder, f"{assay}_{element_type}_{k0}_{k2}_TN_enhancers.bed")
            expected_contrast_file = os.path.join(tf_tre_folder, f"{assay}_{element_type}_{k0}_{k1}_{k2}_contrast.csv")
            if not os.path.exists(expected_true_tre_file) or not os.path.exists(
                    expected_false_tre_file) or not os.path.exists(expected_contrast_file):
                logger.warning(
                    f"Expecting the existences for both {expected_true_tre_file} and {expected_false_tre_file}")
                seed = int(datetime.now().strftime("%m%d%I%M"))
                logger.warning(f"Now generating the two dependent files, using {seed} as seed")

                _build_assay_specific_true_false_sets(assay=assay, random_seed=seed,
                                                      pl_cov_bw=coverage_tracks[assay].pl,
                                                      mn_cov_bw=coverage_tracks[assay].mn, output_dir=tf_tre_folder,
                                                      true_tre_set=positive_set, false_tre_set=negative_set,
                                                      gencode_annotation=gencode_annotation, genome_size=genome_size,
                                                      total_RNA=total_RNA, polyA_RNA=polyA_RNA, CAGE_pl=CAGE_pl,
                                                      CAGE_mn=CAGE_mn, overlapping_gene=overlapping_gene,
                                                      fpkm_threshold=fpkm_threshold, condition=element_type,
                                                      track_sig=k0)
            ref_contrasts.append(pd.read_csv(expected_contrast_file))
            for nr, row in sdf.iterrows():
                incode_name = tool_name_mapping[row['Tool']]
                model = f"Profiler{incode_name}"
                if model in globals():
                    k3 = ""  # fingerprint for peaks
                    if os.path.isfile(row["File"]):
                        k3 = get_file_hash(row["File"])[:7]
                    else:
                        target_path = ""
                        if os.path.exists(row["File"]):
                            target_path = row["File"]
                        else:
                            parent_dir, _ = os.path.split(row["File"])
                            if os.path.exists(parent_dir):
                                target_path = parent_dir
                        if target_path != "":
                            stdout, stderr, rc = run_command(f"ls -alR --full-time {target_path} | sha1sum")
                            if rc == 0:
                                k3 = stdout.split()[0][:7]
                    if comment != "":
                        tmp_result_fn = os.path.join(tmp_dir, f"{row['Tool']}_{assay}_{k0}_{k1}_{k2}_{comment}")
                    else:
                        tmp_result_fn = os.path.join(tmp_dir, f"{row['Tool']}_{assay}_{k0}_{k1}_{k2}")

                    if k3 != "":
                        tmp_result_fn += f"_{k3}"

                    model_obj = globals()[model](output_prefix=tmp_result_fn,
                                                 positive_set=expected_true_tre_file,
                                                 negative_set=expected_false_tre_file,
                                                 cache_dir=tmp_dir, sample_points=101,
                                                 pairing_distance=300)

                    logger.info(f"Evaluating {row['Tool']} on {assay}")
                    try:
                        fn = tmp_result_fn + "_ROC.dat"
                        if not os.path.exists(fn):
                            model_obj.generate_scored_pairs(row['File'])
                            model_obj.roc()
                        else:
                            logger.info(f"Cache for {incode_name} exists ({fn}), loading results from cache...")
                        with open(fn, "rb") as fh:
                            per_tool_per_assay_result = pickle.load(fh)

                        tdf = pd.DataFrame({
                            "Recall": per_tool_per_assay_result["recall"],
                            "Specificity": per_tool_per_assay_result["specificity"],
                            "Precision": per_tool_per_assay_result["precision"],
                            "Fscore": per_tool_per_assay_result["fscore"],
                            "Fscore_alt": per_tool_per_assay_result["fscore_alt"],
                            "Cutoff": per_tool_per_assay_result["cutoff"]
                        })
                        tdf["Assay"] = assay
                        tdf["Tool"] = row["Tool"]
                        pre_merge_sdfs.append(tdf)
                    except Exception as e:
                        logger.error(f"Cannot generate ROC for {row['Tool']} on {assay}")
                        logger.exception(e)
                else:
                    logger.error(f"No model available for {row['Tool']}")

        result_df = pd.concat(pre_merge_sdfs)
        result_df.to_csv(final_result_file1)
        pd.concat(ref_contrasts).to_csv(final_result_file2)
    return final_result_file1, final_result_file2


def main(roc_calls, data_save_to, true_tres, false_tres, cov_tracks, gencode_annotation,
         genome_size, tool_name_mapping, total_RNA, polyA_RNA, CAGE_pl, CAGE_mn, overlapping_gene,
         local_registry="", data_prefix="", **kwargs):
    analysis_summaries = {
        "roc": [],
        "consumption": [],
    }
    roc_calls = kwargs.get("roc_calls")  # peak call dict
    true_tres = kwargs.get("true_tres")  # positive set
    false_tres = kwargs.get("false_tres")  # negative set
    cov_tracks = kwargs.get("cov_tracks")  # coverage tracks for generating assay-specific references
    gencode_annotation = kwargs.get("gencode_annotation")  # gencode annotation in gtf format
    genome_size = kwargs.get("genome_size")  # chromosome size tab file
    total_RNA = kwargs.get("total_RNA")  # path to expression level from total RNA-seq
    polyA_RNA = kwargs.get("polyA_RNA")  # path to expressino level from polyA+ RNA-seq
    CAGE_pl = kwargs.get("CAGE_pl")  # CAGE + bigwig for promoter selection
    CAGE_mn = kwargs.get("CAGE_mn")  # CAGE - bigwig for promoter selection
    overlapping_gene = kwargs.get("overlapping_gene")  # db of overlapping genes

    peak_calls_pre_df = []
    for k, v in roc_calls.items():
        # TSSCall_hg38_K562_GROcap_1_1
        tool, genome_release, cell_line, assay, bio_rep, tech_rep = k.split("_")
        peak_calls_pre_df.append((tool, genome_release, cell_line, assay, bio_rep, tech_rep, v))

    peak_calls_df = pd.DataFrame(peak_calls_pre_df, columns=("Tool", "Genome release", "Cell line", "Assay",
                                                                "Biorep", "Techrep", "File"))
    # ROCs
    hg38_calls = peak_calls_df.loc[
                    (peak_calls_df["Genome release"] == "hg38") & (peak_calls_df["Cell line"] == "K562"), :]

    rocs, ref_contrast = generate_rocs(roc_dedicated_jobs=hg38_calls, save_to=data_save_to, positive_set=true_tres,
                                        negative_set=false_tres, gencode_annotation=gencode_annotation,
                                        genome_size=genome_size, tool_name_mapping=tool_name_mapping, 
                                        total_RNA=total_RNA, polyA_RNA=polyA_RNA,
                                        CAGE_pl=CAGE_pl, CAGE_mn=CAGE_mn, overlapping_gene=overlapping_gene,
                                        local_registry=(f"{local_registry}_{data_prefix}_ROC",
                                                        f"{local_registry}_{data_prefix}_RefContrast"),
                                        coverage_tracks=cov_tracks, fpkm_threshold=1,
                                        element_type="bidirectional", comment="")
    analysis_summaries["roc"].append(rocs)
    analysis_summaries["roc"].append(ref_contrast)

    with open(os.path.join(data_save_to, f"{data_prefix}_summary.json"), "w") as fh:
        json.dump(analysis_summaries, fh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("method", default="all", choices=("generate_data", "plot_figures", "all"),
                        help="Generate data file, plot figures, or both (default: %(default)s)")
    parser.add_argument("--data-save-to", required=True, help="Save output data to")
    parser.add_argument("--data-prefix", required=False, default="tool_roc_consumption", help="Prefix for all outputs")
    parser.add_argument("--global-registry", type=int, help="Global registry", default=4)
    parser.add_argument("--config-file", default="config.txt", help="Configuration file")
    parser.add_argument("--tmp-dir", default="./tmp", help="All temporary files will be written into this folder")
    parser.add_argument("--python-path", required=False, help="Additional Python PATH", default=None)

    args = parser.parse_args()

    assert os.path.exists(args.config_file)

    if not os.path.exists(args.tmp_dir):
        os.mkdir(args.tmp_dir)

    logger.addHandler(logging.FileHandler(os.path.join(os.getcwd(), "ToolSensitivityResources.log")))

    global cfg
    cfg = ConfigParser()
    cfg.optionxform = str
    cfg.read(args.config_file)
    collapsed_assays = cfg.get("assays", "plot_order_simplified").split("|")

    global tmp_dir, global_registry

    tmp_dir = args.tmp_dir
    pybedtools.set_tempdir(tmp_dir)
    global_registry = args.global_registry
    full_assays = cfg.get("assays", "plot_order_full").split("|")
    highlight_assays = cfg.get("assays", "plot_order_simplified").split("|")
    layouts = dict()
    for k, v in cfg["dataset_layouts"].items():
        layouts[k] = v
    unified_color_map = dict()
    plot_order = cfg.get("tools", "plot_order").split("|")
    plot_color = cfg.get("tools", "plot_color").split("|")
    unified_color_map = dict()
    for k, v in enumerate(plot_order):
        unified_color_map[v] = plot_color[k]

    roc_calls = dict()
    ds_roc_calls = dict()
    ds_roc_calls_15m = dict()
    ds_roc_calls_10m = dict()
    cov_tracks = defaultdict(bidirectional_coverage_tracks)
    ds_cov_tracks = defaultdict(bidirectional_coverage_tracks)
    ds_cov_tracks_15m = defaultdict(bidirectional_coverage_tracks)
    ds_cov_tracks_10m = defaultdict(bidirectional_coverage_tracks)
    tool_name_mapping = dict()
    
    import socket
    import hashlib

    server1_home_dir = "/fs/cbsuhy01/storage/ly349/" if socket.gethostname().find(
        "cbsuhy01") == -1 else "/local/storage/ly349/"
    server2_home_dir = "/fs/cbsuhy02/storage/ly349/" if socket.gethostname().find(
        "cbsuhy02") == -1 else "/local/storage/ly349/"
    bioq_dir = os.path.join(server2_home_dir, "BioQueue/workspace")
    from .utils import load_bioq_datasets

    # key: cell_line_assay_tool
    # value: server,user,job_id;file_name
    roc_calls = load_bioq_datasets("peak_call_roc_pr", bioq_dir, cfg_file=args.config_file)
    ds_roc_calls = load_bioq_datasets("peak_calls_ds_19m_roc", bioq_dir, cfg_file=args.config_file)
    ds_roc_calls_15m = load_bioq_datasets("peak_calls_ds_15m_roc", bioq_dir, cfg_file=args.config_file)
    ds_roc_calls_10m = load_bioq_datasets("peak_calls_ds_10m_roc", bioq_dir, cfg_file=args.config_file)
        
    pl_covs = load_bioq_datasets("unique_pl_bigwig_merged", bioq_dir, cfg_file=args.config_file)
    mn_covs = load_bioq_datasets("unique_mn_bigwig_merged", bioq_dir, cfg_file=args.config_file)
    for k, v in pl_covs.items():
        if k in mn_covs:
            cl, assay = k.split("_")
            if cl == "K562":
                cov_tracks[assay] = bidirectional_coverage_tracks._make(
                        (pl_covs[k], mn_covs[k], hashlib.md5(
                            f"{get_file_hash(pl_covs[k])}{get_file_hash(mn_covs[k])}".encode("utf-8")).hexdigest()[:7]
                         )
                    )
    
        dspl_covs = load_bioq_datasets("downsampled_clean_pl_bw", bioq_dir, cfg_file=args.config_file)
        dsmn_covs = load_bioq_datasets("downsampled_clean_mn_bw", bioq_dir, cfg_file=args.config_file)
        for k, v in dspl_covs.items():
            if k in dsmn_covs:
                cl, assay, dsn = k.split("_")
                if cl == "K562":
                    ds_cov_tracks[f"{assay}{dsn}"] = bidirectional_coverage_tracks._make(
                        (dspl_covs[k], dsmn_covs[k], hashlib.md5(
                            f"{get_file_hash(dspl_covs[k])}{get_file_hash(dsmn_covs[k])}".encode("utf-8")).hexdigest()[:7]
                         )
                    )

        dspl_covs_15m = load_bioq_datasets("downsampled_clean_pl_bw_15m", bioq_dir, cfg_file=args.config_file)
        dsmn_covs_15m = load_bioq_datasets("downsampled_clean_mn_bw_15m", bioq_dir, cfg_file=args.config_file)
        for k, v in dspl_covs_15m.items():
            if k in dsmn_covs_15m:
                cl, assay, dsn = k.split("_")
                if cl == "K562":
                    ds_cov_tracks_15m[f"{assay}{dsn}"] = bidirectional_coverage_tracks._make(
                        (dspl_covs[k], dsmn_covs[k], hashlib.md5(
                            f"{get_file_hash(dspl_covs[k])}{get_file_hash(dsmn_covs[k])}".encode("utf-8")).hexdigest()[:7]
                         )
                    )

        dspl_covs_10m = load_bioq_datasets("downsampled_clean_pl_bw_10m", bioq_dir, cfg_file=args.config_file)
        dsmn_covs_10m = load_bioq_datasets("downsampled_clean_mn_bw_10m", bioq_dir, cfg_file=args.config_file)
        for k, v in dspl_covs_10m.items():
            if k in dsmn_covs_10m:
                cl, assay, dsn = k.split("_")
                if cl == "K562":
                    ds_cov_tracks_10m[f"{assay}{dsn}"] = bidirectional_coverage_tracks._make(
                        (dspl_covs[k], dsmn_covs[k],
                         hashlib.md5(f"{get_file_hash(dspl_covs[k])}{get_file_hash(dsmn_covs[k])}".encode(
                             "utf-8")).hexdigest()[:7])
                    )


    formal_names = cfg.get("tools", "formal_names").split("|")
    compatible_names = cfg.get("tools", "compatible_names").split("|")
    for k, v in enumerate(formal_names):
        tool_name_mapping[v] = compatible_names[k]

    if args.python_path is not None and os.path.exists(args.python_path):
        sys.path.append(args.python_path)

    # check peak call files and keys make sure their descriptions match
    # peak_call_roc_pr:
    #   keys: tool_genome_cellline_assay_biorep_techrep
    for k, v in roc_calls.items():
        if not os.path.exists(v):
            logger.error(f"Cannot locate peak calls for {k} {v}")
        tool, genome, cellline, assay, biorep, techrep = k.split("_")
        tool = tool.replace(".", "")
        if not all([v.find(j) != -1 for j in (tool, genome, assay)]):
            logger.error(f"{k} and {v} don't seem to match with each other")
    for k, v in ds_roc_calls.items():
        if not os.path.exists(v):
            logger.error(f"Cannot locate peak calls for {k} {v}")
        tool, genome, cellline, assay, biorep, techrep = k.split("_")
        assay = assay[:-1]
        tool = tool.replace(".", "")
        if not all([v.find(j) != -1 for j in (tool, genome, assay)]):
            logger.error(f"{k} and {v} don't seem to match with each other")
    for k, v in ds_roc_calls_15m.items():
        if not os.path.exists(v):
            logger.error(f"Cannot locate peak calls for {k} {v}")
        tool, genome, cellline, assay, biorep, techrep = k.split("_")
        assay = assay[:-1]
        tool = tool.replace(".", "")
        if not all([v.find(j) != -1 for j in (tool, genome, assay)]):
            logger.error(f"{k} and {v} don't seem to match with each other")
    for k, v in ds_roc_calls_10m.items():
        if not os.path.exists(v):
            logger.error(f"Cannot locate peak calls for {k} {v}")
        tool, genome, cellline, assay, biorep, techrep = k.split("_")
        assay = assay[:-1]
        tool = tool.replace(".", "")
        if not all([v.find(j) != -1 for j in (tool, genome, assay)]):
            logger.error(f"{k} and {v} don't seem to match with each other")

    main(method=args.method, roc_calls=roc_calls, data_save_to=args.data_save_to,
         true_tres=cfg.get("references", "true_enhancers"), false_tres=cfg.get("references", "non_enhancers"),
         cov_tracks=cov_tracks, gencode_annotation=cfg.get("references", "gencode_comprehensive_gtf"),
         genome_size=cfg.get("references", "hg38_chromsize_genome"), tool_name_mapping=tool_name_mapping,
         total_RNA=cfg.get("references", "K562_expression_totalRNA"),
         polyA_RNA=cfg.get("references", "K562_expression_polyA"), CAGE_pl=cfg.get("references", "K562_CAGE_pl"),
         CAGE_mn=cfg.get("references", "K562_CAGE_mn"), overlapping_gene=cfg.get("references", "hs_overlapping_genes"),
         local_registry="SensRes_ROC", data_prefix=args.data_prefix, 
         plot_assays=("GROcap", "CoPRO", "csRNAseq", "NETCAGE", "RAMPAGE", "CAGE", "STRIPEseq"))

    # downsampled libs
    main(method=args.method, roc_calls=ds_roc_calls, data_save_to=args.data_save_to, fig_save_to=args.fig_save_to,
         true_tres=nfs_mapping(cfg.get("references", "true_enhancers")),
         false_tres=nfs_mapping(cfg.get("references", "non_enhancers")),
         cov_tracks=ds_cov_tracks,
         gencode_annotation=nfs_mapping(cfg.get("references", "gencode_comprehensive_gtf")),
         genome_size=nfs_mapping(cfg.get("references", "hg38_chromsize_genome")),
         tool_name_mapping=tool_name_mapping,
         total_RNA=nfs_mapping(cfg.get("references", "K562_expression_totalRNA")),
         polyA_RNA=nfs_mapping(cfg.get("references", "K562_expression_polyA")),
         CAGE_pl=nfs_mapping(cfg.get("references", "K562_CAGE_pl")),
         CAGE_mn=nfs_mapping(cfg.get("references", "K562_CAGE_mn")),
         overlapping_gene=nfs_mapping(cfg.get("references", "hs_overlapping_genes")),
         local_registry="SensRes_DS19M", data_prefix=args.data_prefix+"_DS19M", is_supplement=True,
         plot_assays=("GROcap1", "GROcap2", "GROcap3", "CoPRO1", "CoPRO2", "CoPRO3",
                      "csRNAseq1", "csRNAseq2", "csRNAseq3", "NETCAGE1", "NETCAGE2", "NETCAGE3",
                      "RAMPAGE1", "RAMPAGE2", "RAMPAGE3", "CAGE1", "CAGE2", "CAGE3",
                      "STRIPEseq1", "STRIPEseq2", "STRIPEseq3")
         )
    main(method=args.method, roc_calls=ds_roc_calls_15m, data_save_to=args.data_save_to, fig_save_to=args.fig_save_to,
         true_tres=nfs_mapping(cfg.get("references", "true_enhancers")),
         false_tres=nfs_mapping(cfg.get("references", "non_enhancers")),
         cov_tracks=ds_cov_tracks_15m,
         gencode_annotation=nfs_mapping(cfg.get("references", "gencode_comprehensive_gtf")),
         genome_size=nfs_mapping(cfg.get("references", "hg38_chromsize_genome")),
         tool_name_mapping=tool_name_mapping,
         total_RNA=nfs_mapping(cfg.get("references", "K562_expression_totalRNA")),
         polyA_RNA=nfs_mapping(cfg.get("references", "K562_expression_polyA")),
         CAGE_pl=nfs_mapping(cfg.get("references", "K562_CAGE_pl")),
         CAGE_mn=nfs_mapping(cfg.get("references", "K562_CAGE_mn")),
         overlapping_gene=nfs_mapping(cfg.get("references", "hs_overlapping_genes")),
         local_registry="SensRes_DS15M", data_prefix=args.data_prefix + "_DS15M", is_supplement=True,
         plot_assays=("GROcap1", "GROcap2", "GROcap3", "CoPRO1", "CoPRO2", "CoPRO3",
                      "csRNAseq1", "csRNAseq2", "csRNAseq3", "NETCAGE1", "NETCAGE2", "NETCAGE3",
                      "RAMPAGE1", "RAMPAGE2", "RAMPAGE3", "CAGE1", "CAGE2", "CAGE3",
                      "STRIPEseq1", "STRIPEseq2", "STRIPEseq3")
         )
    main(method=args.method, roc_calls=ds_roc_calls_10m, data_save_to=args.data_save_to, fig_save_to=args.fig_save_to,
         true_tres=nfs_mapping(cfg.get("references", "true_enhancers")),
         false_tres=nfs_mapping(cfg.get("references", "non_enhancers")),
         cov_tracks=ds_cov_tracks_10m,
         gencode_annotation=nfs_mapping(cfg.get("references", "gencode_comprehensive_gtf")),
         genome_size=nfs_mapping(cfg.get("references", "hg38_chromsize_genome")),
         tool_name_mapping=tool_name_mapping,
         total_RNA=nfs_mapping(cfg.get("references", "K562_expression_totalRNA")),
         polyA_RNA=nfs_mapping(cfg.get("references", "K562_expression_polyA")),
         CAGE_pl=nfs_mapping(cfg.get("references", "K562_CAGE_pl")),
         CAGE_mn=nfs_mapping(cfg.get("references", "K562_CAGE_mn")),
         overlapping_gene=nfs_mapping(cfg.get("references", "hs_overlapping_genes")),
         local_registry="SensRes_DS10M", data_prefix=args.data_prefix + "_DS10M", is_supplement=True,
         plot_assays=("GROcap1", "GROcap2", "GROcap3", "CoPRO1", "CoPRO2", "CoPRO3",
                      "csRNAseq1", "csRNAseq2", "csRNAseq3", "NETCAGE1", "NETCAGE2", "NETCAGE3",
                      "RAMPAGE1", "RAMPAGE2", "RAMPAGE3", "CAGE1", "CAGE2", "CAGE3",
                      "STRIPEseq1", "STRIPEseq2", "STRIPEseq3")
         )

