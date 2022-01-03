#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 8/16/20
import argparse
import json
import logging
import os
from configparser import ConfigParser
from cornerstone.file_parser_mem.bed import parse_bed
from collections import defaultdict, namedtuple
import numpy as np
import pandas as pd
import pyBigWig
import pybedtools
from pybedtools import BedTool
from .utils import bin_scores, read_bed, midpoint_generator, run_command, nfs_mapping, get_file_hash

struct_ppc = namedtuple("PeakPairedComparison", field_names=("a_unique", "b_unique",
                                                             "a_shared", "b_shared"))

logging.basicConfig(format='%(name)s - %(asctime)s - %(levelname)s: %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO,
                    handlers=[
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("PINTS - Resolution and robustness")


def tu_benchmark(annotated_TUs, gencode_ptt_bed, save_to, local_registry):
    """
    Evaluate TU annotation tools

    Parameters
    ----------
    annotated_TUs : dict
        key : {tool}_{assay}
        value : path to a bed file
    gencode_ptt_bed : str
        Path to a bed file which contains gencode protein-coding transcripts
        In this study, for transcripts from the same gene, only the ones with the
        highest expression level were kept
    save_to : str
        Path to write output
    local_registry : str
        Local registration (order) for outputs

    Returns
    -------
    final_result_file : str
        csv file contains the final result
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        def per_gene_overlap(x):
            if x["a_overlap"] == 0:
                return np.nan
            else:
                return x["a_overlap"] / (x["g_end"] - x["g_start"])

        def per_tu_overlap(x):
            try:
                return x["a_overlap"] / (x["a_end"] - x["a_start"])
            except:
                return np.nan

        def per_gene_jaccard(x):
            if x["a_overlap"] == 0:
                return np.nan
            else:
                min_coord = min(x["g_start"], x["g_end"], x["a_start"], x["a_end"])
                max_coord = max(x["g_start"], x["g_end"], x["a_start"], x["a_end"])
                return x["a_overlap"] / (max_coord - min_coord)

        gencode_ptt_bed = BedTool(gencode_ptt_bed)
        ga_tu_sdfs = []
        for tool_assay in annotated_TUs:
            tool, assay = tool_assay.split("_")
            gencode_tool_df = gencode_ptt_bed.intersect(annotated_TUs[f"{tool}_{assay}"], wao=True, s=True).to_dataframe(
                names=(
                    "g_chrom", "g_start", "g_end", "g_gene", "g_exp", "g_strand", "a_chrom", "a_start", "a_end", "a_id",
                    "a_score", "a_strand", "a_overlap")
                )
            gencode_tool_df["overlapG"] = gencode_tool_df.apply(per_gene_overlap, axis=1)
            gencode_tool_df["overlapT"] = gencode_tool_df.apply(per_tu_overlap, axis=1)
            gencode_tool_df["Jaccard"] = gencode_tool_df.apply(per_gene_jaccard, axis=1)
            ga_tu_sdfs.append(pd.DataFrame(
                {"Jaccard": gencode_tool_df.Jaccard.values,
                 "overlapG": gencode_tool_df.overlapG.values,
                 "overlapT": gencode_tool_df.overlapT.values,
                    "Assay": assay,
                 "Tool": tool, }
            ))
        ga_tu_df = pd.concat(ga_tu_sdfs)
        ga_tu_mdf = ga_tu_df.melt(value_vars=("Jaccard", "overlapG", "overlapT"), id_vars=("Assay", "Tool"),
                                  value_name="Consistency")
        ga_tu_mdf.to_csv(final_result_file)
    return final_result_file


def evaluate_systematic_bias(peak_call_jobs_df, score_bws, ref_regions, chromosome_size,
                             save_to, local_registry, ref_region_extension=500):
    """
    Evaluate systematic bias when using different peak calling tools for a set of assays

    Parameters
    ----------
    peak_call_jobs_df : pd.DataFrame
        DataFrame containing Tool, Genome release, Cell line, Assay, Biorep, Techrep, File
    score_bws : dict
        Dictionary of bigwig score files. Keys are in this format: cellLine_scoreName
        Values are full paths to the bigwig files
    ref_regions : BedTool
        reference regions as BedTool
    chromosome_size : str
        Path to a tab-separated file containing chromosome sizes (chr    size)
    save_to : str
        Path to write output
    local_registry : str
        Local registration (order) for outputs
    ref_region_extension : int
        Bps to extend from the center of regions. By default, 500bp

    Returns
    -------
    final_result_file : str
        csv file contains the final result
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        ref_region_bed = BedTool(ref_regions)
        mids_ref = BedTool(midpoint_generator(ref_region_bed))
        mids_ref_extended = mids_ref.slop(l=ref_region_extension, r=ref_region_extension, g=chromosome_size)
        results_pre_df = {
            "tool": [],
            "marker": [],
            "loc": [],
            "mean": [],
            "assay": []
        }
        sdf = peak_call_jobs_df.loc[
              (peak_call_jobs_df["Genome release"] == "hg38") & (peak_call_jobs_df["Cell line"] == "K562"), :]
        for corroborative_name, corroborative_path in score_bws.items():
            with pyBigWig.open(corroborative_path) as ref_bw:
                # score_mat, means, stds
                sm_tpe, m_tpe, s_tpe = bin_scores(mids_ref_extended, ref_bw, bins=100)
                results_pre_df["tool"].extend(["Ref"] * m_tpe.shape[0])
                results_pre_df["marker"].extend([corroborative_name] * m_tpe.shape[0])
                results_pre_df["loc"].extend(np.arange(m_tpe.shape[0]).tolist())
                results_pre_df["mean"].extend(m_tpe.tolist())
                results_pre_df["assay"].extend(["Generic"] * m_tpe.shape[0])

        for assay, peak_calls_sdf in sdf.groupby(by=["Assay", ]):
            for _, row in peak_calls_sdf.iterrows():
                try:
                    # tmp = pd.read_csv(row["File"], sep="\t", header=None, comment="#")
                    tmp = read_bed(row["File"])
                    tmp = tmp.loc[tmp[0].str.startswith("chr"), :]
                    if tmp.shape[1] > 6:
                        tmp = tmp.loc[:, tmp.columns[:6]]
                except Exception as e:
                    logger.error(row["File"])
                    logger.exception(e)
                    continue
                if tmp.shape[0] == 0:
                    logger.warning(f"Empty peak calls: {assay} {row['Tool']} {row['File']}")
                    continue
                tmp = tmp.loc[tmp[2] - tmp[1] < 2000, :]
                tmp[1] = tmp[1].astype(int)
                tmp[2] = tmp[2].astype(int)
                tmp_bed = BedTool.from_dataframe(tmp)
                x = tmp_bed.intersect(ref_region_bed, u=True).saveas()
                if x.__len__() == 0:
                    logger.warning(f"No peak from {row['File']} locates in reference regions.")
                    continue

                mids = BedTool(midpoint_generator(x)).saveas()
                try:
                    mids_extended = mids.slop(l=ref_region_extension, r=ref_region_extension, g=chromosome_size)
                except Exception as e:
                    logger.error(row["File"])
                    logger.exception(e)
                    continue

                for corroborative_name, corroborative_path in score_bws.items():
                    with pyBigWig.open(corroborative_path) as ref_bw:
                        sm, m, s = bin_scores(mids_extended, ref_bw, bins=100)
                        results_pre_df["tool"].extend([row["Tool"]] * m.shape[0])
                        results_pre_df["marker"].extend([corroborative_name] * m.shape[0])
                        results_pre_df["loc"].extend(np.arange(m.shape[0]).tolist())
                        results_pre_df["mean"].extend(m.tolist())
                        results_pre_df["assay"].extend([assay] * m_tpe.shape[0])
        pd.DataFrame(results_pre_df).to_csv(final_result_file)
    return final_result_file


def retrieve_tre_size_dist(peak_call_jobs_df, save_to, local_registry):
    """
    Retrieve all sizes of TREs identified by different tools among different assays
    Parameters
    ----------
    peak_call_jobs_df : pd.DataFrame
        DataFrame containing Tool, Genome release, Cell line, Assay, Biorep, Techrep, File
    save_to : str
        Path to write output
    local_registry : str
        Local registration (order) for outputs

    Returns
    -------
    final_result_file : str
        csv file contains the final result
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        peak_calls_df = peak_call_jobs_df.loc[
                        (peak_call_jobs_df["Genome release"] == "hg38") & (peak_call_jobs_df["Cell line"] == "K562"), :]
        pre_results_df = {
            "Assay": [],
            "Tool": [],
            "Size": []
        }
        for assay, peak_calls_sdf in peak_calls_df.groupby(by=["Assay", ]):
            for _, row in peak_calls_sdf.iterrows():
                try:
                    tmp = read_bed(row["File"])
                except Exception as e:
                    logger.info(row["File"])
                    logger.exception(e)
                    continue
                tmp = tmp.loc[tmp[2] - tmp[1] < 2000, :]
                l = tmp[2] - tmp[1]
                pre_results_df["Assay"].extend([assay] * l.shape[0])
                pre_results_df["Tool"].extend([row["Tool"]] * l.shape[0])
                pre_results_df["Size"].extend(l.tolist())
        pd.DataFrame(pre_results_df).to_csv(final_result_file)
    return final_result_file


def compile_distributions_of_gc(gc_files, save_to, local_registry):
    """
    Compile distributions of genomic coverage

    Parameters
    ----------
    gc_files : dict
        Key: "cell line"_"assay"
        Value: str to coverage tables (contains %d)
    save_to : str
        Save file to this folder
    local_registry : str
        Local registry

    Returns
    -------
    final_result_file : str
        csv file contains the final result
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        pre = []
        for hg in (19, 38):
            for k, v in gc_files.items():
                cell_line, assay = k.split("_")
                hg_file = v % hg
                tdf = pd.read_csv(hg_file, sep="\t", names=("Contig", "Read counts", "RC loci", "Total loci", "Freq"))
                tdf["Assay"] = assay
                tdf["Release"] = hg
                pre.append(tdf)

        pd.concat(pre).to_csv(final_result_file)
    return final_result_file


def _calc_reproducibility_ub(raw_hg19, raw_hg38, tools, assay, cell_line, hg19_promoters, hg38_promoters,
                             hg19_to_38_chain_file="~/refs/liftover_chain/hg19ToHg38.over.chain.gz",
                             hg38_to_19_chain_file="~/refs/liftover_chain/hg38ToHg19.over.chain.gz"):
    """
    Calculate upper bound of reproducibility

    Parameters
    ----------
    raw_hg19 : list or tuple
        list of peak calls from hg19
    raw_hg38 : list or tuple
        list of peak calls from hg38
    tools : list or tuple
    assay : list or tuple
    cell_line : list or tuple
    hg19_to_38_chain_file
    hg38_to_19_chain_file

    Returns
    -------

    """
    result = []

    def _liftover(to_lifts, chain_file, mark=19):
        lifted_files = []
        cached_files = []
        for raw_file in to_lifts:
            _, fname = os.path.split(raw_file)
            fn, ext = os.path.splitext(fname)
            try:
                raw_df = parse_bed(raw_file)
                reformatted_raw = os.path.join(tmp_dir, fn + ".bed")
                raw_df.loc[:, ("chrom", "chromStart", "chromEnd")].to_csv(reformatted_raw, sep="\t",
                                                                          header=False, index=False)

                lifted_file = os.path.join(tmp_dir, "%s_lifted_to_hg%d.bed" % (fn, mark))
                log_file = os.path.join(tmp_dir, "%s_lifted_to_hg%d.log" % (fn, mark))
                stdout, stderr, rc = run_command(
                    "liftOver %s %s %s %s" % (reformatted_raw, chain_file, lifted_file, log_file))
                if rc != 0:
                    logger.warning(stderr + " %s" % raw_file)
                    cached_files.append(None)
                    lifted_files.append(None)
                    continue
                cached_files.append(reformatted_raw)
                lifted_files.append(lifted_file)
            except pd.errors.ParserError:
                cached_files.append(None)
                lifted_files.append(None)
                logger.exception("Failed to parse file %s" % raw_file)

        return lifted_files, cached_files

    lifted_hg19, cached_raw_hg19 = _liftover(raw_hg19, chain_file=hg19_to_38_chain_file, mark=38)
    lifted_hg38, cached_raw_hg38 = _liftover(raw_hg38, chain_file=hg38_to_19_chain_file, mark=19)

    hg38_promoter_bed = BedTool(hg38_promoters)
    hg19_promoter_bed = BedTool(hg19_promoters)

    for k, raw in enumerate(cached_raw_hg19):
        if raw is not None and lifted_hg38[k] is not None:
            a = BedTool(raw).sort().intersect(hg38_promoter_bed, v=True)
            b = BedTool(lifted_hg38[k]).sort().intersect(hg38_promoter_bed, v=True)
            jaccard = a.jaccard(b)["jaccard"]
            result.append((assay[k], cell_line[k], raw, lifted_hg38[k], jaccard, tools[k]))

    for k, raw in enumerate(cached_raw_hg38):
        if raw is not None and lifted_hg19[k] is not None:
            a = BedTool(raw).sort().intersect(hg19_promoter_bed, v=True)
            b = BedTool(lifted_hg19[k]).sort().intersect(hg19_promoter_bed, v=True)
            jaccard = a.jaccard(b)["jaccard"]
            result.append((assay[k], cell_line[k], raw, lifted_hg19[k], jaccard, tools[k]))
    df = pd.DataFrame(result, columns=["Assay", "Cell line", "Raw", "Lifted", "Jaccard", "Software"])

    return df


def evaluate_upper_bound_robustness(peak_calls, hg19_to_38_chain_file, hg38_to_19_chain_file,
                                    hg19_promoters, hg38_promoters, save_to, local_registry):
    """
    Evaluate upper bounds of robustness

    Parameters
    ----------
    peak_calls : dict
        Dict of peak calls, keys must be in this format: tool_genomeRelease_cellLine_assay_biologicalRep_techRep.
        Values must be a full path to the peak calls
    hg19_to_38_chain_file : str
        Liftover chain from hg19 to hg38
    hg38_to_19_chain_file : str
        Liftover chain from hg38 to hg19
    hg19_promoters : str
        Promoter bed file for hg19
    hg38_promoters : str
        Promoter bed file for hg38
    save_to : str
        Path to write output
    local_registry : str
        Local registration (order) for outputs

    Returns
    -------
    final_result_file : str
            csv file contains the final result
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        pk_jobs_pre = []
        for k, v in peak_calls.items():
            t = k.split("_")
            if os.path.exists(v):
                t.append(v)
                pk_jobs_pre.append(t)
            else:
                logger.error(f"Cannot find peak calls: {v}")
        pk_jobs = pd.DataFrame(pk_jobs_pre, columns=("tool", "genome_release",
                                                     "cell_line", "assay",
                                                     "biological_rep", "tech_rep",
                                                     "file"))
        hg19_jobs = []
        hg38_jobs = []
        tools = []
        assays = []
        cell_lines = []
        for m, sdf in pk_jobs.groupby(by=["tool", "cell_line", "assay", "biological_rep", "tech_rep"]):
            assay = set(sdf["assay"].values).pop()
            cell_line = set(sdf["cell_line"].values).pop()
            grs = set(sdf["genome_release"])
            if len(grs) != 2:
                logger.warning(f"Skipping {m}, because it doesn't have peak calls on both hg19 and hg38.")
                continue
            hg19_file = sdf.loc[sdf["genome_release"] == "hg19", "file"].values[0]
            hg38_file = sdf.loc[sdf["genome_release"] == "hg38", "file"].values[0]
            if not os.path.exists(hg19_file) or not os.path.exists(hg38_file):
                logger.warning(f"Skipping {m}, because file doesn't exist.")
                continue
            hg19_jobs.append(hg19_file)
            hg38_jobs.append(hg38_file)
            assays.append(assay)
            cell_lines.append(cell_line)
            tools.append(m[0])
        reproducibility = _calc_reproducibility_ub(raw_hg19=hg19_jobs, raw_hg38=hg38_jobs, tools=tools,
                                                   assay=assays, cell_line=cell_lines,
                                                   hg19_promoters=hg19_promoters,
                                                   hg38_promoters=hg38_promoters,
                                                   hg19_to_38_chain_file=hg19_to_38_chain_file,
                                                   hg38_to_19_chain_file=hg38_to_19_chain_file)
        reproducibility.to_csv(final_result_file)
    return final_result_file


def evaluate_mapping_similarity_same_lib(u, v):
    from collections import Counter
    u_counter = Counter(u)
    v_counter = Counter(v)

    u_set = set(u)
    v_set = set(v)
    common = u_set.intersection(v_set)

    n_u = len(u)
    n_v = len(v)

    shared = 0

    for k in common:
        shared += min(u_counter[k], v_counter[k])
    return np.mean((shared / n_u, shared / n_v)), (shared, n_u, n_v)


def evaluate_real_case_robustness(peak_calls_tech, peak_calls_bio, correlations, save_to, local_registry):
    """

    Parameters
    ----------
    peak_calls_tech : dict
        Dict of peak calls from technical replicates
        Keys must be in this format: cellLine_assay_condition_rep_tool.
        Values must be a full path to the peak calls
    peak_calls_bio : dict
        Dict of peak calls from biological replicates
        Keys must be in this format: cellLine_assay_rep_tool.
        Values must be a full path to the peak calls
    correlations : str

    save_to : str
        Path to write output
    local_registry : str
        Local registration for outputs

    Returns
    -------
    final_result_file : str
            csv file contains the final result
    """
    cor_df = pd.read_csv(correlations, sep="\t", skiprows=1, index_col=0, quotechar="'")

    def _atom(peak_calls):
        rep_datasource_pre = []
        for k, v in peak_calls.items():
            cell_line, assay, condition, rep, tool = k.split("_")
            rep_datasource_pre.append((cell_line, assay, condition, rep, tool,
                                       cell_line + "_" + assay + "_" + condition + "_" + rep, v))
        rep_datasource = pd.DataFrame(rep_datasource_pre,
                                      columns=("Cell line", "Assay", "Condition", "Rep", "Tool", "PrimaryKey", "Peak"))
        grouped = rep_datasource.groupby(by=["Assay", "Condition", "Tool"])
        if grouped.ngroups == rep_datasource.shape[0]:
            grouped = rep_datasource.groupby(by=["Assay", "Rep", "Tool"])
        rep_result_pre = []
        for m, sdf in grouped:
            # avged bam correlation
            from itertools import combinations
            for file_a, file_b in combinations(sdf.Peak.values, 2):
                # if sdf.shape[0] == 2:
                try:
                    cor = np.mean([cor_df.loc[a, b] for a, b in combinations(sdf.PrimaryKey.values, 2)])
                    a = read_bed(file_a)
                    a.sort_values(by=[0, 1], inplace=True)
                    a = a.loc[a[1] >= 0, :]
                    a[1] = a[1].astype(int)
                    a[2] = a[2].astype(int)
                    b = read_bed(file_b)
                    b.sort_values(by=[0, 1], inplace=True)
                    b = b.loc[b[1] >= 0, :]
                    b[1] = b[1].astype(int)
                    b[2] = b[2].astype(int)
                    c = BedTool.from_dataframe(a).jaccard(BedTool.from_dataframe(b))
                    j = c["jaccard"]
                    # print(sdf["Peak"].values[0], sdf["Peak"].values[1], j)
                    # rep_result_pre.append(
                    #     (sdf["Tool"].values[0], sdf["Assay"].values[0], sdf["Condition"].values[0], j, cor))
                except Exception as e:
                    logger.error(file_a + ";" + file_b)
                    logger.exception(e)
                    j = np.nan
                rep_result_pre.append((m[2], m[0], m[1], j, cor))
        robustness = pd.DataFrame(rep_result_pre, columns=["Software", "Assay", "Condition", "Jaccard", "Correlation"])
        robustness_agged = robustness.groupby(by=["Software", ]).agg({"Jaccard": ["mean", "std"],
                                                                      "Correlation": ["mean", ]})
        robustness_agged.reset_index(inplace=True)
        robustness_agged.columns = [" ".join(col).strip() for col in robustness_agged.columns.values]
        robustness_agged_for_pivot = robustness.pivot_table(index="Assay", columns="Software", values="Jaccard",
                                                            aggfunc="mean")
        robustness_agged_for_pivot["Correlation"] = 0
        for nr, _ in robustness_agged_for_pivot.iterrows():
            cor = robustness.loc[robustness.Assay == nr, "Correlation"].mean()
            robustness_agged_for_pivot.loc[nr, "Correlation"] = cor
        return robustness_agged_for_pivot

    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        try:
            tech_results = _atom(peak_calls_tech)
            bio_results = _atom(peak_calls_bio)

            tech_bio_agged = pd.concat([tech_results, bio_results], sort=False)
            tuples = list(zip(tech_bio_agged.index,
                              ["Technical"] * tech_results.shape[0] + ["Biological"] * bio_results.shape[0]))
            index = pd.MultiIndex.from_tuples(tuples, names=["Replicate type", "Assay"])
            tech_bio_agged_dual_index = tech_bio_agged.copy()
            tech_bio_agged_dual_index.index = index
            tech_bio_agged_dual_index.to_csv(final_result_file)
        except:
            # logger.error(peak_calls_tech)
            pass
    return final_result_file


def _uniq_ele_histone_profile_atom(tool_name, shared, pints_unique, other_unique, score_bws, chromosome_size,
                                   region_extension):
    """
    Atom function for profiling histone modification signal

    Parameters
    ----------
    tool_name : str
        ss
    shared : pybedtools.BedTool
        Shared elements
    pints_unique : pybedtools.BedTool
        PINTS unique elements
    other_unique : pybedtools.BedTool
        Unique elements reported by other tools
    score_bws : dict
        key: Name of the corroborative mark (like H3K4me1, H3K27ac)
        value: Path to the bigwig file
    chromosome_size : str
        Path to a tab-separated file which defines the size of each chromosome
    region_extension : int
        The final elements will be [original_mid-ref_region_extension, original_mid+ref_region_extension]

    Returns
    -------
    results_df : pd.DataFrame
        A dataframe with the following columns:
            tool
            label
            marker
            loc
            mean
            mean_u
            mean_l
            assay
    """
    results_pre_df = {
        "tool": [],
        "label": [],
        "marker": [],
        "loc": [],
        "mean": [],
        "mean_u": [],
        "mean_l": [],
        "assay": []
    }
    pints_s_mids = BedTool(midpoint_generator(shared)).saveas().slop(l=region_extension,
                                                                     r=region_extension,
                                                                     g=chromosome_size).saveas()
    pints_u_mids = BedTool(midpoint_generator(pints_unique)).saveas().slop(l=region_extension,
                                                                           r=region_extension,
                                                                           g=chromosome_size).saveas()
    other_u_mids = BedTool(midpoint_generator(other_unique)).saveas().slop(l=region_extension,
                                                                           r=region_extension,
                                                                           g=chromosome_size).saveas()

    for corroborative_name, corroborative_path in score_bws.items():
        with pyBigWig.open(corroborative_path) as ref_bw:
            for bed_obj, label in zip((pints_s_mids, pints_u_mids, other_u_mids), ("Shared",
                                                                                   "PINTS",
                                                                                   tool_name)):
                sm, m, s = bin_scores(bed_obj, ref_bw, bins=100)
                results_pre_df["tool"].extend([tool_name] * m.shape[0])
                results_pre_df["label"].extend([label] * m.shape[0])
                results_pre_df["marker"].extend([corroborative_name] * m.shape[0])
                results_pre_df["loc"].extend(np.arange(m.shape[0]).tolist())
                results_pre_df["mean"].extend(m.tolist())
                results_pre_df["mean_u"].extend(s[1].tolist())
                results_pre_df["mean_l"].extend(s[0].tolist())
                results_pre_df["assay"].extend([assay] * m.shape[0])
    return pd.DataFrame(results_pre_df)


def _uniq_ele_histone_profile(dict_ppc, score_bws, file_save_to, ref_region_extension=300, chromosome_size="hg38"):
    """
    Profile distribution of histone modification signal for unique elements

    Parameters
    ----------
    dict_ppc : dict
        key: name of peak calls used as control
        value: struct_ppc (PeakPairedComparison)
    score_bws : dict
        key: Name of the corroborative mark (like H3K4me1, H3K27ac)
        value: Path to the bigwig file
    file_save_to : str
        Write result to this file
    ref_region_extension : int
        The final elements will be [original_mid-ref_region_extension, original_mid+ref_region_extension]
    chromosome_size : str
        Path to a tab-separated file which defines the size of each chromosome

    Returns
    -------
    file_save_to : str
        Path to the result file
    """
    from multiprocessing import Pool

    jobs = []
    for tool, ppc_obj in dict_ppc.items():
        # shared, pints_unique, other_unique, score_bws, chromosome_size, region_extension
        jobs.append((tool, ppc_obj.a_shared, ppc_obj.a_unique, ppc_obj.b_unique,
                     score_bws, chromosome_size, ref_region_extension))

    with Pool(16) as pool:
        sub_dfs = pool.starmap(_uniq_ele_histone_profile_atom, jobs)
        pd.DataFrame(pd.concat(sub_dfs)).to_csv(file_save_to)
    return file_save_to


def _uniq_ele_motif_profile(dict_ppc, chromosome_size, genome_fa, output_prefix,
                            extension=200, remove_longer_than=1000, n_motif=746, tmp_dir=".",
                            motif_meme="data/JASPAR_CORE_2020_vertebrates.meme",
                            ame_exec="export PATH=/programs/meme/bin:$PATH && ame"):
    """
    Profile motif enrichment for unique elements

    Parameters
    ----------
    dict_ppc : dict
        key: name of peak calls used as control
        value: struct_ppc (PeakPairedComparison)
    chromosome_size : str
        Path to a tab-separated file containing chromosome sizes (chr    size)
    genome_fa : str
        Path to a FASTA file for sequence extraction
    output_prefix : str
        Output prefix
    extension : int
        Basepairs to be extended from the original elements/peaks. This is necessary, because several tools
        generate very narrow elements that are even shorter than the k-mer in `motif_meme`. By default, 200.
    remove_longer_than : int
        Only peaks shorter than this cutoff will be kept for downstream analysis
    n_motif : int
        Number of motifs included in `motif_meme`
    tmp_dir : str
        Path to write tmp files
    motif_meme : str
        Path to a motif database in MEME format
    ame_exec : str
        Command for running AME

    Returns
    -------
    (fn_lor, fn_loru, fn_lorl, fn_lorp) : (str, str, str, str)
        Output files for log odds ratio
        upper bounds of LORs
        lower bounds of LORs
        p-values for each comparison from z test (p-q=0?)
    """
    from multiprocessing import Pool
    ame_lor_fn = output_prefix + "_AME_LOR.csv"
    ame_loru_fn = output_prefix + "_AME_LORu.csv"
    ame_lorl_fn = output_prefix + "_AME_LORl.csv"
    ame_lorp_fn = output_prefix + "_AME_LORp.csv"
    expected_files = (ame_lor_fn, ame_loru_fn, ame_lorl_fn, ame_lorp_fn)

    if all([os.path.exists(f) for f in expected_files]):
        logger.info(f"Final output files {expected_files} exists, skip...")
    else:
        jobs = []
        ame_result_dirs = dict()
        # first, use AME to scan unique elements (both PINTS and the other tool)
        for tool, obj_ppc in dict_ppc.items():
            pints_uniq = obj_ppc.a_unique
            other_uniq = obj_ppc.b_unique
            if extension > 0:
                pints_uniq = obj_ppc.a_unique.slop(b=extension, g=chromosome_size).saveas()
                other_uniq = obj_ppc.b_unique.slop(b=extension, g=chromosome_size).saveas()
            pints_uniq = pints_uniq.filter(lambda x: len(x) < remove_longer_than).saveas()
            other_uniq = other_uniq.filter(lambda x: len(x) < remove_longer_than).saveas()
            pusf = pints_uniq.sequence(fi=genome_fa).seqfn
            ousf = other_uniq.sequence(fi=genome_fa).seqfn
            ame_save_to = os.path.join(tmp_dir, f"AME_PINTS_{tool}")
            jobs.append(
                f"{ame_exec} --evalue-report-threshold {n_motif} -o {ame_save_to} --control {ousf} {pusf} {motif_meme}"
            )
            ame_result_dirs[tool] = ame_save_to

        with Pool(16) as pool:
            ame_exec_results = pool.map(run_command, jobs)
            for stdout, stderr, rc in ame_exec_results:
                if rc != 0:
                    logger.exception(stderr)

        # then aggregate AME results
        from statsmodels.stats.contingency_tables import Table2x2
        AME_LOR_subdfs = []
        AME_LORu_subdfs = []
        AME_LORl_subdfs = []
        AME_logpval_subdfs = []
        for tool, ame_result_dir in ame_result_dirs.items():
            expected_ame_output = os.path.join(ame_result_dir, "ame.tsv")
            assert os.path.exists(expected_ame_output)
            lor_vec = []
            lor_u_vec = []
            lor_l_vec = []
            lor_pval_vec = []

            known_results = pd.read_csv(expected_ame_output, sep="\t", comment="#")
            for _, row in known_results.iterrows():
                ctable = Table2x2(np.asarray(((row["TP"], row["pos"] - row["TP"]),
                                              (row["FP"], row["neg"] - row["FP"]))))
                ci_l, ci_u = ctable.log_oddsratio_confint(alpha=0.1)
                lor_vec.append(ctable.log_oddsratio)
                lor_u_vec.append(ci_u)  # upper CI
                lor_l_vec.append(ci_l)  # lower CI
                lor_pval_vec.append(ctable.log_oddsratio_pvalue())
            known_results["log_odds_ratio"] = lor_vec
            known_results["log_odds_ratio_upper"] = lor_u_vec
            known_results["log_odds_ratio_lower"] = lor_l_vec
            known_results["log_pval"] = lor_pval_vec
            AME_LOR_subdfs.append(
                known_results.loc[:, ("motif_alt_ID", "log_odds_ratio")].set_index("motif_alt_ID").rename(
                    {"log_odds_ratio": tool}, axis=1)
            )
            AME_LORu_subdfs.append(
                known_results.loc[:, ("motif_alt_ID", "log_odds_ratio_upper")].set_index("motif_alt_ID").rename(
                    {"log_odds_ratio_upper": tool}, axis=1)
            )
            AME_LORl_subdfs.append(
                known_results.loc[:, ("motif_alt_ID", "log_odds_ratio_lower")].set_index("motif_alt_ID").rename(
                    {"log_odds_ratio_lower": tool}, axis=1)
            )
            AME_logpval_subdfs.append(
                known_results.loc[:, ("motif_alt_ID", "log_pval")].set_index("motif_alt_ID").rename(
                    {"log_pval": tool}, axis=1)
            )

        pd.concat(AME_LOR_subdfs, axis=1).replace([np.inf, -np.inf], 0).to_csv(ame_lor_fn)
        pd.concat(AME_LORu_subdfs, axis=1).replace([np.inf, -np.inf], 0).to_csv(ame_loru_fn)
        pd.concat(AME_LORl_subdfs, axis=1).replace([np.inf, -np.inf], 0).to_csv(ame_lorl_fn)
        pd.concat(AME_logpval_subdfs, axis=1).replace([np.inf, -np.inf], 0).to_csv(ame_lorp_fn)
    return expected_files


def profile_unique_set(peak_calls_df, genome_fa, dict_signal_track_bws, promoter, motif_meme, chromosome_size, save_to,
                       local_registry, tmp_dir="."):
    """
    Profile unique peaks

    Parameters
    ----------
    peak_calls_df : pd.DataFrame
        DataFrame containing Tool, Genome release, Cell line, Assay, Biorep, Techrep, File
    genome_fa : str
        Path to a FASTA file for sequence extraction
    dict_signal_track_bws : dict
        key: Name of the corroborative mark (like H3K4me1, H3K27ac)
        value: Path to the bigwig file
    promoter : str
        Path to a bed file, which defines all promoter regions
    motif_meme : str
        Path to a motif database in MEME format
    chromosome_size : str
        Path to a tab-separated file which defines the size of each chromosome
    save_to : str
        Folder to store the figure
    local_registry : str
        Local registration (order) for outputs
    tmp_dir : str
        Path to write tmp files

    Returns
    -------
    histone_result_file : str
        Output file contains results for histone modifications
    (fn_lor, fn_loru, fn_lorl, fn_lorp) : (str, str, str, str)
        Output files for log odds ratio
        upper bounds of LORs
        lower bounds of LORs
        p-values for each comparison from z test (p-q=0?)
    """
    promoter_bed = BedTool(promoter)
    sdf = peak_calls_df.loc[
          (peak_calls_df["Genome release"] == "hg38") & (peak_calls_df["Cell line"] == "K562") & (
                      peak_calls_df["Assay"] == "GROcap"), :]
    # get unique peaks comparing to PINTS (distal only)
    pints_bed = BedTool(sdf.loc[sdf.Tool == "PINTS", "File"].values[0])
    sdf = sdf.loc[sdf.Tool != "PINTS", :]
    assert all([c == 1 for c in sdf.Tool.value_counts()])

    uniq_element_dict = dict()
    for _, row in sdf.iterrows():
        tmp = read_bed(row["File"])
        tmp = tmp.loc[tmp[0].str.startswith("chr"), :]
        if tmp.shape[1] > 6:
            tmp = tmp.loc[:, tmp.columns[:6]]

        tmp[1] = tmp[1].astype(int)
        tmp[2] = tmp[2].astype(int)
        tool_bed = BedTool.from_dataframe(tmp)
        pints_uniq = pints_bed.intersect(tool_bed, v=True).saveas().intersect(promoter_bed, v=True).saveas()
        tool_uniq = tool_bed.intersect(pints_bed, v=True).saveas().intersect(promoter_bed, v=True).saveas()
        pints_shared = pints_bed.intersect(tool_bed, u=True).saveas().intersect(promoter_bed, v=True).saveas()
        tool_shared = tool_bed.intersect(pints_bed, u=True).saveas().intersect(promoter_bed, v=True).saveas()
        logger.info(f"PINTS vs {row['Tool']}: Pu\t{pints_uniq.count()}\tTu\t{tool_uniq.count()}")
        uniq_element_dict[row["Tool"]] = struct_ppc._make((pints_uniq, tool_uniq, pints_shared, tool_shared))

    # histone modifications
    histone_result_file = os.path.join(save_to,
                                       "{global_registry}_{local_registry}_histone.csv".format(
                                           global_registry=global_registry,
                                           local_registry=local_registry))
    if os.path.exists(histone_result_file):
        logger.info(f"Final output file {histone_result_file} exists, skip...")
    else:
        logger.info("Check histone profiles among unique elements...")
        _uniq_ele_histone_profile(dict_ppc=uniq_element_dict, score_bws=dict_signal_track_bws,
                                  file_save_to=histone_result_file,
                                  ref_region_extension=300, chromosome_size=chromosome_size)
        logger.info(f"Result saved to {histone_result_file} (digest: {get_file_hash(histone_result_file)[:7]})")
    # motifs
    logger.info("Profiling motif distribution among unique elements")
    fn_lor, fn_loru, fn_lorl, fn_lorp = _uniq_ele_motif_profile(dict_ppc=uniq_element_dict,
                                                                chromosome_size=chromosome_size, genome_fa=genome_fa,
                                                                output_prefix=os.path.join(save_to,
                                                                                           "{global_registry}_{local_registry}".format(
                                                                                               global_registry=global_registry,
                                                                                               local_registry=local_registry)),
                                                                extension=200, remove_longer_than=1000,
                                                                n_motif=746, tmp_dir=tmp_dir,
                                                                motif_meme=motif_meme,
                                                                ame_exec="export PATH=/programs/meme/bin:$PATH && ame")
    logger.info(f"Result (part 1) saved to {fn_lor} (digest: {get_file_hash(fn_lor)[:7]})")
    logger.info(f"Result (part 2) saved to {fn_loru} (digest: {get_file_hash(fn_loru)[:7]})")
    logger.info(f"Result (part 3) saved to {fn_lorl} (digest: {get_file_hash(fn_lorl)[:7]})")
    logger.info(f"Result (part 4) saved to {fn_lorp} (digest: {get_file_hash(fn_lorp)[:7]})")
    return histone_result_file, fn_lor, fn_loru, fn_lorl, fn_lorp




def main(normal_calls, peak_calls_tech, peak_calls_bio, peak_calls_uniq, bam_cors, data_save_to, 
         gc_files, hg19_to_38_chain_file, hg38_to_19_chain_file, hg19_promoters, 
         hg38_promoters, score_bws, ref_regions, chromosome_size, genome_fa, jaspar_meme, 
         assay_groups=dict(), data_prefix="", **kwargs):
    analysis_summaries = {
        "tu_annotation": [],
        "resolution": [],
        "robustness": [],
        "unique_set": [],
    }

    annotated_TUs = kwargs.get("annotated_TUs", None)
    gencode_ptt_bed = kwargs.get("gencode_ptt_bed", None)
    if annotated_TUs is not None and gencode_ptt_bed is not None:
        tu_results = tu_benchmark(annotated_TUs, gencode_ptt_bed, save_to=data_save_to, local_registry="JaccardTU")
        analysis_summaries["tu_annotation"].append(tu_results)

    peak_calls_pre_df = []
    for k, v in normal_calls.items():
        tool, genome_release, cell_line, assay, bio_rep, tech_rep = k.split("_")
        peak_calls_pre_df.append((tool, genome_release, cell_line, assay, bio_rep, tech_rep, v))

    peak_calls_df = pd.DataFrame(peak_calls_pre_df, columns=("Tool", "Genome release", "Cell line", "Assay",
                                                                "Biorep", "Techrep", "File"))
    # systematic bias
    systematic_bias = evaluate_systematic_bias(peak_call_jobs_df=peak_calls_df, score_bws=score_bws,
                                                ref_regions=ref_regions, chromosome_size=chromosome_size,
                                                save_to=data_save_to,
                                                local_registry="SystematicBias", ref_region_extension=500)
    analysis_summaries["resolution"].append(systematic_bias)

    # TRE sizes
    size_distribution = retrieve_tre_size_dist(peak_call_jobs_df=peak_calls_df, save_to=data_save_to,
                                                local_registry="TRE_sizes")
    analysis_summaries["resolution"].append(size_distribution)

    # GR similarity
    cov_distribution = compile_distributions_of_gc(gc_files, save_to=data_save_to, local_registry="DistGenomeCov")
    analysis_summaries["robustness"].append(cov_distribution)

    # upper bound of robustness
    ub_robustness = evaluate_upper_bound_robustness(peak_calls=normal_calls,
                                                    hg19_to_38_chain_file=hg19_to_38_chain_file,
                                                    hg38_to_19_chain_file=hg38_to_19_chain_file,
                                                    hg19_promoters=hg19_promoters, hg38_promoters=hg38_promoters,
                                                    save_to=data_save_to, local_registry="UpperBoundRobustness")
    analysis_summaries["robustness"].append(ub_robustness)

    # robustness from real case (tech rep/bio rep)
    tech_bio_robustness = evaluate_real_case_robustness(peak_calls_tech=peak_calls_tech,
                                                        peak_calls_bio=peak_calls_bio,
                                                        correlations=bam_cors,
                                                        save_to=data_save_to,
                                                        local_registry="RobustnessRealCases")
    analysis_summaries["robustness"].append(tech_bio_robustness)

    histone_bws = dict()
    for k, v in score_bws.items():
        cell_line, mark = k.split("_")
        if cell_line == "K562" and mark in ("H3K27ac", "H3K4me1"):
            histone_bws[mark] = nfs_mapping(v)
    # for TSScall, use only divergent calls
    peak_calls_pre_df = []
    for k, v in peak_calls_uniq.items():
        tool, genome_release, cell_line, assay, bio_rep, tech_rep = k.split("_")
        peak_calls_pre_df.append((tool, genome_release, cell_line, assay, bio_rep, tech_rep, v))
        
    peak_calls_df = pd.DataFrame(peak_calls_pre_df, columns=("Tool", "Genome release", "Cell line", "Assay",
                                                                 "Biorep", "Techrep", "File"))
    unique_outputs = profile_unique_set(peak_calls_df, genome_fa,
                                        dict_signal_track_bws=histone_bws,
                                        promoter=hg38_promoters,
                                        chromosome_size=chromosome_size,
                                        motif_meme=jaspar_meme,
                                        save_to=data_save_to,
                                        tmp_dir=tmp_dir,
                                        local_registry="UniqueElements")
    analysis_summaries["unique_set"].extend(unique_outputs)

    with open(os.path.join(data_save_to, f"{data_prefix}_summary.json"), "w") as fh:
        json.dump(analysis_summaries, fh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-save-to", required=True, help="Save output data to")
    parser.add_argument("--data-prefix", required=False, default="tools_resolution_robustness",
                        help="Prefix for all outputs")
    parser.add_argument("--global-registry", type=int, help="Global registry", default=3)
    parser.add_argument("--config-file", default="config.txt", help="Configuration file")
    parser.add_argument("--bam-correlation-mat", help="Configuration file",
                        default="/local/storage/ly349/projects/peakcalling/data/best_tool/tech_bio_reps_bam_summaries.cor.txt")
    parser.add_argument("--tmp-dir", default="./tmp", help="All temporary files will be written into this folder")

    args = parser.parse_args()

    assert os.path.exists(args.config_file)

    if not os.path.exists(args.tmp_dir):
        os.mkdir(args.tmp_dir)

    logger.addHandler(logging.FileHandler(os.path.join(os.getcwd(), 'ToolResolutionRobustness.log')))

    cfg = ConfigParser()
    cfg.optionxform = str
    cfg.read(args.config_file)
    
    global tmp_dir, global_registry, full_assays, layouts, unified_color_map, ucm_assay, official_name_map
    tmp_dir = args.tmp_dir
    pybedtools.set_tempdir(tmp_dir)
    global_registry = args.global_registry
    full_assays = cfg.get("assays", "plot_order_full").split("|")
    assay_offical_names = cfg.get("assays", "assay_full_names").split("|")
    assay_plot_color = cfg.get("assays", "plot_colors").split("|")

    layouts = dict()
    for k, v in cfg["dataset_layouts"].items():
        layouts[k] = v

    ucm_assay = dict()
    unified_color_map = dict()
    official_name_map = dict()
    plot_order = cfg.get("tools", "plot_order").split("|")
    plot_color = cfg.get("tools", "plot_color").split("|")
    unified_color_map = dict()
    for k, v in enumerate(full_assays):
        ucm_assay[v] = assay_plot_color[k]
    for k, v in enumerate(plot_order):
        unified_color_map[v] = plot_color[k]
    for k, v in enumerate(cfg.get("assays", "plot_order_full").split("|")):
        official_name_map[v] = assay_offical_names[k]
    official_name_map["GROcap*"] = "GRO-cap"

    tu_calls = {}
    normal_calls = {}
    unique_calls = {}
    tech_calls = {}
    bio_calls = {}
    bi_ref_covs = {}
    if args.method == "generate_data":
        import socket

        server1_home_dir = "/fs/cbsuhy01/storage/ly349/" if socket.gethostname().find(
            "cbsuhy01") == -1 else "/local/storage/ly349/"
        server2_home_dir = "/fs/cbsuhy02/storage/ly349/" if socket.gethostname().find(
            "cbsuhy02") == -1 else "/local/storage/ly349/"
        bioq_dir = os.path.join(server2_home_dir, "BioQueue/workspace")
        from .utils import load_bioq_datasets

        tu_calls = load_bioq_datasets("tu_jobs", bioq_dir, cfg_file=args.config_file)
        normal_calls = load_bioq_datasets("peak_call_job_ids", bioq_dir, cfg_file=args.config_file)
        unique_calls = load_bioq_datasets("peak_call_ue_ids", bioq_dir, cfg_file=args.config_file)
        tech_calls = load_bioq_datasets("peak_call_tech_rep", bioq_dir, cfg_file=args.config_file)
        bio_calls = load_bioq_datasets("peak_call_bio_rep", bioq_dir, cfg_file=args.config_file)
        bi_ref_covs = load_bioq_datasets("genomewide_hist_coverage", bioq_dir, cfg_file=args.config_file)

    # check peak call files and keys make sure their descriptions match
    # peak_call_job_ids:
    #   keys: tool_assay
    for k, v in tu_calls.items():
        if not os.path.exists(v):
            logger.error(f"Cannot locate peak calls for {k} {v}")
        tool, assay = k.split("_")
        tool = tool.replace(".", "")
        if not all([v.find(j) != -1 for j in (tool, assay)]):
            logger.error(f"{k} and {v} don't seem to match with each other")

    # check peak call files and keys make sure their descriptions match
    # peak_call_job_ids:
    #   keys: tool_genome_cellline_assay_biorep_techrep
    for k, v in normal_calls.items():
        if not os.path.exists(v):
            logger.error(f"Cannot locate peak calls for {k} {v}")
        tool, genome, cellline, assay, biorep, techrep = k.split("_")
        tool = tool.replace(".", "")
        if not all([v.find(j) != -1 for j in (tool, genome, assay)]):
            logger.error(f"{k} and {v} don't seem to match with each other")

    # peak_call_tech_rep and peak_call_bio_rep
    #   keys: cellline_assay_biorep_techrep_tool
    for dic in (tech_calls, bio_calls):
        for k, v in dic.items():
            if not os.path.exists(v):
                logger.error(f"Cannot locate peak calls for {k} {v}")
            cellline, assay, biorep, techrep, tool = k.split("_")
            # tool = tool.replace(".", "")
            assay = assay.replace("*", "")
            if v.find(assay) == -1 or (v.find(tool) == -1 and v.find(tool.replace(".", "")) == -1):
                logger.error(f"{k} and {v} don't seem to match with each other")
    end_groups = cfg["dataset_precise_ends"]
    assay_groups = defaultdict(list)
    for k, v in end_groups.items():
        assay_groups[v[:1]].append(k)

    try:
        main(normal_calls=normal_calls, peak_calls_tech=tech_calls,
             peak_calls_bio=bio_calls, peak_calls_uniq=unique_calls, bam_cors=args.bam_correlation_mat,
             data_save_to=args.data_save_to, 
             gc_files=bi_ref_covs, assay_groups=assay_groups,
             hg19_to_38_chain_file=nfs_mapping(cfg.get("references", "liftover_chain_19_38")),
             hg38_to_19_chain_file=nfs_mapping(cfg.get("references", "liftover_chain_38_19")),
             hg19_promoters=nfs_mapping(cfg.get("references", "promoter_hg19")),
             hg38_promoters=nfs_mapping(cfg.get("references", "promoter_hg38")),
             score_bws=cfg["corroborative_bws"], ref_regions=nfs_mapping(cfg.get("references", "true_enhancers")),
             chromosome_size=nfs_mapping(cfg.get("references", "hg38_chromsize_genome")),
             genome_fa=nfs_mapping(cfg.get("references", "hg38_genome")),
             jaspar_meme=nfs_mapping(cfg.get("references", "jaspar_meme")), annotated_TUs=tu_calls,
             gencode_ptt_bed=nfs_mapping(cfg.get("references", "gencode_highest_pt_transcripts")),
             data_prefix=args.data_prefix)
    except Exception as e:
        logger.exception(e)
    finally:
        pybedtools.cleanup()
