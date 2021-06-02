#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 8/16/20
import argparse
import json
import logging
import os
import sys
from configparser import ConfigParser
from file_parser_mem.bed import parse_bed
import numpy as np
import pandas as pd
import pyBigWig
import pybedtools
from pybedtools import BedTool
from .utils import bin_scores, read_bed, midpoint_generator, run_command

logging.basicConfig(format='%(name)s - %(asctime)s - %(levelname)s: %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO,
                    handlers=[
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("PINTS - Resolution and robustness")


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

    lifted_hg19, cached_raw_hg19 = _liftover(raw_hg19, chain_file=hg19_to_38_chain_file, mark=19)
    lifted_hg38, cached_raw_hg38 = _liftover(raw_hg38, chain_file=hg38_to_19_chain_file, mark=38)

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


def main(normal_calls, peak_calls_tech, peak_calls_bio, bam_cors, data_save_to, 
         gc_files, hg19_to_38_chain_file, hg38_to_19_chain_file, hg19_promoters, 
         hg38_promoters, score_bws, ref_regions, chromosome_size, data_prefix=""):
    analysis_summaries = {
        "resolution": [],
        "robustness": [],
    }

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
    collapsed_assays = cfg.get("assays", "plot_order_simplified").split("|")

    global tmp_dir, n_threads, n_samples, global_registry, full_assays, highlight_assays, layouts, unified_color_map, official_name_map
    tmp_dir = args.tmp_dir
    pybedtools.set_tempdir(tmp_dir)
    n_threads = int(cfg.get("global", "n_threads"))
    global_registry = args.global_registry
    full_assays = cfg.get("assays", "plot_order_full").split("|")
    highlight_assays = cfg.get("assays", "plot_order_simplified").split("|")
    assay_offical_names = cfg.get("assays", "assay_full_names").split("|")
    layouts = dict()
    for k, v in cfg["dataset_layouts"].items():
        layouts[k] = v

    unified_color_map = dict()
    official_name_map = dict()
    plot_order = cfg.get("tools", "plot_order").split("|")
    plot_color = cfg.get("tools", "plot_color").split("|")
    unified_color_map = dict()
    for k, v in enumerate(plot_order):
        unified_color_map[v] = plot_color[k]
    for k, v in enumerate(cfg.get("assays", "plot_order_full").split("|")):
        official_name_map[v] = assay_offical_names[k]
    official_name_map["GROcap*"] = "GRO-cap"

    normal_calls = {}
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

        normal_calls = load_bioq_datasets("peak_call_job_ids", bioq_dir, cfg_file=args.config_file)
        tech_calls = load_bioq_datasets("peak_call_tech_rep", bioq_dir, cfg_file=args.config_file)
        bio_calls = load_bioq_datasets("peak_call_bio_rep", bioq_dir, cfg_file=args.config_file)
        bi_ref_covs = load_bioq_datasets("genomewide_hist_coverage", bioq_dir, cfg_file=args.config_file)

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
    try:
        main(normal_calls=normal_calls, peak_calls_tech=tech_calls,
             peak_calls_bio=bio_calls, bam_cors=args.bam_correlation_mat,
             data_save_to=args.data_save_to, 
             gc_files=bi_ref_covs,
             hg19_to_38_chain_file=cfg.get("references", "liftover_chain_19_38"),
             hg38_to_19_chain_file=cfg.get("references", "liftover_chain_38_19"),
             hg19_promoters=cfg.get("references", "promoter_hg19"),
             hg38_promoters=cfg.get("references", "promoter_hg38"),
             score_bws=cfg["corroborative_bws"], ref_regions=cfg.get("references", "true_enhancers"),
             chromosome_size=cfg.get("references", "hg38_chromsize_genome"),
             data_prefix=args.data_prefix)
    except Exception as e:
        logger.exception(e)
    finally:
        pybedtools.cleanup()
