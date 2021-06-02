#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 8/9/20
import logging
import argparse
import os
import pyBigWig
import gzip
import json
import pybedtools
import numpy as np
import pandas as pd
from pybedtools import BedTool
from multiprocessing import Pool
from configparser import ConfigParser
from .utils import bed_bed_strand_specific_coverage, bin_scores, contrast_regions, run_command

logging.basicConfig(format="%(name)s - %(asctime)s - %(levelname)s: %(message)s",
                    datefmt="%d-%b-%y",
                    level=logging.INFO,
                    handlers=[
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("Assay comparison")


def parse_read_distribution(exp_name, input_bed, genome_segmentation):
    """
    Parse distribution of reads

    Parameters
    ----------
    exp_name : str
        Name of the assay
    input_bed : str
        Path to the aligned bed file
    genome_segmentation : BedTool object
        BedTool object for the segmentation
    Returns
    -------
    exp_name : str
        Name of the assay
    x : pd.Series
        Series containing read counts in promoters, introns, exons and intergenic regions
    """
    x = pd.Series({"Promoter": 0, "Intron": 0, "Exon": 0, "Intergenic region": 0})
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    try:
        exp_bed = BedTool(input_bed)
        exp_gencode = exp_bed.intersect(genome_segmentation, wb=True)
        pre_sort = os.path.join(tmp_dir, "%s_presort.bed" % exp_name)
        after_sort = os.path.join(tmp_dir, "%s_sorted.bed" % exp_name)
        exp_gencode.moveto(pre_sort)
        _, stderr, rc = run_command("sort -k4,4 %s > %s" % (pre_sort, after_sort))
        assert rc == 0, "sort failed %s" % stderr
        os.remove(pre_sort)
        prev_read = None
        categories = set()
        with open(after_sort) as fh:
            for nr, row in enumerate(fh):
                items = row.strip().split("\t")
                if prev_read is None:
                    prev_read = items[3]
                if prev_read == items[3]:
                    categories.add(items[10])
                else:
                    if "promoter" in categories or "promoter(NP)" in categories:
                        x["Promoter"] += 1
                    elif "intron" in categories:
                        x["Intron"] += 1
                    elif "5_UTR" in categories:
                        x["Exon"] += 1
                    elif "3_UTR" in categories:
                        x["Exon"] += 1
                    elif "exon" in categories:
                        x["Exon"] += 1
                    else:
                        x["Intergenic region"] += 1
                    categories = set()
                    categories.add(items[10])
                    prev_read = items[3]
        return exp_name, x
    except ValueError as e:
        logger.error(e, exp_name, input_bed)


def read_distribution_atom(bed_files, segmentation, save_to, local_registry):
    """
    Atom for parsing read distribution

    Parameters
    ----------
    bed_files : list or tuple
        List of bed files, each bed file contains a wildcard
    segmentation : BedTool
        BedTool object for the genomic segmentation
    save_to : str
        Path to write outputs
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
        for i in range(1, n_samples + 1):
            args = []
            for k, bed in bed_files.items():
                cell_line, exp = k.split("_")
                if cell_line.startswith("K562"):
                    args.append((exp, bed % i, segmentation))
            with Pool(n_threads) as pool:
                read_distributions = pool.starmap(parse_read_distribution, args)
            rd_cols = []
            rd_colnames = []
            for exp, summary in read_distributions:
                rd_cols.append(summary)
                rd_colnames.append(exp)

            read_distribution_df = pd.concat(rd_cols, axis=1, sort=False)
            read_distribution_df.columns = rd_colnames
            read_distribution_df = read_distribution_df.T
            read_distribution_df.to_csv(os.path.join(save_to, "{global_registry}_{local_registry}_rep{rep}.csv".format(
                global_registry=global_registry, local_registry=local_registry, rep=i)))

        rd_rep1 = pd.read_csv(os.path.join(save_to,
                                           "{global_registry}_{local_registry}_rep1.csv".format(
                                               global_registry=global_registry,
                                               local_registry=local_registry)), index_col=0)
        rd_rep2 = pd.read_csv(os.path.join(save_to, "{global_registry}_{local_registry}_rep2.csv".format(
            global_registry=global_registry,
            local_registry=local_registry)), index_col=0)
        rd_rep3 = pd.read_csv(os.path.join(save_to, "{global_registry}_{local_registry}_rep3.csv".format(
            global_registry=global_registry,
            local_registry=local_registry)), index_col=0)
        read_distribution_df = (rd_rep1 + rd_rep2 + rd_rep3) / 3

        read_distribution_df.to_csv(final_result_file)
    return final_result_file


def get_read_counts_from_bed(bed_files, save_to, global_registry, local_registry):
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        results = []
        # for bed, exp in zip(bed_files, labels):
        for k, bed in bed_files.items():
            cell_line, exp = k.split("_")
            if cell_line.startswith("K562"):
                fn_open = gzip.open if bed.endswith('.gz') else open
                mode = "rt" if bed.endswith('.gz') else "r"
                results.append((exp, sum([1 for t in fn_open(bed, mode)])))
        pd.DataFrame(results, columns=("Assay", "Reads")).to_csv(final_result_file)
    return final_result_file


def _parse_read_distribution(exp_name, input_bed, abundant_rna_bed, tmp_dir="."):
    """
    Parse read distribution and infer proportion of reads from abundant transcripts

    Parameters
    ----------
    exp_name : str
        Source of the input bed
    input_bed : str
        path to a bed file (describing reads)
    abundant_rna_bed : str
        bed annotations for abundant transcripts
    tmp_dir : str
        tmp dir for saving intermediate files

    Returns
    -------
    dist : pd.Series
        # of reads mapped to snRNA, tRNA, snoRNA, scRNA or miRNA loci
    """
    x = pd.Series({
        "snRNA": 0, "rRNA": 0, "tRNA": 0, "snoRNA": 0, "scRNA": 0, "miRNA": 0})
    try:
        exp_bed = BedTool(input_bed)
        exp_arb = exp_bed.intersect(abundant_rna_bed, s=False, wb=True)
        pre_sort = os.path.join(tmp_dir, f"{exp_name}_presort.bed")
        after_sort = os.path.join(tmp_dir, f"{exp_name}_sorted.bed")
        exp_arb.moveto(pre_sort)
        # sort by reads names
        _, err, rc = run_command("sort -k4,4 %s > %s" % (pre_sort, after_sort))
        assert rc == 0, err
        os.remove(pre_sort)
        prev_read = None
        categories = set()
        with open(after_sort) as fh:
            for nr, row in enumerate(fh):
                items = row.strip().split("\t")
                if prev_read is None:
                    prev_read = items[3]
                if prev_read == items[3]:
                    categories.add(items[10])
                else:
                    if "snRNA" in categories:
                        x["snRNA"] += 1
                    elif "rRNA" in categories:
                        x["rRNA"] += 1
                    elif "tRNA" in categories:
                        x["tRNA"] += 1
                    elif "snoRNA" in categories:
                        x["snoRNA"] += 1
                    elif "scRNA" in categories:
                        x["scRNA"] += 1
                    elif "miRNA" in categories:
                        x["miRNA"] += 1
                    categories = set()
                    categories.add(items[10])
                    prev_read = items[3]
        os.remove(after_sort)
        stdout, stderr, ec = run_command("grep '^chrrRNA|^U13369' %s | wc -l" % input_bed)
        rrna = int(stdout.strip()[0])
        x["rRNA"] += rrna
        return exp_name, x
    except ValueError as e:
        logger.error(f"{e} {exp_name} {input_bed}")


def get_reads_from_abundant_transcripts(bed_files, abundant_rna_file, save_to, global_registry, local_registry,
                                        n_threads=16, tmp_dir="."):
    """
    Factory function for getting reads from abundant transcripts

    Parameters
    ----------
    bed_files : dict
        dict of bed files, each bed file represents all mapped reads from a corresponding bam file
    abundant_rna_file : str
        bed annotations for abundant transcripts
    save_to : str
        Path to write outputs
    global_registry : int
        Global registration for outputs
    local_registry : str
        Local registration (order) for outputs
    n_threads : int
        Number of threads that this function can create, by default, 16
    tmp_dir : str


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
        args = []
        abundant_rna_bed = BedTool(abundant_rna_file)
        # for exp, sample in zip(labels, bed_files):
        for k, sample in bed_files.items():
            cell_line, exp = k.split("_")
            if cell_line.startswith("K562"):
                args.append((exp, sample, abundant_rna_bed, tmp_dir))
            else:
                logger.warning(f"Abundant transcript not evaluated {k}, as it's not K562-derived.")

        with Pool(n_threads) as pool:
            rep_result = pool.starmap(_parse_read_distribution, args)
        pybedtools.cleanup()

        rd_cols = []
        rd_colnames = []
        for exp, summary in rep_result:
            rd_cols.append(summary)
            rd_colnames.append(exp)

        abundant_ele_distribution = pd.concat(rd_cols, axis=1, sort=False)
        abundant_ele_distribution.columns = rd_colnames
        abundant_ele_distribution.to_csv(final_result_file)
    return final_result_file


def _atom_ss_agg(bw_path, bed_5ss_flanked, bed_3ss_flanked, save_to, name, bins=400):
    """
    Atom function for generate stats around splicing sites

    Parameters
    ----------
    bw_path : str
        Path to a bw file which provides coverage information
    bed_5ss_flanked : BedTool
        BedTool objects describing 5' splicing sites
    bed_3ss_flanked : BedTool
        BedTool objects describing 3' splicing sites
    save_to : str
        Path to write outputs
    name : str
        Name of the assay, which will affect the output (name_5ss.npz and name_3ss.npz)
    bins : int
        optional, by default: 400.

    Returns
    -------

    """
    with pyBigWig.open(bw_path) as bw:
        sm, m_5, (l_5, u_5) = bin_scores(bed_5ss_flanked, bw, bins=bins, n_boot=1000)
        # np.savez(f"{save_to}/{name}_5ss", sm)
        sm, m_3, (l_3, u_3) = bin_scores(bed_3ss_flanked, bw, bins=bins, n_boot=1000)
        # np.savez(f"{save_to}/{name}_3ss", sm)
    return m_5, (l_5, u_5), m_3, (l_3, u_3)


def get_reads_around_splicing_sites(generic_bws, bed_5ss, bed_3ss, save_to, local_registry, bins=400):
    """

    Parameters
    ----------
    generic_bws : dict
        A dict contains bw files which provides coverage information
    bed_5ss : BedTool
        Path to a bed file describing 5' splicing sites (and flanking regions)
    bed_3ss : BedTool
        Path to a bed file describing 3' splicing sites (and flanking regions)
    save_to : str
        Path to write outputs
    local_registry : str
        Local registration (order) for outputs
    bins : int
        optional, by default: 400.

    Returns
    -------

    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    labels = []
    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        with Pool(n_threads) as pool:
            jobs = []
            for k, bw in generic_bws.items():
                cell_line, exp = k.split("_")
                if cell_line.startswith("K562"):
                    labels.append(exp)
                    jobs.append((bw, bed_5ss, bed_3ss, save_to, exp, bins))
            direct_results = pool.starmap(_atom_ss_agg, jobs)

        results = []

        for i, assay in enumerate(labels):
            m_5, (l_5, u_5), m_3, (l_3, u_3) = direct_results[i]
            results.append(pd.DataFrame({"Mean": m_5, "LB": l_5, "UB": u_5, "Assay": assay, "End": 5}))
            results.append(pd.DataFrame({"Mean": m_3, "LB": l_3, "UB": u_3, "Assay": assay, "End": 3}))

        pd.concat(results).to_csv(final_result_file)
    return final_result_file


def get_signal_to_noise_ratios(bed_files, save_to, pos_region_file, neg_region_file, promoter_file,
                               global_registry, local_registry):
    """
    Get signal to noise ratios

    Parameters
    ----------
    bed_files : dict
        keys: cell_line+"_"+assay_name
    save_to : str
        Save result to this place
    pos_region_file : str
        Path to bed file defining positive regions
    neg_region_file : str
        Path to bed file defining negative regions
    promoter_file : str
        Path to bed file defining promoter regions
    global_registry : int
        Global prefix
    local_registry : str
        Local prefix

    Returns
    -------
    final_result_file : str
        Full path to the result file
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))

    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        pos_regions = pd.read_csv(pos_region_file, sep="\t", header=None)
        pos_regions[4] = "."
        pos_regions[5] = "+"
        neg_regions = pd.read_csv(neg_region_file, sep="\t", header=None)
        neg_regions[4] = "."
        neg_regions[5] = "+"
        distribution_dict = dict()
        pos_ref_bed = BedTool.from_dataframe(pos_regions)
        neg_ref_bed = BedTool.from_dataframe(neg_regions).intersect(BedTool(promoter_file), v=True)
        exps = []
        for ref, distribution_key in zip((pos_ref_bed, neg_ref_bed), ("te_distribution", "ne_np_distribution")):
            tmp_distributions = []
            for rep in range(n_samples):
                args = []

                for k, sample in bed_files.items():
                    cell_line, exp = k.split("_")
                    if cell_line.startswith("K562"):
                        exps.append(exp)
                        args.append((ref, sample % (rep + 1), exp))

                with Pool(n_threads) as pool:
                    rep_result = pool.starmap(bed_bed_strand_specific_coverage, args)
                tmp_distributions.append(pd.concat(rep_result, axis=1))
            distribution_dict[distribution_key] = sum(tmp_distributions) / n_samples

        pos_neg_cov = []
        for exp in exps:
            for nr, row in distribution_dict["te_distribution"].iterrows():
                pos_neg_cov.append((np.log1p(row[f"{exp}_fwd"] + row[f"{exp}_rev"]), exp, "True enhancers"))
            for nr, row in distribution_dict["ne_np_distribution"].iterrows():
                pos_neg_cov.append((np.log1p(row[f"{exp}_fwd"] + row[f"{exp}_rev"]), exp, "Non-enhancers"))
        pos_neg_cov_df = pd.DataFrame(pos_neg_cov, columns=("Coverage", "Assay", "Category"))
        pos_neg_cov_df.to_csv(final_result_file)
        pybedtools.cleanup()
    return final_result_file


def _atom_assay_confusion_mat(pos_ref, neg_ref, bin_df, top_bin_cut, assay):
    X = bin_df[assay]
    X = X[X > 0]
    ranked = X.sort_values(ascending=False)
    top_n_candidate = ranked.iloc[:top_bin_cut, ]
    candidate_bed = BedTool(list(map(lambda x: x.replace("-", ":").split(":"), top_n_candidate.index.tolist())))
    tp = len(pos_ref.intersect(candidate_bed, u=True))
    fp = len(neg_ref.intersect(candidate_bed, u=True))
    fn = len(pos_ref.intersect(candidate_bed, v=True))
    tn = len(neg_ref.intersect(candidate_bed, v=True))
    return tp, fp, fn, tn


def calculate_fdr_per_assay(bins_df_file, pos_ref, neg_ref, save_to, global_registry, local_registry):
    """
    Calculate FDR per assay

    Parameters
    ----------
    bins_df_file : str
        Path to a tab file which contains read coverages in genomic bins
    pos_ref : str
        Path to bed file defining positive regions
    neg_ref : str
        Path to bed file defining negative regions
    save_to : str
        Save result to this place
    global_registry : int
        Global prefix
    local_registry : str
        Local prefix

    Returns
    -------
    final_result_file : str
        Full path to the result file
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))

    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        bins_df = pd.read_csv(bins_df_file, sep="\t")
        bins_df.index = bins_df["chr"] + ":" + bins_df["start"].map(str) + "-" + bins_df["end"].map(str)
        bins_df.drop(columns=["chr", "start", "end"], inplace=True)
        bins_df = bins_df.loc[bins_df.sum(axis=1) > 0, :]
        fdr_results = []
        pos_ref_bed = BedTool(pos_ref)
        neg_ref_bed = BedTool(neg_ref)
        job_pars = []
        job_info = []
        for tbc in (5_000, 10_000, 20_000, 100_000):
            for col in bins_df.columns:
                assay, rep = col.split("_")
                job_pars.append((pos_ref_bed, neg_ref_bed, bins_df, tbc, col))
                job_info.append([assay, tbc])
        with Pool(16) as pool:
            job_results = pool.starmap(_atom_assay_confusion_mat, job_pars)
        for ji, jr in zip(job_info, job_results):
            # jr:  0  1  2  3
            #     tp fp fn tn
            # fdr = fp/(tp+fp)
            fdr_results.append(ji + [jr[1] / (jr[0] + jr[1]), ])
        fdr_df = pd.DataFrame(fdr_results, columns=("Assay", "Cutoff", "FDR"))
        fdr_df.to_csv(final_result_file)
        pybedtools.cleanup()
    return final_result_file


def contrast_te_ne(te, ne, score_bw_dict, chromosome_size, save_to, global_registry, local_registry,
                   region_extension=1000, n_bins=100):
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))

    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        regions = [te, ne]
        labels = ("True enhancer", "Non-enhancer")
        df = contrast_regions(regions, labels, score_bw_dict, region_extension, chromosome_size, n_bins=n_bins)
        df.to_csv(final_result_file)
        pybedtools.cleanup()
    return final_result_file


def main(data_save_to, data_prefix="", **kwargs):
    """
    Factory

    Parameters
    ----------
    data_save_to : str
        generated data will be exported to this path
    data_prefix : str
        All outputs from `generate_data` will have this prefix
    kwargs :
        bed_files : dict
            Alignment in bed format
        dsc_beds : dict
            Bed files from downsampled and cleaned alignments
        precise_end_bw_files : dict
            Dict of bigwig files, which contains signal density from the precisely/approximately mapped ends
        segmentations : list or list-like
            bed files defining genomic segmentations
        te_bed : str
            Path to the true enhancer bed file
        ne_bed : str
            Path to the non-enhancer bed file
        promoter_bed : str
            Path to the promoter bed file
        gb_cov_tab : str
            Path to a tab-separated file which contains coverage info for genomic bins (multiple samples)
        cor_bws : dict
            keys : name of the bw in the form of `CL_MARK`
            values : paths to the bw files
        chromosome_size : str
            Path to a tab-separated file which defines the sizes of chromosomes
        abundant_rna : str
            Path to a bed file, which defines genomic regions for abundant RNAs
        file_5ss : str
            Path to a bed file, which defines the 5' Splicing Sites
        file_3ss : str
            Path to a bed file, which defines the 3' SS

    Returns
    -------

    """
    analysis_summaries = {
        "read_distribution": [],
        "read_counts": [],
        "abundant_transcripts": [],
        "splicing_intermediates": [],
        "signal_to_noise": []
    }

    dsc_beds = kwargs.pop("dsc_beds")
    bed_files = kwargs.pop("bed_files")
    precise_end_bw_files = kwargs.pop("precise_end_bw_files")
    segmentations = kwargs.pop("segmentations")
    abundant_rna = kwargs.pop("abundant_rna")
    file_5ss = kwargs.pop("file_5ss")
    file_3ss = kwargs.pop("file_3ss")
    te_bed = kwargs.pop("te_bed")
    ne_bed = kwargs.pop("ne_bed")
    promoter_bed = kwargs.pop("promoter_bed")
    gb_cov_tab = kwargs.pop("gb_cov_tab")
    cor_bws = kwargs.pop("cor_bws")
    chromosome_size = kwargs.pop("chromosome_size")

    # read distribution
    i = 0
    for segmentation in segmentations:
        seg_bed = BedTool(segmentation)
        distribution_dat = read_distribution_atom(bed_files=dsc_beds, segmentation=seg_bed,
                                                    save_to=data_save_to, local_registry=f"ReadDistribution{i}")
        analysis_summaries["read_distribution"].append(distribution_dat)
        i += 1
    # total read counts
    read_counts = get_read_counts_from_bed(bed_files=bed_files, save_to=data_save_to,
                                            global_registry=global_registry,
                                            local_registry="TotalReadCounts")
    analysis_summaries["read_counts"].append(read_counts)

    # reads from abundant transcripts
    abundant_transcripts_dat = get_reads_from_abundant_transcripts(bed_files=bed_files,
                                                                    abundant_rna_file=abundant_rna,
                                                                    save_to=data_save_to,
                                                                    global_registry=global_registry,
                                                                    local_registry="AbundantTranscript",
                                                                    tmp_dir=tmp_dir)
    analysis_summaries["abundant_transcripts"].append(abundant_transcripts_dat)

    # splicing intermediates
    splicing_dat = get_reads_around_splicing_sites(generic_bws=precise_end_bw_files,
                                                    bed_5ss=BedTool(file_5ss),
                                                    bed_3ss=BedTool(file_3ss), save_to=data_save_to,
                                                    local_registry="SplicingIntermediates")
    analysis_summaries["splicing_intermediates"].append(splicing_dat)

    # FC
    te_ne_fc = get_signal_to_noise_ratios(bed_files=dsc_beds, save_to=data_save_to,
                                            pos_region_file=te_bed, neg_region_file=ne_bed,
                                            promoter_file=promoter_bed, global_registry=global_registry,
                                            local_registry="FC")
    analysis_summaries["signal_to_noise"].append(te_ne_fc)

    # FDR
    fdr_res = calculate_fdr_per_assay(bins_df_file=gb_cov_tab, pos_ref=te_bed, neg_ref=ne_bed,
                                        save_to=data_save_to, global_registry=global_registry,
                                        local_registry="FDR")
    analysis_summaries["signal_to_noise"].append(fdr_res)

    # TE vs. NE
    score_bw_dict = {
        "K562_DHS": cor_bws["K562_DHS"],
        "K562_H3K27ac": cor_bws["K562_H3K27ac"],
        "K562_CTCF": cor_bws["K562_CTCF"],
    }
    r = contrast_te_ne(te=te_bed, ne=ne_bed, score_bw_dict=score_bw_dict, chromosome_size=chromosome_size,
                        save_to=data_save_to, global_registry=global_registry, local_registry="TrueEvsFalseE",
                        region_extension=1000, n_bins=100)
    analysis_summaries["signal_to_noise"].append(r)

    with open(os.path.join(data_save_to, f"{data_prefix}_summary.json"), "w") as fh:
        json.dump(analysis_summaries, fh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-save-to", required=True, help="Save output data to")
    parser.add_argument("--abundant-rna", required=False,
                        help="Abundant RNA loci defined in a bed file")
    parser.add_argument("--te-bed", required=False, help="True enhancers defined in a bed file")
    parser.add_argument("--ne-bed", required=False,
                        default="/local/storage/ly349/projects/peakcalling/data/refs/compiled/Nonenhancers_nonepls.bed",
                        help="Bed file which includes definitions for non-enhancers")
    parser.add_argument("--bed-promoter", required=False,
                        default="/local/storage/ly349/refs/annotations/human/hg38/segmentation/promoters_1kb_tss_centered.bed",
                        help="Bed file which includes definitions for promoters")
    parser.add_argument("--bed-5ss", required=False,
                        help="5' splicing sites and flanking regions (400bp) defined in a bed file")
    parser.add_argument("--bed-3ss", required=False,
                        help="3' splicing sites and flanking regions (400bp) defined in a bed file")
    parser.add_argument("--data-prefix", required=False, default="specificity_reads", help="Prefix for all outputs")
    parser.add_argument("--global-registry", type=int, help="Global registry", default=2)
    parser.add_argument("--genomic-segmentations", nargs="*", required=False,
                        help="Genomic segmentations, each file will be parsed separately")
    parser.add_argument("--config-file", default="config.conf", help="Configuration file")
    parser.add_argument("--tmp-dir", default="./tmp", help="All temporary files will be written into this folder")

    args = parser.parse_args()

    assert os.path.exists(args.config_file)

    if not os.path.exists(args.tmp_dir):
        os.mkdir(args.tmp_dir)

    logger.addHandler(logging.FileHandler(os.path.join(os.getcwd(), 'AssaySpecificityReads.log')))

    cfg = ConfigParser()
    cfg.optionxform = str
    cfg.read(args.config_file)
    collapsed_assays = cfg.get("assays", "plot_order_simplified").split("|")
    end_groups = cfg["dataset_precise_ends"]

    global tmp_dir, n_threads, n_samples, global_registry, full_assays, highlight_assays
    global adapters, layouts, unified_color_map, unified_color_map_s
    global ci_upper, ci_lower
    tmp_dir = args.tmp_dir
    pybedtools.set_tempdir(tmp_dir)
    n_threads = int(cfg.get("global", "n_threads"))
    ci_upper = float(cfg.get("global", "ci_upper"))
    ci_lower = float(cfg.get("global", "ci_lower"))
    n_samples = int(cfg.get("assays", "n_downsamples"))
    global_registry = args.global_registry
    full_assays = cfg.get("assays", "plot_order_full").split("|")
    highlight_assays = cfg.get("assays", "plot_order_simplified").split("|")
    assay_offical_names = cfg.get("assays", "assay_full_names").split("|")
    layouts = dict()
    for k, v in cfg["dataset_layout_conversion"].items():
        layouts[k] = v
    adapters = dict()
    for k, v in cfg["RT_primer"].items():
        adapters[k] = v
    unified_color_map = dict()
    official_name_map = dict()
    plot_order = cfg.get("assays", "plot_order_full").split("|")
    plot_order_simplified = cfg.get("assays", "plot_order_simplified").split("|")
    plot_color = cfg.get("assays", "plot_colors").split("|")
    cor_bws = cfg["corroborative_bws"]
    for k, v in enumerate(plot_order):
        unified_color_map[v] = plot_color[k]
        official_name_map[v] = assay_offical_names[k]
    unified_color_map_s = dict()
    for k, v in enumerate(plot_order_simplified):
        unified_color_map_s[v] = unified_color_map[v]

    p3_bws = dict()
    pe_bws = dict()
    beds = dict()
    dsc_beds = dict()
    gb_cov_tab = None

    import socket

    server1_home_dir = "/fs/cbsuhy01/storage/ly349/" if socket.gethostname().find(
        "cbsuhy01") == -1 else "/local/storage/ly349/"
    server2_home_dir = "/fs/cbsuhy02/storage/ly349/" if socket.gethostname().find(
        "cbsuhy02") == -1 else "/local/storage/ly349/"
    bioq_dir = os.path.join(server2_home_dir, "BioQueue/workspace")

    from .utils import load_bioq_datasets

    pe_bws = load_bioq_datasets("unique_plmn_bigwig_merged", bioq_dir, cfg_file=args.config_file)
    beds = load_bioq_datasets("all_alignments_merged_beds", bioq_dir, cfg_file=args.config_file)
    dsc_beds = load_bioq_datasets("downsampled_clean_bed", bioq_dir, cfg_file=args.config_file)
    p3_bws = load_bioq_datasets("k562_3p_bw_prefixs", bioq_dir, cfg_file=args.config_file)
    if gb_cov_tab is None:
        gb_covs = load_bioq_datasets("genomewide_bin_coverage", bioq_dir, cfg_file=args.config_file)
        gb_cov_tab = gb_covs["500"]

    main(data_save_to=args.data_save_to, data_prefix=args.data_prefix,
         bed_files=beds, dsc_beds=dsc_beds, precise_end_bw_files=pe_bws,
         segmentations=args.genomic_segmentations, te_bed=args.te_bed, ne_bed=args.ne_bed,
         promoter_bed=args.bed_promoter, gb_cov_tab=gb_cov_tab, cor_bws=cor_bws,
         chromosome_size=cfg.get("references", "hg38_chromsize_genome"), abundant_rna=args.abundant_rna,
         file_5ss=args.bed_5ss, file_3ss=args.bed_3ss)
