#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 7/25/20
import logging
import argparse
import os
import json
import pyBigWig
import pysam
import pybedtools
import numpy as np
import pandas as pd
from configparser import ConfigParser
from pybedtools import BedTool
from multiprocessing import Pool
from pyfaidx import Fasta
from statsmodels.stats.contingency_tables import Table2x2
from .utils import bed_bed_strand_specific_coverage

logging.basicConfig(format="%(name)s - %(asctime)s - %(levelname)s: %(message)s",
                    datefmt="%d-%b-%y",
                    level=logging.INFO,
                    handlers=[
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("Assay comparison")


# background 3mer counts from GRCh38
BG_3MER_RATES = {'AAA': 112389107, 'AAC': 43364499, 'AAG': 58394303, 'AAT': 72679287, 'ACA': 59215433,
                 'ACC': 33810569, 'ACG': 7528138, 'ACT': 47385805, 'AGA': 65342927, 'AGC': 41013796,
                 'AGG': 51724957, 'AGT': 47335671, 'ATA': 60258895, 'ATC': 39068740, 'ATG': 53637448,
                 'ATT': 73269819, 'CAA': 55232450, 'CAC': 43977022, 'CAG': 59642556, 'CAT': 53874768,
                 'CCA': 53393897, 'CCC': 38076564, 'CCG': 8035900, 'CCT': 51864644, 'CGA': 6518716,
                 'CGC': 7007139, 'CGG': 8205096, 'CGT': 7572995, 'CTA': 37654024, 'CTC': 49410743,
                 'CTG': 59057983, 'CTT': 59158393, 'GAA': 58838074, 'GAC': 27691396, 'GAG': 49489551,
                 'GAT': 39482704, 'GCA': 42412118, 'GCC': 34519459, 'GCG': 7054930, 'GCT': 40689419,
                 'GGA': 45941459, 'GGC': 34505028, 'GGG': 38179438, 'GGT': 33798585, 'GTA': 33222715,
                 'GTC': 27479150, 'GTG': 44457727, 'GTT': 43092672, 'TAA': 60367538, 'TAC': 32907039,
                 'TAG': 37890853, 'TAT': 60198123, 'TCA': 57705337, 'TCC': 44964401, 'TCG': 6684985,
                 'TCT': 65341272, 'TGA': 57698560, 'TGC': 42149982, 'TGG': 54315027, 'TGT': 59545028,
                 'TTA': 60227918, 'TTC': 58737524, 'TTG': 56555437, 'TTT': 113628363}
N_BG_3MERS = sum(BG_3MER_RATES.values())


def bidirectional_coverage_among_ref_regions(bed_files, ref_bed, save_to, global_registry, local_registry,
                                             detectable_threshold=(5,), n_samples=3, n_threads=16):
    """
    Get read coverages among reference regions defined in a bed file

    Parameters
    ----------
    bed_files : list or tuple
        List of bed files, each bed file contains a wildcard
    ref_bed : str
        Path to the reference bed file
    save_to : str
        Path to write outputs
    global_registry : int
        Global registration for outputs
    local_registry : tuple or list
        Local registrations (order) for outputs
    detectable_threshold : tuple of ints
        Lower bound of read counts among a region to be considered as detected
    n_samples : int
        Number of downsamples that each key has in `bed_files`
    n_threads : int
        Max number of threads this function can create

    Returns
    -------
    aggregate_result : str
        csv file to aggregated result
    pairwise_result : str
        csv file to pairwise result
    """
    final_result_file1 = os.path.join(save_to,
                                      "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                      local_registry=local_registry[0]))
    final_result_file2 = os.path.join(save_to,
                                      "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                      local_registry=local_registry[1]))
    if os.path.exists(final_result_file1) and os.path.exists(final_result_file2):
        logger.info(f"Final output file {final_result_file1} and {final_result_file2} exist, skip...")
    else:
        # init counts df
        pos_regions = pd.read_csv(ref_bed, sep="\t", header=None)
        pos_regions[4] = "."
        pos_regions[5] = "+"
        distribution_per_rep = []
        labels = []
        for rep in range(n_samples):
            args = []
            pos_regions_pl_bed = BedTool.from_dataframe(pos_regions)
            for k, sample in bed_files.items():
                cell_line, exp = k.split("_")
                if cell_line.startswith("K562"):
                    labels.append(exp)
                    args.append((pos_regions_pl_bed, sample % (rep + 1), exp))

            with Pool(n_threads) as pool:
                rep_result = pool.starmap(bed_bed_strand_specific_coverage, args)
            distribution_per_rep.append(pd.concat(rep_result, axis=1))
        distribution = sum(distribution_per_rep) / n_samples

        sub_dfs = []
        for threshold in detectable_threshold:
            d_result = pd.DataFrame(0,
                                    index=[exp.replace("K562_", "") if i == 0 else exp.replace("K562_", "") + f"_{i}"
                                           for i
                                           in range(n_samples) for exp in bed_files],
                                    columns=("Bidirectional", "Unidirectional"))
            for exp in labels:
                d_result.loc[exp, "Bidirectional"] = sum(np.logical_and(distribution["%s_fwd" % exp] > threshold,
                                                                        distribution[
                                                                            "%s_rev" % exp] > threshold))
                d_result.loc[exp, "Unidirectional"] = sum(np.logical_or(distribution["%s_fwd" % exp] > threshold,
                                                                        distribution[
                                                                            "%s_rev" % exp] > threshold)) - \
                                                      d_result.loc[exp, "Bidirectional"]
                for i in range(n_samples):
                    d_result.loc[exp + f"_{i}", "Bidirectional"] = sum(
                        np.logical_and(distribution_per_rep[i]["%s_fwd" % exp] > threshold,
                                       distribution_per_rep[i][
                                           "%s_rev" % exp] > threshold))
                    d_result.loc[exp + f"_{i}", "Unidirectional"] = sum(
                        np.logical_or(distribution_per_rep[i]["%s_fwd" % exp] > threshold,
                                      distribution_per_rep[i][
                                          "%s_rev" % exp] > threshold)) - \
                                                                    d_result.loc[exp + f"_{i}", "Bidirectional"]

            d_result /= pos_regions.shape[0]
            d_result["Threshold"] = threshold
            sub_dfs.append(d_result)
        d_result = pd.concat(sub_dfs)
        counts_comparison = distribution.copy()
        for exp in labels:
            counts_comparison["%s" % exp] = np.log1p(distribution["%s_fwd" % exp] + distribution["%s_rev" % exp])

        d_result.to_csv(final_result_file1)
        counts_comparison.to_csv(final_result_file2)

    return final_result_file1, final_result_file2


def ttseq(bed_files, ttseq_csv, all_assays, save_to, global_registry, local_registries):
    """
    Estimate bias toward stable/unstable transcripts

    Parameters
    ----------
    bed_files : dict
        Dictionary of bed files, each bed file contains a wildcard
    ttseq_csv : str
        Path to ttseq file, which must contains these columns: seqname, start, end, synthesis rate,
        decay rate, strand, feature, GENCODE, group
    all_assays : list or tuple
        Assays of interest
    save_to : str
        Path to write outputs
    global_registry : int
        Global registration for outputs
    local_registries : list of str
        Local registration (order) for outputs

    Returns
    -------

    """
    final_result_file1 = os.path.join(save_to,
                                      "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                      local_registry=local_registries[0]))
    final_result_file2 = os.path.join(save_to,
                                      "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                      local_registry=local_registries[1]))
    if os.path.exists(final_result_file1) and os.path.exists(final_result_file2):
        logger.info(f"Final output file {final_result_file1} and {final_result_file2} exist, skip...")
    else:
        ttseq_fixed_len = pd.read_csv(ttseq_csv)
        ttseq_distribution_per_rep = []
        for rep in range(n_samples):
            args = []
            ttseq_bed = BedTool().from_dataframe(
                ttseq_fixed_len.loc[:, ("seqname", "start", "end", "synthesis rate", "decay rate", "strand")])
            # for exp, sample in zip(labels, bed_files):
            for k, sample in bed_files.items():
                cell_line, exp = k.split("_")
                if cell_line.startswith("K562"):
                    args.append((ttseq_bed, sample % (rep + 1), exp))
            with Pool(n_threads) as pool:
                rep_result = pool.starmap(bed_bed_strand_specific_coverage, args)
            ttseq_distribution_per_rep.append(pd.concat(rep_result, axis=1))
        ttseq_distribution = sum(ttseq_distribution_per_rep) / n_samples

        ttseq_counts_comparison = ttseq_distribution.copy()
        ttseq_counts_comparison["group"] = ttseq_fixed_len["group"]

        for exp in all_assays:
            ttseq_counts_comparison["%s" % exp] = np.log1p(
                ttseq_counts_comparison["%s_fwd" % exp] + ttseq_counts_comparison["%s_rev" % exp])

        ttseq_counts_comparison_melted = ttseq_counts_comparison.melt(id_vars=("group",),
                                                                      value_vars=all_assays,
                                                                      value_name="Read counts",
                                                                      var_name="Experiment")

        ttseq_counts_comparison_melted.to_csv(final_result_file1)
        ttseq_counts_comparison.to_csv(final_result_file2)
    return final_result_file1, final_result_file2


def evaluate_ss_control(paired_bams, ref_bed, save_to, global_registry, local_registry):
    """
    Evaluate strand specificity (taking upstream anti-sense transcription into consideration)

    Parameters
    ----------
    paired_bams : list
        List of bam files
    ref_bed : str
        Path to the reference bed file (bidirectional transcripts free zone)
    save_to : str
        Path to write outputs
    global_registry : int
        Global registry for outputs
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
        args = []
        for (ub, sb) in paired_bams:
            args.append((ub, ref_bed))
            args.append((sb, ref_bed))

        with Pool(n_threads) as pool:
            results = pool.starmap(_estimate_strand_specificity, args)
        res = []

        for i, ss in enumerate(results):
            if i % 2 == 0:
                lab = "Unstranded"
            else:
                lab = "Stranded"
            res.append((lab, ss[0]))
        df = pd.DataFrame(res, columns=("Group", "Strand specificity"))
        df.to_csv(final_result_file)
    return final_result_file


def _identify_mispriming_sites(pl_bw, mn_bw, assay, adapter, genomic_fasta_file, splicing_junctions_file, tmp_dir,
                               fold_change=0.2, flankings=10):
    """
    Identify mispriming sites

    Parameters
    ----------
    pl_bw : str
        Signal on the forward strand
    mn_bw : str
        Signal on the reverse strand
    assay : str
        Name of the assay
    adapter : str
        Sequence of the primer used
    genomic_fasta_file : str
        Path to a fasta file, which contains sequences for each chromosome
    splicing_junctions_file : str
        Path to a bed file, which defines splicing junction loci
    tmp_dir :  str
        Temporary files will be written into this path
    fold_change : float
        Max of fold change for a site to be considered as flush end site
    flankings : int
        range of exclusion zone (no other k-sites are allowed in x bps flanking the current candidate site)

    Returns
    -------
    mispriming_rate : float
        Mispriming rate
    mispriming_sites.shape[0] : int
        # of mispriming sites
    ctable.log_oddsratio : float
        Log odds ratio
    ci_l : float
        Lower bound of CI
    ci_u : float
        Upper bound of CI
    assay : str
        Name of the assay
    """
    genomic_fasta = Fasta(genomic_fasta_file)
    splicing_junctions = BedTool(splicing_junctions_file)

    def _strand_atom(bw_obj, strand_direction="+", rpm_threshold=1):
        blunt_k_sites = []
        counts = 0
        for chrom, size in bw_obj.chroms().items():
            # CPM (RPM) normalization
            signal = np.nan_to_num(bw_obj.values(chrom, 0, size))
            signal[signal < 0] *= -1
            chrom_directed_counts = signal.sum()
            counts += chrom_directed_counts
            real_thred = chrom_directed_counts / 1_000_000 * rpm_threshold
            if real_thred < 10:
                real_thred = 10
            if strand_direction == "+":
                # current/prev, the first element will bf FC for loc_1 instead of loc_0
                with np.errstate(divide="ignore", invalid="ignore"):
                    fcs = signal[1:] / signal[:-1]
                    fc_sites = np.where(fcs < fold_change)[0]
                if fc_sites.shape[0] > 0:
                    candidate_k_sites = set(fc_sites)
                else:
                    candidate_k_sites = {}
            elif strand_direction == "-":
                # next/current, the first element will be FC for
                with np.errstate(divide="ignore", invalid="ignore"):
                    fcs = signal[:-1] / signal[1:]
                    fc_sites = np.where(fcs < fold_change)[0]
                fc_sites += 1
                if fc_sites.shape[0] > 0:
                    candidate_k_sites = set(fc_sites)
                else:
                    candidate_k_sites = {}
            # make FCs 0-based coordinates
            for cks in candidate_k_sites:
                searching_range = range(cks - flankings, cks + flankings, 1)
                flag = 1
                for sr in searching_range:
                    if sr != cks:
                        if sr in candidate_k_sites and signal[sr] > 1:
                            flag = 0
                        elif sum(signal[cks - flankings:cks] > real_thred) > 0 or sum(
                                signal[cks + 1:cks + flankings] > real_thred) > 0:
                            flag = 0
                if flag:
                    try:
                        if strand_direction == "+":
                            context = genomic_fasta[chrom][cks + 1:cks + n_adapter + 1].seq
                        else:
                            context = genomic_fasta[chrom][cks - n_adapter - 1:cks].reverse.complement.seq
                        if True:  # context == adapter:
                            blunt_k_sites.append((chrom, cks, cks + 1, context, signal[cks]))
                    except Exception as e:
                        logger.error(e)
                        pass
        blunt_df = pd.DataFrame(blunt_k_sites, columns=("chrom", "start", "end", "context", "reads"))
        blunt_df["strand"] = strand_direction
        return blunt_df, counts

    n_adapter = len(adapter)
    with pyBigWig.open(pl_bw) as pfh:
        plk, plc = _strand_atom(pfh)
    with pyBigWig.open(mn_bw) as mfh:
        mnk, mnc = _strand_atom(mfh, strand_direction="-")
    # removing candidate sites around splicing junctions
    raw_k_sites = pd.concat([plk, mnk])
    k_sites_bed = BedTool.from_dataframe(raw_k_sites).intersect(splicing_junctions, v=True)
    mispriming_sites = k_sites_bed.to_dataframe(names=("chrom", "start", "end", "context", "reads", "strand"))
    mispriming_sites.columns = ("chrom", "start", "end", "context", "reads", "strand")
    mispriming_sites.to_csv(os.path.join(tmp_dir, f"{assay}_blunt_sites.csv"))
    total_signals = plc + mnc
    mispriming_rate = mispriming_sites.loc[mispriming_sites["context"] == adapter, "reads"].sum() / total_signals
    obs = mispriming_sites.loc[mispriming_sites.context == adapter, "reads"].sum()
    total_blunts = mispriming_sites.reads.sum()
    ctable = Table2x2(np.asarray(((obs, total_blunts - obs),
                                  (BG_3MER_RATES[adapter], N_BG_3MERS - BG_3MER_RATES[adapter]))))
    ci_l, ci_u = ctable.log_oddsratio_confint()
    return mispriming_rate, mispriming_sites.shape[0], ctable.log_oddsratio, ci_l, ci_u, assay


def get_mispriming_loci(bw_prefixes, save_to, genome_fasta, splicing_junctions, adapters, global_registry,
                        local_registry, tmp_dir, n_threads=16, adapter_cut=3, max_fold_change=0.2, non_k_zone=10):
    """
    Get mispriming loci

    Parameters
    ----------
    bw_prefixes : dict
        Dict of bigwig files, which contains signal density from the precisely/approximately mapped ends
    save_to : str
        Save output files to
    genome_fasta : str
        Path to a fasta file, which contains sequences for each chromosome
    splicing_junctions : str
        Path to a bed file, which defines splicing junction loci
    adapters : dict
        RT primer info for each assay
    global_registry : int
        Global registry for outputs
    local_registry : str
        Local registration (order) for outputs
    tmp_dir : str
        Temporary files will be written into this path
    n_threads : int
        # of threads this function can create
    adapter_cut : int
        # of bps to be counted for the adapters/RT primers
    max_fold_change : float
        Max of fold change for a site to be considered as flush end site
    non_k_zone : int
        range of exclusion zone (no other k-sites are allowed in x bps flanking the current candidate site)

    Returns
    -------

    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))
    if adapter_cut != 3:
        raise NotImplemented("adapter_cut: only support 3 nts.")
    prev_jobs = set()
    if os.path.exists(final_result_file):
        prev_job_df = pd.read_csv(final_result_file, index_col=0)
        for assay in prev_job_df.Assay:
            prev_jobs.add(assay)
        logger.info(f"Final output file {final_result_file} exists, {', '.join(prev_jobs)} will not be evaluated...")

    args = []
    for assay, bw_prefix in bw_prefixes.items():
        pl_bw = bw_prefix + "_pl.bw"
        mn_bw = bw_prefix + "_mn.bw"
        assert os.path.exists(pl_bw) and os.path.exists(mn_bw)
        if assay not in prev_jobs and assay in adapters:
            args.append((pl_bw, mn_bw, assay, adapters[assay][:adapter_cut],
                         genome_fasta, splicing_junctions, tmp_dir, max_fold_change, non_k_zone))

    if len(args) > 0:
        with Pool(n_threads) as pool:
            results = pool.starmap(_identify_mispriming_sites, args)

        result_df = pd.DataFrame(results,
                                 columns=("Mispriming rate", "Mispriming sites", "LOR", "CI_l", "CI_u", "Assay"))
        if len(prev_jobs) > 0:
            result_df = pd.concat([prev_job_df, result_df])
        result_df.to_csv(final_result_file)
    return final_result_file


def _estimate_strand_specificity(input_bam, reference_bed="../data/refs/compiled/first_exons_els_trimmed.bed",
                                 first_round=True, debug=False):
    """
    Atom function for estimating strand specificity

    Parameters
    ----------
    input_bam : str
        Full path to the bam file to be tested
    reference_bed : str
        Full path to the bed file which contains reference regions to be compared
    debug : bool
        debug mode? (False by default)
    Returns
    -------
    strand_specificity : float
        Strand specificity
    result_df : pd.DataFrame or None
        If debug == True, returns a detailed result
    """
    bed = pd.read_csv(reference_bed, sep="\t")
    bam = pysam.AlignmentFile(input_bam)
    results = []
    pre_result_dict = {"1++": 0, "1--": 0, "2+-": 0, "2-+": 0,
                       "1+-": 0, "1-+": 0, "2++": 0, "2--": 0}
    for _, line in bed.iterrows():
        for hit in bam.fetch(line[0], line[1], line[2]):
            # just in case there are spliced/fusion reads from another gene
            blocks = hit.get_blocks()
            is_partial_in = 0
            for (bs, be) in blocks:
                if bs <= line[2] and be >= line[1]:
                    is_partial_in = 1
            if not is_partial_in:
                continue
            read_id = "1" if hit.is_read1 else "2"
            read_strand = "-" if hit.is_reverse else "+"
            if debug:
                results.append((line[0], line[1], line[2], line[3], line[5], hit.reference_start, hit.reference_end,
                                read_strand, read_id))
            else:
                key = "{0}{1}{2}".format(read_id, read_strand, line[5])
                pre_result_dict[key] += 1
    if debug:
        df = pd.DataFrame(results)
        pre_result_dict = (df[8] + df[7] + df[4]).value_counts().to_dict()
    else:
        df = None
    # 1++,1--,2+-,2-+
    pos_1 = pre_result_dict["1++"] + pre_result_dict["1--"] + pre_result_dict["2+-"] + pre_result_dict["2-+"]
    # 1+-,1-+,2++,2--
    pos_2 = pre_result_dict["1+-"] + pre_result_dict["1-+"] + pre_result_dict["2++"] + pre_result_dict["2--"]
    dominant_layout = max(pos_1, pos_2)
    # if there isn't enough data points to evaluate strand specificity
    # then loosen the constraint for no ELS-like parts
    if first_round and dominant_layout < 30000:
        logger.info(f"There's only {pos_1 + pos_2} reads useful for evaluating SS on {input_bam}, "
                    f"using raw definition for FEMD...")
        return _estimate_strand_specificity(input_bam, reference_bed=reference_bed.replace("_els_trimmed", ""),
                                            first_round=False, debug=False)
    else:
        return dominant_layout / sum(pre_result_dict.values()), df


def evaluate_strand_specificity(bams, ref_bed, save_to, global_registry, local_registry, n_threads=16):
    """
    Evaluate strand specificity (taking upstream anti-sense transcription into consideration)

    Parameters
    ----------
    bams : list
        List of bam files
    ref_bed : str
        Path to the reference bed file (bidirectional transcripts free zone)
    save_to : str
        Path to write outputs
    global_registry : int
        Global registry for outputs
    local_registry : str
        Local registration (order) for outputs
    n_threads : int
        # of threads this function can create

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
        labels = []
        for l, b in bams.items():
            if l.startswith("K562_") or l.startswith("GM12878_"):
                args.append((b, ref_bed))
                labels.append(l)
        with Pool(n_threads) as pool:
            results = pool.starmap(_estimate_strand_specificity, args)
        res = []
        for ss, lab in zip(results, labels):
            cl, assay, biorep, techrep = lab.split("_")
            res.append((assay, biorep, techrep, ss[0]))
        df = pd.DataFrame(res, columns=("Assay", "Bio rep", "Tech rep", "Strand specificity"))
        df.to_csv(final_result_file)
    return final_result_file


def justify_lor_mispriming(selected_assays, save_to, global_registry, local_registry, tmp_dir, adapter_cut=3):
    """
    Analyze k-mer distribution for flush-end sites

    Parameters
    ----------
    selected_assays : list or tuple
        Assays of interest
    save_to : str
        Save output files to
    global_registry : int
        Global registry for outputs
    local_registry : str
        Local registration (order) for outputs
    tmp_dir : str
        Temporary files will be written into this path
    adapter_cut : int
        # of bps to be counted for the adapters/RT primers

    Returns
    -------
    final_result_file : str
        The plot is exported to this file
    """
    final_result_file = os.path.join(save_to,
                                     "{global_registry}_{local_registry}.csv".format(global_registry=global_registry,
                                                                                     local_registry=local_registry))

    if os.path.exists(final_result_file):
        logger.info(f"Final output file {final_result_file} exists, skip...")
    else:
        from itertools import product
        expected_kmers = set()
        for it in product(("A", "T", "C", "G"), repeat=adapter_cut):
            expected_kmers.add("".join(it))

        results = {
            "kmer": [],
            "freq": [],
            "assay": []
        }

        selected_assays_copy = selected_assays.copy()
        selected_assays_copy.extend(["RIPseq", "TGIRTseq"])

        for assay in selected_assays_copy:
            expected_path = os.path.join(tmp_dir, f"{assay}_blunt_sites.csv")
            if os.path.exists(expected_path) and os.path.isfile(expected_path):
                tdf = pd.read_csv(expected_path, index_col=0)
                blunts = tdf.loc[tdf.context.isin(expected_kmers), "context"].value_counts()
                blunts /= blunts.sum()
                for kmer in expected_kmers:
                    if kmer not in blunts.index:
                        blunts[kmer] = 0
                results["kmer"].extend(blunts.index.tolist())
                results["freq"].extend(blunts.values.tolist())
                results["assay"].extend([assay, ] * blunts.shape[0])
            else:
                logger.warning(f"Blunt sites for {assay} cannot be located.")
        pd.DataFrame(results).to_csv(final_result_file)
    return final_result_file


def main(bed_files, bam_files, data_save_to, data_prefix="", **kwargs):
    """
    Factory

    Parameters
    ----------
    bed_files : list or list-like
        bed files
    bam_files : list or list-like
        bam files
    data_save_to : str
        generated data will be exported to this path
    name_mapping : dict
        Name mapping
    data_prefix : str
        All outputs from `generate_data` will have this prefix
    **kwargs : dict
        se_bed : str
            Path to a bed file, which defines super enhancer elements
        te_bed : str
            Path to a bed file, which defines true enhancer elements
        ttseq_csv : str
            Path to a csv file, which defines decay rates estimated from TT-seq
        prime_3_bw_files : dict
            Dict of bigwig files, which contains signal density from the 3' mapped ends
        strands_control_bams : Tuple or list
            Nested tuples/lists of bam files generated from libraries where the strand specificity info is known
            ((unstranded_lib_1, stranded_lib_1), ..., (unstranded_lib_n, stranded_lib_n))
        btfz_bed : str
            Path to a bed file, which defines bidirectional-transcript free zones
        splicing_junctions : str
            Path to a bed file, which defines splicing junction loci
        genome_fasta : str
            Path to a fasta file, which contains sequences for each chromosome
        adapters : dict
            RT primer info for each assay

    Returns
    -------

    """
    se_bed = kwargs.pop("se_bed")
    te_bed = kwargs.pop("te_bed")
    ttseq_csv = kwargs.pop("ttseq_csv")
    prime_3_bw_files = kwargs.pop("prime_3_bw_files")
    strands_control_bams = kwargs.pop("strands_control_bams")
    btfz_bed = kwargs.pop("btfz_bed")
    splicing_junctions = kwargs.pop("splicing_junctions")
    genome_fasta = kwargs.pop("genome_fasta")
    adapters = kwargs.pop("adapters")

    assert se_bed is not None
    assert te_bed is not None
    assert ttseq_csv is not None

    analysis_summaries = {
        "se_coverage": [],
        "se_pairwise": [],
        "te_coverage": [],
        "te_pairwise": [],
        "stable_unstable": [],
        "mispriming": [],
        "strand_specificity": []
    }
    # se coverage
    se_coverage, se_pairwise = bidirectional_coverage_among_ref_regions(bed_files=bed_files,
                                                                        ref_bed=se_bed,
                                                                        save_to=data_save_to,
                                                                        global_registry=global_registry,
                                                                        local_registry=("SEPercentCoverage",
                                                                                        "SECount"))
    analysis_summaries["se_coverage"].append(se_coverage)
    analysis_summaries["se_pairwise"].append(se_pairwise)
    # te coverage
    te_coverage, te_pairwise = bidirectional_coverage_among_ref_regions(bed_files=bed_files,
                                                                        ref_bed=te_bed,
                                                                        save_to=data_save_to,
                                                                        global_registry=global_registry,
                                                                        local_registry=("TEPercentCoverage",
                                                                                        "TECount"),
                                                                        detectable_threshold=(1, 3, 5, 10, 15, 20))
    analysis_summaries["te_coverage"].append(te_coverage)
    analysis_summaries["te_pairwise"].append(te_pairwise)

    # decay rate from TTseq
    ttseq_dat_melted, ttseq_dat_raw = ttseq(bed_files=bed_files, ttseq_csv=ttseq_csv,
                                            all_assays=full_assays, save_to=data_save_to,
                                            global_registry=global_registry,
                                            local_registries=("DecayRatesMelted", "DecayRates"))
    analysis_summaries["stable_unstable"].append(ttseq_dat_melted)
    analysis_summaries["stable_unstable"].append(ttseq_dat_raw)

    # strand specificity (positive & negative control)
    sspn = evaluate_ss_control(paired_bams=strands_control_bams, ref_bed=btfz_bed,
                                save_to=data_save_to, global_registry=global_registry,
                                local_registry="StrandSpecificityControl")
    analysis_summaries["strand_specificity"].append(sspn)

    # strand specificity (real evaluation)
    ss = evaluate_strand_specificity(bams=bam_files, ref_bed=btfz_bed, save_to=data_save_to,
                                        global_registry=global_registry, local_registry="StrandSpecificity")
    analysis_summaries["strand_specificity"].append(ss)

    # mispriming
    adapter_cut = kwargs.get("adapter_cut", 3)
    max_fold_change = kwargs.get("max_fold_change", 0.2)
    non_k_zone = kwargs.get("non_k_zone", 10)

    mispriming_dat = get_mispriming_loci(bw_prefixes=prime_3_bw_files, save_to=data_save_to,
                                            genome_fasta=genome_fasta, splicing_junctions=splicing_junctions,
                                            adapters=adapters, global_registry=global_registry,
                                            local_registry="MisprimingSites", tmp_dir=tmp_dir, adapter_cut=adapter_cut,
                                            max_fold_change=max_fold_change, non_k_zone=non_k_zone)
    analysis_summaries["mispriming"].append(mispriming_dat)

    blunt_sites_dist = justify_lor_mispriming(selected_assays=full_assays, save_to=data_save_to,
                                                global_registry=global_registry, local_registry="KMerDist",
                                                tmp_dir=tmp_dir, adapter_cut=adapter_cut)
    analysis_summaries["mispriming"].append(blunt_sites_dist)

    with open(os.path.join(data_save_to, f"{data_prefix}_summary.json"), "w") as fh:
        json.dump(analysis_summaries, fh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, help="Number of downsamples", default=3)
    parser.add_argument("--data-save-to", required=True, help="Save output data to")
    parser.add_argument("--se-bed", required=False, help="Super enhancers defined in a bed file")
    parser.add_argument("--te-bed", required=False, help="True enhancers defined in a bed file")
    parser.add_argument("--ttseq-csv", required=False, help="Decay rate from TTseq in a csv file")
    parser.add_argument("--genomic-fasta", required=False,
                        help="Reference sequences")
    parser.add_argument("--btfz-bed", required=False,
                        help="Path to the bed file which defines Bidirectional-transcript free zones")
    parser.add_argument("--data-prefix", required=False, default="sensitivity_reads", help="Prefix for all outputs")
    parser.add_argument("--global-registry", type=int, help="Global registry", default=1)
    parser.add_argument("--config-file", default="config.txt", help="Configuration file")
    parser.add_argument("--tmp-dir", default="./tmp", help="All temporary files will be written into this folder")

    args = parser.parse_args()

    assert os.path.exists(args.config_file)

    if not os.path.exists(args.tmp_dir):
        os.mkdir(args.tmp_dir)

    logger.addHandler(logging.FileHandler(os.path.join(os.getcwd(), "AssaySensitivityReads.log")))

    cfg = ConfigParser(interpolation=None)
    cfg.optionxform = str
    cfg.read(args.config_file)
    collapsed_assays = cfg.get("assays", "plot_order_simplified").split("|")
    end_groups = cfg["dataset_precise_ends"]
    cor_bws = cfg["corroborative_bws"]

    global tmp_dir, n_threads, n_samples, global_registry, full_assays, highlight_assays, layouts
    tmp_dir = args.tmp_dir
    pybedtools.set_tempdir(tmp_dir)
    n_threads = int(cfg.get("global", "n_threads"))
    n_samples = int(cfg.get("assays", "n_downsamples"))
    global_registry = args.global_registry
    full_assays = cfg.get("assays", "plot_order_full").split("|")
    assay_offical_names = cfg.get("assays", "assay_full_names").split("|")
    highlight_assays = cfg.get("assays", "plot_order_simplified").split("|")

    layouts = dict()
    for k, v in cfg["dataset_layouts"].items():
        layouts[k] = v
    official_name_map = dict()
    plot_order = cfg.get("assays", "plot_order_full").split("|")
    plot_order_simplified = cfg.get("assays", "plot_order_simplified").split("|")
    
    for k, v in enumerate(plot_order):
        official_name_map[v] = assay_offical_names[k]

    adapters = dict()
    for k, v in cfg["RT_primer"].items():
        adapters[k] = v

    beds = dict()
    ubams = dict()
    p3_bws = dict()
    known_ss_libraries = []

    import socket

    server1_home_dir = "/fs/cbsuhy01/storage/ly349/" if socket.gethostname().find(
        "cbsuhy01") == -1 else "/local/storage/ly349/"
    server2_home_dir = "/fs/cbsuhy02/storage/ly349/" if socket.gethostname().find(
        "cbsuhy02") == -1 else "/local/storage/ly349/"
    bioq_dir = os.path.join(server2_home_dir, "BioQueue/workspace")

    from .utils import load_bioq_datasets, load_bioq_dataset

    p3_bws = load_bioq_datasets("k562_3p_bw_prefixs", bioq_dir, cfg_file=args.config_file)
    known_ss_libraries = (
        # unstranded, stranded
        (load_bioq_dataset(cfg.get("auxiliary_alignment_per_rep", "K562_uRNAseq_1_1"), bioq_dir),
         load_bioq_dataset(cfg.get("unique_alignments_per_rep", "K562_totalRNAseq_1_1"), bioq_dir)),
        (load_bioq_dataset(cfg.get("auxiliary_alignment_per_rep", "K562_uRNAseq_2_1"), bioq_dir),
         load_bioq_dataset(cfg.get("unique_alignments_per_rep", "K562_totalRNAseq_2_1"), bioq_dir)),
        (load_bioq_dataset(cfg.get("auxiliary_alignment_per_rep", "K562_uRNAseq_3_1"), bioq_dir),
         load_bioq_dataset(cfg.get("unique_alignments_per_rep", "K562_totalRNAseq_3_1"), bioq_dir))
    )
    beds = load_bioq_datasets("downsampled_clean_bed", bioq_dir, cfg_file=args.config_file)
    ubams = load_bioq_datasets("unique_alignments_per_rep", bioq_dir, cfg_file=args.config_file)

    main(bed_files=beds, data_save_to=args.data_save_to, 
         se_bed=args.se_bed, te_bed=args.te_bed, ttseq_csv=args.ttseq_csv,
         prime_3_bw_files=p3_bws, strands_control_bams=known_ss_libraries, btfz_bed=args.btfz_bed,
         bam_files=ubams, splicing_junctions=cfg.get("references", "gencode_splicing_junctions"),
         genome_fasta=args.genomic_fasta, adapters=adapters, data_prefix=args.data_prefix)
