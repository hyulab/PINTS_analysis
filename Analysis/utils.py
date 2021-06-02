#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 4/20/20


def load_bioq_dataset(file_path, bioq_path,
                      server="http://localhost:987%s/ui/api?model=Job&obj_id=%s"):
    import os
    import requests
    from requests.adapters import HTTPAdapter
    from requests.packages.urllib3.util.retry import Retry

    retry_strategy = Retry(
        total=10,
        status_forcelist=[429, 500, 502, 503, 504],
        method_whitelist=["HEAD", "GET", "OPTIONS"],
        backoff_factor=3
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    with requests.Session() as http:
        http.mount("https://", adapter)
        http.mount("http://", adapter)

        server_id, user_id, job_id, file_name = file_path.split(";")
        try:
            req = http.get(url=server % (server_id, job_id))
            rj = req.json()
            assert rj["status"] == 1

            full_path = os.path.join(bioq_path, user_id, rj["info"])
            full_path2 = os.path.join("/fs/cbsuhyfs1/storage/bioqueue/workspace", user_id, rj["info"])
            if os.path.exists(full_path) and os.path.isdir(full_path):
                full_path = os.path.join(full_path, file_name)
            elif os.path.exists(full_path2) and os.path.isdir(full_path2):
                full_path = os.path.join(full_path2, file_name)
            else:
                raise Exception(f"Cannot locate result for {job_id} ({rj['info']})")
        except Exception as e:
            assert 0 == 1, str(e) + server % (server_id, job_id)
        return full_path


def load_bioq_datasets(cfg_section, bioq_path, cfg_file="config.txt"):
    from configparser import ConfigParser
    cfg = ConfigParser(interpolation=None)
    cfg.optionxform = str
    cfg.read(cfg_file)
    res = dict()
    for k, v in cfg[cfg_section].items():
        if not k.startswith("#"):
            res[k] = load_bioq_dataset(v, bioq_path)
    return res


def bin_scores(regions, score_bw, bins=100, lower_bound=2.5, upper_bound=97.5, n_boot=10000):
    """

    Parameters
    ----------
    regions : BedTool
        BedTool object for regions to be analyzed
    score_bw : pyBigWig.bigWigFile
        bigWigFile object for scores
    bins : int
        Number of bins, by default: 100.
    lower_bound : float
        Lower bound of statistic from bootstrap to be reported, by default 2.5
    upper_bound : float
        Upper bound of statistic from bootstrap to be reported, by default 97.5
    n_boot : int
        Number of bootstraps to run, by default 10000

    Returns
    -------
    score_mat : numpy.ndarray
        Score matrix for each region (binned)
    means : numpy.ndarray
        Mean of scores among each column
    (lower, upper) : (numpy.ndarray, numpy.ndarray)
        Lower and upper bound of stats
    """
    import numpy as np
    from scipy.stats import binned_statistic
    from seaborn.algorithms import bootstrap
    score_per_window = []
    for row in regions:
        assert row.end - row.start > bins, "Number of bins must be smaller than the region (%s)" % str(row).strip()
        values = np.nan_to_num(score_bw.values(row.chrom, row.start, row.end))
        y, x, z = binned_statistic(np.arange(values.shape[0]), values, statistic="mean", bins=bins)
        score_per_window.append(y)
    score_mat = np.vstack(score_per_window)
    means = np.mean(score_mat, axis=0)
    # stds = np.std(score_mat, axis=0)
    boots = bootstrap(score_mat, axis=0, n_boot=n_boot)
    lower = np.percentile(boots, lower_bound, axis=0)
    upper = np.percentile(boots, upper_bound, axis=0)
    return score_mat, means, (lower, upper)


def bed_bed_strand_specific_coverage(bed_needle, bed_haystack, exp):
    """
    Strand specific coverage computation

    Parameters
    ----------
    bed_needle: pybedtools.BedTool
        Needle/query bed
    bed_haystack: pybedtools.BedTool
        Haystack/dest bed
    exp: str
        name of the experiment

    Returns
    -------
    res: pd.DataFrame
        coverage for needles (forward and reverse separately)
    """
    import pandas as pd
    pl = bed_needle.intersect(bed_haystack, c=True, s=True)  # signals from pl strand
    mn = bed_needle.intersect(bed_haystack, c=True, S=True)  # signals from mn strand
    pl_df = pl.to_dataframe(names=("seqname", "start", "end", "evidence", "score", "strand", "fwd"))
    mn_df = mn.to_dataframe(names=("seqname", "start", "end", "evidence", "score", "strand", "rev"))
    res = pd.concat([pl_df["fwd"], mn_df["rev"]], axis=1)
    res.columns = ["%s_fwd" % exp, "%s_rev" % exp]
    return res


def evaluation_core_element(pred, threshold, true_set="../data/criteria/v3_20200614/GROcap_bidirectional_positives.bed",
                            false_set="../data/criteria/v3_20200614/GROcap_bidirectional_TN_enhancers.bed"):
    import numpy as np
    from pybedtools import BedTool
    # predictions = predicted positive
    predictions = pred.copy()
    predictions.columns = ("chrom", "start", "end", "score")
    predictions = predictions.loc[predictions["score"] >= threshold, :]
    pred_obj = BedTool.from_dataframe(predictions)
    condition_positive_bed = BedTool(true_set)
    condition_negative_bed = BedTool(false_set)
    TP = len(condition_positive_bed.intersect(pred_obj.intersect(condition_negative_bed, v=True), u=True))
    FP = len(condition_negative_bed.intersect(pred_obj, u=True))
    FN = len(condition_positive_bed.intersect(pred_obj, v=True))
    TN = len(condition_negative_bed.intersect(pred_obj, v=True))
    if TP > 0:
        recall_sensitivity = TP / (TP + FN)
        specificity = TN / (FP + TN)
        precision = TP / (TP + FP)
    else:
        recall_sensitivity = 0
        specificity = 1
        precision = 0
    prob_fx_true = (TP + FP) / (TP + FP + FN + TN)
    if prob_fx_true > 0:
        f_pu = recall_sensitivity ** 2 / prob_fx_true
    else:
        f_pu = np.nan
    if precision + recall_sensitivity > 0:
        f = 2 * (precision * recall_sensitivity) / (precision + recall_sensitivity)
    else:
        f = np.nan
    return recall_sensitivity, specificity, precision, f, f_pu


def cohens_d(x, y):
    """
    Calculate Cohen's d
    """
    import numpy as np
    m1 = np.mean(x)
    m2 = np.mean(y)
    v1 = np.std(x, ddof=1) ** 2
    v2 = np.std(y, ddof=1) ** 2
    n1 = len(x)
    n2 = len(y)
    s = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    d = (m1 - m2) / s
    return d


def run_command(cmd):
    """
    Run command
    :param cmd:
    :return: (stdout, stderr, return_code)
    """
    from subprocess import Popen, PIPE
    p = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate()
    stderr = stderr.decode("utf-8")
    stdout = stdout.decode("utf-8")
    return stdout, stderr, p.returncode


def run_command_must_success(cmd):
    stdout, stderr, rc = run_command(cmd)
    if rc != 0:
        raise RuntimeError(stderr)
    else:
        return stdout, stderr, rc


def read_bed(file, remove_mt=True):
    import pandas as pd
    try:
        df = pd.read_csv(file, sep="\t", header=None, comment="#")
    except pd.errors.ParserError:
        df = pd.read_csv(file, sep="\t", header=None, comment="#", skiprows=1)
    if remove_mt:
        df = df.loc[df[0] != "chrM", :]
    return df


def midpoint_generator(bed_regions):
    from pybedtools.featurefuncs import midpoint
    try:
        for region in bed_regions:
            yield midpoint(region)
    except Exception as e:
        print(e)


def contrast_regions(regions, labels, score_bw_dict, region_extension, chromosome_size, n_bins=100):
    """
    Contrast regions

    Parameters
    ----------
    regions : list of str
        List of paths to region bed files
    labels : list of str
        List of corresponding labels for each region bed file
    score_bw_dict : dict
        keys : name of the bw in the form of `CL_MARK`
        values : paths to the bw files
    region_extension : int
        Number of bps to be extended from the mid point
    chromosome_size : str
        Path to a tab file describing length of each chromosome
    n_bins : int
        Number of bins to generate in each region (`mid`-`region_extension`, `mid`+`region_extension`)

    Returns
    -------
    result_df : pd.DataFrame
        Mean :
        Upper :
        Lower :
        Marker : Marker name (bw)
        Set : Region name (bed)
    """
    from pybedtools import BedTool
    import pybedtools
    import pyBigWig
    import pandas as pd
    pybedtools.set_tempdir(".")

    ref_beds = []
    for region in regions:
        ref_beds.append(BedTool(midpoint_generator(BedTool(region))).slop(l=region_extension, r=region_extension,
                                                                          g=chromosome_size))

    result_df = None
    results = []
    for k, v in score_bw_dict.items():
        with pyBigWig.open(v) as bw_obj:
            for label, region in zip(labels, ref_beds):
                mat, m, (u, l) = bin_scores(regions=region, score_bw=bw_obj, bins=n_bins)
                marker_label = k.split("_")[1]
                results.append(pd.DataFrame({"mean": m,
                                             "upper": u,
                                             "lower": l,
                                             "marker": marker_label,
                                             "set": label}))

    if len(results) > 0:
        result_df = pd.concat(results)
    return result_df


LOGGING_CFG = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {
            "class": "logging.Formatter",
            "format": "%(name)s - %(asctime)s - %(levelname)s: %(message)s",
            "datefmt": "%d-%b-%y %H:%M:%S",
        },
    },
    "handlers": {
        "file": {
            "level": "INFO",
            "class": "logging.FileHandler",
            "formatter": "default",
            "filename": "PINTS_analysis.log",
            "mode": "a",
            "encoding": "utf-8"
        },
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "default",
            "level": "INFO",
        }
    },
    "root": {
        "handlers": ["file", "console"],
        "level": "INFO"
    }
}


class DependencyError(Exception):
    pass
