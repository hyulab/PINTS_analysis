#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 8/16/20
import logging.config
import os
import sys
import pickle
import numpy as np
import pandas as pd
import pybedtools
from abc import abstractmethod
from pybedtools import BedTool
from decimal import Decimal, getcontext
from .utils import LOGGING_CFG, DependencyError, run_command

_x = LOGGING_CFG.copy()
_x["handlers"]["file"]["filename"] = "ROC.log"
logging.config.dictConfig(_x)
logger = logging.getLogger("roc")


class ROCProfiler(object):
    def __init__(self, output_prefix, positive_set, negative_set, cache_dir=".",
                 sample_points=101, pairing_distance=300):
        """
        Init

        Parameters
        ----------
        output_prefix : str
            Prefix of the output file. All outputs will share this prefix. For example, `output_prefix`.bed
        cache_dir : str
            Path to write tmp file. By default, current folder
        sample_points : int

        """
        output_prefix = output_prefix if output_prefix[-1] != "_" else output_prefix[:-1]
        self.output_prefix = output_prefix
        self._cache_dir = cache_dir
        self._sample_points = sample_points
        self._positive_set = positive_set
        self._negative_set = negative_set
        self._pairing_distance = pairing_distance
        self.roc_ready_peak = None
        self.recalls = np.ones(sample_points) * -1
        self.specificities = np.ones(sample_points) * -1
        self.precisions = np.ones(sample_points) * -1
        self.fscores = np.ones(sample_points) * -1
        self.fscores_alt = np.ones(sample_points) * -1
        self.cutoffs = np.empty(sample_points, dtype=object)
        pybedtools.set_tempdir(self._cache_dir)

    def _data_sanity_check(self):
        # empty_probe = np.zeros(self._sample_points)
        if self.roc_ready_peak is not None:
            n = len(self.roc_ready_peak)
            self.recalls = np.zeros(n)
            self.specificities = np.zeros(n)
            self.precisions = np.zeros(n)
            self.fscores = np.zeros(n)
            self.fscores_alt = np.zeros(n)
        else:
            raise NotImplemented

    @staticmethod
    def naive_pairing(merged_peaks, pairing_distance=300):
        """
        Generate naive bidirectional pairs

        Parameters
        ----------
        merged_peaks : pd.DataFrame
            Bed file in pd.DataFrame
        pairing_distance : int
            Max distance between two peaks. By default, 300.
        Returns
        -------
        paired_peaks : pd.DataFrame
            Paired peaks
        """
        merged_peak_bed = BedTool.from_dataframe(merged_peaks)
        pl_peaks = merged_peak_bed.filter(lambda x: x.strand == "+").saveas()
        mn_peaks = merged_peak_bed.filter(lambda x: x.strand == "-").saveas()
        pl_clear_peaks = pl_peaks.intersect(mn_peaks, v=True).sort()
        mn_clear_peaks = mn_peaks.intersect(pl_peaks, v=True).sort()
        clear_peaks = BedTool.cat(*[pl_clear_peaks, mn_clear_peaks], postmerge=False).sort()
        # Sm: Only report hits in B that overlap A on the opposite strand.
        # By default, overlaps are reported without respect to strand.
        paired_peaks = clear_peaks.window(clear_peaks, w=pairing_distance, Sm=True)
        if len(paired_peaks) > 0:
            return paired_peaks.to_dataframe()
        else:
            return None

    @abstractmethod
    def generate_scored_pairs(self, inputs, **kwargs):
        pass

    def evaluation_core_element(self, pred, threshold):
        predictions = pred.copy()
        if predictions.shape[0] == 0:
            return np.nan, np.nan, np.nan, np.nan, np.nan
        predictions.columns = ("chrom", "start", "end", "score")
        predictions = predictions.loc[predictions["score"] >= threshold, :]
        pred_obj = BedTool.from_dataframe(predictions)
        condition_positive_bed = BedTool(self._positive_set)
        condition_negative_bed = BedTool(self._negative_set)
        TP = len(condition_positive_bed.intersect(pred_obj.intersect(condition_negative_bed, v=True), u=True, F=0.5))
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

    def evaluation_core_strict(self, pred, threshold):
        predictions = pred.copy()
        predictions.columns = ("chrom", "start", "end", "score")
        predictions["predict label"] = 1
        predictions.loc[predictions["score"] <= threshold, "predict label"] = 0
        predictions["true label"] = np.nan
        predictions["k"] = predictions["chrom"] + ":" + predictions["start"].map(str) + "-" + predictions["end"].map(
            str)
        condition_positive_bed = BedTool(self._positive_set)
        condition_negative_bed = BedTool(self._negative_set)
        # time consuming
        pred_obj = BedTool.from_dataframe(predictions)
        pred_cond_pos = pred_obj.intersect(condition_positive_bed, u=True).to_dataframe(
            names=("chrom", "start", "end", "score", "predict label", "true label", "k"))["k"]
        pred_cond_neg = pred_obj.intersect(condition_negative_bed, u=True).to_dataframe(
            names=("chrom", "start", "end", "score", "predict label", "true label", "k"))["k"]
        predictions.loc[predictions["k"].isin(pred_cond_pos), "true label"] = 1
        predictions.loc[predictions["k"].isin(pred_cond_neg), "true label"] = 0
        pred_result = predictions.dropna().copy()
        pred_result.loc[:, "true label"] = pred_result["true label"].astype(int)

        confusion_mat = pred_result.loc[:, ("predict label", "true label", "k")].groupby(
            by=["predict label", "true label"]).count()
        # predict, true
        try:
            tp = confusion_mat.loc[1, 1]["k"]
        except:
            tp = 0
        try:
            fp = confusion_mat.loc[1, 0]["k"]
        except:
            fp = 0
        try:
            tn = confusion_mat.loc[0, 0]["k"]
        except:
            tn = 0
        try:
            fn = confusion_mat.loc[0, 1]["k"]
        except:
            fn = 0

        prob_fx_true = (len(pred_cond_pos) + len(pred_cond_neg)) / (
                len(condition_positive_bed) + len(condition_negative_bed))
        recall_sensitivity = tp / max(tp + fn, 1)
        specificity = tn / max(fp + tn, 1)
        precision = tp / max(tp + fp, 1)
        f_pu = recall_sensitivity ** 2 / prob_fx_true
        if precision + recall_sensitivity > 0:
            f = 2 * (precision * recall_sensitivity) / (precision + recall_sensitivity)
        else:
            f = np.nan
        return recall_sensitivity, specificity, precision, f, f_pu

    def _roc_formatter(self):
        specificity_order = np.argsort(self.specificities)
        specificity_sorted = self.specificities[specificity_order]
        precision_sorted = self.precisions[specificity_order]
        recall_sorted = self.recalls[specificity_order]
        cutoff_sorted = self.cutoffs[specificity_order]
        fscore_sorted = self.fscores[specificity_order]
        fscore_alt_sorted = self.fscores_alt[specificity_order]
        if specificity_sorted[-1] != 1:
            specificity_sorted = np.append(specificity_sorted, 1)
            recall_sorted = np.append(recall_sorted, 0)
            precision_sorted = np.append(precision_sorted, np.nan)
            fscore_sorted = np.append(fscore_sorted, np.nan)
            fscore_alt_sorted = np.append(fscore_alt_sorted, np.nan)
            cutoff_sorted = np.append(cutoff_sorted, "")
        if specificity_sorted[0] != 0:
            specificity_sorted = np.insert(specificity_sorted, 0, 0)
            recall_sorted = np.insert(recall_sorted, 0, recall_sorted[0])
            precision_sorted = np.insert(precision_sorted, 0, np.nan)
            fscore_sorted = np.insert(fscore_sorted, 0, np.nan)
            fscore_alt_sorted = np.insert(fscore_alt_sorted, 0, np.nan)
            cutoff_sorted = np.insert(cutoff_sorted, 0, "")

        result = {
            "recall": recall_sorted,
            "specificity": specificity_sorted,
            "precision": precision_sorted,
            "fscore": fscore_sorted,
            "fscore_alt": fscore_alt_sorted,
            "cutoff": cutoff_sorted
        }
        with open(self.output_prefix + "_ROC.dat", "wb") as fh:
            pickle.dump(result, fh)
        return result

    def _unbalanced_range_for_pqs(self, local_min=0):
        n_small_pvals = int(self._sample_points * 0.25)
        n_exp_pvals = int(self._sample_points * 0.25)
        if n_exp_pvals > 51:
            n_exp_pvals = 51
        small_pvals = np.linspace(local_min, 0.1, n_small_pvals).tolist()
        exp_pvals = [10 ** i for i in range(-1 * n_exp_pvals, -1)]
        # make sure enough precision
        getcontext().prec = n_exp_pvals + 1
        actual_values = len(set(small_pvals + exp_pvals))
        n_general_pvals = self._sample_points - actual_values + 1
        general_pvals = np.linspace(0.1, 1, n_general_pvals).tolist()
        return set(small_pvals + exp_pvals + general_pvals)

    def roc(self, definition="element"):
        """
        ROC core

        Parameters
        ----------
        definition : str
            `element`, if the ROC is for evaluating tool's ability in terms of discovering TREs
            `canonical`, if the ROC is for the score

        Returns
        -------
        result : dict
            recall : np.array
                Sorted recalls
            specificity : np.array
                Sorted specificities
            precisions : np.array
                Sorted precisions
            fscores : np.array
                Sorted F-scores
            fscores_alt : np.array
                Sorted F-scores (alternative, optimized for Positive and Unlabeled learning)
        """
        if self.roc_ready_peak is None or not os.path.exists(self.roc_ready_peak) or not os.path.isfile(
                self.roc_ready_peak):
            raise DependencyError("Scored divergent pairs not found")
        df = pd.read_csv(self.roc_ready_peak, sep="\t", header=None)
        assert df.shape[1] == 4, "File must have 4 columns, chr\tstart\tend\tscore"
        score_max = df[3].max()
        score_min = df[3].min()

        pos_bed = BedTool(self._positive_set)
        neg_bed = BedTool(self._negative_set)
        scope_bed = BedTool.cat(*[pos_bed, neg_bed], postmerge=False).sort()
        # scope_bed = pos_bed.cat(neg_bed).sort()
        prediction_bedobj = BedTool.from_dataframe(df)
        # prediction_bedobj.merge(c=(4,), o=("max",))
        prediction_bedobj = prediction_bedobj.intersect(scope_bed, u=True)
        df = prediction_bedobj.to_dataframe(names=("chr", "start", "end", "score"))

        for i, threshold in enumerate(np.linspace(score_min, score_max, self._sample_points)):
            if definition == "element":
                res = self.evaluation_core_element(df, threshold)
            else:
                res = self.evaluation_core_strict(df, threshold)
            self.recalls[i] = res[0]
            self.specificities[i] = res[1]
            self.precisions[i] = res[2]
            self.fscores[i] = res[3]
            self.fscores_alt[i] = res[4]
            self.cutoffs[i] = threshold

        return self._roc_formatter()

    def _pq_evaluation_core_element(self, pred, threshold, F_restraint=0.5):
        predictions = pred.copy()
        pred_obj = BedTool.from_dataframe(predictions)
        condition_positive_bed = BedTool(self._positive_set)
        condition_negative_bed = BedTool(self._negative_set)
        TP = len(
            condition_positive_bed.intersect(pred_obj.intersect(condition_negative_bed, v=True), u=True, F=F_restraint))
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

    def _pq_evaluation_core_strict(self, pred, threshold):
        raise NotImplemented("not implemented")

    def _pq_roc(self, definition="element", F_restraint=0.5):
        if self.roc_ready_peak is None:
            raise DependencyError

        pos_bed = BedTool(self._positive_set)
        neg_bed = BedTool(self._negative_set)
        scope_bed = BedTool.cat(*[pos_bed, neg_bed], postmerge=False).sort()
        if self.recalls.shape[0] != len(self.roc_ready_peak):
            logger.warning("# of peak calls doesn't match the settings when initiating ROCProfiler")
            self._data_sanity_check()
        logger.info(f"{len(self.roc_ready_peak)} sample points will be drawn")
        i = 0

        for score, bidirectional_peak_file in self.roc_ready_peak.items():
            try:
                df = pd.read_csv(bidirectional_peak_file, sep="\t", header=None)
                df = df.loc[:, (0, 1, 2)]
                prediction_bedobj = BedTool.from_dataframe(df)
                prediction_bedobj = prediction_bedobj.intersect(scope_bed, u=True)

                df = prediction_bedobj.to_dataframe(names=("chr", "start", "end"))

                if definition == "element":
                    res = self._pq_evaluation_core_element(df, score, F_restraint=F_restraint)
                else:
                    res = self._pqevaluation_core_strict(df, score)
                self.recalls[i] = res[0]
                self.specificities[i] = res[1]
                self.precisions[i] = res[2]
                self.fscores[i] = res[3]
                self.fscores_alt[i] = res[4]
                self.cutoffs[i] = score
            except Exception as e:
                self.recalls[i] = np.nan
                self.specificities[i] = np.nan
                self.precisions[i] = np.nan
                self.fscores[i] = np.nan
                self.fscores_alt[i] = np.nan
                self.cutoffs[i] = score
                logger.error(f"{score} {bidirectional_peak_file}")
                logger.error(e)
            i += 1

        return self._roc_formatter()


class ComplexScoreProfiler(ROCProfiler):
    def roc(self, definition="element"):
        return self._pq_roc(definition=definition)

    def generate_scored_pairs(self, inputs, **kwargs):
        from file_parser_mem.bed import parse_bed
        b = BedTool.from_dataframe(parse_bed(inputs))
        tidy_file = inputs + ".roctidy.bed"
        b.saveas(tidy_file)
        self.roc_ready_peak = {0: tidy_file}


class ProfilerMACS2(ROCProfiler):
    def generate_scored_pairs(self, inputs, **kwargs):
        self.roc_ready_peak = inputs

    def roc(self, definition="element"):
        """
        ROC core

        Parameters
        ----------
        definition : str
            `element`, if the ROC is for evaluating tool's ability in terms of discovering TREs
            `canonical`, if the ROC is for the score

        Returns
        -------
        result : dict
            recall : np.array
                Sorted recalls
            specificity : np.array
                Sorted specificities
            precisions : np.array
                Sorted precisions
            fscores : np.array
                Sorted F-scores
            fscores_alt : np.array
                Sorted F-scores (alternative, optimized for Positive and Unlabeled learning)
        """
        if self.roc_ready_peak is None or not os.path.exists(self.roc_ready_peak) or not os.path.isfile(
                self.roc_ready_peak):
            raise DependencyError("Scored divergent pairs not found")
        df = pd.read_csv(self.roc_ready_peak, sep="\t", header=None)
        assert df.shape[1] == 4, f"File must have 4 columns, chr\tstart\tend\tscore ({self.roc_ready_peak}, {df.shape})"
        score_min = df[3].min()

        pos_bed = BedTool(self._positive_set)
        neg_bed = BedTool(self._negative_set)
        scope_bed = BedTool.cat(*[pos_bed, neg_bed], postmerge=False).sort()
        prediction_bedobj = BedTool.from_dataframe(df)
        # prediction_bedobj.merge(c=(4,), o=("max",))
        prediction_bedobj = prediction_bedobj.intersect(scope_bed, u=True)
        df = prediction_bedobj.to_dataframe(names=("chr", "start", "end", "score"))
        exp_pvals = np.linspace(1, 1000, 30).tolist()
        general_pvals = (-1 * np.log10(np.linspace(0.105, 1, self._sample_points - len(exp_pvals)))).tolist()
        c = set(exp_pvals + general_pvals)
        for i, threshold in enumerate(c):
            if definition == "element":
                res = self.evaluation_core_element(df, threshold)
            else:
                res = self.evaluation_core_strict(df, threshold)
            self.recalls[i] = res[0]
            self.specificities[i] = res[1]
            self.precisions[i] = res[2]
            self.fscores[i] = res[3]
            self.fscores_alt[i] = res[4]
            self.cutoffs[i] = threshold

        return self._roc_formatter()


class ProfilerHOMERTSS(ROCProfiler):
    def __init__(self, output_prefix, positive_set, negative_set, cache_dir=".",
                 sample_points=101, pairing_distance=300):
        super(ProfilerHOMERTSS, self).__init__(output_prefix=output_prefix, positive_set=positive_set,
                                                negative_set=negative_set, cache_dir=cache_dir,
                                                sample_points=sample_points, pairing_distance=pairing_distance)
        self.bid_pairs_with_score = dict()
    def roc(self, definition="element"):
        return self._pq_roc(definition=definition)

    def generate_scored_pairs(self, inputs, **kwargs):
        """
        Generate scored pairs for HOMER (TSS mode)

        This function will pair all bidirectional peaks using kwargs["pairing_distance"] as the upper bound
        For each pair, the min score from the two peaks will be served as the score for the pair

        Parameters
        ----------
        inputs : str
            reformatted bed file for TSS peaks identified by HOMER

        Returns
        -------

        """
        # tss_file = inputs
        from pints.io_engine import index_bed_file
        pl_file = inputs + "pl.bed.gz"
        mn_file = inputs + "mn.bed.gz"
        if os.path.exists(pl_file) and os.path.exists(mn_file):
            pass
        elif os.path.exists(inputs):
            peak_bed = pd.read_csv(inputs, sep="\t", names=("chrom", "start", "end", "name",
                                                            "score", "strand"))
            pl_peaks = peak_bed.loc[peak_bed.strand == "+", :]
            pl_peaks.to_csv(f"{inputs}pl.bed", sep="\t", header=False, index=False)
            index_bed_file(f"{inputs}pl.bed")
            mn_peaks = peak_bed.loc[peak_bed.strand == "-", :]
            mn_peaks.to_csv(f"{inputs}mn.bed", sep="\t", header=False, index=False)
            index_bed_file(f"{inputs}mn.bed")
        else:
            raise DependencyError("Peak calls for ROC not found.")
        max_bidirectional_distance = self._pairing_distance
        pe_extended_bed = pd.read_csv(self._positive_set, sep="\t", header=None)
        ne_extended_bed = pd.read_csv(self._negative_set, sep="\t", header=None)
        known_regions = pd.concat([pe_extended_bed, ne_extended_bed])
        known_regions[1] -= self._pairing_distance
        known_regions[2] += self._pairing_distance
        known_regions = BedTool.from_dataframe(known_regions)
        logger.info("Selecting peaks only in known regions (SPOIKR)")
        pl_df = pd.read_csv(pl_file, sep="\t", header=None)
        all_pl_peaks_df = BedTool.from_dataframe(pl_df).intersect(known_regions, u=True).to_dataframe(names=(
            "chrom", "start", "end", "name", "score", "strand"))
        del pl_df
        logger.info("SPOIKR pl done")
        mn_df = pd.read_csv(mn_file, sep="\t", header=None)
        all_mn_peaks_df = BedTool.from_dataframe(mn_df).intersect(known_regions, u=True).to_dataframe(names=(
            "chrom", "start", "end", "name", "score", "strand"))
        del mn_df
        logger.info("SPOIKR mn done")
        all_pl_peaks_df.columns = np.arange(all_pl_peaks_df.shape[1])
        all_mn_peaks_df.columns = np.arange(all_mn_peaks_df.shape[1])

        ls = np.log1p(all_pl_peaks_df[4])
        scope = np.linspace(ls.min(), ls.max(), 100)
        for transformed_score in scope:
            score = np.expm1(transformed_score)
            # real_k = score
            sig_pl_df = all_pl_peaks_df.loc[all_pl_peaks_df[4] >= score, :]
            sig_mn_df = all_mn_peaks_df.loc[all_mn_peaks_df[4] >= score, :]
            if sig_pl_df.shape[0] == 0 or sig_mn_df.shape[0] == 0:
                logger.warning(f"No peak passed this threshold: {score}")
                continue

            fn_pl_sig = os.path.join(self._cache_dir, "sig_pl.bed")
            fn_mn_sig = os.path.join(self._cache_dir, "sig_mn.bed")
            sig_pl_df.to_csv(fn_pl_sig, sep="\t", header=None, index=False)
            sig_mn_df.to_csv(fn_mn_sig, sep="\t", header=None, index=False)

            pl_merged_peaks = pd.concat([sig_pl_df, all_mn_peaks_df])
            mn_merged_peaks = pd.concat([sig_mn_df, all_pl_peaks_df])
            pl_np_df = self.naive_pairing(pl_merged_peaks, max_bidirectional_distance)
            mn_np_df = self.naive_pairing(mn_merged_peaks, max_bidirectional_distance)
            merged_bidirectional_file = os.path.join(self._cache_dir,
                                                     "HOMER_TSS_%s_bidirectional_peaks.bed" % transformed_score)
            if pl_np_df is not None and mn_np_df is not None:
                bid_pl = BedTool.from_dataframe(pl_np_df)
                bid_mn = BedTool.from_dataframe(mn_np_df)
                bid = BedTool.cat(*[bid_pl, bid_mn])
                bid.saveas(merged_bidirectional_file)

                self.bid_pairs_with_score[score] = merged_bidirectional_file
            elif pl_np_df is not None:
                bid_pl = BedTool.from_dataframe(pl_np_df)
                bid_pl.saveas(merged_bidirectional_file)
                self.bid_pairs_with_score[score] = merged_bidirectional_file
            elif mn_np_df is not None:
                bid_mn = BedTool.from_dataframe(mn_np_df)
                bid_mn.saveas(merged_bidirectional_file)
                self.bid_pairs_with_score[score] = merged_bidirectional_file
            else:
                logger.warning(
                    f"Candidate peaks with score higher than {score} are not available: pl {pl_merged_peaks.shape}; mn {mn_merged_peaks.shape} ")
        self.roc_ready_peak = self.bid_pairs_with_score


class ProfilerHOMERgroseq(ROCProfiler):
    def roc(self, definition="element"):
        return self._pq_roc(definition=definition, F_restraint=0.01)

    def generate_scored_pairs(self, inputs, **kwargs):
        self.roc_ready_peak = {0: inputs}


class ProfilerdREG(ROCProfiler):
    def generate_scored_pairs(self, inputs, **kwargs):
        """
        Generate scored pairs for dREG

        Original dREG scores will be served as scores for pairs

        Parameters
        ----------
        inputs : str
            riginal peak file from dREG (.dREG.peak.full.bed.gz)

        Returns
        -------

        """
        peak_file = inputs
        peak_df = pd.read_csv(peak_file, sep="\t", header=None)
        peak_df = peak_df.loc[:, (0, 1, 2, 3)]
        peak_df.to_csv(self.output_prefix + ".bed", sep="\t", header=False, index=False)
        self.roc_ready_peak = self.output_prefix + ".bed"


class ProfilerTfit(ROCProfiler):
    def generate_scored_pairs(self, inputs, **kwargs):
        """
        Generate scored pairs for Tfit

        Parameters
        ----------
        inputs : str
            Output from Tfit, miscs in the fourth column will be parsed, then scores will be extracted from this column.
            Tfit release newer than 20eb2aa no longer output score, so they cannot be used to calc ROC.

        Returns
        -------

        """
        tre_file = inputs
        df = pd.read_csv(tre_file, sep="\t", comment="#", header=None, names=("chr", "start", "end", "misc"))
        df["score"] = df["misc"].apply(lambda x: x.split("|")[1].split(",")[0])
        df = df.loc[:, ("chr", "start", "end", "score")]
        df.to_csv(self.output_prefix + ".bed", sep="\t", header=False, index=False)
        self.roc_ready_peak = self.output_prefix + ".bed"


class ProfilerTSScall(ROCProfiler):
    def roc(self, definition="element"):
        return self._pq_roc(definition=definition)

    def generate_scored_pairs(self, inputs, **kwargs):
        """
        Generate scored pairs for TSScall

        Parameters
        ----------
        inputs : str
            Path to the folder contains ROC data for TSScall.
            Peak calls can be generated by running `helpers/TSScall/roc_pr_pre.py`

        Returns
        -------

        """
        roc_folder = inputs
        summary_file = os.path.join(roc_folder, "TSScall.csv")
        tsscall_summary = pd.read_csv(summary_file)
        bid_pairs_with_score = dict()

        # # make sure enough precision
        getcontext().prec = 52
        for nr, row in tsscall_summary.iterrows():
            file_path = os.path.join(roc_folder, row["Output_prefix"] + ".bid.bed")
            real_k = Decimal(1) - Decimal(row["Threshold"])
            bid_pairs_with_score[real_k] = file_path
        self.roc_ready_peak = bid_pairs_with_score
        sample_points = len(self.roc_ready_peak)
        self.recalls = np.zeros(sample_points)
        self.specificities = np.zeros(sample_points)
        self.precisions = np.zeros(sample_points)
        self.fscores = np.zeros(sample_points)
        self.fscores_alt = np.zeros(sample_points)
        self.cutoffs = np.empty(sample_points, dtype=object)


class ProfilerdREGHD(ComplexScoreProfiler):
    pass


class ProfilerGROcapTSSHMM(ComplexScoreProfiler):
    pass


class ProfilerFivePrime(ComplexScoreProfiler):
    pass


class ProfilerPINTS(ROCProfiler):
    def __init__(self, output_prefix, positive_set, negative_set, cache_dir=".",
                 sample_points=101, pairing_distance=300):
        super(ProfilerPINTS, self).__init__(output_prefix=output_prefix, positive_set=positive_set,
                                            negative_set=negative_set, cache_dir=cache_dir,
                                            sample_points=sample_points, pairing_distance=pairing_distance)
        self.bid_pairs_with_score = dict()

    def roc(self, definition="element"):
        return self._pq_roc(definition=definition)

    def generate_scored_pairs(self, inputs, **kwargs):
        import shutil
        pints_executable = shutil.which("pints_caller")
        assert pints_executable is not None, "PINTS is not properly installed"
        if not os.path.exists(pints_executable+".py"):
            os.symlink(pints_executable, pints_executable+".py")
        pints_folder, _ = os.path.split(pints_executable)
        if pints_folder not in sys.path:
            sys.path.insert(0, pints_folder)
        from pints_caller import merge_opposite_peaks
        from pints.io_engine import index_bed_file, parse_gtf
        stringent_pairs_only = kwargs.get("stringent_only", False)
        min_len_opposite_peaks = kwargs.get("min_len_opposite_peaks", 0)
        pl_file = inputs + "pl.bed.gz"
        mn_file = inputs + "mn.bed.gz"
        gtf_file = inputs + "peaks" + ".gtf"
        # be compatible to old versions of PINTS
        if os.path.exists(pl_file) and os.path.exists(mn_file):
            pass
        elif os.path.exists(gtf_file):
            peak_gtf = parse_gtf(gtf_file)
            # peak_gtf.score: qvalue
            # peak_gtf.pval: pvalue
            COMMON_HEADER = ('chromosome', 'start', 'end', 'name', 'padj', 'strand', 'reads',
                             'pval', 'mu_0', 'pi_0', 'mu_1', 'pi_1', 'ler_1', 'ler_2', 'ler_3', 'summit')
            col_name_mapping = {
                "seqname": "chromosome", "start": "start", "end": "end", "score": "padj", "strand": "strand",
                "reads": "reads", "pval": "pval", "mu_bg": "mu_0", "pi_bg": "pi_0", "mu_peak": "mu_1",
                "pi_peak": "pi_1", "ler1": "ler_1", "ler2": "ler_2", "ler3": "ler_3", "summit": "summit",
                "peak_name": "name"
            }
            peak_gtf.columns = peak_gtf.columns.map(col_name_mapping)
            peak_gtf.padj = peak_gtf.padj.astype(float)
            pl_peaks = peak_gtf.loc[peak_gtf.strand == "+", COMMON_HEADER]
            pl_peaks.to_csv(f"{inputs}pl.bed", sep="\t", header=False, index=False)
            index_bed_file(f"{inputs}pl.bed")
            mn_peaks = peak_gtf.loc[peak_gtf.strand == "-", COMMON_HEADER]
            mn_peaks.to_csv(f"{inputs}mn.bed", sep="\t", header=False, index=False)
            index_bed_file(f"{inputs}mn.bed")
        else:
            raise DependencyError("Peak calls for ROC not found.")
        max_bidirectional_distance = self._pairing_distance
        pe_extended_bed = pd.read_csv(self._positive_set, sep="\t", header=None)
        ne_extended_bed = pd.read_csv(self._negative_set, sep="\t", header=None)
        known_regions = pd.concat([pe_extended_bed, ne_extended_bed])
        known_regions[1] -= self._pairing_distance
        known_regions[2] += self._pairing_distance
        known_regions = BedTool.from_dataframe(known_regions)
        logger.info("Selecting peaks only in known regions (SPOIKR)")
        pl_df = pd.read_csv(pl_file, sep="\t", header=None)
        pints_all_pl_peaks_df = BedTool.from_dataframe(pl_df).intersect(known_regions, u=True).to_dataframe(names=(
            'chromosome', 'start', 'end', 'name', 'padj', 'strand', 'reads', 'pval', 'mu_0', 'pi_0', 'mu_1', 'pi_1',
            'ler_1', 'ler_2', 'ler_3', 'summit'))
        del pl_df
        logger.info("SPOIKR pl done")
        mn_df = pd.read_csv(mn_file, sep="\t", header=None)
        pints_all_mn_peaks_df = BedTool.from_dataframe(mn_df).intersect(known_regions, u=True).to_dataframe(names=(
            'chromosome', 'start', 'end', 'name', 'padj', 'strand', 'reads', 'pval', 'mu_0', 'pi_0', 'mu_1', 'pi_1',
            'ler_1', 'ler_2', 'ler_3', 'summit'))
        del mn_df
        logger.info("SPOIKR mn done")
        pints_all_pl_peaks_df.columns = np.arange(pints_all_pl_peaks_df.shape[1])
        pints_all_mn_peaks_df.columns = np.arange(pints_all_mn_peaks_df.shape[1])

        # small_pvals = [10**(-1*e) for e in range(100, 1, -1)]
        # whether to transform p-values to log scale
        # if converted to LS, the minimum scores should be assigned to TREs
        # pints_all_pl_peaks_df[7] = -1*np.log1p(pints_all_pl_peaks_df[7])
        # threshold = -1*np.log(score)
        # keep_score_criteria = "min"
        # if use original p-vals, then the max p-val should be kept
        keep_score_criteria = "max"
        for score in sorted(self._unbalanced_range_for_pqs(local_min=pints_all_pl_peaks_df[7].min())):
            # exp_pvals = np.linspace(1, 1000, 30).tolist()
            # general_pvals = (-1 * np.log10(np.linspace(0.105, 1, self._sample_points - len(exp_pvals)))).tolist()
            # c = set(exp_pvals + general_pvals)
            # for i, score in enumerate(c):
            real_k = Decimal(1) - Decimal(score)
            # real_k = score
            sig_pl_df = pints_all_pl_peaks_df.loc[pints_all_pl_peaks_df[7] <= score, :]
            sig_mn_df = pints_all_mn_peaks_df.loc[pints_all_mn_peaks_df[7] <= score, :]
            if sig_pl_df.shape[0] == 0 or sig_mn_df.shape[0] == 0:
                logger.warning(f"No peak passed this threshold: {real_k}")
                continue

            fn_pl_sig = os.path.join(self._cache_dir, "sig_pl.bed")
            fn_mn_sig = os.path.join(self._cache_dir, "sig_mn.bed")
            fn_pl_div_peak = os.path.join(self._cache_dir, "sig_pl_divergent_peaks.bed")
            fn_mn_div_peak = os.path.join(self._cache_dir, "sig_mn_divergent_peaks.bed")
            fn_pl_bid_peak = os.path.join(self._cache_dir, "sig_pl_bidirectional_peaks.bed")
            fn_mn_bid_peak = os.path.join(self._cache_dir, "sig_mn_bidirectional_peaks.bed")
            fn_pl_single_peak = os.path.join(self._cache_dir, "sig_pl_singletons.bed")
            fn_mn_single_peak = os.path.join(self._cache_dir, "sig_mn_singletons.bed")
            sig_pl_df.to_csv(fn_pl_sig, sep="\t", header=None, index=False)
            sig_mn_df.to_csv(fn_mn_sig, sep="\t", header=None, index=False)
            merge_opposite_peaks(fn_pl_sig, mn_file,
                                 divergent_output_bed=fn_pl_div_peak,
                                 bidirectional_output_bed=fn_pl_bid_peak,
                                 singleton_bed=fn_pl_single_peak,
                                 fdr_target=score, close_threshold=max_bidirectional_distance,
                                 div_size_min=0, summit_dist_min=0,
                                 min_len_opposite_peaks=min_len_opposite_peaks,
                                 stringent_only=stringent_pairs_only)
            merge_opposite_peaks(fn_mn_sig, pl_file,
                                 divergent_output_bed=fn_mn_div_peak,
                                 bidirectional_output_bed=fn_mn_bid_peak,
                                 singleton_bed=fn_mn_single_peak,
                                 fdr_target=score, close_threshold=max_bidirectional_distance,
                                 div_size_min=0, summit_dist_min=0,
                                 min_len_opposite_peaks=min_len_opposite_peaks,
                                 stringent_only=stringent_pairs_only)

            merged_bidirectional_file = os.path.join(self._cache_dir, "PINTS_%s_bidirectional_peaks.bed" % (1 - score))
            if os.path.exists(merged_bidirectional_file):
                try:
                    os.remove(merged_bidirectional_file)
                except Exception as e:
                    logger.error(e)

            command = "cat %s %s | sort -k1,1 -k2,2n | awk 'BEGIN{OFS=\"\t\"}{print $1,$2,$3}' >> %s" % (
                fn_pl_bid_peak, fn_mn_bid_peak, merged_bidirectional_file)
            self.bid_pairs_with_score[real_k] = merged_bidirectional_file
            run_command(command)
        self.roc_ready_peak = self.bid_pairs_with_score


def interpolate_roc(specificity, recall, n_points=100, kind="linear"):
    """
    Interpolate ROCs

    Parameters
    ----------
    specificity : array-like
        Specificity values
    recall : array-like
        Recall / sensitivity values
    n_points : int
        Number of points between 0 and 1 to be drawn
    kind : str
        'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic',
        'previous', 'next', where 'zero', 'slinear', 'quadratic' and 'cubic'

    Returns
    -------
    X: array-like
        1-specificity
    Y: array-like
        Interpolated recall
    """
    from scipy.interpolate import interp1d
    x = 1 - specificity
    y = recall
    f = interp1d(x, y, kind=kind)
    standard_x = np.linspace(1, 0, n_points)
    interpolated_y = f(standard_x)
    return standard_x, interpolated_y
