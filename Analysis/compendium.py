#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Li Yao
# @Date: 3/13/21
import logging
import argparse
import os
import json
import pybedtools
import pandas as pd
from pybedtools import BedTool
from configparser import ConfigParser
from file_parser_mem.bed import parse_bed

logging.basicConfig(format='%(name)s - %(asctime)s - %(levelname)s: %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO,
                    handlers=[
                        logging.FileHandler(os.path.join(os.getcwd(), 'Compendium.log')),
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("PINTS - Compendium")


def percent_peaks_covered_by_other_assays(peak_calls, promoters_bed, read_table_file, save_to, global_registry, local_registry):
    """

    Parameters
    ----------
    peak_calls : dict
        key: name of assays
        value: peak calls generated for this assay
    promoters_bed : str
        Path to bed file defining promoter regions
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
        all_assays_in_scope = set(peak_calls.keys())
        covered_by_other_assays = []
        reads_dict = pd.read_csv(read_table_file, index_col=0).set_index("Assay").to_dict()["Reads"]
        for assay in all_assays_in_scope:
            # distal peaks from one assay
            pkf = BedTool.from_dataframe(parse_bed(peak_calls[assay])).intersect(promoters_bed, v=True)

            # other assays to be compared
            other_assays = all_assays_in_scope.difference({assay})

            # merge all peaks from other assays
            other_peak_calls = [BedTool.from_dataframe(parse_bed(peak_calls[oa])).sort() for oa in other_assays]
            print(len(other_peak_calls))
            other_unidirectional_peak_calls = [BedTool.from_dataframe(parse_bed(peak_calls[oa].replace("bidirectional", "unidirectional"))).sort() for oa in other_assays]
            other_peak_calls.extend(other_unidirectional_peak_calls)
            print(len(other_peak_calls))
            peak_calls_other_assay = BedTool.cat(*other_peak_calls).intersect(promoters_bed, v=True)

            # get shared and unique peaks
            shared_peaks = pkf.intersect(peak_calls_other_assay, u=True)
            covered_by_other_assays.append((assay, len(pkf), len(shared_peaks), reads_dict[assay]))
        covered_by_other_assays_df = pd.DataFrame(covered_by_other_assays)
        covered_by_other_assays_df.to_csv(final_result_file)
    return final_result_file


def main(peak_calls, reads_table, data_save_to, promoters_bed, data_prefix):
    analysis_summaries = {
        "n_tres": [],
    }
    analysis_summaries["n_tres"].append(percent_peaks_covered_by_other_assays(peak_calls,
                                                                                promoters_bed,
                                                                                reads_table,
                                                                                data_save_to,
                                                                                global_registry,
                                                                                "nTREs"))

    with open(os.path.join(data_save_to, f"{data_prefix}_summary.json"), "w") as fh:
        json.dump(analysis_summaries, fh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-save-to", required=True, help="Save output data to")
    parser.add_argument("--bed-promoter", required=False,
                        default="/local/storage/ly349/refs/annotations/human/hg38/segmentation/promoters_1kb_tss_centered.bed",
                        help="Bed file which includes definitions for promoters")
    parser.add_argument("--data-prefix", required=False, default="compendium", help="Prefix for all outputs")
    parser.add_argument("--reads-table", required=False, help="Reads table file",
                        default="/local/storage/ly349/projects/peakcalling/data/data4plot/2_TotalReadCounts.csv")
    parser.add_argument("--global-registry", type=int, help="Global registry", default=5)
    parser.add_argument("--config-file", default="config.conf", help="Configuration file")
    parser.add_argument("--tmp-dir", default="./tmp", help="All temporary files will be written into this folder")

    args = parser.parse_args()

    assert os.path.exists(args.config_file)

    if not os.path.exists(args.tmp_dir):
        os.mkdir(args.tmp_dir)

    global unified_color_map, unified_color_map_s
    cfg = ConfigParser()
    cfg.optionxform = str
    cfg.read(args.config_file)
    collapsed_assays = cfg.get("assays", "plot_order_simplified").split("|")
    unified_color_map = dict()
    official_name_map = dict()
    assay_offical_names = dict()
    plot_order = cfg.get("assays", "plot_order_full").split("|")
    plot_order_simplified = cfg.get("assays", "plot_order_simplified").split("|")
    plot_color = cfg.get("assays", "plot_colors").split("|")
    assay_offical_names = cfg.get("assays", "assay_full_names").split("|")
    for k, v in enumerate(plot_order):
        unified_color_map[v] = plot_color[k]
        official_name_map[v] = assay_offical_names[k]

    global tmp_dir, n_threads, n_samples, global_registry

    tmp_dir = args.tmp_dir
    pybedtools.set_tempdir(tmp_dir)
    global_registry = args.global_registry

    peak_calls_from_raw_data = dict()
    if args.method == "generate_data":
        import socket

        server1_home_dir = "/fs/cbsuhy01/storage/ly349/" if socket.gethostname().find(
            "cbsuhy01") == -1 else "/local/storage/ly349/"
        server2_home_dir = "/fs/cbsuhy02/storage/ly349/" if socket.gethostname().find(
            "cbsuhy02") == -1 else "/local/storage/ly349/"
        bioq_dir = os.path.join(server2_home_dir, "BioQueue/workspace")

        from .utils import load_bioq_datasets
        tmp = load_bioq_datasets("peak_call_job_ids", bioq_path=bioq_dir, cfg_file=args.config_file)

        _5p_assays = {"GROcap", "CoPRO", "csRNAseq", "NETCAGE", "RAMPAGE", "CAGE", "STRIPEseq"}
        _assay_of_interests = {"GROcap", "csRNAseq", "NETCAGE", "RAMPAGE", "CAGE", "STRIPEseq"}
        for k, v in tmp.items():
            caller, gr, cl, assay, br, tr = k.split("_")
            if caller == "PINTS" and gr == "hg38" and assay in _5p_assays:
                peak_calls_from_raw_data[assay] = v

    main(method=args.method, peak_calls=peak_calls_from_raw_data,
         data_save_to=args.data_save_to, fig_save_to=args.fig_save_to,
         promoters_bed=args.bed_promoter, data_prefix=args.data_prefix,
         reads_table=args.reads_table, name_mapping=official_name_map)
