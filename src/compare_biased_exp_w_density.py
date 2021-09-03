#!/usr/bin/env/python

"""
Compare TE Density values between syntelogs
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import pandas as pd
import numpy as np
import re

from TE_Density.transposon.gene_data import GeneData
from TE_Density.transposon.density_data import DensityData
from TE_Density.transposon.import_filtered_genes import import_filtered_genes

from src.bias_utils import verify_te_type_sequence
from src.bias_utils import verify_window_sequence
from src.bias_utils import verify_te_index_dictionary_sequence
from src.bias_utils import return_list_density_data_objs
from src.bias_utils import add_TE_values_to_table
from src.bias_utils import find_missing_genes
from src.bias_utils import read_processed_syntelog_table
from src.bias_utils import plot_reg_plot_bias

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    data_dir = os.path.abspath(os.path.join(dir_main, "../../", "Domestication_Data"))
    parser = argparse.ArgumentParser(
        description="compare density values between syntelogs"
    )

    parser.add_argument(
        "processed_syntelog_input_file",
        type=str,
        help="parent path to syntelog file",
    )

    parser.add_argument(
        "density_562_data_folder",
        type=str,
        help="Parent path of folders containing TE Density results",
    )

    parser.add_argument(
        "density_1008_data_folder",
        type=str,
        help="Parent path of folders containing TE Density results",
    )

    parser.add_argument(
        "gene_data_562",
        type=str,
        help="Path to 562's filtered gene data file",
    )

    parser.add_argument(
        "gene_data_1008",
        type=str,
        help="Path to 1008's filtered gene data file",
    )

    parser.add_argument(
        "biased_genes_file",
        type=str,
        help="path to a file that contains syntelog pairs",
    )

    parser.add_argument("genome_id", type=str, help="genome string")

    parser.add_argument(
        "level_of_bias",
        type=int,
        help="integer value representing the strength of bias in the user supplied files",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=data_dir,
        help="parent directory to output results",
    )

    # Redefine args to abspath
    args = parser.parse_args()
    args.processed_syntelog_input_file = os.path.abspath(
        args.processed_syntelog_input_file
    )
    args.density_562_data_folder = os.path.abspath(args.density_562_data_folder)
    args.density_1008_data_folder = os.path.abspath(args.density_1008_data_folder)
    args.gene_data_562 = os.path.abspath(args.gene_data_562)
    args.gene_data_1008 = os.path.abspath(args.gene_data_1008)
    args.biased_genes_file = os.path.abspath(args.biased_genes_file)
    args.output_dir = os.path.abspath(args.output_dir)

    # Get logging set up
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Read the preprocessed syntelog/syntelog file into the pandas dataframe structure
    syntelog_table_w_indices = read_processed_syntelog_table(
        args.processed_syntelog_input_file
    )

    # Read the gene files into the pandas dataframe structure
    all_genes_562 = import_filtered_genes(args.gene_data_562, logger)
    all_genes_1008 = import_filtered_genes(args.gene_data_1008, logger)

    # Break the pandaframe representing the complete genome gene annotation
    # into chromosome chunks, wrap the chunks in GeneData class and use each
    # chunk to initialize DensityData for that chromosome
    # MAGIC 'Chromosome' corresponds to column in pandaframe
    chromosomes_562_panda_list = [
        dataframe for k, dataframe in all_genes_562.groupby("Chromosome")
    ]
    chromosomes_1008_panda_list = [
        dataframe for k, dataframe in all_genes_1008.groupby("Chromosome")
    ]

    # MAGIC 0 for the chromosome ID
    gene_data_562_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in chromosomes_562_panda_list
    ]
    gene_data_1008_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in chromosomes_1008_panda_list
    ]

    # Initialize all DensityData objects for 562
    processed_562_density_data = return_list_density_data_objs(
        gene_data_562_list,
        args.density_562_data_folder,
        "Fragaria_562_(.*?).h5",
        logger,
    )
    # Initialize all DensityData objects for 1008
    processed_1008_density_data = return_list_density_data_objs(
        gene_data_1008_list,
        args.density_1008_data_folder,
        "Fragaria_1008_2339_(.*?).h5",
        logger,
    )

    # NOTE reordering to same chromosome order for easier processing later
    processed_562_density_data = sorted(
        processed_562_density_data, key=lambda x: x.genome_id, reverse=False
    )
    processed_1008_density_data = sorted(
        processed_1008_density_data, key=lambda x: x.genome_id, reverse=False
    )

    biased_genes_pandas = pd.read_csv(
        args.biased_genes_file,
        header=None,
        sep="\t",
        names=["1008_Biased", "562_Biased"],
    )

    biased_genes_list = biased_genes_pandas[str(args.genome_id + "_Biased")].tolist()

    #############################################
    # TODO refactor
    # Try a different way of getting the value pairs
    directions = ["LEFT", "RIGHT"]
    major_te_groupings = ["ORDERS", "SUPERFAMILIES"]

    for direction in directions:
        for major_te_grouping in major_te_groupings:
            for window_idx, window_val in enumerate(
                processed_1008_density_data[0].window_list
            ):
                if major_te_grouping == "ORDERS":
                    my_iterable = processed_1008_density_data[0].order_list
                    my_index_dict = "order_index_dict"
                if major_te_grouping == "SUPERFAMILIES":
                    my_iterable = processed_1008_density_data[0].super_list
                    my_index_dict = "super_index_dict"
                for te_grouping in my_iterable:
                    chromosome_subset_w_te_values_1008 = []
                    chromosome_subset_w_te_values_562 = []
                    for fragaria_562_dd_obj, fragaria_1008_dd_obj in zip(
                        processed_562_density_data, processed_1008_density_data
                    ):
                        verify_te_type_sequence(
                            fragaria_562_dd_obj, fragaria_1008_dd_obj
                        )  # NB maybe move because performance?
                        verify_window_sequence(
                            fragaria_562_dd_obj, fragaria_1008_dd_obj
                        )  # NB maybe move because performance?
                        verify_te_index_dictionary_sequence(
                            fragaria_562_dd_obj, fragaria_1008_dd_obj
                        )  # NB maybe move because performance?

                        # Copy this for the other genome
                        chromosome_subset_w_te_values_1008.append(
                            add_TE_values_to_table(
                                syntelog_table_w_indices,
                                fragaria_1008_dd_obj,
                                major_te_grouping,
                                te_grouping,
                                my_index_dict,
                                window_idx,
                                direction,
                                "1008_Chromosome",
                                "1008_Indices",
                                "1008_TE_Values",
                            )
                        )
                        chromosome_subset_w_te_values_562.append(
                            add_TE_values_to_table(
                                syntelog_table_w_indices,
                                fragaria_562_dd_obj,
                                major_te_grouping,
                                te_grouping,
                                my_index_dict,
                                window_idx,
                                direction,
                                "562_Chromosome",
                                "562_Indices",
                                "562_TE_Values",
                            )
                        )

                    chromosome_subset_w_te_values_1008 = pd.concat(
                        chromosome_subset_w_te_values_1008
                    )
                    chromosome_subset_w_te_values_562 = pd.concat(
                        chromosome_subset_w_te_values_562
                    )
                    graphing_table = chromosome_subset_w_te_values_1008.merge(
                        chromosome_subset_w_te_values_562, how="inner"
                    )

                    logger.warning(
                        find_missing_genes(
                            biased_genes_list,
                            graphing_table[args.genome_id + "_Gene"].tolist(),
                        )
                    )

                    graphing_table.loc[
                        graphing_table[str(args.genome_id + "_Gene")].isin(
                            biased_genes_list
                        ),
                        "Bias",
                    ] = "Present"

                    graphing_table["Bias"].fillna("Not_Present", inplace=True)

                    # Redefine for plot
                    to_graph = graphing_table.loc[graphing_table["Bias"] == "Present"]
                    bias_grouping = str(
                        str(args.level_of_bias) + "X_Towards_" + args.genome_id
                    )

                    plot_reg_plot_bias(
                        to_graph,
                        bias_grouping,
                        te_grouping,
                        window_val,
                        direction,
                        logger,
                        args.output_dir,
                    )
