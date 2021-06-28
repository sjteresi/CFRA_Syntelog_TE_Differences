#!/usr/bin/env/python

"""
Take the previously filtered syntelog file and add columns to it representing
the gene's index in the TE Density (HDF5) data.

Later, when comparing the TE density values of interesting genes (biased or
unbiased genes) we use this table as a lookup table to get the TE values
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
from src.bias_utils import edit_syntelog_table
from src.bias_utils import return_list_density_data_objs
from src.bias_utils import find_missing_genes
from src.bias_utils import read_syntelog_table
from src.bias_utils import add_indices_of_genes_in_HDF5_to_syntelog_table
from src.bias_utils import save_syntelog_table


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    data_dir = os.path.abspath(os.path.join(dir_main, "../../", "Domestication_Data"))
    parser = argparse.ArgumentParser(
        description="compare density values between syntelogs"
    )

    parser.add_argument(
        "syntelog_input_file",
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
    args.syntelog_input_file = os.path.abspath(args.syntelog_input_file)
    args.density_562_data_folder = os.path.abspath(args.density_562_data_folder)
    args.density_1008_data_folder = os.path.abspath(args.density_1008_data_folder)
    args.gene_data_562 = os.path.abspath(args.gene_data_562)
    args.gene_data_1008 = os.path.abspath(args.gene_data_1008)
    args.output_dir = os.path.abspath(args.output_dir)

    # Get logging set up
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Read the syntelog/syntelog file into the pandas dataframe structure
    syntelogs = read_syntelog_table(args.syntelog_input_file)

    # Read the gene files into the pandas dataframe structure
    all_genes_562 = import_filtered_genes(args.gene_data_562, logger)
    all_genes_1008 = import_filtered_genes(args.gene_data_1008, logger)

    # Add chromosome of each gene to the syntelog table
    syntelogs_w_chromosomes = edit_syntelog_table(
        syntelogs, all_genes_562, all_genes_1008
    )

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

    # Create lists of the genes for each genome
    fragaria_562_orth_genes = syntelogs_w_chromosomes["562_Gene"].tolist()
    fragaria_1008_orth_genes = syntelogs_w_chromosomes["1008_2339_Gene"].tolist()

    syntelogs_w_chromosomes = add_indices_of_genes_in_HDF5_to_syntelog_table(
        processed_562_density_data,
        fragaria_562_orth_genes,
        syntelogs_w_chromosomes,
        "562_Indices",
    )
    syntelogs_w_chromosomes = add_indices_of_genes_in_HDF5_to_syntelog_table(
        processed_1008_density_data,
        fragaria_1008_orth_genes,
        syntelogs_w_chromosomes,
        "1008_Indices",
    )

    # Only get the syntelog table where we have a corresponding syntelog pair
    syntelog_table_w_indices = syntelogs_w_chromosomes[
        syntelogs_w_chromosomes[["562_Indices", "1008_Indices"]].notnull().all(1)
    ]  # the all command wants them both to be true

    # MAGIC name for file
    save_syntelog_table(
        syntelog_table_w_indices, "Syntelog_Matches.tsv", args.output_dir, logger
    )
