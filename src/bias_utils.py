#!/usr/bin/env/python

"""
Helper functions to compare TE Density between genes that have a bias in
expression
"""

__author__ = "Scott Teresi"

# General imports
import argparse
import os
import logging
import coloredlogs
import pandas as pd
import numpy as np
import re

# Function imports
from TE_Density.transposon.gene_data import GeneData
from TE_Density.transposon.density_data import DensityData
from TE_Density.transposon.import_filtered_genes import import_filtered_genes

# Graphing imports
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


def find_missing_genes(list_of_interesting_genes, gene_list_from_HDF5):
    """
    Reports the genes provided in our interesting genes file (the file of genes
    with an expression bias), that are not in our dataframe containing TE
    values (thus they were not assessed when TE Density was run).

    This could be due to issues of gene name formatting, or the fact that the
    genes were never truly present during the running of TE Density

    Args:
        list_of_interesting_genes (list): list of genes derived from biased
        genes file, len will be less than the second argument
        gene_list_from_HDF5 (list): list of genes derived from HDF5 file
    Returns: a list of genes that were not present in the HDF5 (TE Density
        data)
    """
    return [
        item for item in list_of_interesting_genes if item not in gene_list_from_HDF5
    ]


def read_processed_syntelog_table(strawberry_syntelog_table):
    """
    Read a pandaframe from disk

    Args:
        strawberry_syntelog_table (str): Path to pandaframe saved on disk which is in
        .tsv format with the columns specified below

    Returns:
        dataframe (pandas.DataFrame): Pandas dataframe of the file that was
        provided
    """
    dataframe = pd.read_csv(
        strawberry_syntelog_table,
        sep="\t",
        header="infer",
    )
    return dataframe


def read_syntelog_table(strawberry_syntelog_table):
    """
    Read a pandaframe from disk

    Args:
        strawberry_syntelog_table (str): Path to pandaframe saved on disk which is in
        .tsv format with the columns specified below

    Returns:
        dataframe (pandas.DataFrame): Pandas dataframe of the file that was
        provided
    """
    dataframe = pd.read_csv(
        strawberry_syntelog_table,
        sep="\t",
        header="infer",
    )
    return dataframe


def supply_density_data_files(path_to_folder):
    """
    Iterate over a folder containing the H5 files of TE Density output and
    return a list of absolute file paths.

    Args:
        path_to_folder (str): path to the folder containing multiple h5 files
        of density data

    Returns:
        raw_file_list (list of str): A list containing the absolute paths to
        each relevant H5 file of density data
    """
    raw_file_list = []  # init empty list to store filenames
    for root, dirs, files in os.walk(path_to_folder):
        for a_file_object in files:
            # N.B very particular usage of abspath and join.
            a_file_object = os.path.abspath(os.path.join(root, a_file_object))
            if a_file_object.endswith(".h5"):  # MAGIC
                raw_file_list.append(a_file_object)

    return raw_file_list


def save_syntelog_table(syntelog_table, filename, output_dir, logger):
    """
    Saves an syntelog table (pandas.DataFrame) to disk

    Args:
        syntelog_table (pandas.DataFrame): A pandas dataframe of syntelog
        relationships between genes

        filename (str): String representing the name of the panda frame the
        user intends to save

        output_dir (str): The path to the output directory

        logger (logging.Logger): Object to pass logging commands through

    Returns:
        None
    """
    logger.info("Saving data to: %s" % (os.path.join(output_dir, filename)))
    syntelog_table.to_csv(
        os.path.join(output_dir, filename), sep="\t", header="True", index=False
    )


def verify_te_type_sequence(density_data_one, density_data_two):
    """
    Simple function to verify that the actual order of the TE types in the
    genomes are the same, as this would cause an error during the comparison
    stage if they had different indices in the HDF5. This could be due to
    different TE annotation sources but should not be an issue for this
    example. I am including this note for future work as it could be an easy
    thing to overlook.

    Args:
        density_data_one (DensityData): A density data instance representing a
        single genome

        density_data_two (DensityData): A density data instance representing a
        single genome

    Returns: None, raises errors if detects errors as described above
    """
    if density_data_one.order_list != density_data_two.order_list:
        raise ValueError(
            """The TE Orders in %s do not match the TE Orders in
                         %s"""
            % (density_data_one, density_data_two)
        )
    if density_data_one.super_list != density_data_two.super_list:
        raise ValueError(
            """The TE Superfamilies in %s do not match the TE Superfamilies in
                         %s"""
            % (density_data_one, density_data_two)
        )


def verify_window_sequence(density_data_one, density_data_two):
    """
    Simple function to verify that the actual order of the window values in the
    genomes are the same, as this would cause an error during the comparison
    stage if they had different indices in the HDF5. This could be due to
    running the program with different window settings. I am including this
    note for future work as it could be an easy thing to overlook.

    Args:
        density_data_one (DensityData): A density data instance representing a
        single genome

        density_data_two (DensityData): A density data instance representing a
        single genome

    Returns: None, raises errors if detects errors as described above
    """
    if density_data_one.window_list != density_data_two.window_list:
        raise ValueError(
            """The windows in %s do not match the windows in
                         %s"""
            % (density_data_one, density_data_two)
        )


def verify_te_index_dictionary_sequence(density_data_one, density_data_two):
    """
    This could be due to
    running the program with different window settings. I am including this
    note for future work as it could be an easy thing to overlook.

    Args:
        density_data_one (DensityData): A density data instance representing a
        single genome

        density_data_two (DensityData): A density data instance representing a
        single genome

    Returns: None, raises errors if detects errors as described above
    """
    if density_data_one.order_index_dict != density_data_two.order_index_dict:
        raise ValueError(
            """The order index dictionary in %s do not match the order index
            dictionary in %s"""
            % (density_data_one, density_data_two)
        )

    if density_data_one.super_index_dict != density_data_two.super_index_dict:
        raise ValueError(
            """The super index dictionary in %s do not match the super index
            dictionary in %s"""
            % (density_data_one, density_data_two)
        )


def edit_syntelog_table(syntelogs, gene_pandaframe_562, gene_pandaframe_1008):
    """
    Have to add chromosome IDs to the syntelog table because it was not given with
    chromosome IDs in order to appropriately subset out the DensityData

    Args:
        syntelogs (pandas.Data.Frame): Pandaframe representing the syntelogs.
            Has format columns=['1008_Gene', '562_Gene']
        gene_pandaframe_562 (pandas.Data.Frame): Pandaframe representing the
            GeneData info for the 562 genome.
        gene_pandaframe_1008 (pandas.Data.Frame): Pandaframe representing the
            GeneData info for the 1008 genome.

    Returns:
        syntelogs_w_chromosomes (pandas.Data.Frame): Has format
        columns=['1008_Gene', '562_Gene', '562_Chromosome',
        '1008_Chromosome']
    """
    syntelogs_w_chromosomes = pd.merge(
        left=syntelogs,
        right=gene_pandaframe_562,
        left_on="562_Gene",
        right_on="Gene_Name",
    )
    syntelogs_w_chromosomes.rename(
        columns={"Chromosome": "562_Chromosome"}, inplace=True
    )
    syntelogs_w_chromosomes.drop(
        columns=["Feature", "Start", "Stop", "Strand", "Length"], inplace=True
    )
    # 1008
    syntelogs_w_chromosomes = pd.merge(
        left=syntelogs_w_chromosomes,
        right=gene_pandaframe_1008,
        left_on="1008_Gene",
        right_on="Gene_Name",
    )
    syntelogs_w_chromosomes.rename(
        columns={"Chromosome": "1008_Chromosome"}, inplace=True
    )
    syntelogs_w_chromosomes.drop(
        columns=["Feature", "Start", "Stop", "Strand", "Length"], inplace=True
    )
    return syntelogs_w_chromosomes


def return_list_density_data_objs(
    list_of_gene_data, HDF5_folder, file_substring, logger
):
    """
    Returns a list of DensityData instances from a list of GeneData and
    HDF5 files that are gathered from an input directory

    Args:
        list_of_gene_data (list): List of GeneData instances
        HDF5_folder (str): Path to HDF5 folder containing results
        file_substring (str, MAGIC): MAGIC substring with which to identify the
            genome and chromosome IDs.

    Returns: processed_density_data (list of DensityData)
    """
    processed_density_data = []
    for raw_hdf5_data_file in supply_density_data_files(HDF5_folder):
        current_hdf5_file_chromosome = re.search(
            file_substring, raw_hdf5_data_file
        ).group(1)
        for gene_data_obj in list_of_gene_data:
            if gene_data_obj.chromosome_unique_id == current_hdf5_file_chromosome:
                dd_data_obj = DensityData.verify_h5_cache(
                    raw_hdf5_data_file, gene_data_obj, logger
                )
                processed_density_data.append(dd_data_obj)
    return processed_density_data


def add_indices_of_genes_in_HDF5_to_syntelog_table(
    list_of_density_data, gene_list, syntelog_table, genome_substring_col_to_create
):
    """
    Take a list of DensityData instances and a syntelog table (with gene IDs as
    columns), and create a new column in the syntelog table which contains that
    gene ID's index in the DensityData data

    Args:
        list_of_density_data (list):
        gene_list (list): list of genes from the syntelog table
        syntelog_table (pandas.Data.Frame): columns=['1008_Gene', '562_Gene', '562_Chromosome',
        '1008_Chromosome']
        genome_substring_col_to_create (str): String identifer for the new
            column to be created which represents the index values of the genes in
            the HDF5
    """
    super_dict = {}  # create empty dictionary
    for density_data_obj in list_of_density_data:  # loop over the DD objs
        to_rename = {
            dd_gene: dd_index
            for dd_index, dd_gene in enumerate(density_data_obj.gene_list)
        }

        super_dict.update(to_rename)  # Update the dictionary
    syntelog_table[genome_substring_col_to_create] = [
        super_dict.get(gene_name, None) for gene_name in gene_list
    ]  # Add the dictionary information to the syntelog_table
    return syntelog_table


def add_TE_values_to_table(
    syntelog_table_w_indices,
    fragaria_dd_obj,
    major_te_grouping,
    te_grouping,
    te_index_dictionary,
    window_idx,
    direction,
    chromosome_id,
    index_id,
    TE_values_id,
):
    """
    Taking the complete syntelog table (which at this point represents ALL
    genes (all chromosomes) and contains the HDF5 index information), for a
    specific TE type, window, and direction combination, add the genes' TE
    density values to the syntelog table

    Args:
        syntelog_table_w_indices ():
        fragaria_dd_obj ():
        major_te_grouping ():
        te_grouping ():
        te_index_dictionary ():
        window_idx ():
        direction ():
        chromosome_id ():
        index_id ():
        TE_values_id ():
    Returns:
    """

    # Create subsets via chromosome
    chromosome_subset = syntelog_table_w_indices.loc[
        syntelog_table_w_indices[chromosome_id] == fragaria_dd_obj.unique_chromosome_id
    ]

    # Need to sort by index value
    chromosome_subset = chromosome_subset.sort_values(by=[index_id])

    # Adding in TE values now
    chromosome_subset[TE_values_id] = fragaria_dd_obj.data_frame[
        str("RHO_" + major_te_grouping + "_" + direction)
    ][
        getattr(fragaria_dd_obj, te_index_dictionary)[te_grouping],
        window_idx,
        chromosome_subset[index_id].tolist(),
    ]

    return chromosome_subset


def plot_reg_plot_bias(
    graphing_matrix,
    bias_descriptor,
    te_grouping,
    window_val,
    direction,
    logger,
    output_dir,
):
    """
    TODO


    Args:

        graphing_matrix ():
        bias_descriptor ():
        te_grouping ():
        window_val ():
        direction ():
        logger ():
        output_dir ():

    Returns:
    """
    # MAGIC
    if direction == "LEFT":
        direction_desc = "Upstream"
    if direction == "RIGHT":
        direction_desc = "Downstream"

    plt.figure(figsize=(8, 6))
    sns.regplot(
        data=graphing_matrix,
        x="562_TE_Values",
        y="1008_TE_Values",
        line_kws={"color": "skyblue"},
    )
    plt.ylim(0.0, 1.0)
    plt.xlim(0.0, 1.0)

    # Remove this if we don't want the 1:1 diagonal line
    diagonal = [0.0, 1.0]
    plt.plot(diagonal, diagonal)

    plt.title(bias_descriptor.replace("_", " "))
    N = mpatches.Patch(
        label="Total Plotted Genes: %s \nTE type: %s \nWindow %s \nDirection: %s"
        % (len(graphing_matrix), te_grouping, window_val, direction_desc)
    )
    plt.legend(handles=[N], loc="upper left")
    plt.tight_layout()
    # plt.show()
    filename_to_save = os.path.join(
        output_dir,
        (
            te_grouping
            + "_"
            + str(window_val)
            + "_BP_"
            + bias_descriptor
            + "_"
            + direction_desc
        ),
    )
    logger.info("Saving figure %s" % filename_to_save)
    plt.savefig(filename_to_save)
    plt.clf()
    plt.close("all")
