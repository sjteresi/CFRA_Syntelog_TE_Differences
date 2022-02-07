#!/usr/bin/env/python

"""
Generate an updates dotplot of avg TE density values for H4
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np

from TE_Density.transposon.gene_data import GeneData
from TE_Density.transposon.density_data import DensityData
from TE_Density.transposon.import_filtered_genes import import_filtered_genes

from src.bias_utils import verify_te_type_sequence
from src.bias_utils import verify_window_sequence
from src.bias_utils import verify_te_index_dictionary_sequence

from collections import defaultdict

import matplotlib.pyplot as plt


def calculate_mean_over_chromosomes(list_of_density_data, te_type, window_idx):
    """
    Calculate a tuple (left, right) of means for a given TE type, and window
    value, when considering all chromosomes (list of density data).


    Args:
        list_of_density_data (list of DensityData)
        te_type (str)
        window_idx (int)

    Returns:
        (left_mean, right_mean): A tuple of two floats that represents the mean
        of a specific te type for a specific window for all genes when
        considering all density data.
    """
    left_container = []
    right_container = []
    for a_density_data in list_of_density_data:
        all_genes_left = a_density_data.left_orders[
            a_density_data.order_index_dict[te_type], window_idx, :
        ]
        all_genes_right = a_density_data.right_orders[
            a_density_data.order_index_dict[te_type], window_idx, :
        ]
        left_container.append(all_genes_left)
        right_container.append(all_genes_right)

    # Now we have a list of lists, where each individual list is all genes for
    # that individual chromosome (density data)
    # We will now flatten, and take a mean
    left_flattened = np.concatenate(left_container).ravel()
    right_flattened = np.concatenate(right_container).ravel()
    left_mean = np.mean(left_flattened)
    right_mean = np.mean(right_flattened)
    return (left_mean, right_mean)


def calulate_intra_mean_over_chromosomes(list_of_density_data, te_type):
    """
    Take all DensityData, extract the arrays of data for all genes for a
    specific te_type and window. Concatenate the information for several
    chromosomes and then take the mean, because we cannot average the means of
    individual chromosomes.

    Args:
        list_of_density_data (list of DensityData)
        te_type (str)

    Returns:
        intra_mean: A float that represents the mean
        of a specific te type for all genes when
        considering all density data.


    """
    intra_container = []
    for a_density_data in list_of_density_data:
        all_genes_intra = a_density_data.intra_orders[
            a_density_data.order_index_dict[te_type], :, :
        ]
        intra_container.append(all_genes_intra.ravel())

    intra_flattened = np.concatenate(intra_container).ravel()
    intra_mean = np.mean(intra_flattened)
    return intra_mean


def gather_avg_for_dds(list_of_density_data, output_dir):
    """
    Iterate over the TE types and windows. Call the functions that calculate
    mean over all DensityData objs for a given window/TE type combo. Store the
    data in an array and call the graphing command when the arrays are complete

    Args:
        list_of_density_data (list of DensityData): Each DensityData represents
        the TE Density data for one chromosome.
        output_dir (str): Path to output directory, not used in this function
        but in the call to plot_density_dotplot.

    Returns: None, calls plot_density_dotplot to create a graph and save to
    disk
    """
    dd_obj = list_of_density_data[0]
    te_index_dict = dd_obj.order_index_dict
    data_left_dict = defaultdict(list)
    data_intra_dict = defaultdict(list)
    data_right_dict = defaultdict(list)

    for te_type, te_type_idx in te_index_dict.items():
        if "Revision" in te_type:
            continue
        for window_idx, window_val in enumerate(dd_obj.window_list):
            left_mean, right_mean = calculate_mean_over_chromosomes(
                list_of_density_data, te_type, window_idx
            )
            data_left_dict[te_type].append(left_mean)
            data_right_dict[te_type].append(right_mean)
        data_intra_dict[te_type].append(
            calulate_intra_mean_over_chromosomes(list_of_density_data, te_type)
        )

    window_list = dd_obj.window_list
    plot_density_dotplot(
        data_left_dict, data_intra_dict, data_right_dict, window_list, output_dir
    )


def plot_density_dotplot(
    data_left_dict, data_intra_dict, data_right_dict, window_list, output_dir
):
    """
    Use Matplotlib to generate a dotplot of avg TE density values segregated by
    TE type and window value genome wide. Does this for upstream, intragenic,
    and downstream datasets. Heavy use of MAGIC numbers to achieve
    good look of plot.

    Args:
        data_left_dict (dictionary of list): A dictionary with TE types as
            keys, and a list as values. This dictionary corresponds to the
            upstream values. The list is len(window_list). The
            values are the avg TE density value for all genes, for a given
            window. The values are ordered sequentially from least greatest
            window to the greatest window. Thus the last value corresponds to
            the greatest window (here 10KB).

        data_intra_dict (dictionary of list): A dictionary with TE types as
        keys, and a list as values. This dictionary corresponds to the
        intragenic values. The list is length 1. The values are the avg TE
        density values for all genes for the intragenic "window".

        data_right_dict (dictionary of list): A dictionary with TE types as
            keys, and a list as values. This dictionary corresponds to the
            downstream values. The list is len(window_list). The
            values are the avg TE density value for all genes, for a given
            window. The values are ordered sequentially from least greatest
            window to the greatest window. Thus the last value corresponds to
            the greatest window (here 10KB).

        window_list (list of int): Used to label the x-axis, corresponds to the
        values in the data dictionaries

        output_dir (str): Path to output directory

    Returns: None, generates a graph.
    """
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey="col")
    fig.set_size_inches(16, 9.5)
    # define colors
    NUM_COLORS = sum(1 for te_type in data_left_dict.items())
    cm = plt.get_cmap("tab20")
    ax1.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])
    ax2.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])
    ax3.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])

    for key, val in data_left_dict.items():
        ax1.plot(
            window_list,
            val,
            label=key,
            linestyle=(0, (3, 1, 1, 1)),
            marker="o",
        )
    ax1.set(
        xlabel="BP Upstream",
        ylabel="TE Density",
        xlim=[max(window_list), min(window_list)],
        xticks=range(min(window_list), (max(window_list) + 1), 1000),
        ylim=[0.0, 0.6],  # MAGIC
    )

    # MAGIC, add panel 'A' label to the left side graph
    ax1.text(
        -0.01,
        1.05,
        "A",
        transform=ax1.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="right",
    )

    for key, val in data_intra_dict.items():
        ax2.scatter(
            0,
            val,
            label=key,
        )

    ax2.set(xlabel="Intragenic TEs", xticks=[])
    ax2.legend(loc="upper right", bbox_to_anchor=(0.76, 0.8))  # MAGIC
    # Heavily tailored location of legend in the intragenic plot, could just use
    # loc='center' but this occasionally obscured a point

    # MAGIC, add panel 'B' label to the center graph
    ax2.text(
        -0.01,
        1.05,
        "B",
        transform=ax2.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="right",
    )

    for key, val in data_right_dict.items():
        ax3.plot(
            window_list,
            val,
            label=key,
            linestyle=(0, (3, 1, 1, 1)),
            marker="o",
        )
    ax3.set(
        xlabel="BP Downstream",
        xlim=[min(window_list), max(window_list)],
        xticks=range(min(window_list), (max(window_list) + 1), 1000),
        ylim=[0.0, 0.6],  # MAGIC
    )
    ax3.yaxis.tick_right()

    # MAGIC, add panel 'C' label to the right side graph
    ax3.text(
        -0.01,
        1.05,
        "C",
        transform=ax3.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="right",
    )

    # MAGIC, name of output file.
    filename_to_save = os.path.join(output_dir, "Orders_H4_Combined_Density_Plot.png")

    fig.suptitle(
        "Average TE Density of All Genes as a Function of Window Size and Location"
    )
    logger.info("Saving graphic to: %s" % filename_to_save)
    plt.savefig(
        filename_to_save,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    data_dir = os.path.abspath(os.path.join(dir_main, "../../", "Domestication_Data"))
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "density_H4_data_folder",
        type=str,
        help="Parent path of folders containing TE Density results",
    )

    parser.add_argument(
        "gene_data_H4",
        type=str,
        help="Path to H4's filtered gene data file",
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
    args.density_H4_data_folder = os.path.abspath(args.density_H4_data_folder)
    args.gene_data_H4 = os.path.abspath(args.gene_data_H4)
    args.output_dir = os.path.abspath(args.output_dir)

    # Get logging set up
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Read the gene files into the pandas dataframe structure
    all_genes_H4 = import_filtered_genes(args.gene_data_H4, logger)

    # Break the pandaframe representing the complete genome gene annotation
    # into chromosome chunks, wrap the chunks in GeneData class and use each
    # chunk to initialize DensityData for that chromosome
    # MAGIC 'Chromosome' corresponds to column in pandaframe
    chromosomes_H4_panda_list = [
        dataframe for k, dataframe in all_genes_H4.groupby("Chromosome")
    ]

    # MAGIC 0 for the chromosome ID
    gene_data_H4_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in chromosomes_H4_panda_list
    ]

    # Initialize all DensityData objects for H4
    # NOTE MAGIC substring rule for matching the chromosome IDs and initializing
    # the DD objects
    
    # NOTE check once more if this code is ultimately used. Otherwise candidate for deletion.
    processed_H4_density_data = DensityData.from_list_gene_data_and_hdf5_dir(
        gene_data_H4_list,
        args.density_H4_data_folder,
        #"Fragaria_562_(562_scaffold_(\d)).h5",
        "H4_(Fvb(\d)).h5",
        logger,
    )

    # NOTE reordering to same chromosome order for easier processing later
    processed_H4_density_data = sorted(
        processed_H4_density_data, key=lambda x: x.genome_id, reverse=False
    )

    for a_dd in processed_H4_density_data:
        for a_second_dd in processed_H4_density_data:
            verify_te_type_sequence(a_dd, a_second_dd)
            verify_te_index_dictionary_sequence(a_dd, a_second_dd)
            verify_window_sequence(a_dd, a_second_dd)

    gather_avg_for_dds(processed_H4_density_data, args.output_dir)
