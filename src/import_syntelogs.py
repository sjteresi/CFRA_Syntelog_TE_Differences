#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import os
import argparse


def import_syntelogs(syntelog_input_file, data_output_path, filename_to_save):
    """
    Import the syntelogs from the raw file and manage data filtration

    Args:


    """

    syntelog_table = pd.read_csv(
        syntelog_input_file,
        sep=",",
        header="infer",
    )
    # Drop anything that isn't a primary transcript
    syntelog_table = syntelog_table[
        syntelog_table["1008_2339_Gene"].str.contains(r"-mRNA-1")
    ]
    syntelog_table = syntelog_table[syntelog_table["562_Gene"].str.contains(r"-mRNA-1")]

    # Drop any row that has a gene that appears more than once
    syntelog_table.drop_duplicates(subset=["562_Gene"], keep=False, inplace=True)
    syntelog_table.drop_duplicates(subset=["1008_2339_Gene"], keep=False, inplace=True)

    # Get the correct name for the genes
    # MAGIC to split the name correctly
    syntelog_table["1008_2339_Gene"] = syntelog_table["1008_2339_Gene"].str.replace(
        "-mRNA-\d", "", regex=True
    )
    syntelog_table["562_Gene"] = syntelog_table["562_Gene"].str.replace(
        "-mRNA-\d", "", regex=True
    )

    file_to_save = os.path.join(data_output_path, filename_to_save)
    syntelog_table.to_csv(file_to_save, sep="\t", header=True, index=False)


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="get syntelog table names fixed")

    parser.add_argument(
        "syntelog_input_file",
        type=str,
        help="parent path to syntelog file",
    )

    parser.add_argument(
        "filename_to_save", type=str, help="name to save the syntelog file as"
    )

    parser.add_argument(
        "output_dir",
        type=str,
        help="parent directory to output results",
    )

    args = parser.parse_args()
    args.syntelog_input_file = os.path.abspath(args.syntelog_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    import_syntelogs(args.syntelog_input_file, args.output_dir, args.filename_to_save)
