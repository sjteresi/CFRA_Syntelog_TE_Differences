"""
Filter a EDTA-created TE annotation to the appropriate format for TE
Density algorithm
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs

from src.H4.replace_names_H4_strawberry import te_annot_renamer


def check_nulls(my_df, logger):
    """Check the TE dataframe for ANY null values in ANY rows

    Args:
        my_df (pandas.core.DataFrame): Pandas dataframe of TE values from TE
            annotation
    """
    Bool = my_df.isnull().values.any()
    if Bool:
        logger.critical("You have null values in your dataframe!")
        logger.critical("Here are the null values in the output:")
        null_columns = my_df.columns[my_df.isnull().any()]
        print((my_df[my_df.isnull().any(axis=1)][null_columns].head()))


def write_cleaned_transposons(te_pandaframe, output_dir, genome_name, logger):
    file_name = os.path.join(output_dir, ("Cleaned_" + genome_name + "_EDTA_TEs.tsv"))

    logger.info("Writing cleaned TE file to: %s" % file_name)
    te_pandaframe.to_csv(file_name, sep="\t", header=True, index=False)


def import_transposons(tes_input_path, te_annot_renamer, genome_name, logger):
    """Import TE file and read as a dataframe in Pandas

    Args:
        tes_input_path (str): string of the file path to the TE annotation

        logger (logging obj): The object to call logs and info
    """
    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Score",
        "Strand",
        "Phase",
        "Attribute",
    ]

    TE_Data = pd.read_csv(
        tes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        comment="#",
        dtype={
            "Chromosome": str,
            "Software": str,
            "Feature": str,
            "Start": "float64",
            "Stop": "float64",
            "Score": str,
            "Strand": str,
            "Phase": str,
            "Attribute": str,
        },
    )

    # Drop extraneous columns
    TE_Data.drop(
        columns=["Score", "Strand", "Software", "Attribute", "Phase"], inplace=True
    )

    # Create Order and SuperFamily column from Attribute column
    # Because that column contains the detailed TE information
    # Then remove old Attribute column
    # TE_Data["Attribute"] = TE_Data["Attribute"].str.extract(r"Classification=(.*?);")

    TE_Data[["Order", "SuperFamily"]] = TE_Data.Feature.str.split("/", expand=True)
    TE_Data.drop(columns=["Feature"], inplace=True)
    TE_Data.Order = TE_Data.Order.astype(str)
    TE_Data.SuperFamily = TE_Data.SuperFamily.astype(str)

    # Call renamer
    TE_Data = te_annot_renamer(TE_Data)

    # Declare data types
    TE_Data["Length"] = TE_Data.Stop - TE_Data.Start + 1
    check_nulls(TE_Data, logger)

    TE_Data.sort_values(by=["Chromosome", "Start"], inplace=True)

    # MAGIC I only want the first 7
    # MAGIC chromosome/scaffold names
    chromosomes_i_want = ["Fvb" + str(i) for i in range(8)]
    TE_Data = TE_Data.loc[TE_Data["Chromosome"].isin(chromosomes_i_want)]

    diagnostic(TE_Data)
    print(TE_Data)

    return TE_Data


def diagnostic(TE_Data):
    print(sorted(TE_Data["Order"].unique().tolist()))
    print(sorted(TE_Data["SuperFamily"].unique().tolist()))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reformat TE annotation file")
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(dir_main, "../results")
    parser.add_argument(
        "TE_input_file", type=str, help="Parent path of TE annotation file"
    )
    parser.add_argument("genome_name", type=str, help="Name for the genome")

    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=output_default,
        help="Parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.TE_input_file = os.path.abspath(args.TE_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Execute
    cleaned_transposons = import_transposons(
        args.TE_input_file, te_annot_renamer, args.genome_name, logger
    )
    write_cleaned_transposons(
        cleaned_transposons, args.output_dir, args.genome_name, logger
    )
