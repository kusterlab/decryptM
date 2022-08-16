# pg_transform.py
#
# This script transforms the proteinGroups.txt file such that it can be fed into the decryptM pipeline.
# the output file will be placed next to the input file with the suffix restructured.
#
# Run this script as follows:
# $ python ./decryptM/pg_transform.py <in_path>
#
# Florian P. Bayer - Aug. 2022
#
# version 1.0
#

# ####### #
# Imports #
# ####### #

import sys
import os
import argparse

import numpy as np
import pandas as pd


# ######### #
# Functions #
# ######### #
def check_path(path, create_folder=False):
    """
    Check the given path whether it exists and is accessible.

    if create_folder is True and folder dose not exist yet it will be created
    """
    if not os.path.exists(path):
        if create_folder:
            os.mkdir(path)
        else:
            raise FileNotFoundError(f'Path "{path}" does not exist.')
    if not os.access(path, os.R_OK):
        raise PermissionError(f'File at "{path}" cannot be opened. Try to close it elsewhere.')


def parse_tmt_channel(s):
    s_lst = s.split(" ")
    return " ".join(s_lst[:4])


def parse_experiment(s):
    s_lst = s.split(" ")
    return s_lst[-1]


def transform_fullproteome(in_path, out_path):
    """
    Is doing the transformation
    """
    # Load the file
    df = pd.read_csv(in_path, sep='\t', low_memory=False)

    # Overwrite protein names with IDs
    df.loc[df['Protein names'].isna(), 'Protein names'] = df['Protein IDs']

    # Define id columns and get reporter columns
    id_cols = ['Protein names', 'Gene names', "Protein IDs"]
    reporter_cols = list(df.columns[df.columns.str.contains(r'Reporter intensity corrected \d+ .+', regex=True)])

    # Restructure 
    tmt_intensity = df[id_cols + reporter_cols].melt(value_vars=reporter_cols, id_vars=id_cols,
                                                     value_name='Intensity', var_name='Channel_Experiment')

    # Parse channel and experiment information
    parsed_cols = ['Experiment', 'Channel']
    tmt_intensity['Channel'] = tmt_intensity['Channel_Experiment'].apply(parse_tmt_channel)
    tmt_cols = sorted(tmt_intensity['Channel'].unique())
    tmt_intensity['Experiment'] = tmt_intensity['Channel_Experiment'].apply(parse_experiment)

    # Sum total intensity
    tmt_intensity = tmt_intensity.set_index(id_cols + parsed_cols)['Intensity'].unstack().reset_index()
    tmt_intensity['Intensity'] = tmt_intensity[tmt_cols].sum(axis=1)

    # Add a mock mod sequence column to be compatible with the evidence based script which requires this
    tmt_intensity['Modified sequence'] = tmt_intensity['Protein IDs']

    # Add extra information from the input file to the restructured file
    keep_cols = ['Protein IDs', 'Reverse', 'Potential contaminant']
    tmt_intensity = pd.merge(left=tmt_intensity, right=df[keep_cols], on='Protein IDs', how='left')

    # export file
    tmt_intensity.to_csv(out_path, sep='\t', index=False)
    print('Writing file to:', out_path, sep='\n')


# ###### #
# Script #
# ###### #

if __name__ == '__main__':
    # Add a command line parser for the config file
    parser = argparse.ArgumentParser(
        description='Restructure the proteinGroups.txt file',
        epilog="This script transforms the proteinGroups.txt file such that it can be fed into the decryptM pipeline.\nthe output file will be placed next to the input file with the suffix restructured.\n\nFPB-2022"
    )
    parser.add_argument(
        dest="in_path",
        metavar="in path",
        type=str,
        help="Path to the proteinGroups.txt file",
    )

    # Check the input and define the output
    args = parser.parse_args()
    check_path(args.in_path)

    # Extract the input
    in_file_folder = os.path.dirname(args.in_path)
    in_file_name = os.path.basename(args.in_path).split(".")
    assert len(in_file_name) == 2

    # Define the output
    out_file_name = f'{in_file_name[0]}_restructured.{in_file_name[1]}'
    out_path = os.path.join(in_file_folder, out_file_name)

    # Run script
    transform_fullproteome(args.in_path, out_path)
