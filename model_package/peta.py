#! /usr/bin/env python

import sys
import os
import pathlib
import argparse
import subprocess

from model_package.DNS_Ratel import simulation_variables_nominal


def peta_copy(source_directory, output_directory):
    '''Copy DNS results from the CU Peta library to the output directory

    :param str source_directory: The source directory of DNS simulation results
    :param str output_directory: The output directory destination
    '''

    source_files = simulation_variables_nominal.I41_02['DNS_files'] + \
                   [simulation_variables_nominal.I41_02['DNS_forces']] + \
                   simulation_variables_nominal.F83['DNS_files'] + \
                   simulation_variables_nominal.additional_files

    # Create source directory if it doesn't exist
    if os.path.exists(output_directory) == False:
        os.makedirs(output_directory)

    # Get user's CU username
    username = input("Enter your CU identikey: ")

    # Build single command to copy all files at once
    sources = []
    for file in source_files:
        dest_file = f"{output_directory}/{file.split('/')[-1]}"
        if os.path.isfile(dest_file) == False:
            sources.append(f'{source_directory}/{file}')
        else:
            print(f'{output_directory}/{dest_file} already exists!')
    file_string=f"{username}@login.rc.colorado.edu:" + '{' + f"{','.join(sources)}" + '}'

    # Transfer files
    full_command = ["scp", "-P", "22", file_string, f'{output_directory}/.']
    p = subprocess.Popen(full_command)
    sts = os.waitpid(p.pid, 0)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)

    cli_description = "Copy DNS results from the CU Peta library to the output directory"
    parser = argparse.ArgumentParser(description=cli_description)
    parser.add_argument('--source-directory', type=str, required=True,
        help="The source directory of DNS simulation results")
    parser.add_argument('--output-directory', type=pathlib.Path, required=True,
        help="The output directory destination.")

    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(peta_copy(source_directory=args.source_directory,
                       output_directory=args.output_directory,
                       ))