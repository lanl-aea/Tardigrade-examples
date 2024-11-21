#! /usr/bin/env python

import sys
import os
import pathlib
import argparse
import subprocess

from model_package.DNS_Ratel import simulation_variables_nominal as Ratel_variables
from model_package.DNS_GEOS import simulation_variables_nominal as GEOS_variables


def transfer_Ratel_files(source_files, source_directory, output_directory, username):
    '''Transfer all Ratel DNS binary VTK results files to the output_directory using a single scp command (multiple sources, single destination)

    :params list source_files: The source files to copy
    :params str source_directory: The common CU Peta Library DNS results path
    :params str output_directory: The destination directory of DNS files
    :params str username: The user's CU identikey

    '''

    # Create output directory if it doesn't exist
    if os.path.exists(output_directory) == False:
        os.makedirs(output_directory)

    # Build single command to copy all files at once
    sources = []
    for file in source_files:
        dest_file = f"{output_directory}/{file.split('/')[-1]}"
        if os.path.isfile(dest_file) == False:
            sources.append(f'{source_directory}/{file}')
        else:
            print(f'{dest_file} already exists!')

    # Only initiate file transfer if there are new files to copy
    if len(sources) == 0:
        print('No files to copy for Ratel')
        return 0
    else:
        # Transfer files
        file_string = f"{username}@login.rc.colorado.edu:" + '{' + f"{','.join(sources)}" + '}'
        full_command = ["scp", "-P", "22", file_string, f'{output_directory}/.']
        p = subprocess.Popen(full_command)
        os.waitpid(p.pid, 0)

    return 0


def transfer_GEOS_files(source_files_dict, source_directory, output_directory, split_key, username):
    '''Transfer GEOS DNS ascii VTK results and subdirectories a series of scp commands (multiple sources, multiple destinations)

    :params dict source_files_dict: The a dictionary containing lists of files GEOS ascii VTK results to copy
    :params str source_directory: The common CU Peta Library DNS results path
    :params str output_directory: The root destination directory of DNS files
    :params str split_key: a unique string to split the full Petalibrary path for copying
    :params str username: The user's CU identikey

    '''

    # Create output directory if it doesn't exist
    if os.path.exists(output_directory) == False:
        os.makedirs(output_directory)

    # Loop over the individual keys to copy files
    for key in source_files_dict.keys():
        source_files = source_files_dict[key]

        # Hand destination directory
        dest_split = f"{source_files[0].split(split_key)[-1]}"
        if len(dest_split.split('/')) == 1:
            dest_dir = f"{output_directory}/."
        else:
            dest_dir = f"{output_directory}/{'/'.join(dest_split.split('/')[0:-1])}/."
        if os.path.exists(dest_dir.split('.')[0]) == False:
            os.makedirs(dest_dir.split('.')[0])

        # Build single command to copy all relevant files at once
        sources = []
        for file in source_files:
            dest_file = f"{output_directory}/{file.split(split_key)[-1]}"
            if os.path.isfile(dest_file) == False:
                sources.append(f'{source_directory}/{file}')
            else:
                print(f'{dest_file} already exists!')

        # Only initiate file transfer if there are new files to copy
        if len(sources) == 0:
            continue
        else:
            # Transfer files
            file_string = f"{username}@login.rc.colorado.edu:" + '{' + f"{','.join(sources)}" + '}'
            full_command = ["scp", "-P", "22", file_string, dest_dir]
            p = subprocess.Popen(full_command)
            os.waitpid(p.pid, 0)


    return 0


def peta_copy(source_directory, output_directory):
    '''Copy DNS results from the CU Peta library to the output directory

    :param str source_directory: The source directory of DNS simulation results
    :param str output_directory: The output directory destination
    '''

    #Get user's CU identikey
    username = input("Enter your CU identikey: ")

    # Ratel transfer
    source_files = Ratel_variables.I41_02['DNS_files'] + \
                   [Ratel_variables.I41_02['DNS_forces']] + \
                   Ratel_variables.F83['DNS_files'] + \
                   Ratel_variables.additional_files
                   #[Ratel_variables.I43_damage['DNS_forces'][0]] + \
                   #Ratel_variables.I43_damage['DNS_files']
    transfer_Ratel_files(source_files, source_directory, output_directory, username)

    # GEOS transfer
    source_files_dicts = (('GEOS_elastic_cylinder', 'UniaxialCylinderCompression/', GEOS_variables.elastic_cylinder['DNS_files']),
                           )
    for top_directory, split_key, source_files_dict in source_files_dicts:
        GEOS_output_directory = f'{output_directory}/{top_directory}'
        print(GEOS_output_directory)
        transfer_GEOS_files(source_files_dict, source_directory, GEOS_output_directory, split_key, username)


    return 0


def get_parser():

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
