#! /usr/bin/env python

import sys
import os
import argparse
import pathlib
import yaml

import subprocess

def config_software(config_file):
    '''Configure software paths in a YAML file

    :param str config_file: The YAML file to write software paths

    :returns: Writes or updates ``config_file``
    '''

    # Check if config file exists
    if os.path.exists(config_file):
        print(f'Config file found {config_file}!')
        stream = open(config_file, 'r')
        program_paths = yaml.load(stream, Loader=yaml.FullLoader)
        stream.close()
    else:
        new_config_file_query = str(input(f'Make new config file at {config_file} (y/n)?: '))
        if new_config_file_query == 'y':
            print('WARNING! To use this new config file, the SConstruct must be modified to use it')
            print('Generating new config file')
            program_paths = {'Abaqus': [''],
                             'Cubit': [''],
                             'Ratel': [''],
                             'filter': [''],
                             'micromorphic': [''],
                             'constraints': [''],
                             'Tardigrade': [''],
                             'LD_PATH': [''],
                             'mpi': [''],
                             'paraview': [''],
                             }

    # Ask for new program paths
    for program in program_paths.keys():
        if os.path.exists(program_paths[program][-1]):
            print(f'{program} program found')
            add_program_query = str(input("\tAdd a new program path (y/n)?: "))
            if add_program_query == 'y':
                new_program_path = str(input(f"\tSpecify additional path for {program}: "))
                if os.path.exists(new_program_path) == False:
                    print('\tProgram not found!')
                else:
                    print('\tProgram found! Adding to {config_file}')
                    program_paths[program].append(new_program_path)
        else:
            print(f'{program} program not found')
            add_program_query = str(input("\tAdd a new program path (y/n)?: "))
            if add_program_query == 'y':
                new_program_path = str(input(f"\tSpecify path for {program}: "))
                if os.path.exists(new_program_path) == False:
                    print('\tProgram not found!')
                else:
                    print('\tProgram found! Adding to {config_file}')
                    if len(program_paths[program][-1]) > 1:
                        program_paths[program].append(new_program_path)
                    else:
                        program_paths[program] = [new_program_path]

    # Write the config file
    with open(config_file, 'w') as f:
        yaml.dump(program_paths, f)

    # Force git to ignore changes to config.yml
    subprocess.run([f'git update-index --assume-unchanged {config_file}'], shell=True)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)

    cli_description = "Configure software paths in a YAML file"
    parser = argparse.ArgumentParser(description=cli_description)
    parser.add_argument('--config-file', type=str, required=True,
        help="The YAML file to write software paths")

    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(config_software(config_file=args.config_file))