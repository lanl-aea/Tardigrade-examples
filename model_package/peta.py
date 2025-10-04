#! /usr/bin/env python

import sys
import os
import pathlib
import argparse
import subprocess

import xml.etree.ElementTree as ET

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

        # Handle destination directory
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


def remove_block_by_name(parent, block_name_to_remove):
    ''' Recursively search through an XML node to find and remove a specified block

    :params xml_node parent: The XML node to start search
    :params string block_name_to_remove: The name attribute of an XML block to remove
    '''

    # Iterate through parent node
    for child in list(parent):
        if child.tag == "Block" and child.attrib.get("name") == block_name_to_remove:
            parent.remove(child)
        else:
            # Recursively search children
            remove_block_by_name(child, block_name_to_remove)


def transfer_GEOS_files_selective(source_files_dict, source_directory, output_directory, split_key, username):
    '''Transfer GEOS DNS ascii VTK results and subdirectories a series of scp commands (multiple sources, multiple destinations)
    for a specific selection of timesteps and modify associated XML data of PVD and VTM files

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

        # Handle destination directory
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

        # Only initiate file transfer if there are new files to copy
        if len(sources) > 0:
            print(sources[0])
            # Transfer files
            if len(sources) == 1:
                file_string = f"{username}@login.rc.colorado.edu:{sources[0]}"
            else:
                file_string = f"{username}@login.rc.colorado.edu:" + '{' + f"{','.join(sources)}" + '}'
            full_command = ["scp", "-P", "22", file_string, dest_dir]
            p = subprocess.Popen(full_command)
            os.waitpid(p.pid, 0)

        # Trim main pvd file to remove unnecessary timesteps
        if key == 'main':
            print(source_files)
            pvd_file = f"{output_directory}/{source_files[0].split('/')[-1]}"
            pvd_output_file = f"{output_directory}/vtkOutput_trimmed.pvd"
            # Always rewrite the trimmed pvd file
            if os.path.isfile(pvd_output_file) == True:
                os.remove(pvd_output_file)
            # Find the target vtm files to keep in pvd file
            target_files = [file.split('/')[-1] for file in source_files_dict['vtms']]
            # parse the XML
            tree = ET.parse(pvd_file)
            root=tree.getroot()
            # Find collection element
            collection=root.find('Collection')
            # Remove all DataSet elements not in target_files
            for dataset in list(collection):
                vtm_file = dataset.attrib['file'].split('/')[-1]
                print(vtm_file)
                if vtm_file not in target_files:
                    collection.remove(dataset)
                else:
                    dataset.set('file', f'vtkOutput/{vtm_file.split('.')[0]}_trimmed.vtm')
            # Write the filtered xml file
            print(f'writing {pvd_output_file}!!!')
            tree.write(pvd_output_file, encoding='utf-8', xml_declaration=True)

        # Trim vtm files to remove backgroundGrid data
        if key == 'vtms':
            print(source_files)
            for vtm_file in source_files:
                vtm_filename = vtm_file.split('/')[-1]
                vtm_output_name = f"{output_directory}/vtkOutput/{vtm_filename.split('.')[0]}_trimmed.vtm"
                # Always rewrite the trimmed vtm file
                if os.path.isfile(vtm_output_name) == True:
                    os.remove(vtm_output_name)
                # parse the XML
                vtm = f"{output_directory}/vtkOutput/{vtm_filename}"
                tree = ET.parse(vtm)
                root = tree.getroot()
                # Remove all backgroundGrid
                for multiblock in root.findall(".//vtkMultiBlockDataSet"):
                    # Do not keep backgroundGrid
                    remove_block_by_name(multiblock, "backgroundGrid")
                    # Do not keep platen particles
                    remove_block_by_name(multiblock, "ParticleRegion1")
                # Write the filtered xml file
                print(f'writing {vtm_output_name}!!!')
                tree.write(vtm_output_name, encoding='utf-8', xml_declaration=True)

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
                   Ratel_variables.additional_files + \
                   Ratel_variables.I43_damage_coarse_finetime['DNS_files'] + \
                   [Ratel_variables.I43_damage_coarse_finetime['DNS_forces']] + \
                   Ratel_variables.I43_damage_CGNS['DNS_files'] +\
                   Ratel_variables.I43_best_DNS['DNS_files'] +\
                   [Ratel_variables.I43_best_DNS['DNS_forces']]

    # Transfer the main chunk of results for Ratel
    transfer_Ratel_files(source_files, source_directory, output_directory, username)

    # Transfer Ratel I43.09 MAP results
    source_files = Ratel_variables.I43_MAP['DNS_files'] + [Ratel_variables.I43_MAP['DNS_forces']] + \
                   [Ratel_variables.I43_MAP['input_file']]
    MAP_directory = output_directory / "Ratel_I43_MAP"
    transfer_Ratel_files(source_files, source_directory, MAP_directory, username)

    # GEOS transfer for elastic cylinder
    source_files_dicts = (('GEOS_elastic_cylinder', 'UniaxialCylinderCompression/', GEOS_variables.elastic_cylinder['DNS_files']),
                           )
    for top_directory, split_key, source_files_dict in source_files_dicts:
        GEOS_output_directory = f'{output_directory}/{top_directory}'
        print(GEOS_output_directory)
        transfer_GEOS_files(source_files_dict, source_directory, GEOS_output_directory, split_key, username)

    # GEOS transfer for larger DNS
    source_files_dicts = (('GEOS_I43_01_sim38', 'lowRes_ceramicPuck_I43p01_sim38/', GEOS_variables.I43_01_sim38['DNS_files']),
                           )
    for top_directory, split_key, source_files_dict in source_files_dicts:
        GEOS_output_directory = f'{output_directory}/{top_directory}'
        print(GEOS_output_directory)
        transfer_GEOS_files_selective(source_files_dict, source_directory, GEOS_output_directory, split_key, username)

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
