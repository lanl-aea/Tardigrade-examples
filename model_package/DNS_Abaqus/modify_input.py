import numpy as np
import os
import sys
import argparse
import inspect


def modify_input(input_file, output_file):
    '''Modify Abaqus input file to output 'COORD' at integration points

    :param str input_file: Relative or absolute path to Abaqus input (.inp) file.
    :param str output_file: Relative or absolute path to modified Abaqus input (.inp) file.
    '''

    # read in data from original input file, replace relevant lines
    data = []
    with open(input_file) as file:
        for line in file:
            if 'COORD, U' in line:
                new_line = line.replace('COORD, ', '')
                data.append(new_line)
            elif ('EVOL, IVOL, S' in line) and ('COORD' not in line):
                new_line = line + ', COORD\n'
                data.append(new_line)
            else:
                data.append(line)

    # Write newly formatted input file
    fin = open(output_file, 'wt')
    for line in data:
        fin.write(line)
    fin.close()

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Modify Abaqus input file to output 'COORD' at integration points"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
                        help="The Abaqus input file created by ``build_model.py``. ")
    parser.add_argument('-o', '--output-file', type=str, 
                        help="The modified Abaqus input file")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(modify_input(input_file=args.input_file,
                          output_file=args.output_file))