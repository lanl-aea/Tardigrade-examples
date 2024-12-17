import sys
import os
import inspect
import argparse
import yaml

import numpy


def parse_sets_from_inp(input_file, output_file):
    """Extract element IDs associated with element sets from Abaqus input file

    """

    # Read input file
    with open(input_file, 'r') as f:
        input_contents = numpy.asarray(f.read().splitlines())

    # tag line numbers
    indices = []
    set_names = []
    for index, text in enumerate(input_contents):
        # store lines including element set information
        if ("*ELEMENT" in text) and ("ELSET" in text):
            set_names.append(text.split("=")[-1])
            indices.append(index)
        # store index of line corresponding to the end of the element set specification
        if "N O D E S E T S" in text:
            end = index

    # construct dictionary of start and end indices to search of elements in each set
    set_search = {}
    for position, set_name in enumerate(set_names):
        if position == len(indices)-1:
            set_search[set_name] = [indices[position]+1, end]
        else:
            set_search[set_name] = [indices[position]+1, indices[position+1]]

    # extract elements, we only need the element id specified at the beginning of each line
    sets = {}
    for key in set_search.keys():
        elements = []
        start, stop = set_search[key][0], set_search[key][1]
        for line in input_contents[start:stop]:
            elements.append(int(line.split(',')[0]))
        sets[key] = elements

    # Output
    with open(output_file, 'w') as f:
        yaml.dump(sets, f)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Extract element IDs associated with element sets from Abaqus input file"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
        help='Specify the Abaqus input file to containing element set information')
    parser.add_argument('-o', '--output-file', type=str,
        help='Specify the output yaml file name containing element set information')


    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(parse_sets_from_inp(input_file=args.input_file,
                                 output_file=args.output_file,
                                 ))