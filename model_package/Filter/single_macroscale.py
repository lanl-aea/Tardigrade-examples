import sys
import os
import inspect
import argparse

import numpy as np
import pandas

import file_io.xdmf

file_path = os.path.dirname(os.path.abspath(__file__))


def single_domain(X1, X2, Y1, Y2, Z1, Z2, output_file):
    ''' Creates an XDMF file containing a single element

    :param float X1: The minimum X value
    :param float X2: The maximum X value
    :param float Y1: The minimum Y value
    :param float Y2: The maximum Y value
    :param float Z1: The minimum Z value
    :param float Z2: The maximum Z value
    :param str output_file: The output filename for the h5 + XDMF file pair

    :returns: ``{output_file}.h5`` and ``{output_file}.xdmf``
    '''

    filter_points = np.array([[X1, Y1, Z1],
                              [X2, Y1, Z1],
                              [X2, Y2, Z1],
                              [X1, Y2, Z1],
                              [X1, Y1, Z2],
                              [X2, Y1, Z2],
                              [X2, Y2, Z2],
                              [X1, Y2, Z2]])

    filter_connectivity = np.array([0, 1, 2, 3, 4, 5, 6, 7,])
    filter_connectivity = filter_connectivity.reshape((1,-1)).astype(int)

    # Write the filter to a file
    xdmf = file_io.xdmf.XDMF(output_filename=output_file)

    grid = xdmf.addGrid(xdmf.output_timegrid, {})
    xdmf.addPoints(grid, filter_points)

    xdmf.addConnectivity(grid, "HEXAHEDRON", filter_connectivity)

    print("Writing single filter domain file!")
    xdmf.write()
    print("filter domain file written!")

    return 0


def write_filter_domain(output_file, single_points=None, csv_file=None):
    '''Write a single macroscale domain file for the Micromorphic Filter

    :param str output_file: The output filename for the h5 + XDMF file pair
    :param list single_points: The X, Y, and Z extens for a single element filter domain
    :param str csv_file: CSV filename containing the bounds of a DNS file

    Calls single_domain function to generate XDMF file
    '''

    if single_points:
        X1, X2, Y1, Y2, Z1, Z2 = (float(item) for item in args.single_points)
        single_domain(X1, X2, Y1, Y2, Z1, Z2, output_file)
    elif csv_file:
        csv_data = pandas.read_csv(csv_file, sep=",")
        X1, X2 = csv_data['xmin'][0], csv_data['xmax'][0]
        Y1, Y2 = csv_data['ymin'][0], csv_data['ymax'][0]
        Z1, Z2 = csv_data['zmin'][0], csv_data['zmax'][0]
        single_domain(X1, X2, Y1, Y2, Z1, Z2, output_file)
    else:
        print('Specify either "single_points" or "csv_file" argument!')
        print('No macro domain file will be written')

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Write a single macroscale domain file for the\
                       Micromorphic Filter"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='Specify the output filename for the h5 + XDMF file pair')
    parser.add_argument('--single-points', nargs=6,
        help='Specify the X, Y, and Z extents for the a single element\
              macro domain')
    parser.add_argument('--csv-file', type=str,
        help='Specify a csv file containing the bounds of a DNS file')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args = parser.parse_args()
    sys.exit(write_filter_domain(output_file=args.output_file,
                                 single_points=args.single_points,
                                 csv_file=args.csv_file,
                                 ))