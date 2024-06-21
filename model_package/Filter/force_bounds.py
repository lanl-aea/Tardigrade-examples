import sys
import os
import argparse
import inspect
import pathlib

import numpy
import pandas

import xdmf_reader_tools as XRT


def force_bounds(output_file, xmin, xmax, ymin, ymax, zmin, zmax):
    '''Create a csv file containing information for a bounding box encompassing all DNS points

    :param str output_file: The name of the output csv file
    :param float xmin: The minimum x-value
    :param float xmax: The maximum x-value
    :param float ymin: The minimum y-value
    :param float ymax: The maximum y-value
    :param float zmin: The maximum z-value
    :param float zmax: The maximum z-value

    :returns: Write ``output_file``
    '''

    df = pandas.DataFrame({'xmin': [xmin],
                           'xmax': [xmax],
                           'ymin': [ymin],
                           'ymax': [ymax],
                           'zmin': [zmin],
                           'zmax': [zmax]})
    df.to_csv(output_file, header=True, sep=',', index=False)

    return 0


def get_parser():

    basename = pathlib.Path(__file__).name
    cli_description = "Create a csv file containing information for a bounding box encompassing all DNS points"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=basename)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='The name of the output csv file of bounding informaiton')
    parser.add_argument('--xmin', type=float, required=True,
        help='The minimum x-value')
    parser.add_argument('--xmax', type=float, required=True,
        help='The maximum x-value')
    parser.add_argument('--ymin', type=float, required=True,
        help='The minimum y-value')
    parser.add_argument('--ymax', type=float, required=True,
        help='The maximum y-value')
    parser.add_argument('--zmin', type=float, required=True,
        help='The minimum z-value')
    parser.add_argument('--zmax', type=float, required=True,
        help='The maximum z-value')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(force_bounds(output_file=args.output_file,
                          xmin=args.xmin,
                          xmax=args.xmax,
                          ymin=args.ymin,
                          ymax=args.ymax,
                          zmin=args.zmin,
                          zmax=args.zmax,
                          ))