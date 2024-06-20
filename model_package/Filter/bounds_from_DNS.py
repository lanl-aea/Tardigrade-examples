import sys
import os
import argparse
import inspect
import pathlib

import numpy
import pandas

import xdmf_reader_tools as XRT


def bounds_from_DNS(DNS_file, output_file, coord='coord'):
    '''Create a csv file containing the extents of a DNS results file

    :param str DNS_file: The name of the input XDMF file containing DNS results
    :param str output_file: The name of the output csv file of bounding information
    :param str coord: The name of the coordinate field

    :returns: ``output_file``
    '''

    # load in xdmf data
    data, geometry, topology = XRT.parse_xdmf_output(DNS_file)

    # Find bounds
    xmin = numpy.min(data[f'{coord}_0'][:,0])
    xmax = numpy.max(data[f'{coord}_0'][:,0])
    ymin = numpy.min(data[f'{coord}_0'][:,1])
    ymax = numpy.max(data[f'{coord}_0'][:,1])
    zmin = numpy.min(data[f'{coord}_0'][:,2])
    zmax = numpy.max(data[f'{coord}_0'][:,2])
    
    # Output to labeled csv file
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
    cli_description = "Create a csv containing the extents of a DNS file"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=basename)
    parser.add_argument('-d', '--dns-file', type=str, required=True,
        help='The name of the input XDMF file containing DNS results')
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='The name of the output csv file of bounding information')
    parser.add_argument('--coord', type=str, default='coord',
        help='The name of the coordinate field')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(bounds_from_DNS(DNS_file=args.dns_file,
                             output_file=args.output_file,
                             coord=args.coord,
                             ))