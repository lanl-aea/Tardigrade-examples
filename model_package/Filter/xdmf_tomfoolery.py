import sys
import os
import inspect
import argparse

import numpy
import meshio

import file_io.xdmf

file_path = os.path.dirname(os.path.abspath(__file__))


def XDMF_tomfoolery(input_file, output_file):
    '''Modify an XDMF file by combining elements from separate 'blocks'

    :param str input_file: The XDMF mesh file to operate on
    :param str output_file: The output filename for the h5 + XDMF file pair

    :returns: ``{output_file}.xdmf``
    '''

    mesh=meshio.read(input_file)
    filter_points = mesh.points
    # manipulate blocks and stack into single connectivity array
    num_blocks = len(mesh.cells)
    connect = []
    for i in range(num_blocks):
        connect.append(mesh.cells[i].data)
    connect = numpy.array(connect)
    tup = numpy.shape(connect)
    filter_connectivity = connect.reshape((tup[0]*tup[1],tup[2]))

    # Write the filter to a file
    xdmf = file_io.xdmf.XDMF(output_filename=output_file)

    grid = xdmf.addGrid(xdmf.output_timegrid, {})
    xdmf.addPoints(grid, filter_points)

    xdmf.addConnectivity(grid, "HEXAHEDRON", filter_connectivity)

    print("Writing single filter domain file!")
    xdmf.write()
    print("filter domain file written!")

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Modify an XDMF file by combining elements from separate 'blocks'"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='Specify the output filename for the h5 + XDMF file pair')
    parser.add_argument('--input-file', type=str, required=True,
        help='Specify the XDMF mesh file to operate on')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args = parser.parse_args()
    sys.exit(XDMF_tomfoolery(input_file=args.input_file,
                             output_file=args.output_file,
                             ))