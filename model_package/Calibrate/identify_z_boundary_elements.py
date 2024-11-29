#!python
import argparse
import pathlib
import sys

import numpy
import pandas

import file_io.xdmf


def identify_z_boundary_elements(macro_file, output_file):
    '''Read in macroscale XDMF file of a cylindrical geometry and identify element found on the z-boundary

    :param str macro_file: The macroscale filter domain XDMF file, less extension
    :param str output_file: Output csv filename containing list of boundary elements

    :returns: ``output_file``
    '''

    # open XDMF file
    macro = file_io.xdmf.XDMF(macro_file)
    macro.open()

    # Get nodes and element connectivity
    nodes = macro.getIncrementReferenceNodePositions(0)[0][0]
    elements = macro.getIncrementConnectivity(0)[0][0]

    # find boundaries from extent of node positions in the z-direcion
    min = numpy.min(nodes,axis=0)[2]
    max = numpy.max(nodes,axis=0)[2]

    # append elements to a list if they have nodes on the z-boundary
    boundary_elements = []
    for i, e in enumerate(elements):
        flag = False
        for n in e:
            x, y, z = nodes[n]
            if (z <= min) or (z >= max):
                flag = True
        if flag:
            boundary_elements.append(i)

    # Write csv file containing boundary elements
    boundary_df = pandas.DataFrame({'boundary_elements': boundary_elements})
    boundary_df.to_csv(output_file, sep=',', index=False)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Read in macroscale XDMF file of a cylindrical geometry and identify element found on the z-boundary"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--macro-file', type=str, required=True,
        help='The macroscale filter domain XDMF file, less extension')
    parser.add_argument('--output-file', type=str, required=True,
        help='Output csv filename containing list of boundary elements')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(identify_z_boundary_elements(macro_file=args.macro_file,
                                          output_file=args.output_file,
                                          ))