#!python
import sys
import os
import argparse
import pathlib

import cubit


def add_element_blocks_to_mesh(input_mesh, output_mesh, elements):
    '''Take an existing exodus mesh, add element blocks for each element, save with new name

    :param str input_mesh: The input exodus mesh file to modify
    :param str output_mesh: The output exodus mesh file with block names defined
    :param int elements: The number of elements in the mesh for which to define a block name

    :returns: Write ``output_mesh``
    '''

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')
    cubit.cmd(f'import mesh geometry "{input_mesh}"')

    # Delete any existing blocks
    cubit.cmd('reset block')

    element_list = list(range(0, elements))
    for i in element_list:
        cubit.cmd(f'block {i+1} add hex {i+1}')
        cubit.cmd(f'block {i+1} name "element_{i}"')

    cubit.cmd(f'export mesh "{output_mesh}"  overwrite')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "

    cli_description = "Take an existing exodus mesh, add element blocks for each element, save with new name"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--input-mesh', type=str, required=True,
        help="The input exodus mesh file to modify")
    parser.add_argument('--output-mesh', type=str, required=True,
        help="The output exodus mesh file with block names defined")
    parser.add_argument('--elements', type=int, required=True,
        help="The number of elements in the mesh for which to define a block name")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(add_element_blocks_to_mesh(input_mesh=args.input_mesh,
                                        output_mesh=args.output_mesh,
                                        elements=args.elements,
                                        ))
