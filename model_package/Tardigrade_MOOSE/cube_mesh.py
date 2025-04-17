#!python
import argparse
import os
import pathlib
import sys

import cubit


def cube_mesh(side_length, seed_size, output_file, x0=0., y0=0., z0=0.):
    '''Create a cube mesh

    :param float side_length: Cube side length
    :param float seed_size: The approximate mesh size
    :param str output_file: The output filename
    :param float x0: The x-distance to translate the cube
    :param float y0: The y-distance to translate the cube
    :param float z0: The z-distance to translate the cube

    :returns: ``{output_file}.e``
    '''

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')

    # Geometry and move
    cubit.cmd(f'brick x {side_length} y {side_length} z {side_length}')
    cubit.cmd(f'move Volume all x {x0} y {y0} z {z0} include_merged')

    # Mesh
    cubit.cmd(f'volume all size {seed_size}')
    cubit.cmd('mesh volume all')

    # Export
    cubit.cmd(f'export mesh "{output_file}.e"  overwrite')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Create a cube mesh"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)

    parser.add_argument('--side-length', type=float, required=True,
        help='Cube side length')
    parser.add_argument('--seed-size', type=float, required=True,
        help='The approximate mesh size')
    parser.add_argument('--output-file', type=str, required=True,
        help="The output filename")
    parser.add_argument('--x0', type=float, required=False, default=0.,
        help='The x-distance to translate the cube')
    parser.add_argument('--y0', type=float, required=False, default=0.,
        help='The y-distance to translate the cube')
    parser.add_argument('--z0', type=float, required=False, default=0.,
        help='The z-distance to translate the cube')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(cube_mesh(side_length=args.side_length,
                       seed_size=args.seed_size,
                       output_file=args.output_file,
                       x0=args.x0,
                       y0=args.y0,
                       z0=args.z0
                       ))

