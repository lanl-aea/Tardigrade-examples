#!python
import argparse
import os
import pathlib
import sys

import cubit
import pandas
import numpy
import subprocess


def mesh(rad, height, annulus_ratio, x0, y0, z0, seed_size, output_file):
    ''' Mesh a cylinder using Cubit

    :param float rad: Cylinder radius
    :param float height: Cylinder height
    :param float annulus_ratio: The fraction of the radius to keep in the final geometry
    :param float x0: The x-distance to translate the cylinder
    :param float y0: The y-distance to translate the cylinder
    :param float z0: The z-distance to translate the cylinder
    :param float seed-size: The approximate mesh size
    :param str output_file: The output filename

    :returns: ``{output_file}.e``
    '''

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')
    cubit.cmd(f'create Cylinder height {height} radius {rad}')
    rad2 = (1. - annulus_ratio)*rad
    cubit.cmd(f'create Cylinder height {height} radius {rad2}')
    cubit.cmd('subtract volume 2 from volume 1')

    # Cut with planes
    cubit.cmd('webcut volume all with plane xplane offset 0')
    cubit.cmd('webcut volume all with plane yplane offset 0')
    cubit.cmd('webcut volume all with plane zplane offset 0')

    # Mesh and move
    cubit.cmd('imprint volume all')
    cubit.cmd('merge volume all')
    cubit.cmd(f'volume all size {seed_size}')
    cubit.cmd('mesh volume all')
    cubit.cmd(f'move Volume all x {x0} y {y0} z {z0} include_merged')

    # Make a new block for all elements to export
    cubit.cmd('block 9 add hex all')
    cubit.cmd('block 9 name "all"')

    # Export
    cubit.cmd(f'export mesh "{output_file}.e" block 9 overwrite')

    return 0


def annulus_from_bounds(output_file, bounds_file, seed_size, annulus_ratio):
    '''Create a cylinder mesh from the bounds of a DNS file

    :param str output_file: The output filename
    :param str bounds_file: The file containing the bounds of the DNS
    :param float seed_size: The approximate mesh size
    :param float annulus_ratio: The fraction of the radius to keep in the final geometry

    Calls "mesh" function
    '''

    # Process the bounds data to calculate cylinder geometry info
    bounds_data = pandas.read_csv(bounds_file, sep=',')

    xmin = bounds_data['xmin'][0]
    xmax = bounds_data['xmax'][0]
    ymin = bounds_data['ymin'][0]
    ymax = bounds_data['ymax'][0]
    zmin = bounds_data['zmin'][0]
    zmax = bounds_data['zmax'][0]

    radx = (xmax - xmin) / 2
    rady = (ymax - ymin) / 2
    rad = numpy.mean([radx, rady])

    height = zmax - zmin

    x0 = xmin + rad
    y0 = ymin + rad
    z0 = zmin + (height / 2)

    # create mesh
    mesh(rad, height, annulus_ratio, x0, y0, z0, seed_size, output_file)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Create a cylinder mesh from the bounds of a DNS file."
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--output-file', type=str, required=True,
        help="The output filename")
    parser.add_argument('--bounds-file', type=str, required=True,
        help='The file containing the bounds of the DNS')
    parser.add_argument('--seed-size', type=float, required=True,
        help='The approximate mesh size')
    parser.add_argument('--annulus-ratio', type=float, required=True,
        help='The fraction of the radius to keep in the final geometry')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(annulus_from_bounds(output_file=args.output_file,
                                 bounds_file=args.bounds_file,
                                 seed_size=args.seed_size,
                                 annulus_ratio=args.annulus_ratio,
                                 ))

