#!python
import argparse
import os
import pathlib
import sys

import cubit
import pandas
import numpy
import subprocess


def mesh(rad, height, x0, y0, z0, seed_size, output_file, platen_rad_factor, platen_depth, gap=None):
    ''' Mesh a cylinder using Cubit with platens

    :param float rad: Cylinder radius
    :param float height: Cylinder height
    :param float x0: The x-distance to translate the cylinder
    :param float y0: The y-distance to translate the cylinder
    :param float z0: The z-distance to translate the cylinder
    :param float seed-size: The approximate mesh size
    :param str output_file: The output filename
    :param float platen_rad_factor: The factor to multiply the specimen radius by to determine the specimen radius
    :param float platen_depth: The thickness of the platens

    :returns: ``{output_file}.e``
    '''

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')
    cubit.cmd(f'create Cylinder height {height} radius {rad}')

    # Cut with planes
    cubit.cmd('webcut volume all with plane xplane offset 0')
    cubit.cmd('webcut volume all with plane yplane offset 0')
    cubit.cmd('webcut volume all with plane zplane offset 0')
    # Side sets
    cubit.cmd('sideset 1 add surface 28 20 26 14')
    cubit.cmd('sideset 1 name "top"')
    cubit.cmd('sideset 2 add surface 24 16 30 18')
    cubit.cmd('sideset 2 name "bottom"')
    cubit.cmd('sideset 3 add surface 41 59 47 61')
    cubit.cmd('sideset 3 name "x_plane"')
    cubit.cmd('sideset 4 add surface 58 50 62 54')
    cubit.cmd('sideset 4 name "y_plane"')

    # Top Platen
    cubit.cmd(f'brick x {2*rad*platen_rad_factor} y {2*rad*platen_rad_factor} z {platen_depth}')
    cubit.cmd('webcut volume 9 with plane xplane offset 0')
    cubit.cmd('webcut volume 9 10 with plane yplane offset 0')
    cubit.cmd(f'move Volume 9 10 11 12 x 0 y 0 z {0.5*(height + platen_depth)} include_merged')

    # Bottom Platen
    cubit.cmd(f'brick x {2*rad*platen_rad_factor} y {2*rad*platen_rad_factor} z {platen_depth}')
    cubit.cmd('webcut volume 13 with plane xplane offset 0')
    cubit.cmd('webcut volume 13 14 with plane yplane offset 0')
    cubit.cmd(f'move Volume 13 14 15 16 x 0 y 0 z {-0.5*(height + platen_depth)} include_merged')

    # Platen Side sets
    cubit.cmd('sideset 5 add surface 84 86 92 98')
    cubit.cmd('sideset 5 name "top_platen_bottom"')
    cubit.cmd('sideset 6 add surface 82 88 94 96')
    cubit.cmd('sideset 6 name "top_platen_top"')
    cubit.cmd('sideset 7 add surface 120 122 128 134')
    cubit.cmd('sideset 7 name "bottom_platen_bottom"')
    cubit.cmd('sideset 8 add surface 118 124 130 132')
    cubit.cmd('sideset 8 name "bottom_platen_top"')

    # Mesh and Move
    ## Specimen
    cubit.cmd('imprint volume 1 2 3 4 5 6 7 8')
    cubit.cmd('merge volume 1 2 3 4 5 6 7 8')
    cubit.cmd(f'volume 1 2 3 4 5 6 7 8 size {seed_size}')
    cubit.cmd('mesh volume 1 2 3 4 5 6 7 8')
    ## Top platen
    cubit.cmd('imprint volume 9 10 11 12')
    cubit.cmd('merge volume 9 10 11 12')
    cubit.cmd(f'volume 9 10 11 12 size {seed_size}')
    cubit.cmd('mesh volume 9 10 11 12')
    ## Bottom platen
    cubit.cmd('imprint volume 13 14 15 16')
    cubit.cmd('merge volume 13 14 15 16')
    cubit.cmd(f'volume 13 14 15 16 size {seed_size}')
    cubit.cmd('mesh volume 13 14 15 16')
    ## Move
    cubit.cmd(f'move Volume all x {x0} y {y0} z {z0} include_merged')
    if gap:
        cubit.cmd(f'move Volume 9 10 11 12 17 18 19 20 x 0 y 0 z {gap} include_merged')
    # Blocks
    cubit.cmd('block 1 add volume 1 2 3 4 5 6 7 8')
    cubit.cmd('block 1 name "specimen"')
    cubit.cmd('block 2 add volume 9 10 11 12')
    cubit.cmd('block 2 name "top_platen"')
    cubit.cmd('block 3 add volume 13 14 15 16')
    cubit.cmd('block 3 name "bottom_platen"')

    # Export
    #cubit.cmd(f'
    cubit.cmd(f'export mesh "{output_file}.e" block 1 2 3 overwrite')


    return 0


def cylinder_from_bounds_with_platens(output_file, bounds_file, seed_size, platen_rad_factor=1.5, platen_depth=1.0, xdmf=True, ascii=False, gap=None):
    '''Create a cylinder mesh from the bounds of a DNS file with platens

    :param str output_file: The output filename
    :param str bounds_file: The file containing the bounds of the DNS
    :param float seed_size: The approximate mesh size
    :param float platen_rad_factor: The factor to multiply the specimen radius by to determine the specimen radius
    :param float platen_depth: The thickness of the platens
    :param bool xdmf: The option to convert default exodus mesh to XDMF (binary)
    :param bool ascii: The option to convert binary XDMF mesh to ascii

    Calls "mesh" function and converts ``{output_file}.e`` to ``{output_file}.xdmf``
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
    mesh(rad, height, x0, y0, z0, seed_size, output_file, platen_rad_factor, platen_depth, gap)

    # convert to XDMF with subprocess
    if xdmf:
        subprocess.run([f'meshio convert {output_file}.e {output_file}.xdmf'], shell=True)
        if ascii:
            subprocess.run([f'meshio ascii {output_file}.xdmf'], shell=True)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Create a cylinder mesh from the bounds of a DNS file with platens."
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--output-file', type=str, required=True,
        help="The output filename")
    parser.add_argument('--bounds-file', type=str, required=True,
        help='The file containing the bounds of the DNS')
    parser.add_argument('--seed-size', type=float, required=True,
        help='The approximate mesh size')
    parser.add_argument('--platen-rad-factor', type=float, required=False, default=1.5,
        help='The factor to multiply the specimen radius by to determine the specimen radius')
    parser.add_argument('--platen-depth', type=float, required=False, default=1.0,
        help='The thickness of the platens')
    parser.add_argument('--xdmf', type=str, required=False,
        help='The option to convert default exodus mesh to XDMF (binary)')
    parser.add_argument('--ascii', type=str, required=False, default=False,
        help='The option to convert binary XDMF mesh to ascii')
    parser.add_argument('--gap', type=str, required=False, default=None,
        help='An initial gap to place between specimen and platens')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(cylinder_from_bounds_with_platens(output_file=args.output_file,
                                               bounds_file=args.bounds_file,
                                               seed_size=args.seed_size,
                                               platen_rad_factor=args.platen_rad_factor,
                                               platen_depth=args.platen_depth,
                                               xdmf=bool(args.xdmf),
                                               ascii=bool(args.ascii),
                                               gap=args.gap,
                                               ))

