#!python
import sys
import os
import argparse
import pathlib

import cubit
import pandas
import numpy
import subprocess


def mesh(rad, height, x0, y0, z0, seed_size, output_file, cut=False):

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')
    cubit.cmd(f'create Cylinder height {height} radius {rad}')
    if cut:
        cubit.cmd('webcut volume all with plane xplane offset 0')
        cubit.cmd('webcut volume all with plane yplane offset 0')
        cubit.cmd('webcut volume all with plane zplane offset 0')
        cubit.cmd('sideset 1 add surface 28 20 26 14')
        cubit.cmd('sideset 1 name "top"')
        cubit.cmd('sideset 2 add surface 24 16 30 18')
        cubit.cmd('sideset 2 name "bottom"')
        cubit.cmd('sideset 3 add surface 41 59 47 61')
        cubit.cmd('sideset 3 name "x_plane"')
        cubit.cmd('sideset 4 add surface 58 50 62 54')
        cubit.cmd('sideset 4 name "y_plane"')
    else:
        cubit.cmd('sideset 1 add surface 3')
        cubit.cmd('sideset 1 name "top"')
        cubit.cmd('sideset 2 add surface 2')
        cubit.cmd('sideset 2 name "bottom"')
    cubit.cmd(f'volume all size {seed_size}')
    cubit.cmd('mesh volume all')
    cubit.cmd(f'move Volume all x {x0} y {y0} z {z0} include_merged')

    cubit.cmd(f'export mesh "{output_file}.e"  overwrite')
    #cubit.cmd(f'export abaqus "{output_file}inp"  overwrite  everything')
    #cubit.cmd(f'export fluent "{output_file}.msh"  overwrite  everything')
    #cubit.cmd(f'export stl ascii "{output_file}.stl" mesh  overwrite')

    return 0


def main(output_file, bounds_file, seed_size, cut=False):

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
    new_rad = rad*1.02

    height = zmax - zmin

    x0 = xmin + rad
    y0 = ymin + rad
    z0 = zmin + (height / 2)

    # create mesh
    mesh(new_rad, height, x0, y0, z0, seed_size, output_file, cut)

    # convert to XDMF with subprocess
    subprocess.run([f'meshio convert {output_file}.e {output_file}.xdmf'], shell=True)
    subprocess.run([f'meshio ascii {output_file}.xdmf'], shell=True)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Create a cylinder mesh from the bounds of a DNS file."
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--output-file', type=str, required=True,
        help="The output file for the Cubit model without extension. Will be appended with the " \
             "required extension, e.g. ``output_file``.cub")
    parser.add_argument('--bounds-file', type=str, required=True,
        help='The file containing the bounds of the DNS')
    parser.add_argument('--seed-size', type=float, required=True,
        help='The approximate mesh size')
    parser.add_argument('--cut', type=str, required=False,
        help='The option to cut geometry into octants, pass string "True" if desired')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(main(output_file=args.output_file,
                  bounds_file=args.bounds_file,
                  seed_size=args.seed_size,
                  cut=args.cut,
                  ))

