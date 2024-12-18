import sys
import os
import argparse
import pathlib

import cubit
import pandas
import numpy
import subprocess


def brazilian_disk_apparatus(output_file, seed_size, height, width, chord, app_rad, app_dep, spec_rad, spec_dep, tol):
    '''Create a Brazilian Disk specimen and loading apparatus

    :param str output_file: The output filename
    :param float seed_size: The approximate mesh size
    :param float height: The height of a single Brazilian disk compression platen
    :param float width: The base width of a Brazilian disk compression platen
    :param float chord: The chord distance of the Brazilian disk compression platen
    :param float app_rad: The radius of curvature of the Brazilian disk compression platen
    :param float depth: The extrusion depth of the Brazilian disk compression platen
    :param float spec_rad: The radius of the Brazilian disk compression specimen
    :param float spec_dep: The extrusion depth of the Brazilian disk compression specimen
    :param float tol: A tolerance / gap distance to insert between Brazilian disk \
                      compression specimen and platens

    :returns: ``output_file``
    '''

    perp_dist = app_rad - numpy.sqrt((app_rad**2 - (0.5*chord)**2))
    print(perp_dist)
    bottom_spec = height - perp_dist
    gap = 2*(spec_rad - perp_dist)

    width_2 = (width - chord)/2.

    # Platen 1
    cubit.cmd('create vertex 0 0 0')
    cubit.cmd(f'create vertex 0 {height} 0')
    cubit.cmd(f'create vertex {width} 0 0') 
    cubit.cmd(f'create vertex {width} {height} 0')
    cubit.cmd(f'create vertex {width_2} {height} 0')
    cubit.cmd(f'create vertex {width - width_2} {height} 0')
    cubit.cmd('create curve vertex 1 2')  
    cubit.cmd('create curve vertex 1 3')
    cubit.cmd('create curve vertex 3 4')
    cubit.cmd('create curve vertex 2 5')
    cubit.cmd('create curve vertex 6 4')
    cubit.cmd(f'create curve arc vertex 5 6 radius {app_rad} normal 0 0 1')
    cubit.cmd('create surface curve all')
    cubit.cmd(f'sweep surface 1  vector 0 0 1  distance {app_dep}')

    # Platen 2
    cubit.cmd('Volume 1 copy rotate 180 about x')
    cubit.cmd(f'move Volume 2 x 0 y {2*height + gap + 2*tol} z {app_dep} include_merged')

    # Specimen
    cubit.cmd(f'create curve arc radius {spec_rad} center location {width/2} {bottom_spec+spec_rad+tol} \
                0 normal 0 0 1 start angle 0 stop angle 360')
    cubit.cmd('create surface curve 37')
    cubit.cmd(f'sweep surface 17  vector 0 0 1  distance {spec_dep}')

    # Partition
    cubit.cmd('webcut volume all with general plane curve 18')
    cubit.cmd('webcut volume 3 6 with general plane curve 75')

    # Blocks and sets
    cubit.cmd('block 1 add volume 1 4')
    cubit.cmd('block 1 name "bottom_platen"')
    cubit.cmd('block 2 add volume 2 5')
    cubit.cmd('block 2 name "top_platen"')
    cubit.cmd('block 3 add volume 3 6 7 8')
    cubit.cmd('block 3 name "specimen"')
    cubit.cmd('sideset 1 add surface 23 27')
    cubit.cmd('sideset 1 name "bottom"')
    cubit.cmd('sideset 2 add surface 33 37')
    cubit.cmd('sideset 2 name "top"')
    cubit.cmd('sideset 3 add surface 53 55 60 66')
    cubit.cmd('sideset 3 name "back"')
    cubit.cmd('sideset 4 add surface 31 38')
    cubit.cmd('sideset 4 name "top_platen_contact"')
    cubit.cmd('sideset 5 add surface 21 28')
    cubit.cmd('sideset 5 name "bottom_platen_contact"')
    cubit.cmd('sideset 6 add surface 47 50')
    cubit.cmd('sideset 6 name "specimen_top"')
    cubit.cmd('sideset 7 add surface 46 57')
    cubit.cmd('sideset 7 name "specimen_bottom"')

    # Mesh
    cubit.cmd('imprint volume all')
    cubit.cmd('merge volume 3 6 7 8')
    cubit.cmd('merge volume 1 4')
    cubit.cmd('merge volume 2 5')
    cubit.cmd(f'volume all size {seed_size}')
    cubit.cmd('mesh volume all')

    # Export
    cubit.cmd(f'save as "{output_file}.cub" overwrite')
    cubit.cmd(f'export mesh "{output_file}.e"  overwrite')

    return


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Create a Brazilian Disk specimen and loading apparatus"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--output-file', type=str, required=True,
        help="The output filename")
    parser.add_argument('--seed-size', type=float, required=True,
        help='The approximate mesh size')
    parser.add_argument('--height', type=float, required=True,
        help='The height of a single Brazilian disk compression platen')
    parser.add_argument('--width', type=float, required=True,
        help='The base width of a Brazilian disk compression platen')
    parser.add_argument('--chord', type=float, required=True,
        help='The chord distance of the Brazilian disk compression platen')
    parser.add_argument('--app-rad', type=float, required=True,
        help='The radius of curvature of the Brazilian disk compression platen')
    parser.add_argument('--app-dep', type=float, required=True,
        help='The extrusion depth of the Brazilian disk compression platen')
    parser.add_argument('--spec-rad', type=float, required=True,
        help='The radius of the Brazilian disk compression specimen')
    parser.add_argument('--spec-dep', type=float, required=True,
        help='The extrusion depth of the Brazilian disk compression specimen')
    parser.add_argument('--tol', type=float, required=False, default=0.001,
        help='A tolerance / gap distance to insert between Brazilian disk \
              compression specimen and platens')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(brazilian_disk_apparatus(output_file=args.output_file,
                                      seed_size=args.seed_size,
                                      height=args.height,
                                      width=args.width,
                                      chord=args.chord,
                                      app_rad=args.app_rad,
                                      app_dep=args.app_dep,
                                      spec_rad=args.spec_rad,
                                      spec_dep=args.spec_dep,
                                      tol=args.tol,
                                      ))