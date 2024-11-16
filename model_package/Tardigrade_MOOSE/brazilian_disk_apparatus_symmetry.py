#!python
import argparse
import os
import pathlib
import sys

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

    # Platen
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

    # Specimen
    cubit.cmd(f'create curve arc radius {spec_rad} center location {width/2} {bottom_spec + spec_rad + tol} \
                0 normal 0 0 1 start angle 0 stop angle 360')
    cubit.cmd('create surface curve 19')
    cubit.cmd(f'sweep surface 9  vector 0 0 1  distance {spec_dep}')


    # Cut for symmetry
    cubit.cmd('webcut volume all with general plane location curve 2 fraction 0.5 from vertex 1 direction on curve 2 location last')
    cubit.cmd('webcut volume all with general plane location curve 48 fraction 0.5 from vertex 35 direction on curve 48 location last')
    cubit.cmd('webcut volume all with general plane location curve 47 fraction 0.5 from vertex 34 direction on curve 47 location last')
    cubit.cmd('delete volume 1 2 4 5 6 7 8 9 11 12')

    # Blocks and sets
    cubit.cmd('block 1 add volume 3')
    cubit.cmd('block 1 name "bottom_platen"')
    cubit.cmd('block 2 add volume 10')
    cubit.cmd('block 2 name "specimen"')
    cubit.cmd('sideset 1 add surface 53')
    cubit.cmd('sideset 1 name "bottom"')
    cubit.cmd('sideset 2 add surface 82')
    cubit.cmd('sideset 2 name "top_sym"')
    cubit.cmd('sideset 3 add surface 51 84')
    cubit.cmd('sideset 3 name "back_sym"')
    cubit.cmd('sideset 4 add surface 54 85')
    cubit.cmd('sideset 4 name "side_sym"')
    cubit.cmd('sideset 5 add surface 56')
    cubit.cmd('sideset 5 name "platen_contact"')
    cubit.cmd('sideset 6 add surface 83')
    cubit.cmd('sideset 6 name "specimen_contact"')
    cubit.cmd('nodeset 1 add surface 53')
    cubit.cmd('nodeset 1 name "bottom"')
    cubit.cmd('nodeset 2 add surface 82')
    cubit.cmd('nodeset 2 name "top_sym"')
    cubit.cmd('nodeset 3 add surface 51 84')
    cubit.cmd('nodeset 3 name "back_sym"')
    cubit.cmd('nodeset 4 add surface 54 85')
    cubit.cmd('nodeset 4 name "side_sym"')
    cubit.cmd('nodeset 5 add surface 56')
    cubit.cmd('nodeset 5 name "platen_contact"')
    cubit.cmd('nodeset 6 add surface 83')
    cubit.cmd('nodeset 6 name "specimen_contact"')

    # Mesh
    cubit.cmd('imprint volume all')
    cubit.cmd(f'volume all size {seed_size}')
    cubit.cmd('mesh volume all')

    # Output
    cubit.cmd(f'save as "{output_file}.cub" overwrite')
    cubit.cmd(f'export mesh "{output_file}.e"  overwrite')
    cubit.cmd(f'export abaqus "{output_file}_bottom_platen.inp" block 1 source_csys 0 target_csys 0 partial dimension 3 overwrite')
    cubit.cmd(f'export abaqus "{output_file}_specimen.inp" block 2 source_csys 0 target_csys 0 partial dimension 3 overwrite')
    #cubit.cmd(f'export abaqus "{output_file}.inp" instance_per_block source_csys 0 target_csys 0 partial dimension 3 overwrite  everything')
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