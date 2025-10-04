import sys
import os
import argparse
import pathlib

import cubit
import pandas
import numpy
import subprocess

from misc_utilities import str2bool


def brazilian_disk_apparatus(output_file, specimen_seed_size, platen_seed_size,
                             height, width, chord, app_rad, app_dep, spec_rad, spec_dep, tol,
                             x0=0., y0=0., z0=0., export_platens=True):
    '''Create a Brazilian Disk specimen and loading apparatus

    :param str output_file: The output filename
    :param float specimen_seed_size: The approximate mesh size for the specimen
    :param float platen_seed_size: The approximate mesh size for the platen
    :param float height: The height of a single Brazilian disk compression platen
    :param float width: The base width of a Brazilian disk compression platen
    :param float chord: The chord distance of the Brazilian disk compression platen
    :param float app_rad: The radius of curvature of the Brazilian disk compression platen
    :param float depth: The extrusion depth of the Brazilian disk compression platen
    :param float spec_rad: The radius of the Brazilian disk compression specimen
    :param float spec_dep: The extrusion depth of the Brazilian disk compression specimen
    :param float tol: A tolerance / gap distance to insert between Brazilian disk compression specimen and platens
    :param float x0: The x-location to move geometry for the center of the Brazil Disk
    :param float y0: The y-location to move geometry for the center of the Brazil Disk
    :param float z0: The z-location to move geometry for the center of the Brazil Dis
    :param bool export_platens: Flag to export platen meshes of the brazilian disk apparatus

    :returns: Write ``{output_file}.cub``, ``{output_file}_specimen.inp``, and optionally ``{output_file}_bottom_platen.inp`` and ``{output_file}_top_platen.inp``
    '''

    perp_dist = app_rad - numpy.sqrt((app_rad**2 - (0.5*chord)**2))
    print(perp_dist)
    bottom_spec = height - perp_dist
    gap = 2*(spec_rad - perp_dist)
    offset = (app_dep - spec_dep)/2

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
                {offset} normal 0 0 1 start angle 0 stop angle 360')
    cubit.cmd('create surface curve 37')
    cubit.cmd(f'sweep surface 17  vector 0 0 1  distance {spec_dep}')

    # Partition
    cubit.cmd('webcut volume all with general plane curve 18')
    cubit.cmd('webcut volume 3 6 with general plane curve 75')

    # Move
    move_x = width/2
    move_y = bottom_spec+spec_rad+tol
    move_z = offset + spec_dep/2
    cubit.cmd(f'move Volume all x {x0 - move_x} y {y0 - move_y} z {z0 - move_z} include_merged')

    # Blocks and sets
    cubit.cmd('block 1 add volume 1 4')
    cubit.cmd('block 1 name "bottom_platen"')
    cubit.cmd('block 2 add volume 2 5')
    cubit.cmd('block 2 name "top_platen"')
    cubit.cmd('block 3 add volume 3 6 7 8')
    cubit.cmd('block 3 name "specimen"')
    ## sidesets
    cubit.cmd('sideset 1 add surface 23 27')
    cubit.cmd('sideset 1 name "bottom_platen_bottom"')
    cubit.cmd('sideset 2 add surface 33 37')
    cubit.cmd('sideset 2 name "top_platen_top"')
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
    cubit.cmd('sideset 8 add surface 2')
    cubit.cmd('sideset 8 name "bottom_platen_side"')
    cubit.cmd('sideset 9 add surface 9')
    cubit.cmd('sideset 9 name "top_platen_side"')
    cubit.cmd('sideset 10 add surface 22 29')
    cubit.cmd('sideset 10 name "bottom_platen_back"')
    cubit.cmd('sideset 11 add surface 34 36')
    cubit.cmd('sideset 11 name "top_platen_back"')
    cubit.cmd('sideset 12 add surface 51 58 62 64')
    cubit.cmd('sideset 12 name "specimen_front"')
    cubit.cmd('sideset 13 add surface 53 55 60 66')
    cubit.cmd('sideset 13 name "specimen_back"')
    ## nodesets
    cubit.cmd('nodeset 1 add surface 23 27')
    cubit.cmd('nodeset 1 name "bottom_platen_bottom"')
    cubit.cmd('nodeset 2 add surface 33 37')
    cubit.cmd('nodeset 2 name "top_platen_top"')
    cubit.cmd('nodeset 3 add surface 53 55 60 66')
    cubit.cmd('nodeset 3 name "back"')
    cubit.cmd('nodeset 4 add surface 31 38')
    cubit.cmd('nodeset 4 name "top_platen_contact"')
    cubit.cmd('nodeset 5 add surface 21 28')
    cubit.cmd('nodeset 5 name "bottom_platen_contact"')
    cubit.cmd('nodeset 6 add surface 47 50')
    cubit.cmd('nodeset 6 name "specimen_top"')
    cubit.cmd('nodeset 7 add surface 46 57')
    cubit.cmd('nodeset 7 name "specimen_bottom"')
    cubit.cmd('nodeset 8 add surface 2')
    cubit.cmd('nodeset 8 name "bottom_platen_side"')
    cubit.cmd('nodeset 9 add surface 9')
    cubit.cmd('nodeset 9 name "top_platen_side"')
    cubit.cmd('nodeset 10 add surface 22 29')
    cubit.cmd('nodeset 10 name "bottom_platen_back"')
    cubit.cmd('nodeset 11 add surface 34 36')
    cubit.cmd('nodeset 11 name "top_platen_back"')
    cubit.cmd('nodeset 12 add surface 51 58 62 64')
    cubit.cmd('nodeset 12 name "specimen_front"')
    cubit.cmd('nodeset 13 add surface 53 55 60 66')
    cubit.cmd('nodeset 13 name "specimen_back"')
    cubit.cmd('nodeset 15 add curve 74')
    cubit.cmd('nodeset 15 name "top_line_load"')
    cubit.cmd('nodeset 16 add curve 72')
    cubit.cmd('nodeset 16 name "bottom_line_load"')

    # Mesh
    cubit.cmd('merge volume 3 6 7 8')
    cubit.cmd('merge volume 1 4')
    cubit.cmd('merge volume 2 5')
    cubit.cmd(f'volume 3 6 7 8 size {specimen_seed_size}')
    cubit.cmd('mesh volume 3 6 7 8')
    cubit.cmd(f'volume 1 4 size {platen_seed_size}')
    cubit.cmd('mesh volume 1 4')
    cubit.cmd(f'volume 2 5 size {platen_seed_size}')
    cubit.cmd('mesh volume 2 5')

    # All nodes
    cubit.cmd('nodeset 14 add node all')
    cubit.cmd('nodeset 14 name "all"')

    # Export
    if export_platens == True:
        cubit.cmd(f'export abaqus "{output_file}_bottom_platen.inp" block 1 source_csys 0 target_csys 0 partial dimension 3 overwrite')
        cubit.cmd(f'export abaqus "{output_file}_top_platen.inp" block 2 source_csys 0 target_csys 0 partial dimension 3 overwrite')
    else:
        cubit.cmd('delete block 1 2')
        cubit.cmd('delete volume 1 2 4 5')
        cubit.cmd('delete nodeset 1 2 4 5 8 9 10 11')
        cubit.cmd('delete sideset 1 2 4 5 8 9 10 11')
    cubit.cmd(f'save as "{output_file}.cub" overwrite')
    cubit.cmd(f'export mesh "{output_file}.e"  overwrite')
    cubit.cmd(f'export abaqus "{output_file}_specimen.inp" block 3 source_csys 0 target_csys 0 partial dimension 3 overwrite')

    return


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Create a Brazilian Disk specimen and loading apparatus"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--output-file', type=str, required=True,
        help="The output filename")
    parser.add_argument('--specimen-seed-size', type=float, required=True,
        help='The approximate mesh size for the specimen')
    parser.add_argument('--platen-seed-size', type=float, required=True,
        help='The approximate mesh size for the platens')
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
    parser.add_argument('--x0', type=float, required=False, default=0.,
        help='The x-location to move geometry for the center of the Brazil Disk')
    parser.add_argument('--y0', type=float, required=False, default=0.,
        help='The y-location to move geometry for the center of the Brazil Disk')
    parser.add_argument('--z0', type=float, required=False, default=0.,
        help='The z-location to move geometry for the center of the Brazil Disk')
    parser.add_argument('--export-platens', type=str, required=False, default='True',
        help='Flag to export platen meshes of the brazilian disk apparatus')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(brazilian_disk_apparatus(output_file=args.output_file,
                                      specimen_seed_size=args.specimen_seed_size,
                                      platen_seed_size=args.platen_seed_size,
                                      height=args.height,
                                      width=args.width,
                                      chord=args.chord,
                                      app_rad=args.app_rad,
                                      app_dep=args.app_dep,
                                      spec_rad=args.spec_rad,
                                      spec_dep=args.spec_dep,
                                      tol=args.tol,
                                      x0=args.x0,
                                      y0=args.y0,
                                      z0=args.z0,
                                      export_platens=str2bool(args.export_platens),
                                      ))