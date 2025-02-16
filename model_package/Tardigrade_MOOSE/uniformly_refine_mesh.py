#!python
import argparse
import pathlib
import sys

import cubit
import pandas


def refine_mesh(input_mesh, output_mesh, refinement_level, calibration_map_in, calibration_map_out=None):
    '''Uniformly refine an exodus mesh and update a calibration map with new element IDs

    :params str input_mesh: The input exodus mesh file to refine
    :params str output_mesh: The output exodus mesh file
    :params int refinement_level: The uniform refinement level, 1: refine by 1 level (each hex element is split into 8), 2: refinement by 2 levels (each hex element is split into 24)
    :params str calibration_map_in: The original calibration map file
    :params str calibration_map_out: The output calibration map file with updated element ids

    :returns: ``{output_mesh}`` and ``{calibration_map_out}``
    '''

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')
    cubit.cmd(f'import mesh geometry "{input_mesh}"')

    # refine the mesh
    cubit.cmd(f'refine volume all numsplit {refinement_level}')

    # Open the calibration map 
    df = pandas.read_csv(calibration_map_in)

    # Collect calibration values for new element IDs
    dict_out = {}
    block_ids = cubit.get_block_id_list()
    for block_id in block_ids:
        block_elements = cubit.get_volume_hexes(block_id)
        for element in block_elements:
            dict_out[element] = df.loc[block_id-1]
    df_out = pandas.DataFrame(dict_out)
    print(df_out)
    df_out = df_out.T
    print(df_out)

    # Write out new calibration map
    if calibration_map_out:
        df_out.to_csv(calibration_map_out, sep=',', index=False)

    # Save mesh
    cubit.cmd(f'export mesh "{output_mesh}"  overwrite')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Uniformly refine an exodus mesh and update a calibration map with new element IDs"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--input-mesh', type=str, required=True,
        help="The input exodus mesh file to refine")
    parser.add_argument('--output-mesh', type=str, required=True,
        help="The output exodus mesh file")
    parser.add_argument('--refinement-level', type=int, required=True,
        help="The uniform refinement level, 1: refine by 1 level (each hex element is split into 8),\
              2: refinement by 2 levels (each hex element is split into 24)")
    parser.add_argument('--calibration-map-in', type=str, required=True,
        help="The original calibration map file")
    parser.add_argument('--calibration-map-out', type=str, required=False, default=None,
        help="The output calibration map file with updated element ids")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(refine_mesh(input_mesh=args.input_mesh,
                         output_mesh=args.output_mesh,
                         refinement_level=args.refinement_level,
                         calibration_map_in=args.calibration_map_in,
                         calibration_map_out=args.calibration_map_out,
                         ))

