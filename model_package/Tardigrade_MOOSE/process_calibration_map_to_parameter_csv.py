#!python
import argparse
import os
import pathlib
import sys

import pandas


def process_calibration_map_to_parameter_csv(output_file, calibration_map):
    '''Process a calibration map file to a parameter csv for Tardigrade-MOOSE

    :param str output_file: Specify the name of the output csv to write
    :param str calibration_map: CSV file containing calibration data

    :returns: ``output_file``
    '''

    # load calibration map
    parameter_df = pandas.read_csv(calibration_map, index_col='element')
    #parameter_df = pandas.read_csv(calibration_map)

    # drop objective function evaluation
    if 'obj_func_value' in parameter_df.columns:
        parameter_df = parameter_df.drop(columns=['obj_func_value'])

    # output
    parameter_df.to_csv(output_file, header=False, sep=',', index=False)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Process a calibration map file to a parameter csv for Tardigrade-MOOSE"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of the output csv to write")
    parser.add_argument('--calibration-map', type=str, required=True,
        help='CSV file containing calibration data')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(process_calibration_map_to_parameter_csv(output_file=args.output_file,
                                                      calibration_map=args.calibration_map,
                                                      ))