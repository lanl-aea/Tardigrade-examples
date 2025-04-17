#!python
import argparse
import pathlib
import sys

import pandas

from model_package.Tardigrade_MOOSE import write_elastic_material_card


def elastic_map_to_material_card(map_file, element_number, output_file):
    '''Unpack a csv file of elastic parameters and call function to write elastic yaml file

    :params str map_file: CSV file containing previously calibrated elastic parameters
    :params int element_number: The id of the element to extract calibration data
    :params str output_file: The name of the yml material card to write
    '''

    df = pandas.read_csv(map_file)
    row = df[df['element'] == element_number]
    
    write_elastic_material_card.write_elastic_material_card(
        output_file,
        lamb=row['lambda'].values[0], mu=row['mu'].values[0],
        eta=row['eta'].values[0], tau=row['tau'].values[0], kappa=row['kappa'].values[0],
        nu=row['nu'].values[0], sigma=row['sigma'].values[0],
        tau1=row['tau1'].values[0], tau2=row['tau2'].values[0], tau3=row['tau3'].values[0], tau4=row['tau4'].values[0],
        tau5=row['tau5'].values[0], tau6=row['tau6'].values[0], tau7=row['tau7'].values[0], tau8=row['tau8'].values[0],
        tau9=row['tau9'].values[0], tau10=row['tau10'].values[0], tau11=row['tau11'].values[0])

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Unpack a csv file of elastic parameters and call function to write elastic yaml file"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('--map-file', type=str, required=True,
        help="CSV file containing previously calibrated elastic parameters")
    parser.add_argument('--element-number', type=int, required=True,
        help="The id of the element to extract calibration data")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the yml material card to write")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(elastic_map_to_material_card(map_file=args.map_file,
                                          element_number=args.element_number,
                                          output_file=args.output_file,
                                          ))