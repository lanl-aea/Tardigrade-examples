import sys
import pathlib
import argparse

import yaml

def build_calibration_map(calibrated_elements, calibrated_files, ignore_boundary_yml, ignore_boundary_summary_file, output_file):
    '''Create a yaml file to map calibration results

    :params list calibrated_elements: A list of elements with associated calibration files
    :params list calibrated_files: A list of files containing calibration results
    :params str ignore_boundary_yml: A yaml file containing the 'best' calibration using the kernel density estimate
    :params str ignore_boundary_summary_file: A csv file containing a summary of calibrated parameters for each element
    :params str output_file: The name of the output yaml file

    :returns: Write ``output_file``
    '''

    calibration_map = {}
    for element, file in zip(calibrated_elements, calibrated_files):
        calibration_map[str(element)] = file
    calibration_map['ignore_boundary_yml'] = ignore_boundary_yml
    calibration_map['ignore_boundary_summary_file'] = ignore_boundary_summary_file
    with open(output_file, 'w') as f:
        yaml.dump(calibration_map, f)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Create a yaml file to map calibration results"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--calibrated-elements', nargs="+", required=True,
        help="A list of elements with associated calibration files")
    parser.add_argument('--calibrated-files', nargs="+", required=True,
        help="A list of files containing calibration results")
    parser.add_argument('--ignore-boundary-yml', type=str, required=True,
        help="A yaml file containing the 'best' calibration using the kernel density estimate")
    parser.add_argument('--ignore-boundary-summary-file', type=str, required=True,
        help="A csv file containing a summary of calibrated parameters for each element")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output yaml file")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(build_calibration_map(calibrated_elements=args.calibrated_elements,
                                   calibrated_files=args.calibrated_files,
                                   ignore_boundary_yml=args.ignore_boundary_yml,
                                   ignore_boundary_summary_file=args.ignore_boundary_summary_file,
                                   output_file=args.output_file,
                                   ))