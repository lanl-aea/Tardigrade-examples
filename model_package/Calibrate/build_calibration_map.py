import sys
import pathlib
import argparse

import yaml
import pandas
import numpy

import calibration_tools


def ignored_elements_calibration_map(output_file, calibrated_elements, calibrated_files, ignore_boundary_yml, ignore_boundary_summary_file)
    '''Create a yaml file to map calibration results for interior and boundary elements

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


def full_csv_calibration_map(output_file, calibrated_elements, calibrated_files, material_type):
    '''Create a file mapping calibration results for each macroscale element

    :params str output_file: The name of the output file
    :params list calibrated_elements: A list of elements with associated calibration files
    :params list calibrated_files: A list of files containing calibration results

    :returns: Write ``output_file``
    '''

    # set up header
    _, parameter_header = calibration_tools.parse_fparams_file(calibrated_files[0], material_type)
    header = ['element'] + parameter_header

    # unpack values
    out_params = []
    for element, file in zip(calibrated_elements, calibrated_files):
        parameters, _ = calibration_tools.parse_fparams_file(file, material_type)
        if len(out_params) == 0:
            out_params = numpy.hstack([float(element), parameters])
        else:
            out_params = numpy.vstack([out_params, numpy.hstack([float(element), parameters])])

    # DataFrame and output
    df = pandas.DataFrame(output_parameters, columns=header)
    df['element'] = df['element']astype(int)
    df.to_csv(output_file, header=True, sep=',', index=False)

    return 0


def build_calibration_map(output_file, calibrated_elements, calibrated_files, map_type='full_csv',
                          material_type=None,ignore_boundary_yml=None, ignore_boundary_summary_file=None, ):
    '''Create a file mapping calibration results for each macroscale element

    :params str output_file: The name of the output file
    :params list calibrated_elements: A list of elements with associated calibration files
    :params list calibrated_files: A list of files containing calibration results
    :params str map_type: The type of calibration map to generate. 'full_csv' (default) to create a csv
                          file containing material parameters mapped for every element.
                          'ignore_boundary_yaml' to create a yaml file containing names of yaml files
                          containing material parameters for every element.
    :param str material_type: The material type: 'elastic', 'plastic', or 'full_plastic'
    :params str ignore_boundary_yml: A yaml file containing the 'best' calibration using the kernel density estimate
    :params str ignore_boundary_summary_file: A csv file containing a summary of calibrated parameters for each element

    :returns: Call full_csv_calibration_map function if map_type='full_csv', or call ignored_elements_calibration_map function if map_type='ignore_boundary_yaml'
    '''

    if map_type == 'full_csv':
        assert len(calibrated_elements) = len(calibrated_files), "Length of calibrated_elements must equal calibrated_files!"
        assert material_type != None, "material_type must be specified!"
    elif map_type == 'ignore_boundary_yaml':
        assert ignore_boundary_yml != None, "ignore_boundary_yml must be provided!"
        assert ignore_boundary_summary_file != None, "ignore_boundary_yml must be provided!"
        assert len(calibrated_elements) = len(calibrated_files), "Length of calibrated_elements must equal calibrated_files!"
        ignored_elements_calibration_map(output_file, calibrated_elements, calibrated_files,
                                         ignore_boundary_yml, ignore_boundary_summary_file)
    else:
        raise NameError('Specify a valid map_type!')

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Create a yaml file to map calibration results"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output yaml file")
    parser.add_argument('--calibrated-elements', nargs="+", required=True,
        help="A list of elements with associated calibration files")
    parser.add_argument('--calibrated-files', nargs="+", required=True,
        help="A list of files containing calibration results")
    parser.add_argument('--map-type', type=str, required=False, default="full_csv",
        help="The type of calibration map to generate. 'full_csv' (default) to create a csv\
              file containing material parameters mapped for every element.\
              'ignore_boundary_yaml' to create a yaml file containing names of yaml files\
              containing material parameters for every element.")
    parser.add_argument('--material-type', type=str, required=False, default=None,
        help="The material type: 'elastic', 'plastic', or 'full_plastic'")
    parser.add_argument('--ignore-boundary-yml', type=str, required=False, default=None,
        help="A yaml file containing the 'best' calibration using the kernel density estimate")
    parser.add_argument('--ignore-boundary-summary-file', type=str, required=False, default=None,
        help="A csv file containing a summary of calibrated parameters for each element")


    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(build_calibration_map(output_file=args.output_file,
                                   calibrated_elements=args.calibrated_elements,
                                   calibrated_files=args.calibrated_files,
                                   map_type=args.map_type,
                                   material_type=args.material_type,
                                   ignore_boundary_yml=args.ignore_boundary_yml,
                                   ignore_boundary_summary_file=args.ignore_boundary_summary_file,
                                   ))