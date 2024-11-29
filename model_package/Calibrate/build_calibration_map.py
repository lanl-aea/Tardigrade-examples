#!python
import argparse
import pathlib
import sys
import yaml

import numpy
import pandas

import calibration_tools


def ignored_elements_calibration_map(output_file, best_parameters_yml_no_BCs, boundary_csv, previous_calibration_map):
    '''Create a yaml file to map calibration results for interior and boundary elements

    :params str output_file: The name of the output yaml file
    :params str best_parameters_yml_no_BCs: A yaml file containing the 'best' calibration using the kernel density estimate for
                                            elements not on the boundary
    :params str boundary_csv: A csv file containing list of boundary elements
    :params str previous_calibration_map: A csv file containing a previous calibration map for all elements to be modified

    :returns: Write ``output_file``
    '''

    previous_parameter_df = pandas.read_csv(previous_calibration_map)
    boundary_elements_df = pandas.read_csv(boundary_csv)
    boundary_elements = list(boundary_elements_df['boundary_elements'].values)

    # if all elements are on boundary, overwrite name of best_parameters_yml_no_BCs to be best_parameters_yml
    if len(list(previous_parameter_df.index)) == len(boundary_elements):
        print('It looks like all elements are on the boundary! Forcing best_parameters_yml_no_BCs to be best_parameters_yml')
        best_parameters_yml_no_BCs = f'{best_parameters_yml_no_BCs.split("_no_BCs")[0]}.yml'

    # Unpack best_parameters_yml_no_BCs
    temp_dict = {}
    best_params, best_parameter_ordering = calibration_tools.parse_fparams_file(best_parameters_yml_no_BCs, material_type='plastic')
    for param, param_name in zip(best_params, best_parameter_ordering):
        temp_dict[param_name] = param
    temp_dict['obj_func_value'] = 0.0

    # Build new parameter dictionary
    headers = list(previous_parameter_df.columns)
    headers.remove('element')
    for index in previous_parameter_df.index:
        if int(previous_parameter_df.loc[index]['element'] in boundary_elements):
            element_id = int(previous_parameter_df.loc[index]['element'])
            for column in headers:
                #previous_parameter_df.loc[index][column] = temp_dict[column]
                previous_parameter_df.at[index, column] = temp_dict[column]

    # DataFrame and output
    previous_parameter_df['element'] = previous_parameter_df['element'].astype(int)
    previous_parameter_df.to_csv(output_file, header=True, sep=',', index=False)

    return 0


def full_csv_calibration_map(output_file, calibrated_elements, calibrated_files, material_type):
    '''Create a file mapping calibration results for each macroscale element

    :params str output_file: The name of the output file
    :params list calibrated_elements: A list of elements with associated calibration files
    :params list calibrated_files: A list of files containing calibration results
    :param str material_type: The material type: 'elastic', 'plastic', or 'full_plastic'

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
    df = pandas.DataFrame(out_params, columns=header)
    df['element'] = df['element'].astype(int)
    df.to_csv(output_file, header=True, sep=',', index=False)

    return 0


def trim_csv_file_for_tardigrade(output_file, input_file):
    '''Remove 'element' and 'obj_func_value' columns from a csv file containing calibration results

    :params str output_file: The name of the output csv file
    :params str input_csv: An input, previously generated csv file to be trimmed for Tardigrade

    :returns: Write ``output_file``
    '''

    df = pandas.read_csv(input_file)
    headers = list(df.columns)
    if 'element' in headers:
        df = df.drop('element', axis=1)
    if 'obj_func_value' in headers:
        df = df.drop('obj_func_value', axis=1)
    df.to_csv(output_file, header=False, sep=',', index=False)

    return 0


def build_calibration_map(output_file, calibrated_elements=None, calibrated_files=None, map_type='full_csv',
                          material_type=None, best_parameters_yml_no_BCs=None, boundary_csv=None, previous_calibration_map=None,
                          input_csv=None):
    '''Create a file mapping calibration results for each macroscale element

    :params str output_file: The name of the output csv file
    :params list calibrated_elements: A list of elements with associated calibration files
    :params list calibrated_files: A list of files containing calibration results
    :params str map_type: The type of calibration map to generate. 'full_csv' (default) to create a csv
                          file containing material parameters mapped for every element.
                          'ignore_boundary_yaml' to create a yaml file containing names of yaml files
                          containing material parameters for every element.
    :param str material_type: The material type: 'elastic', 'plastic', or 'full_plastic'
    :params str best_parameters_yml_no_BCs: A yaml file containing the 'best' calibration using the kernel density estimate for
                                            elements not on the boundary
    :params str boundary_csv: A csv file containing list of boundary elements
    :params str previous_calibration_map: A csv file containing a previous calibration map for all elements to be modified
    :params str input_csv: An input, previously generated csv file using 'map_type=full_csv' to be trimmed for Tardigrade

    :returns: Call full_csv_calibration_map function if map_type='full_csv', call ignored_elements_calibration_map function if map_type='ignore_boundary_yaml', or call trim_csv_file_for_tardigrade if map_type='trim_for_tardigrade'
    '''

    if map_type == 'full_csv':
        assert calibrated_elements != None, "List of calibrated elements must be provided!"
        assert calibrated_files != None, "List of calibration files must be provided!"
        assert len(calibrated_elements) == len(calibrated_files), "Length of calibrated_elements must equal calibrated_files!"
        assert material_type != None, "material_type must be specified!"
        full_csv_calibration_map(output_file, calibrated_elements, calibrated_files, material_type)
    elif map_type == 'ignore_boundary_yaml':
        assert best_parameters_yml_no_BCs != None, "best_parameters_yml_no_BCs must be provided!"
        assert boundary_csv != None, "boundary_csv must be provided!"
        assert previous_calibration_map != None, "previous_calibration_map must be provided!"
        ignored_elements_calibration_map(output_file, best_parameters_yml_no_BCs, boundary_csv, previous_calibration_map)
    elif map_type == 'trim_for_tardigrade':
        assert input_csv != None, "input_csv must be specified!"
        trim_csv_file_for_tardigrade(output_file, input_csv)
    else:
        raise NameError('Specify a valid map_type!')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Create a yaml file to map calibration results"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output csv file")
    parser.add_argument('--calibrated-elements', nargs="+", required=False,
        help="A list of elements with associated calibration files")
    parser.add_argument('--calibrated-files', nargs="+", required=False,
        help="A list of files containing calibration results")
    parser.add_argument('--map-type', type=str, required=False, default="full_csv",
        help="The type of calibration map to generate. 'full_csv' (default) to create a csv\
              file containing material parameters mapped for every element.\
              'ignore_boundary_yaml' to create a new calibration map with boundary element parameters\
              swapped with best_parameters_yml_no_BCs.\
              'trim_for_tardigrade' to modify a previously generated csv file for Tardigrade")
    parser.add_argument('--material-type', type=str, required=False, default=None,
        help="The material type: 'elastic', 'plastic', or 'full_plastic'")
    parser.add_argument('--best-parameters-yml-no-BCs', type=str, required=False, default=None,
        help="A yaml file containing the 'best' calibration using the kernel density estimate for elements not on the boundary")
    parser.add_argument('--boundary-csv', type=str, required=False, default=None,
        help="A csv file containing list of boundary elements")
    parser.add_argument('--previous-calibration-map', type=str, required=False, default=None,
        help="A csv file containing a previous calibration map for all elements to be modified")
    parser.add_argument('--input-csv', type=str, required=False, default=None,
        help="An input, previously generated csv file using 'map_type=full_csv' to be trimmed for Tardigrade")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(build_calibration_map(output_file=args.output_file,
                                   calibrated_elements=args.calibrated_elements,
                                   calibrated_files=args.calibrated_files,
                                   map_type=args.map_type,
                                   material_type=args.material_type,
                                   best_parameters_yml_no_BCs=args.best_parameters_yml_no_BCs,
                                   boundary_csv=args.boundary_csv,
                                   previous_calibration_map=args.previous_calibration_map,
                                   input_csv=args.input_csv,
                                   ))