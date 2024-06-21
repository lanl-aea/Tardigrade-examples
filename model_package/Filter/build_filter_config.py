#!python
import os
import inspect
import sys
import yaml
import argparse


def write_filter_config(output_file, job_name, dns_file, macro_file, volume, density,  displacement, cauchy_stress, velocity=None, acceleration=None, max_parallel=None):
    '''Write the configuration file for the Micromorphic Filter

    :param str output_file: The output filename for filter configuration
    :param str job_name: The name of the job for the Micromorphic Filter
    :param str dns_file: The name of the XDMF file containing DNS data
    :param str macro_file: The name of the macroscale filter domain file
    :param str volume: The string identifying volume quantities located in "dns_file"
    :param str density: The string identifying density quantities located in "dns-file"
    :param str cauchy_stress: The string identifying stress quantities located in "dns-file"
    :param str displacement: The string identifying displacement quantities located in "dns-file"
    :param str velocity: Optional string identifying velocity quantities located in "dns-file"
    :param str acceleration:  Optional string identifying acceleration quantities located in "dns-file"
    :param int max_parallel: Optional parameter defining the number of parallel processes for the Micromorphic Filter

    returns ``output_file``
    '''

    quantity_dict = {}
    # required quantities
    quantity_dict["volume"] = volume
    quantity_dict["density"] = density
    quantity_dict["displacement"] = displacement
    quantity_dict["cauchy_stress"] = cauchy_stress
    # optional parameters
    if velocity:
        quantity_dict["velocity"] = velocity
    if acceleration:
        quantity_dict["acceleration"] = acceleration

    data = {
        "files":{"output": job_name,\
                 "data": dns_file,\
                 "filter": macro_file},\
        "quantity_names":quantity_dict,\
        "filter":{},}
    # max parallel
    if max_parallel:
        data['filter'].update({'max_parallel':max_parallel})
    # velocity gradient terms for acceleration
    if acceleration:
        data['filter'].update({'add_velocity_gradient_terms':True})

    # dump file
    with open(output_file, 'w') as file:
        yaml.dump(data, file)
    print('configuration file written!')

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Write the configuration file for the Micromorphic Filter"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='Specify the output filename for filter configuration')
    parser.add_argument('--job-name', type=str, required=True,
        help='Specify the name of the job for the Micromorphic Filter')
    parser.add_argument('--dns-file', type=str, required=True,
        help='Specify the name of the XDMF file containing DNS data')
    parser.add_argument('--macro-file', type=str, required=True,
        help='Specify the name of the macroscale filter domain file')
    parser.add_argument('--volume', type=str, required=True,
        help='Specify the string identifying volume quantities located in "dns-file"')
    parser.add_argument('--density', type=str, required=True,
        help='Specify the string identifying density quantities located in "dns-file"')
    parser.add_argument('--cauchy-stress', type=str, required=True,
        help='Specify the string identifying stress quantities located in "dns-file"')
    parser.add_argument('--displacement', type=str, required=True,
        help='Specify the string identifying displacement quantities located in "dns-file"')
    parser.add_argument('--velocity', type=str, required=False, default=None,
        help='Optional string identifying velocity quantities located in "dns-file"')
    parser.add_argument('--acceleration', type=str, required=False, default=None,
        help='Optional string identifying acceleration quantities located in "dns-file"')
    parser.add_argument('--max-parallel', type=int, required=False, default=None,
        help='Optional parameter defining the number of parallel processes for the\
              Micromorphic Filter')
    # TODO: add non-required arguments for optional quantities
    return parser


if __name__ == '__main__':
    parser = get_parser()

    args = parser.parse_args()
    sys.exit(write_filter_config(output_file=args.output_file, 
                                 job_name=args.job_name, 
                                 dns_file=args.dns_file, 
                                 macro_file=args.macro_file, 
                                 volume=args.volume, 
                                 density=args.density, 
                                 displacement=args.displacement, 
                                 cauchy_stress=args.cauchy_stress, 
                                 velocity=args.velocity, 
                                 acceleration=args.acceleration, 
                                 max_parallel=args.max_parallel))