#!python
import argparse
import pathlib
import sys
import yaml


def write_filter_config(output_file, job_name, dns_file, macro_file,
                        volume, density, displacement, cauchy_stress,
                        velocity=None, acceleration=None, damage=None,
                        max_parallel=None, sets_file=None, spectral=None,
                        update_filter_domains=False):
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
    :param str damage: Optional string identifying damage quantities located in "dns-file"
    :param int max_parallel: Optional parameter defining the number of parallel processes for the Micromorphic Filter
    :param str sets_file: Optional yaml file containing prescribed micro-averaging domains
    :param bool update_file_domains: Option ot update filter and microaveraging domains for each time step

    returns ``output_file``
    '''

    quantity_dict = {}
    filter_dict = {}

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
    if damage:
        quantity_dict["damage"] = damage

    # Prescribed or spectral domains
    if (sets_file is not None) and (spectral is None):
        filter_dict["micro_averaging_domains"] = 'prescribed'
        stream = open(sets_file, 'r')
        sets = yaml.load(stream, Loader=yaml.FullLoader)
        stream.close()
        quantity_dict["prescribed_micro_averaging_domains"] = list(sets.keys())
    elif (sets_file is None) and (spectral is not None):
        filter_dict["micro_averaging_domains"] = 'spectral'
    elif (sets_file is None) and (spectral is None):
        print('Using the default "nearest qpt" micro averaging domain method')
    else:
        print('Invalid configuration for "sets_file" and "spectral"')

    # Option ot update filter and microaveraging domains for each timestep
    if bool(update_filter_domains) == True:
        filter_dict["update_filter_domains"] = True

    # assemble main dictionary
    data = {
        "files":{"output": job_name,\
                 "data": dns_file,\
                 "filter": macro_file},\
        "quantity_names":quantity_dict,\
        "filter":filter_dict,}

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

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Write the configuration file for the Micromorphic Filter"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
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
    parser.add_argument('--damage', type=str, required=False, default=None,
        help='Optional string identifying damage quantities located in "dns-file"')
    parser.add_argument('--max-parallel', type=int, required=False, default=None,
        help='Optional parameter defining the number of parallel processes for the\
              Micromorphic Filter')
    parser.add_argument('--sets-file', type=str, required=False, default=None,
        help='Optional yaml file containing prescribed micro-averaging domains')
    parser.add_argument('--update-filter-domains', type=str, required=False, default=None,
        help='Option to update filter and microaveraging domains for each time step')

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
                                 damage=args.damage,
                                 max_parallel=args.max_parallel,
                                 sets_file=args.sets_file,
                                 update_filter_domains=args.update_filter_domains,
                                 ))