#!python
import argparse
import pathlib
import sys
import yaml

import micromorphic_filter.filter_dns


def run_filter(config_file):
    '''Run the Micromorphic Filter

    :param str config_file: The filter configuration file

    Runs the Micromorphic Filter
    '''

    f = micromorphic_filter.filter_dns.FilterMicroDomain(config_file)
    f.filterIncrements()

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Run the Micromorphic Filter"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('--config-file', type=str, required=True,
        help='Specify the filter configuration file')
    return parser


if __name__ == '__main__':
    parser = get_parser()

    args = parser.parse_args()
    sys.exit(run_filter(config_file=args.config_file))