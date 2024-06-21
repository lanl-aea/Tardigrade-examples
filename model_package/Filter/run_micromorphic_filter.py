import os
import inspect
import sys
import yaml
import argparse

import micromorphic_filter.filter_dns


def run_filter(config_file):
    '''Run the Micromorphic Filter

    :param str config_file: The filter configuration file

    Runs the Micromorphic Filter
    '''

    f = micromorphic_filter.filter_dns.Filter(config_file)
    f.filterIncrements()

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Run the Micromorphic Filter"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('--config-file', type=str, required=True,
        help='Specify the filter configuration file')
    return parser


if __name__ == '__main__':
    parser = get_parser()

    args = parser.parse_args()
    sys.exit(run_filter(config_file=args.config_file))