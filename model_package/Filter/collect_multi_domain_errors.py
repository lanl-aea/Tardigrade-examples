import pathlib
import argparse
import sys

import pandas
import numpy as np
import matplotlib.pyplot


def collect_errors(csv_files, num_domains, output_file):
    '''Collect balance equation errors across filter domain studies

    :param list csv_file: A list of csv files containing balance equation errors
    :param list num_domains: A list of integers corresponding to the number of filtering domains associated with results contained in each csv file
    :param str output_file: The name of the output file of collected results

    :returns: ``output_file``
    '''

    fig = matplotlib.pyplot.figure(figsize=(9,6))

    # Create axes
    ax1 = fig.add_subplot(221)   #top left
    ax2 = fig.add_subplot(222)   #top right
    ax3 = fig.add_subplot(223)   #bottom left
    ax4 = fig.add_subplot(224)   #bottom right

    for (csv_file, domain) in zip(csv_files, num_domains):

        df = pandas.read_csv(csv_file, sep=",")

        ax1.plot(df['inc'], df['total_error'], '-o')
        ax1.set_xlabel('Increment')
        ax1.set_ylabel('Total Error')

        ax2.plot(df['inc'], df['rel_error'], '-o')
        ax2.set_xlabel('Increment')
        ax2.set_ylabel('Relative Error')

        ax3.plot(df['inc'], df['max_error'], '-o')
        ax3.set_xlabel('Increment')
        ax3.set_ylabel('Maximum Individual Error')

        ax4.plot(df['inc'], df['max_rel_error'], '-o', label=f"{domain} domain(s)")
        ax4.set_xlabel('Increment')
        ax4.set_ylabel('Maximum Relative Individual Error')

    # Legend
    ax4.legend(bbox_to_anchor=(1.2,0.25))


    # Export figure
    fig.suptitle('Balance Equation Errors')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_file, dpi=300)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Collect balance equation errors across filter domain studies"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--csv-files', nargs="+", required=True,
        help="A list of csv files containing balance equation errors")
    parser.add_argument('--num-domains', nargs="+", required=True,
        help="A list of integers corresponding to the number of filtering domains\
              associated with results contained in each csv file.")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output file of collected results")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(collect_errors(csv_files=args.csv_files,
                            num_domains=args.num_domains,
                            output_file=args.output_file,
                            ))