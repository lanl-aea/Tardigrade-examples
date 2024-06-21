import pathlib
import argparse
import sys

import pandas
import numpy as np
import matplotlib.pyplot as plt


def collect_results(csv_files, num_domains, output_file=None, box_plot=None, narrow=False):
    '''Collect statistics of a homogenized micromorphic quantity across filter domain studies

    :param list csv_files: A list of csv files containing information to collect
    :param list num_domains: A list of integers corresponding to the number of filtering domains associated with results contained in each csv file
    :param str output_file: The name of the output file of collected results
    :param str box_plot: The name of an optional box and whisker plot
    :param str narrow: Optional flag to make a narrow box plot

    :returns: Summary csv file named ``output_file`` and optionally a violin plot named ``box_plot``
    '''

    dfs = []
    for (csv_file, domain) in zip(csv_files, num_domains):
        csv_data = pandas.read_csv(csv_file, sep=",")
        csv_data = csv_data.assign(num_domains=[domain for i in range(len(csv_data.index))])
        dfs.append(csv_data)

    if output_file:
        new_df = pandas.concat(dfs).set_index(['component','num_domains']).unstack()
        new_df.to_csv(output_file)

    if box_plot:
        if narrow == False:
            fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(6,6))
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3,9))
        data = []
        for df in dfs:
            data.append(df['quantity'])
        ax.violinplot(data,showmeans=True, widths=0.7)
        ax.set_xticks([y+1 for y in range(len(data))], labels=num_domains)
        ax.set_xlabel('Number of filter domains')
        ax.set_ylabel(r'$\sigma_{33}$ (MPa)')
        fig.tight_layout()
        fig.savefig(box_plot)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Collect statistics of a homogenized micromorphic quantity across filter domain studies"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--csv-files', nargs="+", required=True,
        help="A list of csv files containing information to collect")
    parser.add_argument('--num-domains', nargs="+", required=True,
        help="A list of integers corresponding to the number of filtering domains\
              associated with results contained in each csv file.")
    parser.add_argument('--output-file', type=str, required=False,
        help="The name of the output file of collected results")
    parser.add_argument('--box-plot', type=str, required=False,
        help="The name of an optional box and whisker plot")
    parser.add_argument('--narrow', type=str, required=False,
        help="Optional flag to make a narrow box plot")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(collect_results(csv_files=args.csv_files,
                             num_domains=args.num_domains,
                             output_file=args.output_file,
                             box_plot=args.box_plot,
                             narrow=bool(args.narrow),
                             ))