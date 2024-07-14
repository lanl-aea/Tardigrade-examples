import pathlib
import argparse
import sys

import pandas
import numpy
import matplotlib.pyplot


def plot_force_displacement(csv_files, plot_labels, output_file, output_csv):
    '''Plot mutliple force displacement plots against each other

    :param list csv_files: The csv files containing force results
    :param list plot_labels: The plot labels, same size as ``csv_files``
    :param str output_file: The name of the output file of collected results
    :param str output_csv: The name of the output csv file

    :returns: Write ``output_file`` and ``output_csv``
    '''

    matplotlib.pyplot.figure()

    dfs = []
    # loop through csvs, plot, and append to output DataFrame
    for csv_file, label in zip(csv_files, plot_labels):
        df = pandas.read_csv(csv_file, sep=",")
        matplotlib.pyplot.plot(df['disp'], df['force'], label=label)
        df = df.assign(id=[label for i in range(len(df.index))])
        dfs.append(df)
    matplotlib.pyplot.xlabel('Displacement (mm)')
    matplotlib.pyplot.ylabel('Force (N)')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.legend()
    matplotlib.pyplot.savefig(output_file, dpi=300)

    # create new dataframe and output csv
    if output_csv:
        out_df = pandas.concat(dfs, axis=1)
        out_df.to_csv(output_csv, header=True, sep=',', index=False)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Plot mutliple force displacement plots against each other"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    # TODO: handle multiple faces to sum forces if needed
    parser.add_argument('--csv-files', nargs="+", required=True,
        help="The csv files containing force results")
    parser.add_argument('--plot-labels', nargs="+", required=True,
        help="The plot labels, same size as '--csv-files'")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output plot")
    parser.add_argument('--output-csv', type=str, required=True,
        help="The name of the output csv file")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(plot_force_displacement(csv_files=args.csv_files,
                                     plot_labels=args.plot_labels,
                                     output_file=args.output_file,
                                     output_csv=args.output_csv,
                                     ))