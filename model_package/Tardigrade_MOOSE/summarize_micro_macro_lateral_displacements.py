import pathlib
import argparse
import sys

import pandas
import numpy
import matplotlib.pyplot

from finite_stVK_calculation import finite_stVK_calculation


def plot_convergence(convergence_plot, y_values, y_label, elements, convergence_value=None):
    '''Plot convergence of a quantity of interest (QoI) vs element count

    :param str convergence_plot: File name for convergence plot
    :param array-like y_values: The QoI values
    :param str y_label: Plot label for QoI
    :param array-lke elements: Array of element counts to plot on x-axis
    :params float convergence_value: Optional "true" value for QoI

    :returns: Write ``convergence_plot``
    '''

    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(elements, y_values, '-o')
    if convergence_value:
        matplotlib.pyplot.plot([elements[0], elements[-1]], [convergence_value, convergence_value], '--', label='solution')
        matplotlib.pyplot.legend()
    matplotlib.pyplot.xlabel('Number of Elements')
    matplotlib.pyplot.ylabel(y_label)
    matplotlib.pyplot.xscale('log')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(convergence_plot, dpi=300)

    return 0


def plot_lateral_displacement(csv_files, plot_labels, output_file, output_csv, convergence_plot=None):
    '''Plot mutliple lateral displacement plots against each other

    :param list csv_files: The csv files containing force results
    :param list plot_labels: The plot labels, same size as ``csv_files``
    :param str output_file: The name of the output file of collected results
    :param str output_csv: The name of the output csv file
    :param str convergence_plot: Optional file name for convergence plot

    :returns: Write ``output_file`` and ``output_csv``, optionally write ``convergence_plot``
    '''

    matplotlib.pyplot.figure()

    dfs = []
    final_disps = []
    # loop through csvs, plot, and append to output DataFrame
    for csv_file, label in zip(csv_files, plot_labels):
        df = pandas.read_csv(csv_file, sep=",")
        matplotlib.pyplot.plot(df['time'], df['lateral_disp'], label=label)
        final_disps.append(df['lateral_disp'].values[-1])
        df = df.assign(id=[label for i in range(len(df.index))])
        dfs.append(df)
    matplotlib.pyplot.xlabel('Time (s)')
    matplotlib.pyplot.ylabel('Lateral Displacement (mm)')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.legend()
    matplotlib.pyplot.savefig(output_file, dpi=300)

    # create new dataframe and output csv
    if output_csv:
        out_df = pandas.concat(dfs, axis=1)
        out_df.to_csv(output_csv, header=True, sep=',', index=False)

    # Convergence plot
    if convergence_plot:
        _, convergence_value = finite_stVK_calculation()
        elements = [int(label.split(' ')[0]) for label in plot_labels]
        plot_convergence(convergence_plot, final_disps, 'lateral disp (mm)', elements, convergence_value)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Plot mutliple lateral displacement plots against each other"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--csv-files', nargs="+", required=True,
        help="The csv files containing force results")
    parser.add_argument('--plot-labels', nargs="+", required=True,
        help="The plot labels, same size as '--csv-files'")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output plot")
    parser.add_argument('--output-csv', type=str, required=True,
        help="The name of the output csv file")
    parser.add_argument('--convergence-plot', type=str, required=False,
        help="Optional file name for convergence plot")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(plot_lateral_displacement(csv_files=args.csv_files,
                                       plot_labels=args.plot_labels,
                                       output_file=args.output_file,
                                       output_csv=args.output_csv,
                                       convergence_plot=args.convergence_plot,
                                       ))