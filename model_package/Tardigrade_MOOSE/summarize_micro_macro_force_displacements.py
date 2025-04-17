#!python
import argparse
import pathlib
import sys

import matplotlib.pyplot
import pandas

from finite_stVK_calculation import finite_stVK_calculation
from summarize_micro_macro_lateral_displacements import plot_convergence


def plot_force_displacement(csv_files, plot_labels, output_file, output_csv, convergence_plot=None,
                            force_field='force', disp_field='disp',
                            x_label='Displacement (mm)', y_label='Force (N)'):
    '''Plot multiple force displacement plots against each other

    :param list csv_files: The csv files containing force results
    :param list plot_labels: The plot labels, same size as ``csv_files``
    :param str output_file: The name of the output file of collected results
    :param str output_csv: The name of the output csv file
    :param str convergence_plot: Optional file name for convergence plot
    :param str force_field: The column label for force values
    :param str disp_field: The column label for displacement values
    :param str x_label: The label (without units) for the x data
    :param str y_label: The label (without units) for the y data

    :returns: Write ``output_file`` and ``output_csv``
    '''

    matplotlib.pyplot.figure()

    dfs = []
    final_forces = []
    # loop through csvs, plot, and append to output DataFrame
    for csv_file, label in zip(csv_files, plot_labels):
        df = pandas.read_csv(csv_file, sep=",")
        matplotlib.pyplot.plot(df[disp_field], df[force_field], label=label)
        final_forces.append(df[force_field].values[-1])
        df = df.assign(id=[label for i in range(len(df.index))])
        dfs.append(df)
    matplotlib.pyplot.xlabel(x_label)
    matplotlib.pyplot.ylabel(y_label)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.legend()
    matplotlib.pyplot.savefig(output_file, dpi=300)

    # create new dataframe and output csv
    if output_csv:
        out_df = pandas.concat(dfs, axis=1)
        out_df.to_csv(output_csv, header=True, sep=',', index=False)

    # Convergence plot
    if convergence_plot:
        convergence_value, _ = finite_stVK_calculation()
        elements = [int(label.split(' ')[0]) for label in plot_labels]
        plot_convergence(convergence_plot, final_forces, '-1*Force (N)', elements, -1*convergence_value)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Plot mutliple force displacement plots against each other"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('--csv-files', nargs="+", required=True,
        help="The csv files containing force results")
    parser.add_argument('--plot-labels', nargs="+", required=True,
        help="The plot labels, same size as '--csv-files'")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output plot")
    parser.add_argument('--output-csv', type=str, required=True,
        help="The name of the output csv file")
    parser.add_argument('--convergence-plot', type=str, required=False, default=None,
        help="Optional file name for convergence plot")
    parser.add_argument('--force-field', type=str, required=False, default='force',
        help="The column label for force values")
    parser.add_argument('--disp-field', type=str, required=False, default='disp',
        help="Optional column label for displacement values")
    parser.add_argument("--x-label", type=str, required=False, default='Displacement (mm)',
        help="The label for the x data")
    parser.add_argument("--y-label", type=str, required=False, default='Force (N)',
        help="The label for the y data")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(plot_force_displacement(csv_files=args.csv_files,
                                     plot_labels=args.plot_labels,
                                     output_file=args.output_file,
                                     output_csv=args.output_csv,
                                     convergence_plot=args.convergence_plot,
                                     force_field=args.force_field,
                                     disp_field=args.disp_field,
                                     x_label=args.x_label,
                                     y_label=args.y_label,
                                     ))