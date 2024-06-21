import inspect
import sys
import os
import argparse

import numpy as np
import pandas
import matplotlib.pyplot as plt


def plot_errors(df, output_plot):
    ''' Plot balance equation errors

    :param dataframe df: Pandas dataframe containing balance equation errors
    :param str output_plot: Output filename of plot

    :returns: ``{output_plot}``
    '''

    # Create new figure
    fig = plt.figure(figsize=(6,6))

    # Create axes
    ax1 = fig.add_subplot(221)   #top left
    ax2 = fig.add_subplot(222)   #top right
    ax3 = fig.add_subplot(223)   #bottom left
    ax4 = fig.add_subplot(224)   #bottom right

    ax1.plot(df['inc'], df['total_error'], '-o')
    ax1.set_xlabel('Increment')
    ax1.set_ylabel('Total Error')

    ax2.plot(df['inc'], df['rel_error'], '-o')
    ax2.set_xlabel('Increment')
    ax2.set_ylabel('Relative Error')

    ax3.plot(df['inc'], df['max_error'], '-o')
    ax3.set_xlabel('Increment')
    ax3.set_ylabel('Maximum Individual Error')

    ax4.plot(df['inc'], df['max_rel_error'], '-o')
    ax4.set_xlabel('Increment')
    ax4.set_ylabel('Maximum Relative Individual Error')

    # Export figure
    fig.suptitle('Balance Equation Errors')
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)

    return 0


def parse_errors(input_file, output_csv, output_plot=None):
    ''' Parse balance equation errors from Micromorphic Filter standard output

    :param str input_file: The standard out file produced when running the Micromorphic Filter
    :param str output_csv: Name of output csv file summarizing output for each timestep
    :param str output_plot: Optional filename to plot balance equation errors

    :returns: ``{output_csv}`` and optionally ``{output_plot}``
    '''

    # Read in specified file as array
    with open(input_file, 'r') as f:
        input_contents = np.asarray(f.read().splitlines())

    # parse contents into dictionary
    results = {}
    results['inc'] = []
    results['total_error'] = []
    results['rel_error'] = []
    results['max_error'] = []
    results['max_rel_error'] = []
    for index, text in enumerate(input_contents):
        # Find beginning of an increment
        if "Processing increment:" in text:
            results['inc'].append(int(text.split(": ")[1]))
        # Find errors
        if "total balance equation error" in text:
            results['total_error'].append(float(text.split(": ")[1]))
        if "total relative percent balance equation error" in text:
            results['rel_error'].append(float(text.split(": ")[1]))
        if "maximum individual error" in text:
            results['max_error'].append(float(text.split(": ")[1]))
        if "maximum individual relative percent error" in text:
            results['max_rel_error'].append(float(text.split(": ")[1]))

    # create DataFrame
    df = pandas.DataFrame(results)

    # Output csv
    df.to_csv(output_csv, header=True, sep=',', index=False)

    # Optional plot
    if output_plot:
        plot_errors(df, output_plot)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Parse balance equation errors from Micromorphic Filter\
                      standard output"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
                        help="The standard out file produced when running the Micromorphic Filter")
    parser.add_argument('--output-csv', type=str, required=True,
                        help="Name of output csv file summarizing output for each timestep")
    parser.add_argument('--output-plot', type=str, required=False, default=None,
                        help="Optional filename to plot balance equation errors")

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(parse_errors(input_file=args.input_file,
                          output_csv=args.output_csv,
                          output_plot=args.output_plot,
                          ))