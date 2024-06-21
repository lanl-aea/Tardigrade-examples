import sys
import pathlib
import argparse

import numpy
import pandas
import matplotlib.pyplot


def plot_objective_evaluation(csv_file, output_file):
    '''Plot calibrated lambda versus mu with transparency based on objective function evaluation

    :param str csv_file: The csv file containing objective function evaluations and parameter sets
    :param str output_file: The name of the output plot

    :returns: Write ``output_file``
    '''

    data = pandas.read_csv(csv_file, sep=',')

    obj = data['obj']
    lambs = data['lamb']
    mus = data['mu']

    # Get rid of some objective functions
    indices = numpy.where(numpy.log10(obj) <= numpy.mean(numpy.log10(obj)))
    good_obj = numpy.log10(obj.iloc[indices])
    good_lambs = lambs.iloc[indices]
    good_mus = mus.iloc[indices]

    # map alphas
    a = numpy.min(good_obj)
    b = numpy.max(good_obj)
    obj_out = (b - good_obj) / (b - a)

    # plot
    matplotlib.pyplot.figure()
    for x, y, c in zip(good_lambs, good_mus, obj_out):
        matplotlib.pyplot.plot(x, y, 'o', color='r', alpha=c)
    matplotlib.pyplot.plot(lambs.iloc[-1], mus.iloc[-1], '*', color='k')
    matplotlib.pyplot.xlabel('lambas (MPa)')
    matplotlib.pyplot.ylabel('mu (MPa)')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_file)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Plot calibrated lambda versus mu with transparency based on objective function evaluation"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--csv-file', type=str, required=True,
        help="The csv file containing objective function evaluations and parameter sets")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output plot")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(plot_objective_evaluation(csv_file=args.csv_file,
                                       output_file=args.output_file,
                                       ))