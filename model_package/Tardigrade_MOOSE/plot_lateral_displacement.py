import pathlib
import argparse
import sys

import pandas
import numpy
import matplotlib.pyplot


def plot_lateral_displacement(csv_file, output_file, output_csv):
    '''Process lateral displacement from Tardigrade-MOOSE results

    :param str csv_file: The csv file containing force results
    :param str output_file: The name of the output file of collected results
    :param str output_csv: The name of the output csv file

    :returns: Write ``output_file`` and ``output_csv``
    '''

    df = pandas.read_csv(csv_file, sep=",")

    # get timesand lateral displacements
    times = numpy.array(df['time'])
    lateral_disp = numpy.array(df['lateral_disp'])

    # plot
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(times, lateral_disp)
    matplotlib.pyplot.xlabel('Time (s)')
    matplotlib.pyplot.ylabel('Lateral Displacement (mm)')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_file)

    # create new dataframe and output csv
    out_df = pandas.DataFrame({'time': times,
                               'lateral_disp': lateral_disp})
    out_df.to_csv(output_csv, header=True, sep=',', index=False)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Process lateral displacement from Tardigrade-MOOSE results"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--csv-file', type=str, required=True,
        help="The csv file containing force results")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output file of collected results")
    parser.add_argument('--output-csv', type=str, required=True,
        help="The name of the output csv file")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(plot_lateral_displacement(csv_file=args.csv_file,
                                       output_file=args.output_file,
                                       output_csv=args.output_csv,
                                       ))