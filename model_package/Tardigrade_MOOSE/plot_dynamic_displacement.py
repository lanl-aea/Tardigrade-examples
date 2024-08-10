import pathlib
import argparse
import sys

import pandas
import numpy
import matplotlib.pyplot


def plot_force_displacement(csv_file, output_file, output_csv, final_disp, disp_factor=1):
    '''Process displacement vs time from Tardigrade-MOOSE results

    :param str csv_file: The csv file containing force results
    :param str output_file: The name of the output file of collected results
    :param str output_csv: The name of the output csv file
    :param float disp_factor: The factor to scale displacement

    :returns: Write ``output_file`` and ``output_csv``
    '''

    df = pandas.read_csv(csv_file, sep=",")

    # get times, forces, and displacements
    times = numpy.array(df['time'])
    disps = disp_factor*numpy.array(df['disp_x_p'])

    # plot
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(times, disps)
    matplotlib.pyplot.ylabel('Time (s)')
    matplotlib.pyplot.ylabel('Displacement (mm)')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_file)

    # create new dataframe and output csv
    out_df = pandas.DataFrame({'time': times,
                               'disp': disps})
    out_df.to_csv(output_csv, header=True, sep=',', index=False)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Process displacement vs time from Tardigrade-MOOSE results"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--csv-file', type=str, required=True,
        help="The csv file containing force results")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output file of collected results")
    parser.add_argument('--output-csv', type=str, required=True,
        help="The name of the output csv file")
    parser.add_argument('--disp-factor', type=float, required=False, default=1,
        help="The factor to scale displacement")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(plot_force_displacement(csv_file=args.csv_file,
                                     output_file=args.output_file,
                                     output_csv=args.output_csv,
                                     disp_factor=args.disp_factor,
                                     ))