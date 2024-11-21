#!python
import argparse
import pathlib
import sys

import matplotlib.pyplot
import numpy
import pandas


def plot_force_displacement(csv_file, output_file, output_csv, force_col=' Rz+', length_col=' length_z', force_factor=1, disp_factor=1, filter_markers=None):
    '''Process force-displacement from GEOS DNS results

    :param str csv_file: The csv file containing force results
    :param str output_file: The name of the output file of collected results
    :param str output_csv: The name of the output csv file
    :param int face_id: The face id (or ids) of forces to process
    :param float final_disp: The final displacement (mm) to linearly ramp over simulation duration
    :param str force_col: The column containing desired force information
    :param int header_col: The row containing the headers
    :param float force_factor: The factor to scale force

    :returns: Write ``output_file`` and ``output_csv``
    '''

    df = pandas.read_csv(csv_file, sep=",")
    print(df)
    print(df.columns)

    # calculate displacement, get initial length and subtract from the rest
    length_0 = df[length_col].values[0]
    lengths = df[length_col].values
    disps = disp_factor*numpy.array(lengths - length_0)[:-2]

    # process forces
    forces = force_factor*numpy.array(df[force_col])[:-2]

    # plot
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(disps, forces)
    if filter_markers:
        filter_markers = [int(i) for i in filter_markers]
        matplotlib.pyplot.plot(disps[filter_markers], forces[filter_markers], 'o')
    matplotlib.pyplot.xlabel('Displacement (mm)')
    matplotlib.pyplot.ylabel('Force (N)')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_file)

    # create new dataframe and output csv
    out_df = pandas.DataFrame({'disp': disps,
                               'force': forces})
    out_df.to_csv(output_csv, header=True, sep=',', index=False)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Process force-displacement from GEOS DNS results"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('--csv-file', type=str, required=True,
        help="The csv file containing force results")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output file of collected results")
    parser.add_argument('--output-csv', type=str, required=True,
        help="The name of the output csv file")
    parser.add_argument('--force-col', type=str, required=False, default='force_z',
        help="The column containing desired force information")
    parser.add_argument('--length-col', type=str, required=False, default='length_z',
        help="The column containing domain length data used for calculating displacement")
    parser.add_argument('--force-factor', type=float, required=False, default=1,
        help="The factor to scale force")
    parser.add_argument('--disp-factor', type=float, required=False, default=1,
        help="The factor to scale displacement")
    parser.add_argument('--filter-markers', nargs="+", required=False, default=None,
        help="Optional list of indices to plot markers on force displacement plot \
              corresponding to frames for upscaling")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(plot_force_displacement(csv_file=args.csv_file,
                                     output_file=args.output_file,
                                     output_csv=args.output_csv,
                                     force_col=args.force_col,
                                     length_col=args.length_col,
                                     force_factor=args.force_factor,
                                     disp_factor=args.disp_factor,
                                     filter_markers=args.filter_markers,
                                     ))