import pathlib
import argparse
import sys

import pandas
import numpy
import matplotlib.pyplot


def plot_force_displacement(csv_file, output_file, output_csv, face_id, final_disp, force_col='force_z', header_row=0, force_factor=1):
    '''Process force-displacement from Ratel DNS results

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

    df = pandas.read_csv(csv_file, sep=",", header=header_row)

    # get only the face_id of interest
    df = df.iloc[numpy.where(df['face_id'] == face_id)[0]]

    # get times, forces, and displacements
    times = numpy.array(df['time'])
    disps = final_disp*times
    if type(df[force_col].iloc[0]) is float:
        forces = force_factor*numpy.array(df[force_col])
    else:
        forces = force_factor*numpy.array([float(str(f).replace('\\','')) for f in df[force_col]])

    # plot
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(disps, forces)
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
    cli_description = "Process force-displacement from Ratel DNS results"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    # TODO: handle multiple faces to sum forces if needed
    parser.add_argument('--csv-file', type=str, required=True,
        help="The csv file containing force results")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output file of collected results")
    parser.add_argument('--output-csv', type=str, required=True,
        help="The name of the output csv file")
    parser.add_argument('--face-id', type=int, required=True,
        help="The face id (or ids) of forces to process")
    parser.add_argument('--final-disp', type=float, required=True,
        help="The final displacement (mm) to linearly ramp over simulation duration")
    parser.add_argument('--force-col', type=str, required=False, default='force_z',
        help="The column containing desired force information")
    parser.add_argument('--header-row', type=int, required=False, default=0,
        help="The row containing the headers")
    parser.add_argument('--force-factor', type=float, required=False, default=1,
        help="The factor to scale force")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(plot_force_displacement(csv_file=args.csv_file,
                                     output_file=args.output_file,
                                     output_csv=args.output_csv,
                                     face_id=args.face_id,
                                     final_disp=args.final_disp,
                                     force_col=args.force_col,
                                     header_row=args.header_row,
                                     force_factor=args.force_factor,
                                     ))