#!python
import argparse
import pathlib
import sys

import matplotlib.pyplot
import numpy
import pandas


def Brazil_disk_normalized_force_vs_displacements(input_file, plot_file, csv_file, radius, thickness, disp_factor=1, force_factor=1):
    '''Process force-displacement from Tardigrade-MOOSE results

    :param str input_file: The csv file containing force vs displacement results
    :param str plot_file: The name of the output file of collected resultss
    :param str csv_file: The name of the output csv file
    :param float radius: The specimen initial radius
    :param float thickness: The specimen initial thickness
    :param float disp_factor: The final displacement (mm) to linearly ramp over simulation duration
    :param float force_factor: The factor to scale force

    :returns: Write ``plot_file`` and ``csv_file``
    '''

    df = pandas.read_csv(input_file, sep=",")

    # get times, forces, and displacements
    disps = 100*disp_factor*numpy.array(df['disp'])/(2*radius)
    forces = force_factor*numpy.array(df['force'])/(radius*thickness)

    # plot
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(disps, forces, 'o-')
    matplotlib.pyplot.xlabel('Displacement, $\delta / D$ (%)')
    matplotlib.pyplot.ylabel(r'Applied Load, $P / \left(BR\right)$ (MPa)')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(plot_file)

    # create new dataframe and output csv
    out_df = pandas.DataFrame({'norm_disp': disps,
                               'norm_force': forces})
    out_df.to_csv(csv_file, header=True, sep=',', index=False)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Process force-displacement from Tardigrade-MOOSE results"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('--input-file', type=str, required=True,
        help="The csv file containing force vs displacement results")
    parser.add_argument('--radius', type=float, required=False, default=1,
        help="The specimen initial radius")
    parser.add_argument('--thickness', type=float, required=False, default=1,
        help="The specimen initial thickness")
    parser.add_argument('--force-factor', type=float, required=False, default=1,
        help="The factor to scale force")
    parser.add_argument('--disp-factor', type=float, required=False, default=1,
        help="The factor to scale displacement")
    parser.add_argument('--plot-file', type=str, required=True,
        help="The name of the output file of collected results")
    parser.add_argument('--csv-file', type=str, required=True,
        help="The name of the output csv file")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(Brazil_disk_normalized_force_vs_displacements(input_file=args.input_file,
                                                           radius=args.radius,
                                                           thickness=args.thickness,
                                                           force_factor=args.force_factor,
                                                           disp_factor=args.disp_factor,
                                                           plot_file=args.plot_file,
                                                           csv_file=args.csv_file,
                                                           ))