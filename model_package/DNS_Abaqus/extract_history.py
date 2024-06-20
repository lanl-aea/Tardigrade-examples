import sys
import argparse
import pathlib
import yaml

import xarray
import pandas
import matplotlib.pyplot
import h5py
import numpy as np


def plot(input_file, output_file, x_path, y_path, x_label, y_label, x_units, y_units, csv_file=None):
    '''Plot Abaqus history output for force versus displacement

    :param str input_file: Relative or absolute path to h5netcdf file containing Xarray
        Datasets of Abaqus simulation results.
    :param str output_file: The plot file name. Relative or absolute path.
    :param str x_units: The independent (x-axis) units
    :param str y_units: The dependent (y-axis) units
    :param str csv_file: path-like or file-like object containing the CSV dataset to compare with the current
        plot data. If the data sets do not match a non-zero exit code is returned.

    :returns: ``output_file``
    '''

    f = h5py.File(input_file[0])
    X = np.array(f[x_path]).flatten()
    Y = np.array(f[y_path]).flatten()

    # Plot
    #combined_data.sel(selection_dict).plot.scatter(x=x_var, y=y_var, hue=concat_coord)
    matplotlib.pyplot.plot(X, Y, 'o-')
    matplotlib.pyplot.xlabel(f"{x_label} ({x_units})")
    matplotlib.pyplot.ylabel(f"{y_label} ({y_units})")
    matplotlib.pyplot.title(None)
    matplotlib.pyplot.savefig(output_file)

    if csv_file:
        headers = ["disp", "force"]
        output = pandas.DataFrame(np.array([-1*X,-1*Y]).T, columns=headers)
        output.to_csv(csv_file, sep=',', index=False)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Plot Abaqus history output for force versus displacement"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("-i", "--input-file", nargs="+", required=True,
                                help="The Xarray Dataset file(s)")
    required_named.add_argument("--x-path", type=str, required=True,
                                help="The HDF5 path to the x data")
    required_named.add_argument("--y-path", type=str, required=True,
                                help="The HDF5 path to the y data")
    required_named.add_argument("--x-label", type=str, required=True,
                                help="The label (without units) for the x data")
    required_named.add_argument("--y-label", type=str, required=True,
                                help="The label (without units) for the y data.")
    required_named.add_argument("--x-units", type=str, required=True,
                                help="The dependent (x-axis) units string.")
    required_named.add_argument("--y-units", type=str, required=True,
                                help="The independent (y-axis) units string.")
    required_named.add_argument("--csv_file", type=str,
                                help="Name of output CSV file.")

    parser.add_argument("-o", "--output-file", type=str, required=True,
                        help="The output file for plotting")
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args, unknown = parser.parse_known_args()

    sys.exit(plot(input_file=args.input_file,
                  output_file=args.output_file,
                  x_path=args.x_path,
                  y_path=args.y_path,
                  x_label=args.x_label,
                  y_label=args.y_label,
                  x_units=args.x_units,
                  y_units=args.y_units,
                  csv_file=args.csv_file))
