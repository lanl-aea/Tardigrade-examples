import pathlib
import argparse
import sys
import math

import pandas
import numpy
import matplotlib.pyplot

from model_package.DNS_Abaqus.dynamic_analytical_comparison import meirovitch
from model_package.Tardigrade_MOOSE import simulation_variables_nominal


def summarize_dynamic_displacements(csv_files, plot_labels, output_file, output_csv, disp_factor=1):
    '''Plot mutliple dynamic displacement plots against each other

    :param list csv_files: The csv files containing force results
    :param list plot_labels: The plot labels, same size as ``csv_files``
    :param str output_file: The name of the output file of collected results
    :param str output_csv: The name of the output csv file
    :param float disp_factor: The factor to scale displacement

    :returns: Write ``output_file`` and ``output_csv``, optionally write ``convergence_plot``
    '''

    

    # Meirovitch solution
    params = simulation_variables_nominal.dynamic_elastic_cylinder
    density = params['material_rho']*1.e12
    d = params['diam']*1.e-3
    L = params['height']*1.e-3
    area = 0.25*math.pi*d*d
    Young = params['material_E']*1.e6
    speed = numpy.sqrt(Young/density)
    total_force = params['pressure']*area*1.e6
    x=L
    param1=8*total_force*L/(Young*area*math.pi*math.pi)
    time = numpy.linspace(0, params['duration'], params['num_steps'])
    u_sum = meirovitch(x, time, speed, 201, L)
    u_displ = 1000*param1*u_sum

    # plot
    fig, (ax1, ax2) = matplotlib.pyplot.subplots(1,2,figsize=[6.5,3.25], gridspec_kw={'width_ratios': [4, 1]})
    ax1.plot(1000*time, disp_factor*u_displ, 'k--', label='Analyical Solution')

    dfs = []
    # loop through csvs, plot, and append to output DataFrame
    for csv_file, label in zip(csv_files, plot_labels):
        df = pandas.read_csv(csv_file, sep=",")
        ax1.plot(1000*df['time'], disp_factor*df['disp'], label=label)
        df = df.assign(id=[label for i in range(len(df.index))])
        dfs.append(df)
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Displacement (mm)')
    ax2.legend(*ax1.get_legend_handles_labels(), loc='center')
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    ax2.set_frame_on(False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300)

    # create new dataframe and output csv
    if output_csv:
        out_df = pandas.concat(dfs, axis=1)
        out_df.to_csv(output_csv, header=True, sep=',', index=False)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Plot mutliple dynamic displacement plots against each other"
    parser=argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--csv-files', nargs="+", required=True,
        help="The csv files containing force results")
    parser.add_argument('--plot-labels', nargs="+", required=True,
        help="The plot labels, same size as '--csv-files'")
    parser.add_argument('--output-file', type=str, required=True,
        help="The name of the output plot")
    parser.add_argument('--output-csv', type=str, required=True,
        help="The name of the output csv file")
    parser.add_argument('--disp-factor', type=float, required=False, default=1,
        help="The factor to scale displacement")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(summarize_dynamic_displacements(csv_files=args.csv_files,
                                             plot_labels=args.plot_labels,
                                             output_file=args.output_file,
                                             output_csv=args.output_csv,
                                             disp_factor=args.disp_factor,
                                             ))