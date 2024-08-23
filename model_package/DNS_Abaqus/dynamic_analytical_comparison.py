import sys
import argparse
import pathlib
import yaml
import math

import xarray
import pandas
import matplotlib.pyplot
import h5py
import numpy


def meirovitch(x, t, c, nmax, L):
    '''Calculate the analytical Meirovitch solution for the dynamic bar problem.

    :param float x: The location on the bar to calculate the solution
    :param linspace t: The discrete time points to calculate the solution
    :param float c: The speed of sound of the bar
    :param int nmax: The number of series expansion terms to calculat the solution
    :param float L: The length of the bar

    :returns: A list of solutions at point `x` for times `t`
    '''

    n_sum = numpy.zeros_like(t)
    for n in range(1, nmax):
        term1 = ((-1)**(n-1))/((2.*n-1)**2)
        trig_fact = ((2.*n-1)*math.pi)/(2.*L)
        sin_term = numpy.sin(trig_fact*x)
        cos_term = numpy.cos(trig_fact*c*t)
        n_sum = n_sum + term1*sin_term*(1. - cos_term)

    return(n_sum)


def plot(input_file, output_file, x_path, y_path, x_label, y_label, x_units, y_units, height, diam, material_E, material_rho, total_force, duration, num_steps, csv_file=None, series_plot=None):
    '''Extracts displacement and reaction force history from Abaqus results. Plots Abaqus results against an analytical solution. Optionally outputs csv_file containing results. Optionally outputs a plot of the converge of the analytical solution.

    :param str input_file: The HDF5 dataset file containing Abaqus results
    :param str output_file: The output file for plotting
    :param str x_path: The HDF5 path to the x data
    :param str y_path: The HDF5 path to the y data
    :param str x_label: The label (without units) for the x data
    :param str y_label: The label (without units) for the y data
    :param str x_units: The independent (x-axis) units string
    :param str y_units: The dependent (y-axis) units string
    :param float height: The height (mm) of the cylinder. This values will be multiplied by 1.e-3 to convert to units of m
    :param float diam: The diameter (mm) of the cylinder. This values will be multiplied by 1.e-3 to convert to units of m
    :param float material_E: The elastic modulus (MPa) of the material. This value will be multiplied by 1.6 to convert to units of Pa
    :param float material_rho: The density (g/cm^3) of the material. This value will be multiplied by 1.00e3 to convert to units of kg/m^3
    :param float total_force: The force (N) applied to cylinder
    :param float duration: The duration of the simulation
    :param int num_steps: The number of fixed time increments
    :param str csv_file: Name of output CSV file
    :param str series_plot: Name of the output series convergence plot for summation terms

    :returns: ``output_file`` and optionally ``csv_file`` and ``series_plot``
    '''

    f = h5py.File(input_file[0])
    X = numpy.array(f[x_path]).flatten()
    Y = numpy.array(f[y_path]).flatten()

    # generate analytical solution
    density = material_rho*1.e3
    d = diam*1.e-3
    L = height*1.e-3
    area = 0.25*math.pi*d*d
    Young = material_E*1.e6
    speed = numpy.sqrt(Young/density)
    x=L
    param1=8*total_force*L/(Young*area*math.pi*math.pi)
    time = numpy.linspace(0, duration, num_steps)
    u_sum = 0
    u_sum = meirovitch(x, time, speed, 201, L)
    u_displ = 1000*param1*u_sum
    # Plot
    matplotlib.pyplot.plot(X, -1*Y, '-', label='numerical')
    matplotlib.pyplot.plot(time, u_displ, '-', label='analytical')
    matplotlib.pyplot.xlabel(f"{x_label} ({x_units})")
    matplotlib.pyplot.ylabel(f"-1*{y_label} ({y_units})")
    matplotlib.pyplot.legend()
    matplotlib.pyplot.title(None)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    matplotlib.pyplot.savefig(output_file)

    if csv_file:
        headers = ["numerical_time, numerical_disp, analytical_time, analytical_disp"]
        output = pandas.DataFrame(numpy.array([X,-1*Y,time,1000*u_displ]).T, columns=headers)
        output.to_csv(csv_file, sep=',', index=False)

    if series_plot:
        f1 = (1/(4*L))*speed
        t = 1/(2*f1)
        plot_series_convergence(L, t, speed, L, series_plot)

    return 0


def plot_series_convergence(x, t, speed, L, output_name):
    '''Plot the convergence of the Meirovitch solution over a range of terms

    :param float x: The location on the bar to calculate the solution
    :param linspace t: The discrete time points to calculate the solution
    :param float speed: The speed of sound of the bar
    :param float L: The length of the bar
    :param str output_name: Name of the output series convergence plot for summation terms

    :returns: ``output_name``
    '''
    ns = list(range(1,1001))
    ns = ns + [2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 50000, 100000]
    #solutions = numpy.zeros_like(ns)
    solutions = []
    #for i, ni in enumerate(ns):
    for ni in ns:
        #solutions[i] = meirovitch(x, t, speed, ni, L)
        solutions.append(meirovitch(x, t, speed, ni, L))

    print(f'converged sum = {solutions[-1]}')
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(ns, solutions)
    ax.set_xscale('log')
    ax.set_xlabel('n')
    ax.set_ylabel(r"Sum, $\bar{N} \left(x = L, t = \frac{T}{2}\right)$")
    matplotlib.pyplot.title(None)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_name)

    return 0


def get_parser():
    script_name = pathlib.Path(__file__)
    default_output_file = f"{script_name.stem}.png"

    prog = f"python {script_name.name} "
    cli_description = "Plot dynamic Abaqus results against an analytical solution"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("-i", "--input-file", nargs="+", required=True,
        help="The HDF5 dataset file containing Abaqus results")
    parser.add_argument("-o", "--output-file", type=str, default=default_output_file,
        help="The output file for plotting")
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
    parser.add_argument('--diam', type=float,
        help='Specify the diameter (mm) of the cylinder. This values will be multiplied by 1.e-3 to convert to units of m')
    parser.add_argument('--height', type=float,
        help='Specify the height (mm) of the cylinder. This values will be multiplied by 1.e-3 to convert to units of m')
    parser.add_argument('--material-E', type=float,
        help='Specify the elastic modulus (MPa) of the material. This value will be multiplied by 1.6 to convert to units of Pa.')
    parser.add_argument('--material-rho', type=float,
        help='Specify the density (g/cm^3) of the material. This value will be multiplied by 1.00e3 to convert to units of kg/m^3')
    parser.add_argument('--total-force', type=float,
        help='Specify the force (N) applied to cylinder.')
    parser.add_argument('--duration', type=float,
        help='Specify the duration of the simulation.')
    parser.add_argument('--num-steps', type=int,
        help='Specify the number of fixed time increments.')
    required_named.add_argument("--csv-file", type=str,
        help="Name of output CSV file.")
    required_named.add_argument("--series-plot", type=str,
        help="Name of the output series convergence plot for summation terms.")

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
                  diam=args.diam,
                  height=args.height,
                  material_E=args.material_E,
                  material_rho=args.material_rho,
                  total_force=args.total_force,
                  duration=args.duration,
                  num_steps=args.num_steps,
                  csv_file=args.csv_file,
                  series_plot=args.series_plot,
                  ))