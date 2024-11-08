import subprocess as sp
import numpy
import os
import sys
import argparse
import time
import glob
import yaml
import inspect

import seaborn
import matplotlib.pyplot
import pandas


def write_plastic_material_card(output_file, input_dict):
    '''Write elastic micromorphic material card

    :param str output_file: The root filename of the output yaml file
    :param dict input_dict: A dictionary containing calibrated parameters

    Returns: Writes ``output_file``.yml
    '''
    multiplier = 10
    defaults = {
        'cu0': 3.192202765, 'Hu': 1e-8,
        'cchi0': 1e-8, 'Hchi': 1e-8,
        'friction': 0,
        'cnablachi0': 1e-8, 'Hnablachi': 1e-8,
        'lambda': 0.0, 'mu': 0.0, 'eta': 0.0, 'tau': 0.0,
        'kappa': 0.0, 'nu': 0.0, 'sigma': 0.0, 'tau1': 0.0,
        'tau2': 0.0, 'tau3': 0.0, 'tau4': 0.0, 'tau5': 0.0,
        'tau6': 0.0, 'tau7': 0.0001, 'tau8': 0.0, 'tau9': 0.0,
        'tau10': 0.0, 'tau11': 0.0,
        'int_params_a': 0.5, 'int_params_b': 1e-9}

    # Use defaults if a parameter is not specified
    for key in defaults.keys():
        if key not in input_dict.keys():
            input_dict[key] = defaults[key]

    output_dict = {}
    # plastic
    output_dict['line 01'] = f"2 {input_dict['cu0']} {input_dict['Hu']}"
    output_dict['line 02'] = f"2 {input_dict['cchi0']} {input_dict['Hchi']}"
    output_dict['line 03'] = f"2 {input_dict['cnablachi0']} {input_dict['Hnablachi']}"
    output_dict['line 04'] = "2 0. 0."
    output_dict['line 05'] = "2 0. 0."
    output_dict['line 06'] = "2 0. 0."
    output_dict['line 07'] = "2 0. 0."
    output_dict['line 08'] = "2 0. 0."
    output_dict['line 09'] = "2 0. 0."
    # elastic
    output_dict['line 10'] = f"2 {input_dict['lambda']} {input_dict['mu']}"
    output_dict['line 11'] = f"5 {input_dict['eta']} {input_dict['tau']} {input_dict['kappa']} {input_dict['nu']} {input_dict['sigma']}"
    output_dict['line 12'] = f"11 {input_dict['tau1']} {input_dict['tau2']} {input_dict['tau3']} {input_dict['tau4']} {input_dict['tau5']} {input_dict['tau6']} {input_dict['tau7']} {input_dict['tau8']} {input_dict['tau9']} {input_dict['tau10']} {input_dict['tau11']}"
    output_dict['line 13'] = f"2 {input_dict['tau']} {input_dict['sigma']}"
    # integration
    output_dict['line 14'] = '0.5 0.5 0.5 1e-9 1e-9'
    with open(f'{output_file}.yml', 'w') as f:
        yaml.dump(output_dict, f)

    return 0


def make_summary_csv(summary_csv, parameter_df):
    '''Make a csv file summarizing the mean, min, max, and standard deviation of calibrated parameters

    :param str summary_csv: Filename to store summary statistics of calibrated parameters
    :param dict results_dict: Results dictionary containing list of parameters with each key corresponding to a parameter name

    :returns: ``summary_csv``
    '''

    params, means, mins, maxs, devs = [], [], [], [], []
    fields = list(parameter_df.columns)
    fields.remove('element')
    for key in fields:
        # get stats
        params.append(key)
        means.append(numpy.mean(parameter_df[key]))
        mins.append(numpy.min(parameter_df[key]))
        maxs.append(numpy.max(parameter_df[key]))
        devs.append(numpy.std(parameter_df[key]))
        # output 
        df = pandas.DataFrame({'param': params,
                               'mean': means,
                               'min': mins,
                               'max': maxs,
                               'dev': devs})
    df.to_csv(summary_csv, header=True, sep=',', index=False)

    return 0


def kde(rootname, parameter_df, type, kde_best_parameters=None):
    '''Create a kernel density estimate (KDE) plot for each calibrated parameter

    :param str rootname: The rootname of the output plot
    :param dict results_dict: Dictionary containing list of parameter results with each key corresponding to a parameter
    :param str type: A string specifying the type of KDE to plot. 'kde' gives a regular KDE plot. 'hist' gives a KDE plot with histograms shown
    :param str kde_best_parameters: Optional root filename to output a yaml file containing the "best" parameters sampled from the kernel density estimate associated with "kde_best"

    :returns: ``{rootname}_{key}_{type}.PNG`` for each key in `results_dict`, write ``{kde_best_parameters}.yml`` if requested
    '''

    output_parameters = {}
    fields = list(parameter_df.columns)
    fields.remove('element')
    for key in fields:
        matplotlib.pyplot.figure()
        if type != 'best':
            if type == 'kde':
                ax = seaborn.displot(parameter_df[key], kind='kde', color='red', fill=True)
            elif type == 'hist':
                ax = seaborn.displot(parameter_df[key], kde=True, color='red', fill=True)
            ax.set(xlabel=f'parameter {key}', ylabel='KDE')
            matplotlib.pyplot.title(key)
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig(f'{rootname}_{key}_{type}.PNG')
        else:
            try:
                ax = seaborn.kdeplot(parameter_df[key], color='red', fill=False)
                xs, ys = ax.lines[-1].get_data()
                ax.fill_between(xs, ys, color='red', alpha=0.1)
                mode_idx = numpy.argmax(ys)
                output_parameters[key] = xs[mode_idx]
                best_value = f'{xs[mode_idx]:.6f}'
                ax.set(xlabel=f'parameter {key}', ylabel='KDE')
                matplotlib.pyplot.title(f'Best {key} = {best_value}')
                matplotlib.pyplot.tight_layout()
                matplotlib.pyplot.savefig(f'{rootname}_{key}.png')
            except:
                pass
        matplotlib.pyplot.close()
    # output
    if kde_best_parameters:
        #write_elastic_material_card(kde_best_parameters, output_parameters)
        write_plastic_material_card(kde_best_parameters, output_parameters)
    return 0


def summarize_calibration_results(parameter_csv,
                                  summary_csv=None,
                                  kde_hist_plot=None,
                                  kde_plot=None,
                                  kde_best=None,
                                  kde_best_parameters=None,):
    '''Main function to drive parameter summary and output

    :param list parameter_sets: List of yaml files containing calibration results
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter,\
              3: 7 parameter plus tau7, 4: all 18 parameters
    :param str results_csv: Optional filename to store all calibrated parameter values
    :param str summary_csv: Optional filename to store summary statistics of calibrated parameters
    :param str kde_hist_plot: Optional root filename to plot kernel density estimate of each calibrated parameter with histogram
    :param str kde_plot: Optional root filename to plot kernel density estimate of each calibrated parameter
    :param str kde_best_parameters: Optional root filename to output a yaml file containing the "best" parameters sampled from the kernel density estimate associated with "kde_best"
    '''

    parameter_df = pandas.read_csv(parameter_csv)

    if summary_csv:
        make_summary_csv(summary_csv, parameter_df)
        
    if kde_hist_plot:
        kde(kde_hist_plot, parameter_df, 'hist')

    if kde_plot:
        kde(kde_plot, parameter_df, 'kde')

    if kde_best:
        kde(kde_best, parameter_df, 'best', kde_best_parameters)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Summarize results of parameter calibration from a calibration map csv"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('--parameter-csv', nargs="+", required=True,
        help='Specify the list of yaml files containing calibration results')
    parser.add_argument('--summary-csv', type=str, required=False,
        help='Optional filename to store summary statistics of calibrated parameters')
    parser.add_argument('--kde-hist-plot', type=str, required=False,
        help='Optional root filename to plot kernel density estimate of each calibrated parameter with histogram')
    parser.add_argument('--kde-plot', type=str, required=False,
        help='Optional root filename to plot kernel density estimate of each calibrated parameter')
    parser.add_argument('--kde-best', type=str, required=False,
        help='Optional root filename to plot kernel density estimate of each calibrated parameter with maximum value in title')
    parser.add_argument('--kde-best-parameters', type=str, required=False, default=None,
        help='Optional root filename to output a yaml file containing the "best" parameters sampled from the kernel density estimate associated with "--kde-best"')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(summarize_calibration_results(parameter_csv=args.parameter_csv,
                                           summary_csv=args.summary_csv,
                                           kde_hist_plot=args.kde_hist_plot,
                                           kde_plot=args.kde_plot,
                                           kde_best=args.kde_best,
                                           kde_best_parameters=args.kde_best_parameters,
                                           ))