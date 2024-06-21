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


def write_elastic_material_card(output_file, input_dict):
    '''Write elastic micromorphic material card

    :param str output_file: The root filename of the output yaml file
    :param dict input_dict: A dictionary containing calibrated parameters

    Returns: Writes ``output_file``.yml
    '''

    defaults = {
        'lamb': 0.0, 'mu': 0.0, 'eta': 0.0, 'tau': 0.0,
        'kappa': 0.0, 'nu': 0.0, 'sigma': 0.0, 'tau1': 0.0,
        'tau2': 0.0, 'tau3': 0.0, 'tau4': 0.0, 'tau5': 0.0,
        'tau6': 0.0, 'tau7': 0.0001, 'tau8': 0.0, 'tau9': 0.0,
        'tau10': 0.0, 'tau11': 0.0}

    # Use defaults if a parameter is not specified
    for key in defaults.keys():
        if key not in input_dict.keys():
            input_dict[key] = defaults[key]

    # elastic
    output_dict = {}
    output_dict['line 1'] = f"2 {input_dict['lamb']} {input_dict['mu']}"
    output_dict['line 2'] = f"5 {input_dict['eta']} {input_dict['tau']} {input_dict['kappa']} {input_dict['nu']} {input_dict['sigma']}"
    output_dict['line 3'] = f"11 {input_dict['tau1']} {input_dict['tau2']} {input_dict['tau3']} {input_dict['tau4']} {input_dict['tau5']} {input_dict['tau6']} {input_dict['tau7']} {input_dict['tau8']} {input_dict['tau9']} {input_dict['tau10']} {input_dict['tau11']}"
    output_dict['line 4'] = f"2 {input_dict['tau']} {input_dict['sigma']}"

    with open(f'{output_file}.yml', 'w') as f:
        yaml.dump(output_dict, f)

    return 0


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
        'cgradchi0': 1e-8, 'Hgradchi': 1e-8,
        'lamb': 0.0, 'mu': 0.0, 'eta': 0.0, 'tau': 0.0,
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
    output_dict['line 01'] = f"2 {input_dict['cu0']} {-1*multiplier*input_dict['mu']}"
    output_dict['line 02'] = f'2 1e8 1e-8'
    output_dict['line 03'] = f'2 1e8 1e-8'
    output_dict['line 04'] = "2 0. 0."
    output_dict['line 05'] = "2 0. 0."
    output_dict['line 06'] = "2 0. 0."
    output_dict['line 07'] = "2 0. 0."
    output_dict['line 08'] = "2 0. 0."
    output_dict['line 09'] = "2 0. 0."
    # elastic
    output_dict['line 10'] = f"2 {input_dict['lamb']} {input_dict['mu']}"
    output_dict['line 11'] = f"5 {input_dict['eta']} {input_dict['tau']} {input_dict['kappa']} {input_dict['nu']} {input_dict['sigma']}"
    output_dict['line 12'] = f"11 {input_dict['tau1']} {input_dict['tau2']} {input_dict['tau3']} {input_dict['tau4']} {input_dict['tau5']} {input_dict['tau6']} {input_dict['tau7']} {input_dict['tau8']} {input_dict['tau9']} {input_dict['tau10']} {input_dict['tau11']}"
    output_dict['line 13'] = f"2 {input_dict['tau']} {input_dict['sigma']}"
    # integration
    output_dict['line 14'] = '0.5 0.5 0.5 1e-9 1e-9'
    with open(f'{output_file}_plastic_{multiplier}.yml', 'w') as f:
        yaml.dump(output_dict, f)

    return 0


def collect_parameters(parameter_sets, case):
    '''Collect calibration results from one or more yaml files

    :param list parameter_sets: List of yaml files containing calibration results
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter,\
              3: 7 parameter plus tau7, 4: all 18 parameters

    :returns: dictionary containing list of parameter results with each key corresponding to a parameter name
    '''

    results_dict = {
        'element':[],
        'lamb':[], 'mu':[], 'eta':[], 'tau':[], 'kappa':[], 'nu':[], 'sigma':[],
        'tau1':[], 'tau2':[], 'tau3':[], 'tau4':[], 'tau5':[], 'tau6':[], 'tau7':[], 
        'tau8':[], 'tau9':[], 'tau10':[], 'tau11':[]}
    # unpack all 18 values even if they're zero
    for i, set in enumerate(parameter_sets):
        # Load yaml file
        stream = open(set, 'r')
        UI = yaml.load(stream, Loader=yaml.FullLoader)
        stream.close()
        mat_line_1 = UI['line 1']
        mat_line_2 = UI['line 2']
        mat_line_3 = UI['line 3']
        mat_line_4 = UI['line 4']
        # Store results into dictionary
        results_dict['element'].append(i)
        results_dict['lamb'].append(float(mat_line_1.split(' ')[1]))
        results_dict['mu'].append(float(mat_line_1.split(' ')[2]))
        results_dict['eta'].append(float(mat_line_2.split(' ')[1]))
        results_dict['tau'].append(float(mat_line_2.split(' ')[2]))
        results_dict['kappa'].append(float(mat_line_2.split(' ')[3]))
        results_dict['nu'].append(float(mat_line_2.split(' ')[4]))
        results_dict['sigma'].append(float(mat_line_2.split(' ')[5]))
        results_dict['tau1'].append(float(mat_line_3.split(' ')[1]))
        results_dict['tau2'].append(float(mat_line_3.split(' ')[2]))
        results_dict['tau3'].append(float(mat_line_3.split(' ')[3]))
        results_dict['tau4'].append(float(mat_line_3.split(' ')[4]))
        results_dict['tau5'].append(float(mat_line_3.split(' ')[5]))
        results_dict['tau6'].append(float(mat_line_3.split(' ')[6]))
        results_dict['tau7'].append(float(mat_line_3.split(' ')[7]))
        results_dict['tau8'].append(float(mat_line_3.split(' ')[8]))
        results_dict['tau9'].append(float(mat_line_3.split(' ')[9]))
        results_dict['tau10'].append(float(mat_line_3.split(' ')[10]))
        results_dict['tau11'].append(float(mat_line_3.split(' ')[11]))

    # remove zero entries depending on case
    remove = []
    if case == 1:
        remove = ['eta','tau','kappa','nu','sigma','tau1','tau2','tau3',
                  'tau4','tau5','tau6','tau7','tau8','tau9','tau10','tau11']
    elif case == 2:
        remove = ['tau1','tau2','tau3','tau4','tau5','tau6','tau7','tau8',
                  'tau9','tau10','tau11']
    elif case == 3:
        remove = ['tau1','tau2','tau3','tau4','tau5','tau6','tau8','tau9',
                  'tau10','tau11']
    for item in remove:
        results_dict.pop(item)

    return(results_dict)


def make_summary_csv(summary_csv, results_dict):
    '''Make a csv file summarizing the mean, min, max, and standard deviation of calibrated parameters

    :param str summary_csv: Filename to store summary statistics of calibrated parameters
    :param dict results_dict: Results dictionary containing list of parameters with each key corresponding to a parameter name

    :returns: ``summary_csv``
    '''

    params, means, mins, maxs, devs = [], [], [], [], []
    for key in set(results_dict) - {'element'}:
        # get stats
        params.append(key)
        means.append(numpy.mean(results_dict[key]))
        mins.append(numpy.min(results_dict[key]))
        maxs.append(numpy.max(results_dict[key]))
        devs.append(numpy.std(results_dict[key]))
        # output 
        df = pandas.DataFrame({'param': params,
                               'mean': means,
                               'min': mins,
                               'max': maxs,
                               'dev': devs})
    df.to_csv(summary_csv, header=True, sep=',', index=False)

    return 0


def kde(rootname, results_dict, type, kde_best_parameters=None):
    '''Create a kernel density estimate (KDE) plot for each calibrated parameter

    :param str rootname: The rootname of the output plot
    :param dict results_dict: Dictionary containing list of parameter results with each key corresponding to a parameter
    :param str type: A string specifying the type of KDE to plot. 'kde' gives a regular KDE plot. 'hist' gives a KDE plot with histograms shown
    :param str kde_best_parameters: Optional root filename to output a yaml file containing the "best" parameters sampled from the kernel density estimate associated with "kde_best"

    :returns: ``{rootname}_{key}_{type}.PNG`` for each key in `results_dict`, write ``{kde_best_parameters}.yml`` if requested
    '''

    output_parameters = {}
    for key in set(results_dict) - {'element'}:
        matplotlib.pyplot.figure()
        if type != 'best':
            if type == 'kde':
                ax = seaborn.displot(results_dict[key], kind='kde', color='red', fill=True)
            elif type == 'hist':
                ax = seaborn.displot(results_dict[key], kde=True, color='red', fill=True)
            ax.set(xlabel=f'parameter {key}', ylabel='KDE')
            matplotlib.pyplot.title(key)
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig(f'{rootname}_{key}_{type}.PNG')
        else:
            ax = seaborn.kdeplot(results_dict[key], color='red', fill=False)
            xs, ys = ax.lines[-1].get_data()
            ax.fill_between(xs, ys, color='red', alpha=0.1)
            mode_idx = numpy.argmax(ys)
            output_parameters[key] = xs[mode_idx]
            best_value = f'{xs[mode_idx]:.6f}'
            ax.set(xlabel=f'parameter {key}', ylabel='KDE')
            matplotlib.pyplot.title(f'Best {key} = {best_value}')
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig(f'{rootname}_{key}.png')

    # output
    if kde_best_parameters:
        write_elastic_material_card(kde_best_parameters, output_parameters)
        write_plastic_material_card(kde_best_parameters, output_parameters)
    return 0


def summarize_calibration_results(parameter_sets, case,
                                  results_csv=None,
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

    results_dict = collect_parameters(parameter_sets, case)

    if results_csv:
        results_df = pandas.DataFrame(results_dict)
        results_df.to_csv(results_csv, sep=',', index=False)

    if summary_csv:
        make_summary_csv(summary_csv, results_dict)
        
    if kde_hist_plot:
        kde(kde_hist_plot, results_dict, 'hist')

    if kde_plot:
        kde(kde_plot, results_dict, 'kde')

    if kde_best:
        kde(kde_best, results_dict, 'best', kde_best_parameters)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Summarize results of parameter calibration"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('--parameter-sets', nargs="+", required=True,
        help='Specify the list of yaml files containing calibration results')
    parser.add_argument('--case', type=int, required=True,
        help='Specify the calibration "case". 1: two parameter, 2: 7 parameter,\
              3: 7 parameter plus tau7, 4: all 18 parameters')
    parser.add_argument('--results-csv', type=str, required=False,
        help='Optional filename to store all calibrated parameter values')
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
    sys.exit(summarize_calibration_results(parameter_sets=args.parameter_sets,
                                           case=args.case,
                                           results_csv=args.results_csv,
                                           summary_csv=args.summary_csv,
                                           kde_hist_plot=args.kde_hist_plot,
                                           kde_plot=args.kde_plot,
                                           kde_best=args.kde_best,
                                           kde_best_parameters=args.kde_best_parameters,
                                           ))