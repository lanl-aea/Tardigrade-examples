import subprocess as sp
import numpy as np
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
import xarray
import numpy

import file_io.xdmf
import summarize_calibration_results


def sort_elements(macro_file):
    '''Read in macroscale XDMF file to identify element not found on the z-boundary

    :param str macro_file: The macroscale filter domain XDMF file, less extension

    :returns: write ``good_elements`` list containing elements not found on the z-boundary
    '''

    # open XDMF file
    macro = file_io.xdmf.XDMF(macro_file)
    macro.open()

    # Get nodes and element connectivity
    nodes = macro.getIncrementReferenceNodePositions(0)[0][0]
    elements = macro.getIncrementConnectivity(0)[0][0]

    # find boundaries from extend of node positions in the z-direcion
    min = numpy.min(nodes,axis=0)[2]
    max = numpy.max(nodes,axis=0)[2]

    # append elements to a list if the do not have nodes on the z-boundary
    good_elements = []
    for i, e in enumerate(elements):
        flag = True
        for n in e:
            x, y, z = nodes[n]
            if (z <= min) or (z >= max):
                flag = False
        if flag:
            good_elements.append(i)

    return good_elements


def collect_parameters_ignore_boundary(parameter_sets, element_sets, good_elements, case):
    '''Collect calibration results from one or more yaml files

    :param list parameter_sets: List of yaml files containing calibration results
    :param list element_sets: List of elements of the macro domain which have been calibrated
    :param list good_elements: List of elements not found on the z-boundary
    :param str parameter_study_file: H5 file with a WAVES parameter study Xarray Dataset
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
    for parameter_set, element in zip(parameter_sets, element_sets):

        if int(element) in good_elements:
            # Load yaml file
            stream = open(parameter_set, 'r')
            UI = yaml.load(stream, Loader=yaml.FullLoader)
            stream.close()
            mat_line_1 = UI['line 1']
            mat_line_2 = UI['line 2']
            mat_line_3 = UI['line 3']
            mat_line_4 = UI['line 4']
            # Store results into dictionary
            results_dict['element'].append(int(element))
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


def summarize_calibration_results_ignore_boundary(parameter_sets,
                                                  element_sets,
                                                  macro_file,
                                                  case,
                                                  results_csv=None,
                                                  summary_csv=None,
                                                  kde_hist_plot=None,
                                                  kde_plot=None,
                                                  kde_best=None,
                                                  kde_best_parameters=None,):
    '''Main function to drive parameter summary and output

    :param list parameter_sets: List of yaml files containing calibration results
    :param list element_sets: List of elements of the macro domain which have been calibrated
    :param str macro_file: The macroscale filter domain XDMF file, less extension
    :param str parameter_study_file: H5 file with a WAVES parameter study Xarray Dataset
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter,\
              3: 7 parameter plus tau7, 4: all 18 parameters
    :param str results_csv: Optional filename to store all calibrated parameter values
    :param str summary_csv: Optional filename to store summary statistics of calibrated parameters
    :param str kde_hist_plot: Optional root filename to plot kernel density estimate of each calibrated parameter with histogram
    :param str kde_plot: Optional root filename to plot kernel density estimate of each calibrated parameter
    :param str kde_best: Optional root filename to plot kernel density estimate of each calibrated parameter with maximum value in title
    :param str kde_best_parameters: Optional root filename to output a yaml file containing the "best" parameters sampled from the kernel density estimate associated with "kde_best"
    '''

    good_elements = sort_elements(macro_file)

    results_dict = collect_parameters_ignore_boundary(parameter_sets, element_sets, good_elements, case)
    if len(results_dict['element']) > 0:

        if results_csv:
            results_df = pandas.DataFrame(results_dict)
            results_df.to_csv(results_csv, sep=',', index=False)

        if summary_csv:
            summarize_calibration_results.make_summary_csv(summary_csv, results_dict)

        if kde_hist_plot:
            summarize_calibration_results.kde(kde_hist_plot, results_dict, 'hist')

        if kde_plot:
            summarize_calibration_results.kde(kde_plot, results_dict, 'kde')

        if kde_best:
            summarize_calibration_results.kde(kde_best, results_dict, 'best', kde_best_parameters)

    else:
        print('All elements touch the boundary! No output generated')

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Summarize results of parameter calibration while ignoring elements on the z-boundary"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('--parameter-sets', nargs="+", required=True,
        help='Specify the list of yaml files containing calibration results')
    parser.add_argument('--element-sets', nargs="+", required=True,
        help='List of elements of the macro domain which have been calibrated')
    parser.add_argument('--macro-file', type=str, required=True,
        help='The macroscale filter domain XDMF file, less extension')
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
    sys.exit(summarize_calibration_results_ignore_boundary(parameter_sets=args.parameter_sets,
                                                           element_sets=args.element_sets,
                                                           macro_file=args.macro_file,
                                                           case=args.case,
                                                           results_csv=args.results_csv,
                                                           summary_csv=args.summary_csv,
                                                           kde_hist_plot=args.kde_hist_plot,
                                                           kde_plot=args.kde_plot,
                                                           kde_best=args.kde_best,
                                                           kde_best_parameters=args.kde_best_parameters,
                                                           ))