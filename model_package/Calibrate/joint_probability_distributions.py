import os
import sys
import argparse
import time
import glob
import yaml
import inspect

import seaborn
import matplotlib
import matplotlib.pyplot
import pandas
import numpy
import scipy


def plot_distributions(dfs, columns, num_domains, distribution_plots):
    ''' Plot normal distributions across all filtering domain studies for each parameter

    :param list dfs: List of Pandas DataFrames containing calibration results for each filtering domain study
    :param dataframe columns: Pandas DataFrame Index containing parameter names to summarize
    :param list num_domains: A list of integers corresponding to the number of filtering domains associated with results contained in each csv file'
    :param str distribution_plots: Root file name for distribution plots of each parameter

    :returns: Write ``{distribution_plots}_{parameter}.png`` plot file for each calibraiton parameter
    '''

    for column in columns:
        parameter = str(column[2:-1])
        matplotlib.pyplot.figure()
        for df, num_domain in zip(dfs, num_domains):
            mean = numpy.mean(df[column])
            std = numpy.mean(df[column])
            x_axis = numpy.linspace(numpy.min(df[column]), numpy.max(df[column]), 100)
            matplotlib.pyplot.plot(x_axis, scipy.stats.norm.pdf(x_axis, mean, std), label=f'{num_domain} domains')
        matplotlib.pyplot.xlabel(column)
        matplotlib.pyplot.legend()
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig(f'{distribution_plots}_{parameter}.png')

    return 0


def full_kde_plot(df, columns, full_kde):
    '''Plot combined KDE across all filtering domain studies for each parameter

    :param dataframe df: Pandas DataFrame containing results from all filtering domains
    :param dataframe columns: Pandas DataFrame Index containing parameter names to summarize
    :param str full_kde: Root file name for KDE of each parameter

    :returns: Write ``{full_kde}_{parameter}.png`` plot file for each calibration parameter
    '''

    df = df.drop('domains', axis=1)
    for column in columns:
        matplotlib.pyplot.figure()
        data = numpy.array(df[column])
        parameter = str(column[2:-1])
        ax = seaborn.kdeplot(data, color='red', fill=False)
        xs, ys = ax.lines[-1].get_data()
        ax.fill_between(xs, ys, color='red', alpha=0.1)
        mode_idx = numpy.argmax(ys)
        #ax.vlines(xs[mode_idx], 0, ys[mode_idx,], ls="--", color='red')
        #ax.text(xs[mode_idx], -0.1, f'{xs[mode_idx]:.6f}', color='k', ha='center', transform=ax.get_xaxis_transform())
        best_value = f'{xs[mode_idx]:.6f}'
        ax.set(xlabel=column)
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.title(f'Best {column} = {best_value}')
        matplotlib.pyplot.savefig(f'{full_kde}_{parameter}.png')

    return 0


def joint_probability_distributions(output_file, csv_files, num_domains, num_params=None, distribution_plots=None, full_kde=None):
    '''Create a joint probability distribution plot to summarize calibration results

    :param str output_file: The output filename for the joint probability distribution plot"
    :param list csv_files: The csv files containing calibration results
    :param list num_domains: A list of integers corresponding to the number of filtering domains associated with results contained in each csv file'
    :param int num_params: The number of parameters to make a joint probability plot with if not all are desired'
    :param str distribution_plots: Optional root file name for distribution plots of each parameter
    :param str full_kde: Optional root file name for KDE of each parameter

    :returns: ``output_file``
    '''

    output_plot = output_file.split('.')[0]
    matplotlib.pyplot.figure(output_plot)
    matplotlib.rcParams["axes.labelsize"] = 20
    seaborn.set(font_scale=1.7)
    
    headers = {'lamb': r'$\lambda$',
               'mu': r'$\mu$',
               'eta': r'$\eta$',
               'tau': r'$\tau$',
               'kappa': r'$\kappa$',
               'nu': r'$\nu$',
               'sigma': r'$\sigma$',
               'tau1': r'$\tau_1$',
               'tau2': r'$\tau_2$',
               'tau3': r'$\tau_3$',
               'tau4': r'$\tau_4$',
               'tau5': r'$\tau_5$',
               'tau6': r'$\tau_6$',
               'tau7': r'$\tau_7$',
               'tau8': r'$\tau_8$',
               'tau9': r'$\tau_9$',
               'tau10': r'$\tau_{10}$',
               'tau11': r'$\tau_{11}$'}

    dfs = []
    for csv_file, num_domain in zip(csv_files, num_domains):
        df = pandas.read_csv(csv_file, sep=',', index_col=0)
        # Reassign headers
        original_headers = list(df.columns)
        total_params = len(original_headers)
        header_map = {}
        for head in original_headers:
            header_map[head] = headers[head]
        df.rename(columns=header_map, inplace=True)

        # assign domain number
        df = df.assign(domains=[num_domain for i in range(0, len(df.index))])
        dfs.append(df)

    # plot normal distributions
    columns = dfs[0].columns[0:-1]
    if distribution_plots:
        plot_distributions(dfs, columns, num_domains, distribution_plots)

    # Concatenate DataFrames
    final_df = pandas.concat(dfs[::-1])
    print(final_df)

    # Full kde
    if full_kde:
        full_kde_plot(final_df, columns, full_kde)

    # Hand the number of parameters to plot
    if num_params == None:
        num_params = total_params
    col_slice = list(range(0, num_params))
    col_slice.append(len(final_df.columns)-1)

    # Set hue order based on number of domains
    hue_order = [str(num_domain) for num_domain in num_domains[::-1]]

    # Create pairplot
    g = seaborn.PairGrid(final_df.iloc[:,col_slice], hue='domains', hue_order=hue_order, corner=True)
    g.map_diag(seaborn.histplot, stat='percent', common_norm=False, multiple="stack", bins=20)
    g.map_lower(seaborn.scatterplot)
    g.add_legend(title='Domains', fontsize='15')
    # Format labels
    for ax in g.axes.flatten():
        if ax:
            ax.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
            if len(ax.get_ylabel()) > 0:
                ax.set_ylabel(ax.get_ylabel(), rotation='horizontal', labelpad=20)
                ax.yaxis.get_label().set_horizontalalignment('right')
            if len(ax.get_xlabel()) > 0:
                ax.set_xlabel(ax.get_xlabel(), labelpad=20)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_file, dpi=300)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Create a joint probability distribution plot to summarize calibration results"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the output filename for the joint probability distribution plot")
    parser.add_argument('--num-params', type=int, required=False,
        help='Optionally specify the number of parameters to make a joint probability plot with if not all are desired')
    parser.add_argument('--csv-files', nargs="+", required=True,
        help='The csv files containing calibration results')
    parser.add_argument('--num-domains', nargs="+", required=True,
        help="A list of integers corresponding to the number of filtering domains\
              associated with results contained in each csv file")
    parser.add_argument('--distribution-plots', type=str, required=False,
        help="Optional root file name for distribution plots of each parameter")
    parser.add_argument('--full-kde', type=str, required=False,
        help="Optional root file name for KDE of each parameter")

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(joint_probability_distributions(output_file=args.output_file,
                                             csv_files=args.csv_files,
                                             num_domains=args.num_domains,
                                             num_params=args.num_params,
                                             distribution_plots=args.distribution_plots,
                                             full_kde=args.full_kde,
                                             ))