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
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas
    

def main(output_file, num_params, csv_1, csv_2, csv_3):

    if num_params > 7:
        fig_size = (num_params, num_params)
    else:
        fig_size = (10, 10)
    #plt.figure(output_file, figsize=fig_size)
    plt.figure(output_file)
    mpl.rcParams["axes.labelsize"] = 20
    seaborn.set(font_scale=1.7)
    
    headers = [r'$\lambda$', r'$\mu$', r'$\eta$', r'$\tau$', r'$\kappa$', r'$\nu$', r'$\sigma$',
               r'$\tau_1$', r'$\tau_2$', r'$\tau_3$', r'$\tau_4$', r'$\tau_5$', r'$\tau_6$', r'$\tau_7$',
               r'$\tau_8$', r'$\tau_9$', r'$\tau_{10}$', r'$\tau_{11}$']
    # Load in csv files
    df_1 = pandas.read_csv(csv_1, sep=',')
    df_2 = pandas.read_csv(csv_2, sep=',')
    df_3 = pandas.read_csv(csv_3, sep=',')
    
    # calculate averages
    print('mesh 1 summary:')
    print(np.mean(df_1, axis=0))
    print('mesh 2 summary:')
    print(np.mean(df_2, axis=0))
    print('mesh 3 summary:')
    print(np.mean(df_3, axis=0))
    
    # manipulate headers
    original_headers = list(df_1.columns)
    header_map = {}
    for og, new in zip(original_headers, headers):
        header_map[og] = new
    df_1.rename(columns=header_map, inplace=True)
    df_2.rename(columns=header_map, inplace=True)
    df_3.rename(columns=header_map, inplace=True)
    
    # assign mesh id
    df_1 = df_1.assign(mesh=['1' for i in range(0, len(df_1.index))])
    df_2 = df_2.assign(mesh=['2' for i in range(0, len(df_2.index))])
    df_3 = df_3.assign(mesh=['3' for i in range(0, len(df_3.index))])
    
    # drop a few outliers

    #df_3 = df_3.drop([36, 538, 711, 825, 826, 827, 828, 830, 832, 896, 898, 899])
    # The following line is the one used for Ratel_F83_summarize_calibration_ALL <-- figure out a good way to parameterize this! 
    #df_3 = df_3.drop([36, 538, 711, 828, 830, 832, 896, 898, 899])
    
    # concatenate dataframes
    final_df = pandas.concat([df_3, df_2, df_1])
    
    # joint_probabilities
    col_slice = list(range(0,num_params))
    col_slice.append(len(final_df.columns)-1)
    print(final_df.iloc[:,col_slice])
    hue_order = ['1', '2', '3']
    # trim data to 99% 
    # tdf = final_df
    # for column in tdf.columns:
        # tdf = tdf[(tdf.column > tdf.column.quantile(0.005)) & (df.column < df.column.quantile(0.995))]
    #colors = ['#00ce00', '#cfcf00', '#cf00cf']
    colors = ['#1b6032', '#DAA520', '#cf00cf']
    seaborn.set_palette(seaborn.color_palette(colors))
    g = seaborn.PairGrid(final_df.iloc[:,col_slice], hue='mesh', hue_order=hue_order, corner=True)
    #g = seaborn.PairGrid(final_df.iloc[:,col_slice], hue='mesh', palette='tab10', corner=True)
    g.map_diag(seaborn.histplot, stat='percent', common_norm=False, multiple="stack", bins=20)
    #g.map_diag(seaborn.histplot, stat='percent', common_norm=False, log_scale=True)
    #g.map_diag(seaborn.histplot, stat='percent')
    #g.map_diag(seaborn.histplot)
    #g.map_diag(plt.hist(density=True))
    g.map_lower(seaborn.scatterplot)
    g.add_legend(title='Mesh', fontsize='15')
    # for ax in g.axes.flatten():
        # ax.set_ylabel(ax.get_ylabel(), rotation = 90)
    plt.savefig(f'{output_file}.PNG', dpi=300)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Write macroscale definition yaml file"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the basename of output files to write")
    parser.add_argument('--num-params', type=int, required=True,
        help='Specify the number of parameters to make a joint probability plot with')
    parser.add_argument('--csv-1', type=str, required=True,
        help='The first csv file to open')
    parser.add_argument('--csv-2', type=str, required=True,
        help='The second csv file to open')
    parser.add_argument('--csv-3', type=str, required=True,
        help='The third csv file to open')
        
    return parser

if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(main(output_file=args.output_file,
                  num_params=args.num_params,
                  csv_1=args.csv_1,
                  csv_2=args.csv_2,
                  csv_3=args.csv_3,
                  ))