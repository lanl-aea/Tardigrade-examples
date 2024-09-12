import sys
import os
import argparse
import inspect
import yaml

import numpy
import matplotlib.pyplot
import scipy
import pandas
from itertools import compress

import micromorphic


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Calibrate micromorphic linear elasticity for averaged output on a single filter domain (i.e. macroscale element)"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
        help="The homogenized XDMF file output by the Micromorphic Filter")
    parser.add_argument('-o', '--output-file', type=str,
        help="The resulting list of parameters stored in a yaml file")
    parser.add_argument('--Emod', type=float,
        help="DNS elastic modulus, used for initial parameter estimation.")
    parser.add_argument('--nu', type=float,
        help="DNS Poisson's ratio, used for initial parameter estimation.")
    parser.add_argument('--L', type=float,
        help="DNS max dimension (width, height, depth, etc.), used for initial parameter estimation.")
    parser.add_argument('--element', type=int, default=0,
        help="The macro (filter) element to calibrate")
    parser.add_argument('--increment', nargs="+", required=False, default=None,
        help="An optional list of one or more increments to perform calibration")
    parser.add_argument('--case', type=int, required=True,
        help="Specify the calibration 'case'. 1: two parameter, 2: 7 parameter,\
              3: 7 parameter plus tau7, 4: all 18 parameters")
    parser.add_argument('--plot-file', type=str, required=False, default=None,
        help="Optional root filename to plot Cauchy and symmetric micro stress\
              comparison between DNS and calibration results")
    parser.add_argument('--average', type=str, required=False, default=True,
        help='Boolean whether or not homogenized DNS results will be averaged')
    parser.add_argument('--UQ-file', type=str, required=False,
        help='Optional csv filename to store function evaluations and parameter sets for UQ')

    return parser


if __name__ == '__main__':
    parser = get_parser()
   
    args, unknown = parser.parse_known_args()
    sys.exit(calibrate(
                    input_file=args.input_file,
                    output_file=args.output_file,
                    Emod=args.Emod,
                    nu=args.nu,
                    L=args.L,
                    element=args.element,
                    increment=args.increment,
                    case=args.case,
                    plot_file=args.plot_file,
                    average=str2bool(args.average),
                    UQ_file=args.UQ_file,
                    ))