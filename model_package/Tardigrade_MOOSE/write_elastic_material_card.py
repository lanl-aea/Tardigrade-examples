import os
import sys
import argparse
import yaml
import inspect


def write_elastic_material_card(output_file,
                                lamb=0.0, mu=0.0, eta=0.0, tau=0.0,
                                kappa=0.0, nu=0.0, sigma=0.0, tau1=0.0,
                                tau2=0.0, tau3=0.0, tau4=0.0, tau5=0.0,
                                tau6=0.0, tau7=0.0001, tau8=0.0, tau9=0.0,
                                tau10=0.0, tau11=0.0):
    '''Write elastic Tardigrade-MOOSE input card (.yml)

    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param float lamb: The lambda parameter
    :param float mu: The mu parameter
    :param float eta: The eta parameter
    :param float tau: The tau parameter
    :param float kappa: The kappa parameter
    :param float nu: The nu parameter
    :param float sigma: The sigma parameter
    :param float tau1: The tau1 parameter
    :param float tau2: The tau2 parameter
    :param float tau3: The tau3 parameter
    :param float tau4: The tau4 parameter
    :param float tau5: The tau5 parameter
    :param float tau6: The tau6 parameter
    :param float tau7: The tau7 parameter
    :param float tau8: The tau8 parameter
    :param float tau9: The tau9 parameter
    :param float tau10: The tau10 parameter
    :param float tau11: The tau11 parameter

    :returns: Write ``output_file``
    '''

    output_dict = {}
    # elastic
    output_dict['line 1'] = f'2 {lamb} {mu}'
    output_dict['line 2'] = f'5 {eta} {tau} {kappa} {nu} {sigma}'
    output_dict['line 3'] = f'11 {tau1} {tau2} {tau3} {tau4} {tau5} {tau6} {tau7} {tau8} {tau9} {tau10} {tau11}'
    output_dict['line 4'] = f'2 {tau} {sigma}'

    with open(output_file, 'w') as f:
        yaml.dump(output_dict, f)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Write elastic Tardigrade-MOOSE input card (.yml)"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of Tardigrade-MOOSE file to write")
    parser.add_argument('--lamb', type=float, required=False, default=0.0,
        help="Specify lambda")
    parser.add_argument('--mu', type=float, required=False, default=0.0,
        help="Specify mu")
    parser.add_argument('--eta', type=float, required=False, default=0.0,
        help="Specify eta")
    parser.add_argument('--tau', type=float, required=False, default=0.0,
        help="Specify tau")
    parser.add_argument('--kappa', required=False, default=0.0,
        help="Specify kappa")
    parser.add_argument('--nu', required=False, default=0.0,
        help="Specify nu")
    parser.add_argument('--sigma', required=False, default=0.0,
        help="Specify sigma")
    parser.add_argument('--tau1', required=False, default=0.0,
        help="Specify tau1")
    parser.add_argument('--tau2', required=False, default=0.0,
        help="Specify tau2")
    parser.add_argument('--tau3', required=False, default=0.0,
        help="Specify tau3")
    parser.add_argument('--tau4', required=False, default=0.0,
        help="Specify tau4")
    parser.add_argument('--tau5', required=False, default=0.0,
        help="Specify tau5")
    parser.add_argument('--tau6', required=False, default=0.0,
        help="Specify tau6")
    parser.add_argument('--tau7', required=False, default=0.001,
        help="Specify tau7")
    parser.add_argument('--tau8', required=False, default=0.0,
        help="Specify tau8")
    parser.add_argument('--tau9', required=False, default=0.0,
        help="Specify tau9")
    parser.add_argument('--tau10', required=False, default=0.0,
        help="Specify tau10")
    parser.add_argument('--tau11', required=False, default=0.0,
        help="Specify tau11")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(write_elastic_material_card(output_file=args.output_file,
                                         lamb=args.lamb,
                                         mu=args.mu,
                                         eta=args.eta,
                                         tau=args.tau,
                                         kappa=args.kappa,
                                         nu=args.nu,
                                         sigma=args.sigma,
                                         tau1=args.tau1,
                                         tau2=args.tau2,
                                         tau3=args.tau3,
                                         tau4=args.tau4,
                                         tau5=args.tau5,
                                         tau6=args.tau6,
                                         tau7=args.tau7,
                                         tau8=args.tau8,
                                         tau9=args.tau9,
                                         tau10=args.tau10,
                                         tau11=args.tau11,
                                         ))