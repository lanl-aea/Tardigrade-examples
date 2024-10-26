import os
import sys
import argparse
import yaml
import inspect


def write_plastic_material_card(output_file,
                                lamb, mu, eta, tau, kappa, nu, sigma,
                                tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11,
                                cu0, fraction):
    output_dict = {}
    # plastic
    output_dict['line 01'] = '2 1e4 1e-8'
    output_dict['line 02'] = f'2 {cu0} {-1*fraction*mu}'
    output_dict['line 03'] = '2 1e4 1e-8'
    output_dict['line 04'] = '2 0. 0.'
    output_dict['line 05'] = '2 0. 0.'
    output_dict['line 06'] = '2 0. 0.'
    output_dict['line 07'] = '2 0. 0.'
    output_dict['line 08'] = '2 0. 0.'
    output_dict['line 09'] = '2 0. 0.'
    # elastic
    output_dict['line 10'] = f'2 {lamb} {mu}'
    output_dict['line 11'] = f'5 {eta} {tau} {kappa} {nu} {sigma}'
    output_dict['line 12'] = f'11 {tau1} {tau2} {tau3} {tau4} {tau5} {tau6} {tau7} {tau8} {tau9} {tau10} {tau11}'
    output_dict['line 13'] = f'2 {tau} {sigma}'
    # integration
    output_dict['line 14'] = '0.5 0.5 0.5 1e-9 1e-9'

    with open(output_file, 'w') as f:
        yaml.dump(output_dict, f)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Write elastoplastic Tardigrade-MOOSE input card (.yml)"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of Tardigrade-MOOSE file to write")
    parser.add_argument('--lamb', type=float, required=True,
        help="Specify lambda")
    parser.add_argument('--mu', type=float, required=True,
        help="Specify mu")
    parser.add_argument('--eta', type=float, required=True,
        help="Specify eta")
    parser.add_argument('--tau', type=float, required=True,
        help="Specify tau")
    parser.add_argument('--kappa', type=float, required=True,
        help="Specify kappa")
    parser.add_argument('--nu', type=float, required=True,
        help="Specify nu")
    parser.add_argument('--sigma', type=float, required=True,
        help="Specify sigma")
    parser.add_argument('--tau1', type=float, required=False, default=0.0,
        help="Specify tau1")
    parser.add_argument('--tau2', type=float, required=False, default=0.0,
        help="Specify tau2")
    parser.add_argument('--tau3', type=float, required=False, default=0.0,
        help="Specify tau3")
    parser.add_argument('--tau4', type=float, required=False, default=0.0,
        help="Specify tau4")
    parser.add_argument('--tau5', type=float, required=False, default=0.0,
        help="Specify tau5")
    parser.add_argument('--tau6', type=float, required=False, default=0.0,
        help="Specify tau6")
    parser.add_argument('--tau7', type=float, required=False, default=0.001,
        help="Specify tau7")
    parser.add_argument('--tau8', type=float, required=False, default=0.0,
        help="Specify tau8")
    parser.add_argument('--tau9', type=float, required=False, default=0.0,
        help="Specify tau9")
    parser.add_argument('--tau10', type=float, required=False, default=0.0,
        help="Specify tau10")
    parser.add_argument('--tau11', type=float, required=False, default=0.0,
        help="Specify tau11")
    parser.add_argument('--cu0', type=float, required=True,
        help="Specify c^{u,0}")
    parser.add_argument('--fraction', type=float, required=True,
        help="Specify fraction of mu for sotening")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(write_plastic_material_card(output_file=args.output_file,
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
                                         cu0=args.cu0,
                                         fraction=args.fraction,
                                         ))
