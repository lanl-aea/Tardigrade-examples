import os
import sys
import argparse
import yaml
import inspect


def write_plastic_material_card(output_file,
                                lamb, mu, eta, tau, kappa, nu, sigma,
                                tau7, cu0, fraction):
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
    output_dict['line 12'] = f'11 0. 0. 0. 0. 0. 0. {tau7} 0. 0. 0. 0.'
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
    parser.add_argument('--tau7', type=float, required=True,
        help="Specify tau7")
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
                                         tau7=args.tau7,
                                         cu0=args.cu0,
                                         fraction=args.fraction,
                                         ))
