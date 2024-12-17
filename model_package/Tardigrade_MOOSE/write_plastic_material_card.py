#!python
import argparse
import os
import pathlib
import sys

import pandas
import yaml


def write_plastic_material_card(output_file, output_type,
                                lamb, mu, eta, tau, kappa, nu, sigma,
                                tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11,
                                cu0, Hu, cchi0, Hchi, cnablachi0, Hnablachi):
    '''Write elastoplastic Tardigrade-MOOSE input card

    :params str output_file: The name of Tardigrade-MOOSE file to write
    :params str output_type: The type of material card to write, either 'yaml' or 'csv'
    :params float lambda: The elastic lambda parameter
    :params float mu: The elastic mu parameter
    :params float eta: The elastic mu parameter
    :params float tau: The elastic tau parameter
    :params float kappa: The elastic kappa parameter
    :params float nu: The elastic nu parameter
    :params float sigma: The elastic sigma parameter
    :params float tau1: The elastic tau1 parameter
    :params float tau2: The elastic tau2 parameter
    :params float tau3: The elastic tau3 parameter
    :params float tau4: The elastic tau4 parameter
    :params float tau5: The elastic tau5 parameter
    :params float tau6: The elastic tau6 parameter
    :params float tau7: The elastic tau7 parameter
    :params float tau8: The elastic tau8 parameter
    :params float tau9: The elastic tau9 parameter
    :params float tau10: The elastic tau10 parameter
    :params float tau11: The elastic tau11 parameter
    :params float cu0: The plastic initial macro cohesion parameter
    :params float Hu: The plastic macro hardening parameter
    :params float cchi0: The plastic initial micro cohesion parameter
    :params float Hchi: The plastic micro hardening parameter
    :params float cnablachi0: The plastic initial micro gradient cohesion parameter
    :params float Hnablachi: The plastic micro gradient hardening parameter

    :returns: ``output_file``
    '''

    # Dump to yaml
    if output_type == 'yaml':
        output_dict = {}
        # plastic
        output_dict['line 01'] = f'2 {cu0} {Hu}'
        output_dict['line 02'] = f'2 {cchi0} {Hchi}'
        output_dict['line 03'] = f'2 {cnablachi0} {Hnablachi}'
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
    # Dump to csv
    elif output_type == 'csv':
        elastic_parameter_ordering = ['lambda', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',
                                      'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',
                                      'tau8', 'tau9', 'tau10', 'tau11']
        plastic_parameter_ordering = ['cu0', 'Hu', 'cchi0', 'Hchi', 'cnablachi0', 'Hnablachi']
        parameter_ordering = ['element'] + plastic_parameter_ordering + elastic_parameter_ordering
        parameter_list = [0, cu0, Hu, cchi0, Hchi, cnablachi0, Hnablachi,
                          lamb, mu, eta, tau, kappa, nu, sigma,
                          tau1, tau2, tau3, tau4, tau5, tau6, tau7,
                          tau8, tau9, tau10, tau11]
        df = pandas.DataFrame([parameter_list], columns=parameter_ordering)
        df.to_csv(output_file, header=True, sep=',', index=False)
    else:
        print('Specify valid output_type!')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "

    cli_description = "Write elastoplastic Tardigrade-MOOSE input card"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="The name of Tardigrade-MOOSE file to write")
    parser.add_argument('--output-type', type=str, required=False, default='yaml',
        help="The type of material card to write, either 'yaml' or 'csv'")
    parser.add_argument('--lamb', type=float, required=True,
        help="The elastic lambda parameter")
    parser.add_argument('--mu', type=float, required=True,
        help="The elastic mu parameter")
    parser.add_argument('--eta', type=float, required=False, default=0.0,
        help="The elastic mu parameter")
    parser.add_argument('--tau', type=float, required=False, default=0.0,
        help="The elastic tau parameter")
    parser.add_argument('--kappa', type=float, required=False, default=0.0,
        help="The elastic kappa parameter")
    parser.add_argument('--nu', type=float, required=False, default=0.0,
        help="The elastic nu parameter")
    parser.add_argument('--sigma', type=float, required=False, default=0.0,
        help="The elastic sigma parameter")
    parser.add_argument('--tau1', type=float, required=False, default=0.0,
        help="The elastic tau1 parameter")
    parser.add_argument('--tau2', type=float, required=False, default=0.0,
        help="The elastic tau2 parameter")
    parser.add_argument('--tau3', type=float, required=False, default=0.0,
        help="The elastic tau3 parameter")
    parser.add_argument('--tau4', type=float, required=False, default=0.0,
        help="The elastic tau4 parameter")
    parser.add_argument('--tau5', type=float, required=False, default=0.0,
        help="The elastic tau5 parameter")
    parser.add_argument('--tau6', type=float, required=False, default=0.0,
        help="The elastic tau6 parameter")
    parser.add_argument('--tau7', type=float, required=False, default=0.001,
        help="The elastic tau7 parameter")
    parser.add_argument('--tau8', type=float, required=False, default=0.0,
        help="The elastic tau8 parameter")
    parser.add_argument('--tau9', type=float, required=False, default=0.0,
        help="The elastic tau9 parameter")
    parser.add_argument('--tau10', type=float, required=False, default=0.0,
        help="The elastic tau10 parameter")
    parser.add_argument('--tau11', type=float, required=False, default=0.0,
        help="The elastic tau11 parameter")
    parser.add_argument('--cu0', type=float, required=False, default=1.e4,
        help="The plastic initial macro cohesion parameter, $c^{u,0}$")
    parser.add_argument('--Hu', type=float, required=False, default=1.e-8,
        help="The plastic macro hardening parameter, $H^u$")
    parser.add_argument('--cchi0', type=float, required=False, default=1.e4,
        help="The plastic initial micro cohesion parameter, $c^{\chi,0}$")
    parser.add_argument('--Hchi', type=float, required=False, default=1.e-8,
        help="The plastic micro hardening parameter, $H^{\chi}$")
    parser.add_argument('--cnablachi0', type=float, required=False, default=1.e4,
        help="The plastic initial micro gradient cohesion parameter, $c^{\nabla\chi,0}$")
    parser.add_argument('--Hnablachi', type=float, required=False, default=1.e-8,
        help="The plastic micro gradient hardening parameter, $H^{\nabla\chi}$")

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(write_plastic_material_card(output_file=args.output_file,
                                         output_type=args.output_type,
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
                                         Hu=args.Hu,
                                         cchi0=args.cchi0,
                                         Hchi=args.Hchi,
                                         cnablachi0=args.cnablachi0,
                                         Hnablachi=args.Hnablachi,
                                         ))
