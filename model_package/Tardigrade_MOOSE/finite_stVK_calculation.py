import sys
import os
import argparse
import pathlib

import numpy


def finite_stVK_calculation(diameter=5.0, material_E=250., material_nu=0.2, eps_z=-0.01):
    '''Solution for uniaxial stress of a cylinder for finite deformation using the St. Venant-Kirchhoff elasticity model

    :params float diameter: The diameter of the cylinder in millimeters
    :params float material_E: The elastic modulus in MPa
    :params float material_nu: The Poisson ratio
    :params float eps_z: The applied nominal strain in the z-direction

    :returns: Print various solution quantities to the terminal
    '''

    # Convert to Lame' parameters
    lamb = (material_E*material_nu)/((1. + material_nu)*(1. - 2.*material_nu))
    mu   = material_E/(2*(1. + material_nu))
    print(f'lambda* = {lamb}')
    print(f'mu* = {mu}\n')

    # Applied stretch
    stretch_z = 1. + eps_z

    # Solve for lateral stretch
    term1 = ((3 - stretch_z**2)*0.5*lamb) + mu
    term2 = lamb + mu
    stretch_x = numpy.sqrt(term1/term2)
    print(f'stretch_x = {stretch_x}\n') 
    stretch_y = stretch_x
    nominal_lateral_strain = stretch_x - 1.
    lateral_displacement = diameter*nominal_lateral_strain/2
    print(f'lateral displacement x = {lateral_displacement}\n')

    # Construct deformation gradient
    F = numpy.array([[stretch_x, 0.0, 0.0],
                     [0.0, stretch_y, 0.0],
                     [0.0, 0.0, stretch_z]])
    print(f'F = {F}\n')

    # Misc tensor calculations
    Ft = numpy.transpose(F)
    Finv = numpy.linalg.inv(F)
    J = numpy.linalg.det(F)

    # Construct Green-Lagrange Strain
    E = 0.5*(numpy.dot(Ft,F) - numpy.eye(3))
    Ekk = numpy.trace(E)
    print(f'E = {E}')
    print(f'Ekk = {Ekk}\n')

    # PK2
    S = lamb*Ekk*numpy.eye(3) + 2*mu*E
    print(f'PK2 = {S}\n')

    # push forward for Cauchy
    cauchy = (1./J)*numpy.einsum("kK,KL,Ll->kl", F, S, F)
    print(f'cauchy = {cauchy}\n')

    # Calculate Forces
    initial_area = 0.25*numpy.pi*diameter**2
    print(f'initial_area area = {initial_area}')
    initial_force = initial_area*cauchy[2,2]
    print(f'Force using initial_area area = {initial_force} N\n')
    current_area = 0.25*numpy.pi*(diameter*stretch_x)**2
    print(f'current area = {current_area}')
    current_force = current_area*cauchy[2,2]
    print(f'Force using original area = {current_force} N\n')

    return current_force, lateral_displacement


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Solution for uniaxial stress of a cylinder for finite deformation using the St. Venant-Kirchhoff elasticity model"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--diameter', type=float, required=False, default=5.0,
        help='The diameter of the cylinder in millimeters')
    parser.add_argument('--material-E', type=float, required=False, default=250.,
        help='The elastic modulus in MPa')
    parser.add_argument('--material-nu', type=float, required=False, default=0.2,
        help='The Poisson ratio')
    parser.add_argument('--eps-z', type=float, required=False, default=-0.01,
        help='The applied nominal strain in the z-direction')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(finite_stVK_calculation(diameter=args.diameter,
                                     material_E=args.material_E,
                                     material_nu=args.material_nu,
                                     eps_z=args.eps_z,
                                     ))