import os
import sys
import argparse
import time
import glob
import yaml
import inspect


def build_options_file(output_file, material_E, material_nu, material_rho, top_id, bottom_id, num_steps, displacement, BCs):
    '''Write Ratel options file

    :param str output_file: The name of the Ratel options file to output
    :param str material_E: The material's elastic modulus
    :param str material_nu: The material's Poisson ratio
    :param str material_rho: The material's density
    :param int top_id: The id of the top surface
    :param int bottom_id: The id of the bottom surface
    :param int num_steps: The number of steps for the simulation
    :param float displacement: The displacement to apply to the top surface
    :param str BCs: The type of boundary conditions, either 'slip' or 'clamp'

    :returns: Write ``output_file``
    '''

    with open(output_file, 'w') as f:
        f.write('#Ratel input file\n')
        f.write('order: 2\n')
        f.write('\n')
        f.write('material: elastic\n')
        f.write('\n')
        f.write('elastic:\n')
        f.write('  model: elasticity-linear\n')
        f.write(f'  E: {material_E}\n')
        f.write(f'  nu: {material_nu}\n')
        f.write(f'  rho: {material_rho}\n')
        f.write('\n')
        f.write('ts:\n')
        f.write('  max_time: 1.0\n')
        f.write(f'  dt: {1.0/num_steps}\n')
        f.write('\n')
        f.write('dm:\n')
        f.write('  view:\n')
        f.write('\n')
        f.write(f'surface_force_faces: {top_id}, {bottom_id}\n')
        f.write('\n')
        f.write('bc:\n')
        if BCs == 'slip':
            f.write(f'  slip: {top_id}, {bottom_id}\n')
            f.write(f'  slip_{top_id}_components: 2\n')
            f.write(f'  slip_{bottom_id}_components: 2\n')
            f.write(f'  slip_{top_id}_translate: {displacement}\n')
        elif BCs == 'clamp':
            f.write(f'  clamp: {top_id}, {bottom_id}\n')
            f.write(f'  clamp_{top_id}_translate: 0, 0, {displacement}\n')
        else:
            print('Specify a valid BC type!')
        f.write('\n')
        f.write('ksp:\n')
        f.write('  type: gmres\n')
        f.write('  norm_type: unpreconditioned\n')
        f.write('\n')

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Write Ratel options file"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="The name of the Ratel options file to output")
    parser.add_argument('--material-E', type=float, required=True,
        help="The material's elastic modulus")
    parser.add_argument('--material-nu', type=float, required=True,
        help="The material's Poisson ratio")
    parser.add_argument('--material-rho', type=float, required=True,
        help="The material's density")
    parser.add_argument('--top-id', type=int, required=True,
        help="The id of the top surface")
    parser.add_argument('--bottom-id', type=int, required=True,
        help="The id of the bottom surface")
    parser.add_argument('--num-steps', type=int, required=True,
        help="The number of steps for the simulation")
    parser.add_argument('--displacement', type=float, required=True,
        help="The displacement to apply to the top surface")
    parser.add_argument('--BCs', type=str, required=True,
        help="The type of boundary conditions, either 'slip' or 'clamp'")

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(build_options_file(output_file=args.output_file,
                                material_E=args.material_E,
                                material_nu=args.material_nu,
                                material_rho=args.material_rho,
                                top_id=args.top_id,
                                bottom_id=args.bottom_id,
                                num_steps=args.num_steps,
                                displacement=args.displacement,
                                BCs=args.BCs,
                                ))