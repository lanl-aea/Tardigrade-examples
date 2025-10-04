#!python
import argparse
import os
import pathlib
import sys
import yaml

import pandas

import MOOSE_input_deck_tools as moose_tools


def unpack_elastic_parameter_csv(parameter_df, i):
    '''Convert a single line of an elastic calibration map into relevant material strings and element number

    :params DataFrame parameter_df: The loaded calibration map
    :params int i: The current DataFrame index

    :returns: mat_line_1, mat_line_2, mat_line_3, mat_line_4, element
    '''


    # A tensor parameters
    lamb, mu = parameter_df.at[i, 'lambda'], parameter_df.at[i, 'mu']
    mat_line_1 = f'2 {lamb} {mu}'

    # B and D tensor parameters
    eta, tau, kappa, = parameter_df.at[i, 'eta'], parameter_df.at[i, 'tau'], parameter_df.at[i, 'kappa']
    nu, sigma = parameter_df.at[i, 'nu'], parameter_df.at[i, 'sigma']
    mat_line_2 = f'5 {eta} {tau} {kappa} {nu} {sigma}'
    mat_line_4 = f'2 {tau} {sigma}'

    # C tensor parameters
    mat_line_3 = '11'
    for j in range(1, 12):
        tau_j = f'tau{j}'
        mat_line_3 = mat_line_3 + f' {parameter_df.at[i, tau_j]}'

    element = parameter_df.at[i, 'element']

    return mat_line_1, mat_line_2, mat_line_3, mat_line_4, element


def build_input(output_file, mesh_file, BCs, disp, duration, disp_point=None, calibration_map=None,
                elastic_material_card=None, phi_BC=None):
    '''Write a Tardigrade-MOOSE input file

    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param str mesh_file: The name of the mesh file
    :param str BCs: The type of boundary conditions, either "slip" or "clamp"
    :param float disp: The compressive displacement to be applied
    :param float duration: The duration of the simulation
    :param str disp_point: Optional string of coordinates to query x-displacement
    :param str calibration_map: CSV file containing calibration data, first method for specifying material parameters
    :param str elastic_material_card: YML file containing elastic material parameters, second method for specifying material parameters
    :param str phi_BC: Optional string specifying nodeset to force micro deformation components to be zero

    :returns: ``output_file``
    '''

    assert os.path.exists(mesh_file), f"Mesh file not found: {mesh_file}"

    # load calibration map
    if calibration_map is not None:
        parameter_df = pandas.read_csv(calibration_map)
        parameter_df = parameter_df.sort_values(by='element')
    else:
        assert elastic_material_card is not None, "Either calibration_map of elastic_material_card must be provided!"

    # Write input file
    with open(output_file, 'w') as f:
        f.write('###############################################################################\n')
        f.write('[Mesh]\n')
        f.write('  type = FileMesh\n')
        f.write('  displacements = "disp_x disp_y disp_z"\n')
        f.write('  dim = 3\n')
        f.write(f'  file = "{mesh_file}"\n')
        f.write('[]\n')
        f.write('\n')

        # Variables
        f.write('# Variables\n')
        moose_tools.write_variables(f)
        f.write('[]\n')
        f.write('\n')

        # Kernels
        f.write('# Kernels\n')
        moose_tools.write_kernels(f)
        f.write('[]\n')
        f.write('\n')

        # Aux variables
        f.write('# Aux variables\n')
        f.write('[AuxVariables]\n')
        moose_tools.write_default_auxvariables(f)
        f.write('[]\n')
        f.write('\n')

        # Aux kernels
        f.write('# Aux kernels\n')
        moose_tools.write_default_auxkernels(f)
        f.write('[]\n')
        f.write('\n')

        # Postprocessor(s)
        f.write('# Get the reaction force\n')
        f.write('[Postprocessors]\n')
        f.write('  [bot_react_z]\n')
        f.write('    type = NodalSum\n')
        f.write('    variable = force_z\n')
        f.write('    boundary = "top"\n')
        f.write('  [../]\n')
        # Custom displacement query for disp_x
        if disp_point:
            f.write('  [lateral_disp]\n')
            f.write('    type = PointValue\n')
            f.write(f'    point = "{disp_point}"\n')
            f.write('    variable = disp_x\n')
            f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # BCs
        if BCs == 'slip':
            f.write('[BCs]\n')
            f.write('  [./x_symm]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "x_plane"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./y_symm]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "y_plane"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./bottom_z]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./top_z]\n')
            f.write('    type = FunctionDirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "top"\n')
            f.write('    preset = true\n')
            f.write('    function = top_bc\n')
            f.write('  [../]\n')
        elif BCs == 'clamp':
            f.write('[BCs]\n')
            f.write('  [./bottom_x]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./bottom_y]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./bottom_z]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./top_x]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "top"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./top_y]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "top"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./top_z]\n')
            f.write('    type = FunctionDirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "top"\n')
            f.write('    preset = true\n')
            f.write('    function = top_bc\n')
            f.write('  [../]\n')
        else:
            print('Specify a valid BC type!')
        # Option to force Phis to be zero
        if phi_BC is not None:
            moose_tools.write_phi_BCs(f, phi_BC)
        f.write('[]\n')
        f.write('\n')
        f.write('[Functions]\n')
        f.write('  [./top_bc]\n')
        f.write('    type  = ParsedFunction\n')
        f.write(f'    expression = -{disp}*t\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')


        # Materials
        f.write('# Materials\n')
        f.write('[Materials]\n')
        # Load in parameter data for each filter domain / element
        if calibration_map is not None:
            if len(list(parameter_df.index)) > 1:
                for index in parameter_df.index:
                    # Unpack parameters
                    mat_line_1, mat_line_2, mat_line_3, mat_line_4, element = unpack_elastic_parameter_csv(parameter_df, index)
                    # Write in material info
                    f.write(f'  [./linear_elastic_{element}]\n')
                    f.write('    type = MicromorphicMaterial\n')
                    f.write(f'    material_fparameters = "{mat_line_1}\n')
                    f.write(f'                            {mat_line_2}\n')
                    f.write(f'                            {mat_line_3}\n')
                    f.write(f'                            {mat_line_4}"\n')
                    f.write(f'    model_name = "LinearElasticity"\n')
                    f.write('\n')
                    f.write('    #Coupled variables\n')
                    f.write('    u1     = "disp_x"\n')
                    f.write('    u2     = "disp_y"\n')
                    f.write('    u3     = "disp_z"\n')
                    f.write('    phi_11 = "phi_xx"\n')
                    f.write('    phi_22 = "phi_yy"\n')
                    f.write('    phi_33 = "phi_zz"\n')
                    f.write('    phi_23 = "phi_yz"\n')
                    f.write('    phi_13 = "phi_xz"\n')
                    f.write('    phi_12 = "phi_xy"\n')
                    f.write('    phi_32 = "phi_zy"\n')
                    f.write('    phi_31 = "phi_zx"\n')
                    f.write('    phi_21 = "phi_yx"\n')
                    f.write(f'    block = "element_{element}"\n')
                    f.write('  [../]\n')
            else:
                # Unpack parameters
                mat_line_1, mat_line_2, mat_line_3, mat_line_4, element = unpack_elastic_parameter_csv(parameter_df, 0)
                # Write in material info
                f.write(f'  [./linear_elastic]\n')
                f.write('    type = MicromorphicMaterial  \n')
                f.write(f'    material_fparameters = "{mat_line_1}\n')
                f.write(f'                            {mat_line_2}\n')
                f.write(f'                            {mat_line_3}\n')
                f.write(f'                            {mat_line_4}"\n')
                f.write(f'    model_name = "LinearElasticity"\n')
                f.write('\n')
                f.write('    #Coupled variables\n')
                f.write('    u1     = "disp_x"\n')
                f.write('    u2     = "disp_y"\n')
                f.write('    u3     = "disp_z"\n')
                f.write('    phi_11 = "phi_xx"\n')
                f.write('    phi_22 = "phi_yy"\n')
                f.write('    phi_33 = "phi_zz"\n')
                f.write('    phi_23 = "phi_yz"\n')
                f.write('    phi_13 = "phi_xz"\n')
                f.write('    phi_12 = "phi_xy"\n')
                f.write('    phi_32 = "phi_zy"\n')
                f.write('    phi_31 = "phi_zx"\n')
                f.write('    phi_21 = "phi_yx"\n')
                f.write('  [../]\n')
        elif elastic_material_card is not None: 
            # Load yaml file
            stream = open(elastic_material_card, 'r')
            UI = yaml.load(stream, Loader=yaml.FullLoader)
            stream.close()
            moose_tools.write_elastic_material_card(f, UI)
        else:
            print('Uh oh! No valid material type has been provided!')
        f.write('[]\n')

        # Preconditioner
        moose_tools.write_preconditioner_block(f)

        # Execution and Timestepping
        f.write('[Executioner]\n')
        f.write('#  type = Steady\n')
        f.write('  type = Transient\n')
        dt = duration / 10
        f.write('  num_steps = 10\n')
        f.write(f'  dt        = {dt}\n')
        f.write('  solve_type = "PJFNK"\n')
        f.write('#  solve_type = "NEWTON"\n')
        f.write('  nl_rel_tol = 1e-8\n')
        f.write('  nl_abs_tol = 1e-8\n')
        f.write('  nl_max_its = 100\n')
        f.write('[]\n')
        f.write('\n')

        # Outputs
        moose_tools.write_outputs_block(f)
        f.write('\n')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Write Tardigrade-MOOSE input file"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="The name of Tardigrade-MOOSE file to write")
    parser.add_argument('--mesh', type=str, required=True,
        help='The mesh file')

    parser.add_argument('--BCs', type=str, required=True,
        help='The type of boundary conditions, either "slip" or "clamp"')
    parser.add_argument('--disp', type=float, required=True,
        help='The compressive displacement to be applied')
    parser.add_argument('--duration', type=float, required=True,
        help='The duration of the simulation')
    parser.add_argument('--calibration-map', type=str, required=False, default=None,
        help='CSV file containing calibration data, first method for specifying material parameters')
    parser.add_argument('--elastic-material-card', type=str, required=False, default=None,
        help='YML file containing elastic material parameters, second method for specifying material parameters')
    parser.add_argument('--disp-point', type=str, required=False, default=None,
        help='Optional string of coordinates to query x-displacement')
    parser.add_argument('--phi-BC', type=str, required=False, default=None,
        help='Optional string specifying nodeset to force micro deformation components to be zero')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(build_input(output_file=args.output_file,
                         mesh_file=args.mesh,
                         BCs=args.BCs,
                         disp=args.disp,
                         duration=args.duration,
                         calibration_map=args.calibration_map,
                         elastic_material_card=args.elastic_material_card,
                         disp_point=args.disp_point,
                         phi_BC=args.phi_BC,
                         ))