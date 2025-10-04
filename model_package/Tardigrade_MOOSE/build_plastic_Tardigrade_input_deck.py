#!python
import argparse
import os
import pathlib
import sys

import pandas

import MOOSE_input_deck_tools as moose_tools


def unpack_plastic_parameter_csv(parameter_df, i):
    '''Convert a single line of a plastic calibration map into relevant material strings and element number

    :params DataFrame parameter_df: The loaded calibration map
    :params int i: The current DataFrame index

    :returns: mat_line_1, mat_line_2, mat_line_3, mat_line_blank, mat_line_10, mat_line_11, mat_line_12, mat_line_14, element
    '''

    # Main platicity lines
    cu0, Hu = parameter_df.at[i, 'cu0'], parameter_df.at[i, 'Hu']
    mat_line_1 = f'2 {cu0} {Hu}'
    cchi0, Hchi = parameter_df.at[i, 'cchi0'], parameter_df.at[i, 'Hchi']
    mat_line_2 = f'2 {cchi0} {Hchi}'
    cnablachi0, Hnablachi = parameter_df.at[i, 'cnablachi0'], parameter_df.at[i, 'Hnablachi']
    mat_line_3 = f'2 {cnablachi0} {Hnablachi}'

    # force lines 4 to 9 to be the same for now
    # TODO: update once calibration includes all 18 plasticity parameters
    mat_line_blank = '2 0. 0.'

    # A tensor parameters
    lamb, mu = parameter_df.at[i, 'lambda'], parameter_df.at[i, 'mu']
    mat_line_10 = f'2 {lamb} {mu}'

    # B and D tensor parameters
    eta, tau, kappa, = parameter_df.at[i, 'eta'], parameter_df.at[i, 'tau'], parameter_df.at[i, 'kappa']
    nu, sigma = parameter_df.at[i, 'nu'], parameter_df.at[i, 'sigma']
    mat_line_11 = f'5 {eta} {tau} {kappa} {nu} {sigma}'
    mat_line_13 = f'2 {tau} {sigma}'
    mat_line_14 = '0.5 0.5 0.5 1e-9 1e-9'

    # C tensor parameters
    mat_line_12 = '11'
    for j in range(1, 12):
        tau_j = f'tau{j}'
        mat_line_12 = mat_line_12 + f' {parameter_df.at[i, tau_j]}'

    element = parameter_df.at[i, 'element']

    return mat_line_1, mat_line_2, mat_line_3, mat_line_blank, mat_line_10, mat_line_11, mat_line_12, mat_line_13, mat_line_14, element


def build_input(output_file, mesh_file, calibration_map, BCs, disp, duration):
    '''Write a Tardigrade-MOOSE input file
    
    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param str mesh_file: The name of the mesh file
    :param str calibration_map: CSV file containing calibration data
    :param str BCs: The type of boundary conditions, either "slip", "slip_plane", "clamp", or "brazil"
    :param float disp: The compressive displacement to be applied
    :param float duration: The duration of the simulation

    :returns: ``output_file``
    '''

    assert os.path.exists(mesh_file), f"Mesh file not found: {mesh_file}"

    # load calibration map
    parameter_df = pandas.read_csv(calibration_map)
    parameter_df = parameter_df.sort_values(by='element')

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
        f.write('## plastic Aux variables\n')
        moose_tools.write_plastic_auxvariables(f)
        f.write('[]\n')
        f.write('\n')

        # Aux kernels
        f.write('# Aux kernels\n')
        moose_tools.write_default_auxkernels(f)
        moose_tools.write_plastic_auxkernels(f)
        f.write('\n')

        # BCs
        if BCs == 'brazil':
            sample_force = 'force_y'
            sample_boundary = 'brazil_load'
        else:
            sample_force = 'force_z'
            sample_boundary = 'top'
        f.write('# Do some cool math to get the reaction force\n')
        f.write('[Postprocessors]\n')
        f.write('  [bot_react_z]\n')
        f.write('    type = NodalSum\n')
        f.write(f'    variable = {sample_force}\n')
        f.write(f'    boundary = "{sample_boundary}"\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')
        if BCs == 'slip':
            f.write('[BCs]\n')
            f.write('  active = "x_symm y_symm bottom_z top_z"\n')
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
            f.write('[]\n')
        elif BCs == 'slip_plane':
            f.write('[BCs]\n')
            f.write('  active = "bottom_z top_z"\n')
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
            f.write('[]\n')
        elif BCs == 'clamp':
            f.write('[BCs]\n')
            f.write('  active = "bottom_x bottom_y bottom_z top_x top_y top_z"\n')
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
            f.write('[]\n')
        elif BCs == 'brazil':
            f.write('[BCs]\n')
            f.write('  active = "bottom_z brazil_fix_x brazil_fix_y brazil_load_x brazil_load_y"\n')
            f.write('  [./bottom_z]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./brazil_fix_x]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "brazil_fix"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./brazil_fix_y]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "brazil_fix"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./brazil_load_x]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "brazil_load"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./brazil_load_y]\n')
            f.write('    type = FunctionDirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "brazil_load"\n')
            f.write('    preset = true\n')
            f.write('    function = top_bc\n')
            f.write('  [../]\n')
            f.write('[]\n')
        else:
            print('Specify a valid BC type!')
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
        if len(list(parameter_df.index)) > 1:
            for index in parameter_df.index:
                # Unpack parameters
                mat_line_1, mat_line_2, mat_line_3, mat_line_blank, mat_line_10, mat_line_11, mat_line_12, mat_line_13, mat_line_14, element = unpack_plastic_parameter_csv(parameter_df, index)
                # Write in material info
                f.write(f'  [./linear_elastic_{element}]\n')
                f.write('    type = MicromorphicMaterial\n')
                f.write(f'    material_fparameters = "{mat_line_1}\n')
                f.write(f'                            {mat_line_2}\n')
                f.write(f'                            {mat_line_3}\n')
                f.write(f'                            {mat_line_blank}\n')
                f.write(f'                            {mat_line_blank}\n')
                f.write(f'                            {mat_line_blank}\n')
                f.write(f'                            {mat_line_blank}\n')
                f.write(f'                            {mat_line_blank}\n')
                f.write(f'                            {mat_line_blank}\n')
                f.write(f'                            {mat_line_10}\n')
                f.write(f'                            {mat_line_11}\n')
                f.write(f'                            {mat_line_12}\n')
                f.write(f'                            {mat_line_13}\n')
                f.write(f'                            {mat_line_14}"\n')
                f.write('    number_SDVS = 55\n')
                f.write(f'    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
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
            mat_line_1, mat_line_2, mat_line_3, mat_line_blank, mat_line_10, mat_line_11, mat_line_12, mat_line_13, mat_line_14, element = unpack_plastic_parameter_csv(parameter_df, 0)
            # Write in material info
            f.write(f'  [./linear_elastic]\n')
            f.write('    type = MicromorphicMaterial\n')
            f.write(f'    material_fparameters = "{mat_line_1}\n')
            f.write(f'                            {mat_line_2}\n')
            f.write(f'                            {mat_line_3}\n')
            f.write(f'                            {mat_line_blank}\n')
            f.write(f'                            {mat_line_blank}\n')
            f.write(f'                            {mat_line_blank}\n')
            f.write(f'                            {mat_line_blank}\n')
            f.write(f'                            {mat_line_blank}\n')
            f.write(f'                            {mat_line_blank}\n')
            f.write(f'                            {mat_line_10}\n')
            f.write(f'                            {mat_line_11}\n')
            f.write(f'                            {mat_line_12}\n')
            f.write(f'                            {mat_line_13}\n')
            f.write(f'                            {mat_line_14}"\n')
            f.write('    number_SDVS = 55\n')
            f.write(f'    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
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
        f.write('[]\n')
        f.write('\n')

        # Preconditioner
        moose_tools.write_preconditioner_block(f)

        # Execution and Timestepping
        dt = duration / 20
        f.write('[Executioner]\n')
        f.write('  type = Transient\n')
        f.write('  solve_type = PJFNK\n')
        f.write('  petsc_options_iname = "-pc_type -pc_factor_mat_solver_package"\n')
        f.write('  petsc_options_value = "lu       superlu_dist                 "\n')
        f.write('  line_search = none\n')
        f.write('  automatic_scaling = true\n')
        f.write('  nl_rel_tol = 1e-5\n')
        f.write('  nl_abs_tol = 1e-7\n')
        f.write('  nl_max_its = 50\n')
        f.write('  start_time = 0.0\n')
        f.write('  end_time = 1.0\n')
        f.write('  dtmin = 1e-12\n')
        f.write('  dtmax= 0.1\n')
        f.write('  [TimeStepper]\n')
        f.write('    type = IterationAdaptiveDT\n')
        f.write('    optimal_iterations = 4\n')
        f.write('    iteration_window = 3\n')
        f.write('    linear_iteration_ratio = 1000\n')
        f.write('    growth_factor=1.2\n')
        f.write('    cutback_factor=0.5\n')
        f.write(f'    dt = {dt}\n')
        f.write('  []\n')
        f.write('[]\n')

        # Outputs
        moose_tools.write_outputs_block(f)
        f.write('\n')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Write Tardigrade-MOOSE input file for a plastic simulation"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of Tardigrade-MOOSE file to write")
    parser.add_argument('--mesh', type=str, required=True,
        help='Specify the mesh file')
    parser.add_argument('--calibration-map', type=str, required=True,
        help='CSV file containing calibration data')
    parser.add_argument('--BCs', type=str, required=True,
        help='Specify the type of boundary conditions, either "slip", "slip_plane", "clamp", or "brazil"')
    parser.add_argument('--disp', type=float, required=True,
        help='Specify the compressive displacement to be applied')
    parser.add_argument('--duration', type=float, required=True,
        help='Specify the duration of the simulation')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(build_input(output_file=args.output_file,
                         mesh_file=args.mesh,
                         calibration_map=args.calibration_map,
                         BCs=args.BCs,
                         disp=args.disp,
                         duration=args.duration,
                         ))
