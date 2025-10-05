#!python
import argparse
import os
import pathlib
import sys

import MOOSE_input_deck_tools as moose_tools


def build_input(output_file, mesh_file, parameter_csv, BCs, disp, duration, damage_parameter=0.095):
    '''Write Tardigrade-MOOSE input file for a gradient-enhanced damage plasticity simulation
    
    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param str mesh_file: The name of the mesh file
    :param list parameter_csv: The csv file containing unique calibrations for each element
    :param str BCs: The type of boundary conditions, either "slip" or "clamp"
    :param float disp: The compressive displacement to be applied
    :param float duration: The duration of the simulation
    :param float damage_parameter: The value of the damage parameter

    :returns: ``output_file``
    '''

    assert os.path.exists(mesh_file), f"Mesh file not found: {mesh_file}"

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
        f.write('[GlobalParams]\n')
        f.write('  displacements = "disp_x disp_y disp_z"\n')
        f.write('  micro_displacement_gradient = "phi_xx phi_xy phi_xz phi_yx phi_yy phi_yz phi_zx phi_zy phi_zz"\n')
        f.write('  nonlocal_damage = "nonlocal_damage"\n')
        f.write('  family = LAGRANGE\n')
        f.write('[]\n')
        f.write('\n')

        # Variables
        f.write('# Variables\n')
        moose_tools.write_variables(f)
        f.write('  [nonlocal_damage]\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # Kernels
        f.write('# Kernels\n')
        moose_tools.write_kernels(f,
                                  internal_force="GradientEnhancedDamagedInternalForce",
                                  internal_couple="GradientEnhancedDamagedInternalCouple")
        f.write('  [GE_Damage]\n')
        f.write('    type = GradientEnhancedMicromorphicDamage\n')
        f.write('    variable = nonlocal_damage\n')
        f.write('    nonlocal_radius = 1.\n')
        f.write('    use_displaced_mesh = false\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')

        # Aux variables
        f.write('# Aux variables\n')
        f.write('[AuxVariables]\n')
        moose_tools.write_default_auxvariables(f)
        f.write('## plastic Aux variables\n')
        moose_tools.write_plastic_auxvariables(f)
        f.write('  [omega]\n')
        f.write('     order = CONSTANT\n')
        f.write('     family = MONOMIAL\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')

        # Aux kernels
        f.write('# Aux kernels\n')
        moose_tools.write_default_auxkernels(f)
        moose_tools.write_plastic_auxkernels(f)
        f.write('[AuxKernels]\n')
        f.write('  [omega]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = ge_damage_statevars\n')
        f.write('    index = 1\n')
        f.write('    variable = omega\n')
        f.write('  []\n')
        f.write('[]\n')
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
        f.write('  [max_omega]\n')
        f.write('    type = ElementExtremeValue\n')
        f.write('    variable = omega\n')
        f.write('  []\n')
        f.write('  [max_nl_dmg]\n')
        f.write('    type = NodalExtremeValue\n')
        f.write('    variable = nonlocal_damage\n')
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

        f.write('[Materials]\n')
        # Write in material info
        f.write('  [./linear_elastic]\n')
        f.write('    type = GradientEnhancedDamagedMicromorphicMaterial\n')
        f.write('    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
        f.write('    material_fparameters = "2 1.0 2.0\n')
        f.write('                            2 4.0 5.0\n')
        f.write('                            2 6.0 7.0\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 28. 29.\n')
        f.write('                            5 31. 32. 33. 34. 35.\n')
        f.write('                            11 37. 38. 39. 40. 41. 42. 43. 44. 45. 46. 47.\n')
        f.write('                            2 32. 35.\n')
        f.write('                            0.5 0.5 0.5 1e-9 1e-9"\n')
        f.write('\n')
        f.write('    user_material_prop_names = "cu0 Hu cchi0 Hchi cnablachi0 Hnablachi lambda mu eta tau kappa nu sigma tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9 tau10 tau11 tauD sigmaD"\n')
        f.write('    user_material_prop_indices = "1 2 4 5 7 8 28 29 31 32 33 34 35 37 38 39 40 41 42 43 44 45 46 47 49 50"\n')
        f.write('\n')
        f.write('    number_SDVS = 55\n')
        #f.write('    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
        f.write('\n')
        f.write(f'    gradient_enhanced_damage_fparameters = "15 {damage_parameter} 1.0 1"\n')
        #f.write('    #Coupled variables\n')
        #f.write('    u1     = "disp_x"\n')
        #f.write('    u2     = "disp_y"\n')
        #f.write('    u3     = "disp_z"\n')
        #f.write('    phi_11 = "phi_xx"\n')
        #f.write('    phi_22 = "phi_yy"\n')
        #f.write('    phi_33 = "phi_zz"\n')
        #f.write('    phi_23 = "phi_yz"\n')
        #f.write('    phi_13 = "phi_xz"\n')
        #f.write('    phi_12 = "phi_xy"\n')
        #f.write('    phi_32 = "phi_zy"\n')
        #f.write('    phi_31 = "phi_zx"\n')
        #f.write('    phi_21 = "phi_yx"\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        # Specify csv file
        f.write('[UserObjects]\n')
        f.write('  [reader_element]\n')
        f.write('    type = PropertyReadFile\n')
        f.write(f'    prop_file_name = "{parameter_csv}"\n')
        f.write('    read_type = "element"\n')
        f.write('    nprop = 24\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')
        # Set up function material map
        f.write('[Materials]\n')
        f.write('  [E_nu]\n')
        f.write('    type = GenericFunctionMaterial\n')
        f.write('    prop_names = "cu0 Hu cchi0 Hchi cnablachi0 Hnablachi lambda mu eta tau kappa nu sigma tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9 tau10 tau11 tauD sigmaD"\n')
        f.write('    prop_values = "func_cu0 func_Hu func_cchi0 func_Hchi func_cnablachi0 func_Hnablachi func_lambda func_mu func_eta func_tau func_kappa func_nu func_sigma func_tau1 func_tau2 func_tau3 func_tau4 func_tau5 func_tau6 func_tau7 func_tau8 func_tau9 func_tau10 func_tau11 func_tau func_sigma"\n')
        f.write('    outputs = exodus\n')
        f.write('    []\n')
        f.write('[]\n')
        f.write('\n')
        # assign parameters
        f.write('[Functions]\n')
        parameters = ['cu0', 'Hu', 'cchi0', 'Hchi', 'cnablachi0', 'Hnablachi',
                      'lambda', 'mu',
                       'eta', 'tau', 'kappa', 'nu', 'sigma',
                       'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7', 'tau8', 'tau9', 'tau10', 'tau11']
        i = 0
        for p in parameters:
            f.write(f'  [func_{p}]\n')
            f.write('    type = PiecewiseConstantFromCSV\n')
            f.write('    read_prop_user_object = "reader_element"\n')
            f.write('    read_type = "element"\n')
            f.write(f'    column_number = "{i}"\n')
            f.write('  []\n')
            i = i + 1
        f.write('[]\n')

        # Preconditioner
        moose_tools.write_preconditioner_block(f)

        # Execution and Timestepping
        dt = duration / 20
        f.write('[Executioner]\n')
        f.write('  type = Transient\n')
        f.write('  solve_type = NEWTON\n')
        f.write('  petsc_options_iname = "-pc_type -pc_factor_mat_solver_package"\n')
        f.write('  petsc_options_value = "lu       superlu_dist                 "\n')
        f.write('  line_search = none\n')
        f.write('  automatic_scaling = true\n')
        f.write('  nl_rel_tol = 1e-8\n')
        f.write('  nl_abs_tol = 1e-8\n')
        f.write('  nl_max_its = 50\n')
        f.write('  start_time = 0.0\n')
        f.write(f'  end_time = {duration}\n')
        f.write('  dtmin = 1e-12\n')
        f.write('  dtmax= 0.1\n')
        f.write('  [TimeStepper]\n')
        f.write('    type = IterationAdaptiveDT\n')
        f.write('    optimal_iterations = 8\n')
        f.write('    iteration_window = 3\n')
        f.write('    linear_iteration_ratio = 1000\n')
        f.write('    growth_factor=1.1\n')
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
    cli_description = "Write Tardigrade-MOOSE input file for a gradient-enhanced damage plasticity simulation"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of Tardigrade-MOOSE file to write")
    parser.add_argument('--mesh', type=str, required=True,
        help='Specify the mesh file')
    parser.add_argument('--parameter-csv', type=str, required=True,
        help='CSV file containing calibration data')
    parser.add_argument('--BCs', type=str, required=True,
        help='Specify the type of boundary conditions, either "slip" or "clamp"')
    parser.add_argument('--disp', type=float, required=True,
        help='Specify the compressive displacement to be applied')
    parser.add_argument('--duration', type=float, required=True,
        help='Specify the duration of the simulation')
    parser.add_argument('--damage-parameter', type=float, required=False, default=0.095,
        help='The value of the damage parameter')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(build_input(output_file=args.output_file,
                         mesh_file=args.mesh,
                         parameter_csv=args.parameter_csv,
                         BCs=args.BCs,
                         disp=args.disp,
                         duration=args.duration,
                         damage_parameter=args.damage_parameter,
                         ))
