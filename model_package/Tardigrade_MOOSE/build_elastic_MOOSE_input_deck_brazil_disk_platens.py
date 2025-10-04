#!python
import argparse
import os
import pathlib
import sys

import MOOSE_input_deck_tools as moose_tools


def build_input(output_file, mesh_file, material_E, material_nu, platen_E, platen_nu,
                disp, duration,
                specimen_top_surface, specimen_bottom_surface, top_platen_contact, bottom_platen_contact,
                top_platen_fixture, top_platen_side, top_platen_back,
                bottom_platen_fixture, bottom_platen_side, bottom_platen_back,
                contact_type='friction', friction_coefficient=None):
    '''Write MOOSE input file for Brazilian disk simulation with platens
    
    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param str mesh_file: The name of the mesh file
    :param float material_E: The elastic modulus of the specimen material
    :param float material_nu: The Poisson ratio of the specimen material
    :param float platen_E: The elastic modulus of the platen material
    :param float platen_nu: The Poisson ratio of the platen material
    :param float disp: The compressive displacement to be applied
    :param float duration: The duration of the simulation
    :param str specimen_top_surface: Specify the name of the specimen top contact surface
    :param str specimen_bottom_surface: Specify the name of the specimen bottom contact surface
    :param str top_platen_contact: Specify the name of the top platen contact surface
    :param str bottom_platen_contact: Specify the name of the bottom platen contact surface
    :param str top_platen_fixture: Specify the name of the top platen fixture surface
    :param str top_platen_side: Specify the name of the top platen side surface
    :param str top_platen_back: Specify the name of the top platen back surface
    :param str bottom_platen_fixture: Specify the name of the bottom platen fixture surface
    :param str bottom_platen_side: Specify the name of the bottom platen side surface
    :param str bottom_platen_back: Specify the name of the bottom platen back surface
    :param str contact_type: The option for specifying contact, either "frictionless" or "friction"
    :param float friction_coefficient: The friction coefficient if contact_type="friction"

    :returns: ``output_file``
    '''

    assert os.path.exists(mesh_file), f"Mesh file not found: {mesh_file}"

    # Write input file
    with open(output_file, 'w') as f:
        f.write('###############################################################################\n')
        f.write('[Mesh]\n')
        f.write('  type = FileMesh\n')
        f.write(f'  file = "{mesh_file}"\n')
        f.write('  patch_update_strategy = iteration\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[GlobalParams]\n')
        f.write('  displacements = "disp_x disp_y disp_z"\n')
        f.write('[]\n')

        # Variables
        f.write('[Variables]\n')
        f.write('  [./disp_x]\n')
        f.write('  [../]\n')
        f.write('  [./disp_y]\n')
        f.write('  [../]\n')
        f.write('  [./disp_z]\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # Kernels
        f.write('[Kernels]\n')
        f.write('  #Define the internal force balance equations\n')
        f.write('  [./rxn]\n')
        f.write('    type = Reaction\n')
        f.write('    variable  = disp_y\n')
        f.write('    save_in = force_y\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # Aux variables
        f.write('[AuxVariables]\n')
        f.write('  [./force_x]\n')
        f.write('  [../]\n')
        f.write('  [./force_y]\n')
        f.write('  [../]\n')
        f.write('  [./force_z]\n')
        f.write('  [../]\n')
        f.write('  [./saved_x]\n')
        f.write('  [../]\n')
        f.write('  [./saved_y]\n')
        f.write('  [../]\n')
        f.write('  [./saved_z]\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # Aux kernels
        f.write('[AuxKernels]\n')
        f.write('[]\n')
        f.write('\n')

        # Postprocessors
        f.write('[Postprocessors]\n')
        f.write('  [./bot_react_y]\n')
        f.write('    type = NodalSum\n')
        f.write('    variable = force_y\n')
        f.write(f'    boundary = "{bottom_platen_fixture}"\n')
        f.write('  [../]\n')
        f.write('  [./stress_yy]\n')
        f.write('    type = ElementAverageValue\n')
        f.write('    variable = stress_yy\n')
        f.write('    block = "specimen"\n')
        f.write('  [../]\n')
        f.write('  [./reaction_top]\n')
        f.write('    type = SidesetReaction\n')
        f.write('    direction = "0 1 0"\n')
        f.write('    stress_tensor = stress\n')
        f.write(f'    boundary = "{top_platen_fixture}"\n')
        f.write('  [../]\n')
        f.write('  [./reaction_bottom]\n')
        f.write('    type = SidesetReaction\n')
        f.write('    direction = "0 1 0"\n')
        f.write('    stress_tensor = stress\n')
        f.write(f'    boundary = "{bottom_platen_fixture}"\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # Finite deformation
        f.write('[Physics/SolidMechanics/QuasiStatic]\n')
        f.write('  [./all]\n')
        f.write('    add_variables = true\n')
        f.write('    strain = FINITE\n')
        f.write('    block = "specimen bottom_platen top_platen"\n')
        f.write('    use_automatic_differentiation = false\n')
        f.write('    generate_output = "stress_xx stress_xy stress_xz stress_yy stress_zz"\n')
        f.write('    save_in = "saved_x saved_y saved_z"\n')
        f.write('    use_finite_deform_jacobian = true\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # BCs
        f.write('[BCs]\n')
        f.write('  active = "bottom_x bottom_y bottom_z top_x top_y top_z spec_z"\n')
        f.write('  [./bottom_x]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_x\n')
        f.write(f'    boundary = "{bottom_platen_side}"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('  [./bottom_y]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_y\n')
        f.write(f'    boundary = "{bottom_platen_fixture}"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('  [./bottom_z]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_z\n')
        f.write(f'    boundary = "{bottom_platen_back}"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('  [./top_x]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_x\n')
        f.write(f'    boundary = "{top_platen_side}"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('  [./top_y]\n')
        f.write('    type = FunctionDirichletBC\n')
        f.write('    variable = disp_y\n')
        f.write(f'    boundary = "{top_platen_fixture}"\n')
        f.write('    #preset = true\n')
        f.write('    function = top_bc\n')
        f.write('  [../]\n')
        f.write('  [./top_z]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_z\n')
        f.write(f'    boundary = "{top_platen_back}"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('  # DELETE THIS ONE! JUST FOR TESTING!!!!\n')
        f.write('  [./spec_z]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_z\n')
        f.write(f'    boundary = "specimen_back"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # Loading function
        f.write('[Functions]\n')
        f.write('  [./top_bc]\n')
        f.write('    type  = ParsedFunction\n')
        f.write(f'    expression = -{disp}*t\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # Materials
        f.write('[Materials]\n')
        f.write('  [./elasticity_tensor_specimen]\n')
        f.write('    type = ComputeIsotropicElasticityTensor\n')
        f.write(f'    youngs_modulus = {material_E}\n')
        f.write(f'    poissons_ratio = {material_nu}\n')
        f.write('    block = "specimen"\n')
        f.write('  [../]\n')
        f.write('  [./linear_elastic]\n')
        f.write('    type = ComputeFiniteStrainElasticStress\n')
        f.write('    block = "specimen"\n')
        f.write('  [../]\n')
        f.write('  [stress]\n')
        f.write('    type=ComputeFiniteStrainElasticStress\n')
        f.write(f'    boundary = "{bottom_platen_fixture}"\n')
        f.write('  []\n')
        f.write('  [./elasticity_tensor_platen]\n')
        f.write('    type = ComputeIsotropicElasticityTensor\n')
        f.write(f'    youngs_modulus = {platen_E}\n')
        f.write(f'    poissons_ratio = {platen_nu}\n')
        f.write('    block = "bottom_platen top_platen"\n')
        f.write('  [../]\n')
        f.write('  [./platen]\n')
        f.write('    type = ComputeFiniteStrainElasticStress\n')
        f.write('    block = "bottom_platen top_platen"\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # Preconditioner
        moose_tools.write_preconditioner_block(f)

        # Execution and Timestepping
        dt = duration / 100
        f.write('[Executioner]\n')
        f.write('  type = Transient\n')
        f.write('  solve_type = PJFNK\n')
        f.write('  petsc_options_iname = "-pc_type -pc_factor_mat_solver_type"\n')
        f.write('  petsc_options_value = "lu    superlu_dist"\n')
        f.write('  line_search = "none"\n')
        f.write('  automatic_scaling = true\n')
        f.write('  nl_rel_tol = 1e-8\n')
        f.write('  nl_abs_tol = 5e-5\n')
        f.write('  l_tol = 1e-3\n')
        f.write('  l_abs_tol = 1e-5\n')
        f.write('  l_max_its = 60\n')
        f.write('  nl_max_its = 50\n')
        f.write('  start_time = 0.0\n')
        f.write(f'  end_time = {duration}\n')
        f.write('  dtmin = 1e-6\n')
        f.write('  dtmax= 0.1\n')
        f.write('\n')
        f.write('  [TimeStepper]\n')
        f.write('    type = IterationAdaptiveDT\n')
        f.write('    growth_factor=1.2\n')
        f.write(f'    dt = {dt}\n')
        f.write('  []\n')
        f.write('[]\n')

        # Outputs
        moose_tools.write_outputs_block(f)
        f.write('\n')

        # Contact
        if contact_type == 'frictionless':
            f.write('[Contact]\n')
            f.write('  [./bottom_center_cont]\n')
            f.write(f'    primary = "{bottom_platen_contact}"\n')
            f.write(f'    secondary = "{specimen_bottom_surface}"\n')
            f.write('    penalty = 1e4\n')
            f.write('    normalize_penalty = true\n')
            f.write('    tangential_tolerance = 1.e-3\n')
            f.write('  [../]\n')
            f.write('  [./top_center_cont]\n')
            f.write(f'    primary = "{top_platen_contact}"\n')
            f.write(f'    secondary = "{specimen_top_surface}"\n')
            f.write('    penalty = 1e4\n')
            f.write('    normalize_penalty = true\n')
            f.write('    tangential_tolerance = 1.e-3\n')
            f.write('  [../]\n')
            f.write('[]\n')
        elif contact_type == 'friction':
            assert (friction_coefficient != None), f"Friction coefficient must be defined!"
            f.write('[Contact]\n')
            f.write('  [./bottom_center_cont]\n')
            #f.write(f'    primary = "{bottom_platen_contact}"\n')
            #f.write(f'    secondary = "{specimen_bottom_surface}"\n')
            f.write(f'    secondary = "{bottom_platen_contact}"\n')
            f.write(f'    primary = "{specimen_bottom_surface}"\n')
            f.write('    model = coulomb\n')
            f.write('    formulation = tangential_penalty\n')
            f.write(f'    friction_coefficient = {friction_coefficient}\n')
            f.write('    penalty = 1e4\n')
            f.write('    normalize_penalty = true\n')
            f.write('    normal_smoothing_distance = 0.1\n')
            f.write('  [../]\n')
            f.write('  [./top_center_cont]\n')
            f.write(f'    primary = "{top_platen_contact}"\n')
            f.write(f'    secondary = "{specimen_top_surface}"\n')
            f.write('    model = coulomb\n')
            f.write('    formulation = tangential_penalty\n')
            f.write(f'    friction_coefficient = {friction_coefficient}\n')
            f.write('    penalty = 1e4\n')
            f.write('    normalize_penalty = true\n')
            f.write('    normal_smoothing_distance = 0.1\n')
            f.write('  [../]\n')
            f.write('[]\n')
            f.write('\n')
            f.write('[Dampers]\n')
            f.write('  [./contact_slip_bottom]\n')
            f.write('    type = ContactSlipDamper\n')
            #f.write(f'    primary = "{bottom_platen_contact}"\n')
            #f.write(f'    secondary = "{specimen_bottom_surface}"\n')
            f.write(f'    secondary = "{bottom_platen_contact}"\n')
            f.write(f'    primary = "{specimen_bottom_surface}"\n')
            f.write('  [../]\n')
            f.write('  [./contact_slip_top]\n')
            f.write('    type = ContactSlipDamper\n')
            f.write(f'    primary = "{top_platen_contact}"\n')
            f.write(f'    secondary = "{specimen_top_surface}"\n')
            f.write('  [../]\n')
            f.write('[]\n')
            f.write('\n')
        else:
            print('Specify a valid contact_type!')
        f.write('\n')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Write MOOSE input file for Brazilian disk simulation with platens"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of Tardigrade-MOOSE file to write")
    parser.add_argument('--mesh', type=str, required=True,
        help='Specify the mesh file')
    parser.add_argument('--material-E', type=float, required=True,
        help='The elastic modulus of the specimen material')
    parser.add_argument('--material-nu', type=float, required=True,
        help='The Poisson ratio of the specimen material')
    parser.add_argument('--platen-E', type=float, required=True,
        help='The elastic modulus of the platen material')
    parser.add_argument('--platen-nu', type=float, required=True,
        help='The Poisson ratio of the platen material')
    parser.add_argument('--disp', type=float, required=True,
        help='Specify the compressive displacement to be applied')
    parser.add_argument('--duration', type=float, required=True,
        help='Specify the duration of the simulation')
    parser.add_argument('--specimen-top-surface', type=str, required=True,
        help='Specify the name of the specimen top contact surface')
    parser.add_argument('--specimen-bottom-surface', type=str, required=True,
        help='Specify the name of the specimen bottom contact surface')
    parser.add_argument('--top-platen-contact', type=str, required=True,
        help='Specify the name of the top platen contact surface')
    parser.add_argument('--bottom-platen-contact', type=str, required=True,
        help='Specify the name of the bottom platen contact surface')
    parser.add_argument('--top-platen-fixture', type=str, required=True,
        help='Specify the name of the top platen fixture surface')
    parser.add_argument('--top-platen-side', type=str, required=True,
        help='Specify the name of the top platen side surface')
    parser.add_argument('--top-platen-back', type=str, required=True,
        help='Specify the name of the top platen back surface')
    parser.add_argument('--bottom-platen-fixture', type=str, required=True,
        help='Specify the name of the bottom platen fixture surface')
    parser.add_argument('--bottom-platen-side', type=str, required=True,
        help='Specify the name of the bottom platen side surface')
    parser.add_argument('--bottom-platen-back', type=str, required=True,
        help='Specify the name of the bottom platen back surface')
    parser.add_argument('--contact-type', type=str, required=False, default='friction',
        help='The option for specifying contact, either "frictionless" or "friction"')
    parser.add_argument('--friction-coefficient', type=float, required=False, default=None,
        help='The fricition coefficient if contact_type="friction"')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(build_input(output_file=args.output_file,
                         mesh_file=args.mesh,
                         material_E=args.material_E,
                         material_nu=args.material_nu,
                         platen_E=args.platen_E,
                         platen_nu=args.platen_nu,
                         disp=args.disp,
                         duration=args.duration,
                         specimen_top_surface=args.specimen_top_surface,
                         specimen_bottom_surface=args.specimen_bottom_surface,
                         top_platen_contact=args.top_platen_contact,
                         bottom_platen_contact=args.bottom_platen_contact,
                         top_platen_fixture=args.top_platen_fixture,
                         top_platen_side=args.top_platen_side,
                         top_platen_back=args.top_platen_back,
                         bottom_platen_fixture=args.bottom_platen_fixture,
                         bottom_platen_side=args.bottom_platen_side,
                         bottom_platen_back=args.bottom_platen_back,
                         contact_type=args.contact_type,
                         friction_coefficient=args.friction_coefficient,
                         ))
