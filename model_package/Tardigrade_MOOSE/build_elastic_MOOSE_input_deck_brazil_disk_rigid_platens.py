#!python
import argparse
import os
import pathlib
import sys

import MOOSE_input_deck_tools as moose_tools


def build_input(output_file, mesh_file, material_E, material_nu, platen_radius,
                disp, duration, specimen_bottom_surface, specimen_top_surface=None, 
                top_symmetry=None, back_symmetry=None, side_symmetry=None,
                xc_bot=0., yc_bot=0., xc_top=0., yc_top=0., geometry='full'):
    '''Write MOOSE input file for Brazilian disk simulation with rigid platens
    
    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param str mesh_file: The name of the mesh file
    :param float material_E: The elastic modulus of the specimen material
    :param float material_nu: The Poisson ratio of the specimen material
    :param float platen_radius: The radius of curvature of the Brazilian disk compression platen
    :param float disp: The compressive displacement to be applied
    :param float duration: The duration of the simulation
    :param str specimen_bottom_surface: The name of the specimen bottom contact surface
    :param str specimen_top_surface: The name of the specimen top contact surface. Required if "geometry" = "full."
    :param str top_symmetry: The name of the top symmetry surface. Required if "geometry" = "quarter" or "eighth."
    :param str back_symmetry: The name of the back symmetry surface. Required if "geometry" = "quarter" or "eighth" or "half"
    :param str side_set: The name of the side symmetry surface. Required if "geometry" = "quarter" or "eighth."
    :param float xc_bot: The x-position of the center of the circular bottom surface arc
    :param float yc_bot: The y-position of the center of the circular bottom surface arc
    :param float xc_top: The x-position of the center of the circular top surface arc
    :param float yc_top: The y-position of the center of the circular top surface arc
    :param str geometry: The geometry/symmetry type: "full," "half," "quarter," or "eighth"

    :returns: ``output_file``
    '''

    assert os.path.exists(mesh_file), f"Mesh file not found: {mesh_file}"

    if geometry == 'full':
        active_BCs = 'bottom_y bottom_x top_y top_x'
        react_surface = specimen_top_surface
        stress_boundary = f'{specimen_bottom_surface} {specimen_top_surface}'
        assert specimen_top_surface != None, "Specimen top surface must be defined if 'geometry' = 'full'!"
    else:
        assert back_symmetry != None, "Specimen back symmetry must be defined if 'geometry' = 'quarter' or 'eighth' or 'half'!"
        stress_boundary = specimen_bottom_surface
        if geometry == 'quarter':
            assert top_symmetry != None, "Specimen top symmetry must be defined if 'geometry' = 'quarter' or 'eighth'!"
            react_surface = top_symmetry
            active_BCs = 'bottom_y bottom_x top_sym back_sym'
        elif geometry == 'eighth':
            assert top_symmetry != None, "Specimen top symmetry must be defined if 'geometry' = 'quarter' or 'eighth'!"
            react_surface = top_symmetry
            active_BCs = 'bottom_y bottom_x top_sym back_sym side_sym'
            assert side_symmetry != None, "Specimen side symmetry must be defined if 'geometry' = 'eighth'!"
        elif geometry == 'half':
            active_BCs = 'bottom_y bottom_x back_sym top_y top_x'
            react_surface = specimen_top_surface
            assert specimen_top_surface != None, "Specimen top surface must be defined if 'geometry' = 'full'!"
        else:
            print('Specify a valid geometry type!')

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
        f.write('  [./react_y_bottom]\n')
        f.write('    type = NodalSum\n')
        f.write('    variable = force_y\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write('  [../]\n')
        f.write('  [./avg_stress_xx]\n')
        f.write('    type = ElementAverageValue\n')
        f.write('    variable = stress_xx\n')
        f.write('    block = "specimen"\n')
        f.write('  [../]\n')
        f.write('  [./reaction_bottom]\n')
        f.write('    type = SidesetReaction\n')
        f.write('    direction = "0 1 0"\n')
        f.write('    stress_tensor = stress\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write('  [../]\n')
        f.write('  [./react_y_top]\n')
        f.write('    type = NodalSum\n')
        f.write('    variable = force_y\n')
        f.write(f'    boundary = "{react_surface}"\n')
        f.write('  [../]\n')
        f.write('  [./reaction_top]\n')
        f.write('    type = SidesetReaction\n')
        f.write('    direction = "0 1 0"\n')
        f.write('    stress_tensor = stress\n')
        f.write(f'    boundary = "{react_surface}"\n')
        f.write('  [../]\n')
        f.write('[]\n')

        # Finite deformation
        f.write('[Physics/SolidMechanics/QuasiStatic]\n')
        f.write('  [./all]\n')
        f.write('    add_variables = true\n')
        f.write('    strain = FINITE\n')
        f.write('    block = "specimen"\n')
        f.write('    use_automatic_differentiation = false\n')
        f.write('    generate_output = "stress_xx stress_xy stress_xz stress_yy stress_zz"\n')
        f.write('    save_in = "saved_x saved_y saved_z"\n')
        f.write('    use_finite_deform_jacobian = true\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

        # BCs
        f.write('[BCs]\n')
        f.write(f'  active = "{active_BCs}"\n')
        f.write('  [./bottom_y]\n')
        f.write('    type = CylindricalSurfaceDirichletBC\n')
        f.write('    variable = disp_y\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write(f'    center = "{xc_bot} {yc_bot} 0"\n')
        f.write(f'    radius = {platen_radius}\n')
        f.write(f'    velocity = {disp}\n')
        f.write('    axis = "0. 0. 1."\n')
        f.write('    normal = "0. 1. 0."\n')
        f.write('    use_sector = true\n')
        f.write('    angle_min = 3.14159265359\n')
        f.write('    angle_max = 6.28318530718\n')
        f.write('    invert_displacement = true\n')
        f.write('  [../]\n')
        f.write('  [./bottom_x]\n')
        f.write('    type = CylindricalSurfaceDirichletBC\n')
        f.write('    variable = disp_x\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write(f'    center = "{xc_bot} {yc_bot} 0"\n')
        f.write(f'    radius = {platen_radius}\n')
        f.write(f'    velocity = {disp}\n')
        f.write('    axis = "0. 0. 1."\n')
        f.write('    normal = "0. 1. 0."\n')
        f.write('    use_sector = true\n')
        f.write('    angle_min = 3.14159265359\n')
        f.write('    angle_max = 6.28318530718\n')
        f.write('    invert_displacement = true\n')
        f.write('  [../]\n')
        if (geometry == 'full') or (geometry == 'half'):
            f.write('  [./top_y]\n')
            f.write('    type = CylindricalSurfaceDirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write(f'    boundary = "{specimen_top_surface}"\n')
            f.write(f'    center = "{xc_top} {yc_top} 0"\n')
            f.write(f'    radius = {platen_radius}\n')
            f.write(f'    velocity = {disp}\n')
            f.write('    axis = "0. 0. 1."\n')
            f.write('    normal = "0. -1. 0."\n')
            f.write('    use_sector = true\n')
            f.write('    angle_min = 0.\n')
            f.write('    angle_max = 3.14159265359\n')
            f.write('    invert_displacement = true\n')
            f.write('  [../]\n')
            f.write('  [./top_x]\n')
            f.write('    type = CylindricalSurfaceDirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write(f'    boundary = "{specimen_top_surface}"\n')
            f.write(f'    center = "{xc_top} {yc_top} 0"\n')
            f.write(f'    radius = {platen_radius}\n')
            f.write(f'    velocity = {disp}\n')
            f.write('    axis = "0. 0. 1."\n')
            f.write('    normal = "0. -1. 0."\n')
            f.write('    use_sector = true\n')
            f.write('    angle_min = 0.\n')
            f.write('    angle_max = 3.14159265359\n')
            f.write('    invert_displacement = true\n')
            f.write('  [../]\n')
        if geometry != 'full':
            f.write('  [./back_sym]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write(f'    boundary = "{back_symmetry}"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            if geometry != 'half':
                f.write('  [./top_sym]\n')
                f.write('    type = DirichletBC\n')
                f.write('    variable = disp_y\n')
                f.write(f'    boundary = "{top_symmetry}"\n')
                f.write('    preset = true\n')
                f.write('    value = 0\n')
                f.write('  [../]\n')
                if geometry == 'eighth':
                    f.write('  [./side_sym]\n')
                    f.write('    type = DirichletBC\n')
                    f.write('    variable = disp_x\n')
                    f.write(f'    boundary = "{side_symmetry}"\n')
                    f.write('    preset = true\n')
                    f.write('    value = 0\n')
                    f.write('  [../]\n')
        f.write('[]\n')

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
        f.write(f'    boundary = "{stress_boundary}"\n')
        f.write('  []\n')
        f.write('[]\n')

        # Preconditioner
        moose_tools.write_preconditioner_block(f)

        # Execution and Timestepping
        dt = duration / 20
        f.write('[Executioner]\n')
        f.write('  type = Transient\n')
        f.write('  solve_type = NEWTON\n')
        f.write('  petsc_options_iname = "-pc_type -pc_factor_mat_solver_type"\n')
        f.write('  petsc_options_value = "lu    superlu_dist"\n')
        f.write('  line_search = "none"\n')
        f.write('  automatic_scaling = true\n')
        f.write('  nl_rel_tol = 1e-8\n')
        f.write('  nl_abs_tol = 1e-8\n')
        #f.write('  l_tol = 1e-3\n')
        #f.write('  l_max_its = 60\n')
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

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Write MOOSE input file for Brazilian disk simulation with rigid platens"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of Tardigrade-MOOSE file to write")
    parser.add_argument('--mesh', type=str, required=True,
        help='Specify the mesh file')
    parser.add_argument('--material-E', type=float, required=True,
        help='The elastic modulus of the specimen material')
    parser.add_argument('--material-nu', type=float, required=True,
        help='The Poisson ratio of the specimen material')
    parser.add_argument('--platen-radius', type=float, required=True,
        help='The radius of curvature of the Brazilian disk compression platen')
    parser.add_argument('--disp', type=float, required=True,
        help='Specify the compressive displacement to be applied')
    parser.add_argument('--duration', type=float, required=True,
        help='Specify the duration of the simulation')
    parser.add_argument('--specimen-bottom-surface', type=str, required=True,
        help='Specify the name of the specimen bottom contact surface')
    parser.add_argument('--specimen-top-surface', type=str, required=False, default=None,
        help='Specify the name of the specimen top contact surface. Required if "geometry" = "full."')
    parser.add_argument('--top-symmetry', type=str, required=False, default=None,
        help='Specify the name of the top symmetry surface. Required if "geometry" = "quarter" or "eighth."')
    parser.add_argument('--back-symmetry', type=str, required=False, default=None,
        help='Specify the name of the back symmetry surface. Required if "geometry" = "quarter" or "eighth" or "half"')
    parser.add_argument('--side-symmetry', type=str, required=False, default=None,
        help='Specify the name of the side symmetry surface. Required if "geometry" = "quarter" or "eighth."')
    parser.add_argument('--xc-bot', type=float, required=False, default=0.,
        help='Specify the x-position of the center of the circular bottom surface arc')
    parser.add_argument('--yc-bot', type=float, required=False, default=0.,
        help='Specify the y-position of the center of the circular bottom surface arc')
    parser.add_argument('--xc-top', type=float, required=False, default=0.,
        help='Specify the x-position of the center of the circular top surface arc')
    parser.add_argument('--yc-top', type=float, required=False, default=0.,
        help='Specify the y-position of the center of the circular top surface arc')
    parser.add_argument('--geometry', type=str, required=False, default='full',
        help='Specify the geometry/symmetry type: "full," "half," "quarter," or "eighth"')

    return parser

if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(build_input(output_file=args.output_file,
                         mesh_file=args.mesh,
                         material_E=args.material_E,
                         material_nu=args.material_nu,
                         platen_radius=args.platen_radius,
                         disp=args.disp,
                         duration=args.duration,
                         specimen_bottom_surface=args.specimen_bottom_surface,
                         specimen_top_surface=args.specimen_top_surface,
                         top_symmetry=args.top_symmetry,
                         back_symmetry=args.back_symmetry,
                         side_symmetry=args.side_symmetry,
                         xc_bot=args.xc_bot,
                         yc_bot=args.yc_bot,
                         xc_top=args.xc_top,
                         yc_top=args.yc_top,
                         geometry=args.geometry,
                         ))
