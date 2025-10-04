import os
import sys
import argparse
import yaml
import inspect

import pandas

import MOOSE_input_deck_tools as moose_tools


def build_input(output_file, mesh_file, BCs, pressure, start, duration, dt, ref_density, height,
                parameter_sets=None, calibration_map=None, phi_BC=None):
    '''Write a Tardigrade-MOOSE input file for dynamic simulation

    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param str mesh_file: The name of the mesh file
    :param str BCs: The type of boundary conditions, either "slip" or "clamp"
    :param float pressure: The pressure to be applied
    :param float duration: The duration of the simulation
    :param float start: The time when heaviside pressure is applied
    :param float dt: The fixed time increment
    :param float ref_density: Density in reference configuration (Mg/mm^3)
    :param float height: Height of the geometry
    :param list parameter_sets: The list of yaml files containing calibration results, required if calibration-map is not provided
    :param str calibration_map: Optional yaml file containing names of calibration files
    :param str phi_BC: Optional string specifying nodeset to force micro deformation components to be zero

    :returns: ``output_file``
    '''

    assert os.path.exists(mesh_file), f"Mesh file not found: {mesh_file}"

    # unpack parameter set files if calibration map is provided
    if calibration_map:
        stream = open(calibration_map, 'r')
        calibrations = yaml.load(stream, Loader=yaml.FullLoader)
        stream.close()
        # Get number of elements and assign the default parameter set
        num_elements = len(calibrations.keys()) - 2
        parameter_sets = [calibrations['ignore_boundary_yml'] for i in range(0, num_elements)]
        # Override defaults for the elements not located on the boundary
        summary_file = calibrations['ignore_boundary_summary_file']
        if os.path.exists(summary_file):
            df = pandas.read_csv(summary_file, sep=',')
            for element in df['element']:
                parameter_sets[element] = calibrations[str(element)]
    else:
        assert parameter_sets is not None

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
        f.write('  #Inertial force balance equations\n')
        moose_tools.write_dynamic_kernels(f, ref_density)
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

        # Reaction Force
        f.write('[Postprocessors]\n')
        f.write('  [bot_react_z]\n')
        f.write('    type = NodalSum\n')
        f.write('    variable = force_z\n')
        f.write('    boundary = "top"\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[Postprocessors]\n')
        f.write('  [disp_x_p]\n')
        f.write('    type = PointValue\n')
        f.write(f'    point = "0 0 {height/2}"\n')
        f.write('    variable = disp_z\n')
        f.write('  []\n')
        f.write('\n')
        f.write('  [middle_s33]\n')
        f.write('    type = PointValue\n')
        f.write(f'    point = "0 0 {height/2}"\n')
        f.write('    variable = sigma_33\n')
        f.write('  []\n')
        f.write('[]\n')

        # BCs
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
            f.write('    type = FunctionNeumannBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "top"\n')
            f.write('    function = heaviside\n')
            f.write('  [../]\n')
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
            f.write('    type = FunctionNeumannBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "top"\n')
            f.write('    function = heaviside\n')
            f.write('  [../]\n')
        else:
            print('Specify a valid BC type!')
        # Option to force Phis to be zero
        if phi_BC is not None:
            moose_tools.write_phi_BCs(f, phi_BC)
        f.write('[]\n')
        f.write('\n')

        # Heaviside load
        f.write('\n')
        f.write('[Functions]\n')
        f.write('  [heaviside]\n')
        f.write('    type  = ADParsedFunction\n')
        if start < dt:
            f.write(f'    expression = "if(t<{start},0.0,{pressure})" # Heaviside function in MPa i.e N/mm2\n')
        else:
            f.write(f'    expression = "if(t<{start},t*{pressure},{pressure})" # Heaviside function in MPa i.e N/mm2\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')

        # Materials
        f.write('# Materials\n')
        f.write('[Materials]\n')
        # Load in parameter data for each filter domain / element
        if len(parameter_sets) > 1:
            for i, set in enumerate(parameter_sets):
                # Load yaml file
                stream = open(set, 'r')
                UI = yaml.load(stream, Loader=yaml.FullLoader)
                stream.close()
                moose_tools.write_elastic_material_card(f, UI, i)
        else:
            # Load yaml file
            set = parameter_sets[0]
            stream = open(set, 'r')
            UI = yaml.load(stream, Loader=yaml.FullLoader)
            stream.close()
            moose_tools.write_elastic_material_card(f, UI)
        f.write('[]\n')
        f.write('\n')

        # Preconditioner
        moose_tools.write_preconditioner_block(f)

        # Execution and Timestepping
        f.write('[Executioner]\n')
        f.write('  type = Transient\n')
        f.write('  scheme = newmark-beta\n')
        f.write(f'  end_time = {duration + start}\n')
        f.write(f'  dt        = {dt}\n')
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

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Write Tardigrade-MOOSE input file for dynamic simulation"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="The name of Tardigrade-MOOSE file to write")
    parser.add_argument('--mesh', type=str, required=True,
        help='The mesh file')
    parser.add_argument('--parameter-sets', nargs="+", required=False, default=None,
        help='List of yaml files containing calibration results, required if calibration-map is not provided')
    parser.add_argument('--calibration-map', type=str, required=False, default=None,
        help='Optional yaml file containing names of calibration files')
    parser.add_argument('--BCs', type=str, required=True,
        help='The type of boundary conditions, either "slip" or "clamp"')
    parser.add_argument('--pressure', type=float, required=True,
        help='The pressure to be applied')
    parser.add_argument('--start', type=float, required=True,
        help='The time when heaviside pressure is applied')
    parser.add_argument('--duration', type=float, required=True,
        help='The duration of the simulation')
    parser.add_argument('--dt', type=float, required=True,
        help='The fixed time increment')
    parser.add_argument('--ref-density', type=float, required=True,
        help='Density in reference configuration (Mg/mm^3)')
    parser.add_argument('--height', type=float, required=True,
        help='Height of the geometry')
    parser.add_argument('--phi-BC', type=str, required=False, default=None,
        help='Optional string specifying nodeset to force micro deformation components to be zero')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(build_input(output_file=args.output_file,
                         mesh_file=args.mesh,
                         BCs=args.BCs,
                         pressure=args.pressure,
                         start=args.start,
                         duration=args.duration,
                         dt=args.dt,
                         parameter_sets=args.parameter_sets,
                         calibration_map=args.calibration_map,
                         ref_density=args.ref_density,
                         height=args.height,
                         phi_BC=args.phi_BC,
                         ))