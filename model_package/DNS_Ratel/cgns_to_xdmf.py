#!python
import argparse
import pathlib
import sys

import CGNS.MAP
import CGNS.PAT.cgnsutils as CGU
import meshio
import numpy
import pandas

import file_io.xdmf
import calibrate_element


def unpack_CGNS_coordinates(tree, dist_factor):
    '''Collect and scale the reference positions of nodes from a CGNS tree

    :params tree tree: A CGNS tree
    :param float dist_factor: Optional argument to scale DNS displacements and coordinates

    '''

    root_path = '/Diagnostic/Zone/GridCoordinates'
    x_path =f'{root_path}/CoordinateX'
    y_path =f'{root_path}/CoordinateY'
    z_path =f'{root_path}/CoordinateZ'

    coord_x = dist_factor*CGU.getValueByPath(tree, x_path)
    coord_y = dist_factor*CGU.getValueByPath(tree, y_path)
    coord_z = dist_factor*CGU.getValueByPath(tree, z_path)

    reference_positions = numpy.array([coord_x, coord_y, coord_z]).T

    return reference_positions


def unpack_diagnostic(tree, diagnostic_path, quantity):
    '''Collect values of a diagnostic quantity from a CGNS tree

    :params tree tree: A CGNS tree
    :params str diagnostic_path: The root CGNS path for all diagnostic quantities in a given timestep
    :params str quantity: The diagnostic quantity to extract

    :returns: Array of diagnostic quantity data
    '''

    array_out = CGU.getValueByPath(tree, f'{diagnostic_path}{quantity}')

    return array_out


def collect_and_convert_to_XDMF(input_files, output_file, dist_factor, stress_factor, density_factor, damage):
    '''Write XDMF file of collected Ratel DNS results for Micromorphic Filter

    :param dict results: dictionary of results
    :param dict node_dict: dictionary of nodal coordinates
    :param output_file: Name for XDMF file pair output for the Micromorphic Filter
    :param float dist_factor: Argument to scale DNS displacements and coordinates
    :param float stress_factor: Argument to scale DNS stresses
    :param float density_factor: Factor to scale current density (if provided in the DNS results\
                                 to Mg/mm^3

    :returns: ``{output_file}.xdmf`` and ``{outptu_file}.h5``
    '''

    # shorthand for field keys
    time_path = '/Diagnostic/TimeIterValues/TimeValues'
    root_diagnostic_path = '/Diagnostic/Zone/FlowSolution'
    disp_x = 'projected.displacement_x'
    disp_y = 'projected.displacement_y'
    disp_z = 'projected.displacement_z'
    sig_xx = 'projected.Cauchy_stress_xx'
    sig_xy = 'projected.Cauchy_stress_xy'
    sig_xz = 'projected.Cauchy_stress_xz'
    sig_yy = 'projected.Cauchy_stress_yy'
    sig_yz = 'projected.Cauchy_stress_yz'
    sig_zz = 'projected.Cauchy_stress_zz'
    n_vol = 'dual.nodal_volume'
    n_dens = 'dual.nodal_density'
    Jdef = 'projected.J'
    damage_field = 'projected.damage'

    data_filename=output_file
    xdmf = file_io.xdmf.XDMF(output_filename=data_filename)

    point_name = 'points'
    conn_name = 'connectivity'

    # assume times are even spaced pseudo-timesteps
    num_steps = len(input_files)
    step_names = [f'timestep_{i}' for i in range(0, num_steps)]

    #for step_name in step_names:
    for step_name, input_file in zip(step_names, input_files):
        print(f'step = {step_name}')

        # open CGNS file
        (tree,links,paths)=CGNS.MAP.load(input_file)

        # time
        t = CGU.getValueByPath(tree, time_path)[0]

        # Collect paths
        all_paths = CGU.getAllPaths(tree)
        diagnostic_paths = [i for i in all_paths if root_diagnostic_path in i]

        # update diagnostic path for specific timestep_
        diagnostic_path = f"{root_diagnostic_path}{diagnostic_paths[0].split('FlowSolution')[-1].split('/')[0]}/"

        # get the reference positions
        if t == 0.:
            print('collecting reference positions')
            reference_positions = unpack_CGNS_coordinates(tree, dist_factor)
            ndata = reference_positions.shape[0]

        ## initialization stuff
        grid = xdmf.addGrid(xdmf.output_timegrid, {})
        xdmf.addTime(grid, t)
        xdmf.addPoints(grid, reference_positions, duplicate=point_name)
        xdmf.addConnectivity(grid, "POLYVERTEX", numpy.array([v for v in range(ndata)]).reshape((-1,1)), duplicate=conn_name)

        # get the unique displacements
        unique_displacements = dist_factor*numpy.array([
            unpack_diagnostic(tree, diagnostic_path, disp_x),
            unpack_diagnostic(tree, diagnostic_path, disp_y),
            unpack_diagnostic(tree, diagnostic_path, disp_z)]).T
        print(f"shape of unique displacements = {numpy.shape(unique_displacements)}")
        xdmf.addData(grid, "disp", unique_displacements, "Node", dtype='d')

        # get the unique positions
        unique_positions = unique_displacements + reference_positions
        print(f"shape of unique positions = {numpy.shape(unique_positions)}")
        xdmf.addData(grid, "coord", unique_positions, "Node", dtype='d')

        # get the velocity <-- fix for dynamic DNS!
        unique_velocities = numpy.zeros(numpy.shape(unique_positions))
        print(f"shape of unique velocities = {numpy.shape(unique_velocities)}")
        xdmf.addData(grid, "vel", unique_velocities, "Node", dtype='d')

        # get the acceleration <-- fix for dynamic DNS!
        unique_accelerations = numpy.zeros(numpy.shape(unique_positions))
        print(f"shape of unique accelerations = {numpy.shape(unique_accelerations)}")
        xdmf.addData(grid, "acc", unique_accelerations, "Node", dtype='d')

        # get the stresses
        unique_stresses = stress_factor*numpy.array([
            unpack_diagnostic(tree, diagnostic_path, sig_xx),
            unpack_diagnostic(tree, diagnostic_path, sig_xy),
            unpack_diagnostic(tree, diagnostic_path, sig_xz),
            unpack_diagnostic(tree, diagnostic_path, sig_xy),
            unpack_diagnostic(tree, diagnostic_path, sig_yy),
            unpack_diagnostic(tree, diagnostic_path, sig_yz),
            unpack_diagnostic(tree, diagnostic_path, sig_xz),
            unpack_diagnostic(tree, diagnostic_path, sig_yz),
            unpack_diagnostic(tree, diagnostic_path, sig_zz)]).T
        print(f"shape of stresses = {numpy.shape(unique_stresses)}")
        xdmf.addData(grid, "stress", unique_stresses, "Node", dtype='d')

        # Get the reference volumes and Jacobian to calculate current volumes
        vol_factor = dist_factor*dist_factor*dist_factor
        reference_volumes = vol_factor*unpack_diagnostic(
            tree, diagnostic_path, n_vol).reshape((-1,1))
        unique_jacobians = unpack_diagnostic(tree, diagnostic_path, Jdef).reshape((-1,1))
        current_volumes = unique_jacobians*reference_volumes
        print(f"shape of vol = {numpy.shape(current_volumes)}")
        print(f"total volume = {numpy.sum(current_volumes)}")
        xdmf.addData(grid, "volume", current_volumes, "Node", dtype='d')

        # Get the densities
        unique_densities = density_factor*unpack_diagnostic(
            tree, diagnostic_path, n_dens).reshape((-1,1))
        print(f"shape of density = {numpy.shape(unique_densities)}")
        xdmf.addData(grid, "density", unique_densities, "Node", dtype='d')

        # Option for damage
        if damage == True:
            unique_damage = unpack_diagnostic(tree, diagnostic_path, damage_field).reshape((-1,1))
            # "clamp" values within [0, 1]
            n_damages_below_zero = numpy.shape(unique_damage[unique_damage < 0.])[0]
            n_damages_above_one = numpy.shape(unique_damage[unique_damage > 1.])[0]
            if n_damages_below_zero > 0:
                print(f'\tnumber of points with damage < 0 = {n_damages_below_zero}. Setting values to zero.')
                unique_damage[unique_damage < 0.] = 0.
            if n_damages_above_one > 0:
                print(f'\tnumber of points with damage > 1 = {n_damages_above_one}. Setting values to one.')
                unique_damage[unique_damage > 1.] = 1.
            print(f"shape of damage = {numpy.shape(unique_damage)}")
            xdmf.addData(grid, "damage", unique_damage, "Node", dtype='d')

    xdmf.write()
    print("XDMF file written!")

    return 0


def convert_CGNS_to_XDMF(input_files, output_file, dist_factor=1, stress_factor=1, density_factor=1, dump_all_33_stresses=None, damage=False):
    '''Driving function to call functions for parsing Ratel CGNS results and writing XDMF output

    :param list input_file: The input VTK files containing Ratel DNS results
    :param str output_file: The output filename for the h5 + XDMF file pair
    :param float dist_factor: Optional argument to scale DNS displacements and coordinates, default=1
    :param float stress_factor: Optional argument to scale DNS stresses, default=1
    :param float density_factor: Optional factor to scale current density (if provided in the DNS results\
                                 to Mg/mm^3, default=1
    '''
    # collect VTU results and convert to XDMF
    collect_and_convert_to_XDMF(input_files, output_file, dist_factor, stress_factor, density_factor, damage)

    # TODO: Dump Cauchy 33 stresses to csv
    # if dump_all_33_stresses:
        # last_step = [key for key in results.keys()][-1]
        # cauchy33 = results[last_step]['Cauchy_stress_zz']
        # df = pandas.DataFrame({'quantity': cauchy33,})
        # df.to_csv(dump_all_33_stresses, header=True, sep=',', index=False)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Convert Ratel DNS results to XDMF format"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-i', '--input-files', nargs="+",
        help='Specify the input VTK files containing Ratel DNS results')
    parser.add_argument('-o', '--output-file', type=str,
        help='Specify the output filename for the h5 + XDMF file pair')
    parser.add_argument('--dist-factor', type=float, required=False, default=1,
        help='Optional argument to scale DNS displacements and coordinates')
    parser.add_argument('--stress-factor', type=float, required=False, default=1,
        help='Optional argument to scale DNS stresses')
    parser.add_argument('--density-factor', type=float, required=False, default=1,
         help='Optional factor to scale current density (if provided in the DNS results\
               to Mg/mm^3')
    parser.add_argument('--dump-all-33-stresses', type=str, required=False, default=None,
        help='Optional filename to dump all 33 stresses from DNS')
    parser.add_argument('--damage', type=str, required=False, default=False,
        help='Optional filename to dump all 33 stresses from DNS')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(convert_CGNS_to_XDMF(input_files=args.input_files,
                                  output_file=args.output_file,
                                  dist_factor=args.dist_factor,
                                  stress_factor=args.stress_factor,
                                  density_factor=args.density_factor,
                                  dump_all_33_stresses=args.dump_all_33_stresses,
                                  damage=calibrate_element.str2bool(args.damage),
                                  ))