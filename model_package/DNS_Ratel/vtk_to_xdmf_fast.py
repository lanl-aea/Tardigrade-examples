#!python
import argparse
import pathlib
import sys

import meshio
import numpy
import pandas

import file_io.xdmf
import calibrate_element

def collect_and_convert_to_XDMF(input_files, output_file, dist_factor, stress_factor, ref_density, density_factor, damage):
    '''Write XDMF file of collected Ratel DNS results for Micromorphic Filter

    :param dict results: dictionary of results
    :param dict node_dict: dictionary of nodal coordinates
    :param output_file: Name for XDMF file pair output for the Micromorphic Filter
    :param float dist_factor: Optional argument to scale DNS displacements and coordinates, default=1
    :param float stress_factor: Optional argument to scale DNS stresses, default=1
    :param float ref_density: Optional argument to specify the reference density to be converted to
                              current density by the Jacobian of deformation if current density is
                              not reported in the DNS results, default=2.e-9
    :param float density_factor: Optional factor to scale current density (if provided in the DNS results\
                                 to Mg/tonne^3, default=1

    :returns: ``{output_file}.xdmf`` and ``{outptu_file}.h5``
    '''

    # shorthand for field keys
    disp_x = 'diagnostic_quantitiesprojected.displacement_x'
    disp_y = 'diagnostic_quantitiesprojected.displacement_y'
    disp_z = 'diagnostic_quantitiesprojected.displacement_z'
    sig_xx = 'diagnostic_quantitiesprojected.Cauchy_stress_xx'
    sig_xy = 'diagnostic_quantitiesprojected.Cauchy_stress_xy'
    sig_xz = 'diagnostic_quantitiesprojected.Cauchy_stress_xz'
    sig_yy = 'diagnostic_quantitiesprojected.Cauchy_stress_yy'
    sig_yz = 'diagnostic_quantitiesprojected.Cauchy_stress_yz'
    sig_zz = 'diagnostic_quantitiesprojected.Cauchy_stress_zz'
    n_vol = 'diagnostic_quantitiesdual.nodal_volume'
    n_dens = 'diagnostic_quantitiesdual.nodal_density'
    Jdef = 'diagnostic_quantitiesprojected.J'
    damage_field = 'diagnostic_quantitiesprojected.damage'

    data_filename=output_file
    xdmf = file_io.xdmf.XDMF(output_filename=data_filename)

    point_name = 'points'
    conn_name = 'connectivity'

    # assume times are even spaced pseudo-timesteps
    num_times = len(input_files)
    times = numpy.linspace(0, 1, num_times)
    step_names = [f'timestep_{i}' for i in range(0, num_times)]

    #for step_name in step_names:
    for t, step_name, input_file in zip(times, step_names, input_files):
        print(f'step = {step_name}')

        # open vtk file
        mesh=meshio.read(input_file)

        # get the reference positions
        if t == 0.:
            print('collecting reference positions')
            reference_positions = dist_factor*(mesh.points)
            ndata = reference_positions.shape[0]

        ## initialization stuff
        grid = xdmf.addGrid(xdmf.output_timegrid, {})
        xdmf.addTime(grid, t)
        xdmf.addPoints(grid, reference_positions, duplicate=point_name)
        xdmf.addConnectivity(grid, "POLYVERTEX", numpy.array([v for v in range(ndata)]).reshape((-1,1)), duplicate=conn_name)

        # get the unique displacements
        unique_displacements = dist_factor*numpy.array([mesh.point_data[disp_x].flatten(),
                                                        mesh.point_data[disp_y].flatten(),
                                                        mesh.point_data[disp_z].flatten()]).T
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
        unique_stresses = stress_factor*numpy.array([mesh.point_data[sig_xx].flatten(),
                                                     mesh.point_data[sig_xy].flatten(),
                                                     mesh.point_data[sig_xz].flatten(),
                                                     mesh.point_data[sig_xy].flatten(),
                                                     mesh.point_data[sig_yy].flatten(),
                                                     mesh.point_data[sig_yz].flatten(),
                                                     mesh.point_data[sig_xz].flatten(),
                                                     mesh.point_data[sig_yz].flatten(),
                                                     mesh.point_data[sig_zz].flatten()]).T
        print(f"shape of stresses = {numpy.shape(unique_stresses)}")
        xdmf.addData(grid, "stress", unique_stresses, "Node", dtype='d')

        # Get the volumes
        vol_factor = dist_factor*dist_factor*dist_factor
        unique_volumes = vol_factor*mesh.point_data[n_vol].flatten().reshape((-1,1))
        print(f"shape of vol = {numpy.shape(unique_volumes)}")
        print(f"total volume = {numpy.sum(unique_volumes)}")
        xdmf.addData(grid, "volume", unique_volumes, "Node", dtype='d')

        # Get the densities
        if n_dens in mesh.point_data.keys():
            unique_densities = density_factor*mesh.point_data[n_dens].flatten().reshape((-1,1))
        else:
            unique_densities = ref_density*numpy.reciprocal(mesh.point_data[Jdef].flatten().reshape((-1,1)))
        print(f"shape of density = {numpy.shape(unique_densities)}")
        xdmf.addData(grid, "density", unique_densities, "Node", dtype='d')

        # Option for damage
        if damage == True:
            unique_damage = density_factor*mesh.point_data[damage_field].flatten().reshape((-1,1))
            print(f"shape of damage = {numpy.shape(unique_damage)}")
            xdmf.addData(grid, "damage", unique_damage, "Node", dtype='d')

    xdmf.write()
    print("XDMF file written!")

    return 0


def convert_VTK_to_XDMF(input_files, output_file, dist_factor=1, stress_factor=1, ref_density=2.e-9, density_factor=1, dump_all_33_stresses=None, damage=False):
    '''Driving function to call functions for parsing Ratel VTK results and writing XDMF output

    :param list input_file: The input VTK files containing Ratel DNS results
    :param str output_file: The output filename for the h5 + XDMF file pair
    :param float dist_factor: Optional argument to scale DNS displacements and coordinates, default=1
    :param float stress_factor: Optional argument to scale DNS stresses, default=1
    :param float ref_density: Optional argument to specify the reference density to be converted to
                              current density by the Jacobian of deformation if current density is
                              not reported in the DNS results, default=2.e-9
    :param float density_factor: Optional factor to scale current density (if provided in the DNS results\
                                 to Mg/tonne^3, default=1
    '''
    # collect VTU results and convert to XDMF
    collect_and_convert_to_XDMF(input_files, output_file, dist_factor, stress_factor, ref_density, density_factor, damage)

    # Dump Cauchy 33 stresses to csv
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
    parser.add_argument('--ref-density', type=float, required=False, default=2.e-9,
        help='Optional argument to specify the reference density to be converted to\
              current density by the Jacobian of deformation if current density is\
              not reported in the DNS results')
    parser.add_argument('--density-factor', type=float, required=False, default=1,
         help='Optional factor to scale current density (if provided in the DNS results\
               to Mg/tonne^3')
    parser.add_argument('--dump-all-33-stresses', type=str, required=False, default=None,
        help='Optional filename to dump all 33 stresses from DNS')
    parser.add_argument('--damage', type=str, required=False, default=False,
        help='Optional filename to dump all 33 stresses from DNS')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(convert_VTK_to_XDMF(input_files=args.input_files,
                                 output_file=args.output_file,
                                 dist_factor=args.dist_factor,
                                 stress_factor=args.stress_factor,
                                 ref_density=args.ref_density,
                                 density_factor=args.density_factor,
                                 dump_all_33_stresses=args.dump_all_33_stresses,
                                 damage=calibrate_element.str2bool(args.damage),
                                 ))