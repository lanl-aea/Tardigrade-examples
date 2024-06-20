# Imports
import sys
import os
import inspect
import argparse

import meshio
import numpy as np
import pandas

import file_io.xdmf


def collect_VTK_output(input_files):
    '''Parse the Ratel DNS VTK output into a results dictionary

    :param list input_file: The input VTK files containing Ratel DNS results

    :returns: dictionary of results, dictionary of nodal coordinates
    '''

    # collect into results dictionary
    results = {}

    # get time
    num_times = len(input_files)
    times = np.linspace(0, 1, num_times)

    # setup node_dict
    node_dict = {}
    # loop over input files:
    for i, input_file in enumerate(input_files):
        # store temporary results into another dictionary
        small_dict = {}

        # store time
        small_dict['time'] = times[i]
        small_dict['frame'] = int(i)
        print(f"TIME = {times[i]}")

        # read data with meshio
        mesh=meshio.read(input_file)

        # get initial node coordinates for i==0
        if i == 0:
            x = list(mesh.points[:,0])
            y = list(mesh.points[:,1])
            z = list(mesh.points[:,2])
            node_dict['CoordinateX'] = x
            node_dict['CoordinateY'] = y
            node_dict['CoordinateZ'] = z

        # get field data
        for key in mesh.point_data.keys():
            row=key.split('.')
            field_name = row[1]

            # store field into small_dict
            field_data = list(mesh.point_data[key].flatten())
            small_dict[field_name] = field_data

        # store small_dict into larger results dictionary
        results[f'timestep_{i}'] = small_dict

    return(results, node_dict)


def convert_to_XDMF(results, node_dict, output_file, dist_factor, stress_factor, ref_density, density_factor):
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

    data_filename=output_file
    xdmf = file_io.xdmf.XDMF(output_filename=data_filename)

    # get the reference positions
    reference_positions = []
    for i, x in enumerate(node_dict['CoordinateX']):
        y = node_dict['CoordinateY'][i]
        z = node_dict['CoordinateZ'][i]
        reference_positions.append([dist_factor*x,dist_factor*y,dist_factor*z])
    reference_positions = np.array(reference_positions)
    ndata = reference_positions.shape[0]

    # # get reference volumes
    # reference_volumes = np.array([vol for vol in results['timestep_0']['nodal_volume']
    # reference_volumes = reference_volumes.reshape((-1,1))

    point_name = 'points'
    conn_name = 'connectivity'

    # get steps names
    step_names = [key for key in results.keys()]

    for step_name in step_names:
        print(f'step = {step_name}')

        # Grab current time
        t = results[step_name]['time']

        ## initialization stuff
        grid = xdmf.addGrid(xdmf.output_timegrid, {})
        xdmf.addTime(grid, t)
        xdmf.addPoints(grid, reference_positions, duplicate=point_name)
        xdmf.addConnectivity(grid, "POLYVERTEX", np.array([v for v in range(ndata)]).reshape((-1,1)), duplicate=conn_name)

        # get the unique displacements
        unique_displacements = []
        for i, x in enumerate(results[step_name]['displacement_x']):
            y = results[step_name]['displacement_y'][i]
            z = results[step_name]['displacement_z'][i]
            unique_displacements.append([dist_factor*x,dist_factor*y,dist_factor*z])
        unique_displacements = np.array(unique_displacements)
        print(f"shape of unique displacements = {np.shape(unique_displacements)}")
        xdmf.addData(grid, "disp", unique_displacements, "Node", dtype='d')

        # get the unique positions
        unique_positions = unique_displacements + reference_positions
        print(f"shape of unique positions = {np.shape(unique_positions)}")
        xdmf.addData(grid, "coord", unique_positions, "Node", dtype='d')

        # get the velocity <-- fix for dynamic DNS!
        unique_velocities = np.zeros(np.shape(unique_positions))
        print(f"shape of unique velocities = {np.shape(unique_velocities)}")
        xdmf.addData(grid, "vel", unique_velocities, "Node", dtype='d')

        # get the acceleration <-- fix for dynamic DNS!
        unique_accelerations = np.zeros(np.shape(unique_positions))
        print(f"shape of unique accelerations = {np.shape(unique_accelerations)}")
        xdmf.addData(grid, "acc", unique_accelerations, "Node", dtype='d')

        # get the stresses
        data = []
        for i, s in enumerate(results[step_name]['Cauchy_stress_xx']):
            Sxx = s*stress_factor
            Syy = results[step_name]['Cauchy_stress_yy'][i]*stress_factor
            Szz = results[step_name]['Cauchy_stress_zz'][i]*stress_factor
            Syz = results[step_name]['Cauchy_stress_yz'][i]*stress_factor
            Sxz = results[step_name]['Cauchy_stress_xz'][i]*stress_factor
            Sxy = results[step_name]['Cauchy_stress_xy'][i]*stress_factor
            data.append([Sxx, Sxy, Sxz,
                         Sxy, Syy, Syz,
                         Sxz, Syz, Szz])
        unique_stresses = np.array(data)
        print(f"shape of stresses = {np.shape(unique_stresses)}")
        xdmf.addData(grid, "stress", unique_stresses, "Node", dtype='d')

        # Get the volumes
        unique_volumes = np.array([vol for vol in results[step_name]['nodal_volume']])
        unique_volumes = unique_volumes.reshape((-1,1))
        print(f"shape of vol = {np.shape(unique_volumes)}")
        print(f"total volume = {np.sum(unique_volumes)}")
        xdmf.addData(grid, "volume", unique_volumes, "Node", dtype='d')

        # Get the densities
        if 'mass_density' in results[step_name].keys():
            unique_densities = np.array([den*density_factor for den in results[step_name]['mass_density']])
        else:
            unique_densities = np.array([(ref_density/J) for J in results[step_name]['J']])
        unique_densities = unique_densities.reshape((-1,1))
        print(f"shape of density = {np.shape(unique_densities)}")
        xdmf.addData(grid, "density", unique_densities, "Node", dtype='d')

    xdmf.write()
    print("XDMF file written!")

    return 0


def convert_VTK_to_XDMF(input_files, output_file, dist_factor=1, stress_factor=1, ref_density=2.e-9, density_factor=1, dump_all_33_stresses=None):
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

    # parse vtk file
    results, node_dict = collect_VTK_output(input_files)

    # convert to XDMF
    convert_to_XDMF(results, node_dict, output_file, dist_factor, stress_factor, ref_density, density_factor)

    # Dump Cauchy 33 stresses to csv
    if dump_all_33_stresses:
        last_step = [key for key in results.keys()][-1]
        cauchy33 = results[last_step]['Cauchy_stress_zz']
        df = pandas.DataFrame({'quantity': cauchy33,})
        df.to_csv(dump_all_33_stresses, header=True, sep=',', index=False)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Convert Ratel DNS results to XDMF format"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
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
                                 ))