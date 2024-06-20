# Imports
import sys
import os
import inspect
import argparse

import h5py
import numpy as np
import pandas

import file_io.xdmf

file_path = os.path.dirname(os.path.abspath(__file__))


def str2bool(v):
    '''Function for converting string to Boolean. Borrowed from: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse

    :param str/bool v: A string or boolean indicating a True or False value

    :returns: True or False
    '''

    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def interpolate_to_ip_c3d8(node_array, mesh):
    '''interpolate a vector field from the nodes to the integration points of a trilinear hexahedral element (C3D8)

    :param array-like node_array: nodal data to be interpolated
    :param array-like mesh: the element connectivity for all elements

    :returns: dictionary of interpolated results
    '''

    numips = 8
    numelem = np.shape(mesh)[0]
    results = np.zeros([numelem, numips, 3])
    # set Gauss point coordinates in xi,eta,zeta space
    const=1/(np.sqrt(3))
    xi_vect=np.array([[-const, -const, -const],
                      [ const, -const, -const],
                      [ const,  const, -const],
                      [-const,  const, -const],
                      [-const, -const,  const],
                      [ const, -const,  const],
                      [ const,  const,  const],
                      [-const,  const,  const]])

    # loop over all elements in mesh
    for e, n in enumerate(mesh):
        # loop over each node of the element
        node_field = []
        for node in n:
            # get the field values for that node
            node_field.append(node_array[node-1])

        node_field = np.array(node_field).flatten(order='C')

        # interpolate from nodes to integration points
        for ip in range(0,numips):
            xi, eta, zeta = xi_vect[ip,0], xi_vect[ip,1], xi_vect[ip,2]
            N1 = (1-xi)*(1-eta)*(1-zeta)/8
            N2 = (1+xi)*(1-eta)*(1-zeta)/8
            N3 = (1+xi)*(1+eta)*(1-zeta)/8
            N4 = (1-xi)*(1+eta)*(1-zeta)/8
            N5 = (1-xi)*(1-eta)*(1+zeta)/8
            N6 = (1+xi)*(1-eta)*(1+zeta)/8
            N7 = (1+xi)*(1+eta)*(1+zeta)/8
            N8 = (1-xi)*(1+eta)*(1+zeta)/8

            Nu = np.array([
                [ N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0, N5, 0, 0, N6, 0, 0, N7, 0, 0, N8, 0, 0],
                [ 0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0, N5, 0, 0, N6, 0, 0, N7, 0, 0, N8, 0],
                [ 0, 0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0, N5, 0, 0, N6, 0, 0, N7, 0, 0, N8]])
            solve = np.matmul(Nu, node_field)

            results[e, ip, :] = solve

    print(f'nodal_field shape = {np.shape(results)}')
    return(results)


def interpolate_to_center_c3d8(node_array, mesh):
    '''Average a vector or tensor field from the nodes to the center of a trilinear hexahedral element (C3D8)

    :param array-like node_array: nodal data to be interpolated
    :param array-like mesh: the element connectivity for all elements

    :returns: dictionary of interpolated results
    '''

    numpts = 1
    numelem = np.shape(mesh)[0]
    results = np.zeros([numelem, numpts, 3])

    for e, n in enumerate(mesh):
        node_field = []
        for node in n:
            # get the field values for that node
            node_field.append(node_array[node-1])

        #node_field = np.array(node_field).flatten(order='C')
        node_field = np.array(node_field)

        # average over node_field
        solve = np.mean(node_field, axis=0)
        results[e, 0, :] = solve

    print(f'nodal_field shape = {np.shape(results)}')
    return(results)


def parse_input(input_file, elem_path, node_path, mesh_path, collocation_option, velocities=False, accelerations=False, specific_frames=None):
    '''Parse the HDF5 file output by ODBextract (WAVES tool)

    :param str input_file: HDF5 file of Abaqus results extracted using the
        ODB_extract module of WAVES.
    :param str elem_path: HDF5 path to element data
    :param str node_path: HDF5 path to node data
    :param str mesh_path: HDF5 path to mesh data
    :param str collocation option: String specifying "center" to collocate to element center or "ip" for integration points
    :param bool velocities: Boolean whether or not to collect DNS velocity data
    :param bool accelerations: Boolean whether or not to collect DNS accelerations data
    :param list specific_frames: An optional list of frame numbers for converting XDMF data

    :returns: dictionary of results, list of frames, list of time increments
    '''

    with h5py.File(input_file, 'r') as file:
        elem_fields = file[elem_path]
        node_fields = file[node_path]
        times = np.array(elem_fields['time'])

        # frames
        if specific_frames:
            num_frames = len(specific_frames)
            frames = [int(f) for f in specific_frames]
        else:
            num_frames = len(times)
            frames = [i for i in range(0, num_frames)]

        # get nodal locations
        node = np.array(file[mesh_path]['node'])
        node_location = np.array(file[mesh_path]['node_location'])
        c3d8_mesh = np.array(file[mesh_path]['C3D8_mesh']) # (1640 elem) x (8 nodes / elem)

        # loop over frames and get indices for each field
        results = {}
        for f in frames:

            # unpack element field results
            IVOL_elem   = np.array(elem_fields['IVOL'][0][f])    # (1) x (frames) x (elem) x (8 ips) --> (elem) x (8 ips)
            EVOL_elem   = np.array(elem_fields['EVOL'][0][f])    # (1) x (frames) x (elem) x (8 ips) --> (elem) x (8 ips), most of ips are nan
            S_elem      = np.array(elem_fields['S'][0][f])       # (1) x (frames) x (elem) x (8 ips) x (6 comp) --> (elem) x (8 ips) x (6 comp)
            COORD_elem  = np.array(elem_fields['COORD'][0][f])   # (1) x (frames) x (elem) x (8 ips) x (3 comp) --> (elem) x (8 ips) x (3 comp)
            if collocation_option == 'ip':
                EVOL = EVOL_elem.flatten(order='C')
                IVOL = IVOL_elem.flatten(order='C')
                COORDSx = COORD_elem[:,:,0].flatten(order='C')
                COORDSy = COORD_elem[:,:,1].flatten(order='C')
                COORDSz = COORD_elem[:,:,2].flatten(order='C')
                S11 = S_elem[:,:,0].flatten(order='C')
                S22 = S_elem[:,:,1].flatten(order='C')
                S33 = S_elem[:,:,2].flatten(order='C')
                S12 = S_elem[:,:,3].flatten(order='C')
                S13 = S_elem[:,:,4].flatten(order='C')
                S23 = S_elem[:,:,5].flatten(order='C')
                null = np.zeros_like(S11)
            elif collocation_option == 'center':
                EVOL = np.nansum(EVOL_elem,axis=1)
                IVOL = np.sum(IVOL_elem, axis=1)
                COORD_elem = np.mean(COORD_elem, axis=1)
                S_elem = np.mean(S_elem, axis=1)
                IVOL_elem = np.sum(IVOL_elem, axis=1)
                COORDSx = COORD_elem[:,0]
                COORDSy = COORD_elem[:,1]
                COORDSz = COORD_elem[:,2]
                S11 = S_elem[:,0]
                S22 = S_elem[:,1]
                S33 = S_elem[:,2]
                S12 = S_elem[:,3]
                S13 = S_elem[:,4]
                S23 = S_elem[:,5]
                null = np.zeros_like(S11)
            else:
                print('Specify valid collocation options')

            # unpack nodal fields
            ## 1. grab relevant nodal fields
            u_nodes = np.array(node_fields['U'][0][f])
            if velocities == True:
                v_nodes = np.array(node_fields['V'][0][f])
            if accelerations == True:
                a_nodes = np.array(node_fields['A'][0][f])
            ## 2. collocate fields
            if collocation_option == 'ip':
                u_elem = interpolate_to_ip_c3d8(u_nodes, c3d8_mesh)
                U1, U2, U3 = u_elem[:,:,0].flatten(order='C'), u_elem[:,:,1].flatten(order='C'), u_elem[:,:,2].flatten(order='C')
                if velocities:
                    v_elem = interpolate_to_ip_c3d8(v_nodes, c3d8_mesh)
                    V1, V2, V3 = v_elem[:,:,0].flatten(order='C'), v_elem[:,:,1].flatten(order='C'), v_elem[:,:,2].flatten(order='C')
                else:
                    V1, V2, V3 = null, null, null
                if accelerations:
                    a_elem = interpolate_to_ip_c3d8(a_nodes, c3d8_mesh)
                    A1, A2, A3 = a_elem[:,:,0].flatten(order='C'), a_elem[:,:,1].flatten(order='C'), a_elem[:,:,2].flatten(order='C')
                else:
                    A1, A2, A3 = null, null, null
            elif collocation_option == 'center':
                u_elem = interpolate_to_center_c3d8(u_nodes, c3d8_mesh)
                U1, U2, U3 = u_elem[:,:,0], u_elem[:,:,1], u_elem[:,:,2]
                if velocities:
                    v_elem = interpolate_to_center_c3d8(u_nodes, c3d8_mesh)
                    V1, V2, V3 = v_elem[:,:,0], v_elem[:,:,1], v_elem[:,:,2]
                else:
                    V1, V2, V3 = null, null, null
                if accelerations:
                    a_elem = interpolate_to_center_c3d8(u_nodes, c3d8_mesh)
                    A1, A2, A3 = a_elem[:,:,0], a_elem[:,:,1], a_elem[:,:,2]
                else:
                    A1, A2, A3 = null, null, null
            else:
                print('Specify valid collocation options')

            if f == 0:
                COORD_ref = COORD_elem

            # Collect results
            time = times[f]
            result = {'time':time, 'COORDx':COORDSx, 'COORDy':COORDSy, 'COORDz':COORDSz,
                      'EVOL':EVOL, 'IVOL':IVOL,
                      'S11':S11, 'S22':S22, 'S33':S33, 'S12':S12, 'S13':S13, 'S23':S23,
                      'U1': U1, 'U2': U2, 'U3':U3,
                      'V1':V1, 'V2':V2, 'V3':V3,
                      'A1':A1, 'A2':A2, 'A3':A3,
                      }
            results[f] = result

    return(results, frames, times)


def new_XDMF_writer(results, output_file, times, ref_density):
    '''Write XDMF file of collected ABaqus DNS results for Micromorphic Filter

    :param dict results: dictionary of results
    :param output_file: Name for XDMF file pair output for the Micromorphic Filter
    :param list times: Time increments of DNS
    :param float ref_density: The reference density of the material in g/cm^3 which is then converted to Mg/mm^3

    :returns: ``{output_file}.xdmf`` and ``{outptu_file}.h5``
    '''

    #data_filename = os.path.join(file_path, output_file)
    data_filename=output_file
    xdmf = file_io.xdmf.XDMF(output_filename=data_filename)

    # get the reference positions
    reference_positions = []
    for i, x in enumerate(results[0]['COORDx']):
        y = results[0]['COORDy'][i]
        z = results[0]['COORDz'][i]
        reference_positions.append([x,y,z])
    reference_positions = np.array(reference_positions)
    ndata = reference_positions.shape[0]

    # get reference volumes
    reference_volumes = np.array([vol for vol in results[0]['IVOL']])
    reference_volumes = reference_volumes.reshape((-1,1))

    point_name = 'points'
    conn_name = 'connectivity'

    # get step names
    step_names = [key for key in results.keys()]

    for j, t in enumerate(times):
        step_name = step_names[j]
        print(f'step = {step_name}')

        # initialization stuff
        grid = xdmf.addGrid(xdmf.output_timegrid, {})
        xdmf.addTime(grid, t)
        xdmf.addPoints(grid, reference_positions, duplicate=point_name)
        xdmf.addConnectivity(grid, "POLYVERTEX", np.array([v for v in range(ndata)]).reshape((-1,1)), duplicate=conn_name)

        # Get the unique positions
        unique_positions = []
        for i, x in enumerate(results[step_name]['COORDx']):
            y = results[step_name]['COORDy'][i]
            z = results[step_name]['COORDz'][i]
            unique_positions.append([x, y, z])
        unique_positions = np.array(unique_positions)
        print('unique positions', np.shape(unique_positions))
        xdmf.addData(grid, "coord", unique_positions, "Node", dtype='d')

        # get the displacement
        other_displacements = unique_positions - reference_positions
        unique_displacements = []
        for i, x, in enumerate(results[step_name]['U1']):
            y = results[step_name]['U2'][i]
            z = results[step_name]['U3'][i]
            unique_displacements.append([x, y, z])
        unique_displacements = np.array(unique_displacements)
        print(f'interpolation error = {np.mean(abs(unique_displacements - other_displacements),axis=0)}')
        print('unique displacements', np.shape(unique_displacements))
        xdmf.addData(grid, "disp", unique_displacements, "Node", dtype='d')

        # get the velocity
        unique_velocities = []
        for i, x, in enumerate(results[step_name]['V1']):
            y = results[step_name]['V2'][i]
            z = results[step_name]['V3'][i]
            unique_velocities.append([x, y, z])
        unique_velocities = np.array(unique_velocities)
        print(f"shape of unique velocities = {np.shape(unique_velocities)}")
        xdmf.addData(grid, "vel", unique_velocities, "Node", dtype='d')

        # get the acceleration
        unique_accelerations = []
        for i, x, in enumerate(results[step_name]['A1']):
            y = results[step_name]['A2'][i]
            z = results[step_name]['A3'][i]
            unique_accelerations.append([x, y, z])
        unique_accelerations = np.array(unique_accelerations)
        print(f"shape of unique accelerations = {np.shape(unique_accelerations)}")
        xdmf.addData(grid, "acc", unique_accelerations, "Node", dtype='d')

        # Stresses
        #grid_data['attributes'].update(attribute_dict)
        attribute_dict = {'S':{'name':'S'}}
        data = []
        for i, s in enumerate(results[step_name]['S11']):
            Sxx = s
            Syy = results[step_name]['S22'][i]
            Szz = results[step_name]['S33'][i]
            Syz = results[step_name]['S23'][i]
            Sxz = results[step_name]['S13'][i]
            Sxy = results[step_name]['S12'][i]
            data.append([Sxx, Sxy, Sxz,
                         Sxy, Syy, Syz,
                         Sxz, Syz, Szz])
        unique_stresses = np.array(data)
        print(f"shape of stresses = {np.shape(unique_stresses)}")
        xdmf.addData(grid, "stress", unique_stresses, "Node", dtype='d')

        # Volumes
        unique_volumes = np.array([vol for vol in results[step_name]['IVOL']])
        unique_volumes = unique_volumes.reshape((-1,1))
        print(f"shape of vol = {np.shape(unique_volumes)}")
        print(f"total volume = {np.sum(unique_volumes)}")
        xdmf.addData(grid, "volume", unique_volumes, "Node", dtype='d')

        # Density
        reference_density = ref_density * 1.e-9
        unique_densities = []
        for ref, cur in zip(reference_volumes, unique_volumes):
            J = cur / ref
            unique_densities.append(reference_density / J)
        unique_densities = np.array(unique_densities)
        unique_densities = unique_densities.reshape((-1,1))
        print(f"shape of density = {np.shape(unique_densities)}")
        xdmf.addData(grid, "density", unique_densities, "Node", dtype='d')

    xdmf.write()
    print("XDMF file written!")

    return 0


def ODBextract_to_XDMF(input_file, output_file, elem_path, node_path, mesh_path, collocation_option, ref_density, velocities=False, accelerations=False, specific_frames=None, dump_all_33_stresses=None):
    '''Convert Abaqus DNS results to XDMF format

    :param str input_file: HDF5 file of Abaqus results extracted using the
        ODB_extract module of WAVES.
    :param str output_file: Name for XDMF file pair output for the Micromorphic
        Filter.
    :param str elem_path: HDF5 path to element data
    :param str node_path: HDF5 path to node data
    :param str mesh_path: HDF5 path to mesh data
    :param str collocation option: String specifying "center" to collocate to element center or "ip" for integration points
    :param float ref_density: The reference density of the material in g/cm^3
    :param bool velocities: Boolean whether or not to collect DNS velocity data
    :param bool accelerations: Boolean whether or not to collect DNS accelerations data
    :param list specific_frames: An optional list of frame numbers for converting XDMF data
    :param str dump_all_33_stresses: Optional filename to dump all 33 stresses from DNS
    '''

    # print input argmuments and values
    print(f'input_file = {input_file}')
    print(f'output_file = {output_file}')
    print(f'collocation_option = {collocation_option}')
    if specific_frames:
        print(f'specific_frames = {specific_frames}')

    # parse field output of frames from hdf5 file
    results, frames, times = parse_input(input_file, elem_path, node_path, mesh_path, collocation_option, velocities=velocities, accelerations=accelerations, specific_frames=specific_frames)

    # output contents to XDMF file pair
    times = [results[f]['time'] for f in results.keys()]
    print(f"times = {times}")
    print(results[0].keys())
    new_XDMF_writer(results, output_file, times, ref_density)

    # Dump Cauchy 33 stresses to csv
    if dump_all_33_stresses:
        cauchy33 = results[frames[-1]]['S33']
        df = pandas.DataFrame({'quantity': cauchy33,})
        df.to_csv(dump_all_33_stresses, header=True, sep=',', index=False)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Convert Abaqus DNS results to XDMF format"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
        help='Specify the input hdf5 file generated from odb_extract')
    parser.add_argument('-o', '--output-file', type=str,
        help='Specify the output filename for the h5 + XDMF file pair')
    parser.add_argument('--elem-path', type=str,
        help='Specify the hdf5 group path to element fields')
    parser.add_argument('--node-path', type=str,
        help='Specify the hdf5 group path to nodal fields')
    parser.add_argument('--mesh-path', type=str,
        help='Specify the hdf5 group path to mesh data')
    parser.add_argument('-c', '--collocation-option', type=str, default="ip",
        help='Specify the method for collocation, either "qp" for quadrature points or "center" for element center.')
    parser.add_argument('--velocities', type=str, required=False, default="False",
        help='String specifying "True" or "False" if velocities are to be extracted')
    parser.add_argument('--accelerations', type=str, required=False, default="False",
        help='String specifying "True" or "False" if accelerations are to be extracted')
    parser.add_argument('--specific-frames', nargs="+", required=False,
        help='A list of floats corresponding to the frames to extract')
    parser.add_argument('--ref-density', type=float, required=False, default=2.00,
        help='The reference density of the material in g/cm^3')
    parser.add_argument('--dump-all-33-stresses', type=str, required=False, default=None,
        help='Optional filename to dump all 33 stresses from DNS')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(ODBextract_to_XDMF(
                input_file=args.input_file,
                output_file=args.output_file,
                elem_path=args.elem_path,
                node_path=args.node_path,
                mesh_path=args.mesh_path,
                collocation_option=args.collocation_option,
                velocities=str2bool(args.velocities),
                accelerations=str2bool(args.accelerations),
                specific_frames=args.specific_frames,
                ref_density=args.ref_density,
                dump_all_33_stresses=args.dump_all_33_stresses,
                ))
