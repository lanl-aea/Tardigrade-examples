# Imports
import sys
import os
import inspect
import argparse
import yaml

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


def interpolate_to_center(node_array, mesh):
    '''Average a vector or tensor field from the nodes to the center of an element

    :param array-like node_array: nodal data to be interpolated
    :param array-like mesh: the element connectivity for all elements

    :returns: dictionary of interpolated results
    '''

    numpts = 1
    numelem = np.shape(mesh)[0]
    results = np.zeros([numelem, 2])

    for e, n in enumerate(mesh):
        node_field = []
        for node in n:
            # get the field values for that node
            node_field.append(node_array[node-1])

        #node_field = np.array(node_field).flatten(order='C')
        node_field = np.array(node_field)

        # average over node_field
        solve = np.mean(node_field, axis=0)
        results[e, :] = solve[0:2]

    print(f'nodal_field shape = {np.shape(results)}')
    return(results)


def parse_input(input_file, elem_path, node_path, mesh_path, collocation_option, element_type,
                velocities=False, accelerations=False, specific_frames=None):
    '''Parse the HDF5 file output by ODBextract (WAVES tool)

    :param str input_file: HDF5 file of Abaqus results extracted using the
        ODB_extract module of WAVES.
    :param str elem_path: HDF5 path to element data
    :param str node_path: HDF5 path to node data
    :param str mesh_path: HDF5 path to mesh data
    :param str collocation option: String specifying "center" to collocate to element center or "ip" for integration points
    :param str element_type: Abaqus element type
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
        mesh = np.array(file[mesh_path][f'{element_type}_mesh'])

        # loop over frames and get indices for each field
        results = {}
        for f in frames:

            # unpack element field results
            IVOL_elem   = np.array(elem_fields['IVOL'][0][f])
            EVOL_elem   = np.array(elem_fields['EVOL'][0][f])
            S_elem      = np.array(elem_fields['S'][0][f])
            COORD_elem  = np.array(elem_fields['COORD'][0][f])
            if collocation_option == 'ip':
                EVOL = EVOL_elem.flatten(order='C')
                IVOL = IVOL_elem.flatten(order='C')
                COORDSx = COORD_elem[:,:,0].flatten(order='C')
                COORDSy = COORD_elem[:,:,1].flatten(order='C')
                S11 = S_elem[:,:,0].flatten(order='C')
                S22 = S_elem[:,:,1].flatten(order='C')
                S33 = S_elem[:,:,2].flatten(order='C')
                S12 = S_elem[:,:,3].flatten(order='C')
                null = np.zeros_like(S11)
            elif collocation_option == 'center':
                EVOL = np.nansum(EVOL_elem,axis=1)
                IVOL = np.nansum(IVOL_elem, axis=1)
                COORD_elem = np.nanmean(COORD_elem, axis=1)
                S_elem = np.nanmean(S_elem, axis=1)
                IVOL_elem = np.nansum(IVOL_elem, axis=1)
                COORDSx = COORD_elem[:,0]
                COORDSy = COORD_elem[:,1]
                S11 = S_elem[:,0]
                S22 = S_elem[:,1]
                S33 = S_elem[:,2]
                S12 = S_elem[:,3]
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
            if collocation_option == 'center':
                u_elem = interpolate_to_center(u_nodes, mesh)
                U1, U2 = u_elem[:,0], u_elem[:,1]
                if velocities:
                    v_elem = interpolate_to_center(u_nodes, mesh)
                    V1, V2 = v_elem[:,0], v_elem[:,1]
                else:
                    V1, V2 = null, null
                if accelerations:
                    a_elem = interpolate_to_center(u_nodes, mesh)
                    A1, A2 = a_elem[:,0], a_elem[:,1]
                else:
                    A1, A2 = null, null
            else:
                print('Specify valid collocation options')

            if f == 0:
                COORD_ref = COORD_elem

            # Collect results
            time = times[f]
            result = {'time':time, 'COORDx':COORDSx, 'COORDy':COORDSy,
                      'EVOL':EVOL, 'IVOL':IVOL,
                      'S11':S11, 'S22':S22, 'S33':S33, 'S12':S12,
                      'U1': U1, 'U2': U2,
                      'V1':V1, 'V2':V2,
                      'A1':A1, 'A2':A2,
                      }
            results[f] = result

    return(results, frames, times)


def c3d4_volume(node_location, mesh):
    '''Calculate volume of tetrahedral elements in a mesh
    https://math.stackexchange.com/questions/3616760/how-to-calculate-the-volume-of-tetrahedron-given-by-4-points
    http://tamivox.org/redbear/tetra_calc/index.html

    :param array-like node_location: location of nodes in reference state`
    :param array-like mesh: the element connectivity for all elements

    :returns: array of calculated tetrahedral volumes
    '''

    output = []
    for element in mesh:
        n1, n2, n3, n4 = node_location[element[0]-1], node_location[element[1]-1], node_location[element[2]-1], node_location[element[3]-1]
        matrix = np.array([[n1[0], n1[1], n1[2], 1],
                           [n2[0], n2[1], n2[2], 1],
                           [n3[0], n3[1], n3[2], 1],
                           [n4[0], n4[1], n4[2], 1]])
        volume = abs((1/6)*np.linalg.det(matrix))
        output.append(volume)
    return np.array(output)


def collect_tri_areas(points, connectivity):
    '''Calculate the area of all triangular elements/cells

    :params array points: nodal coordinates
    :params list connectivity: list of arrays containing element connectivity

    :returns: array of triangular areas
    '''

    def tri_area(nodes):
        r'''
        A = \frac{1}{2} det \left(
            \begin{bmatrix}
                 1 &  1 & 1\\
                x1 & x2 & x3\\
                y1 & y2 & y3
            \end{bmatrix} \right)
        '''
        # Add one to the bottom row for the determinant calculation
        nodes[:,-1] += 1
        return np.linalg.det(nodes.T)

    volumes = np.array([tri_area(points[i-1]) for i in connectivity])

    return volumes


def initialize_reference(old_results, old_times, input_file, mesh_path, collocation_option):
    '''Manually initialize reference increment

    :param dict old_results: dictionary of results without reference state
    :param list old_times: Time increments of DNS without reference state
    :param str input_file: HDF5 file of Abaqus results extracted using the
        ODB_extract module of WAVES.
    :param str mesh_path: HDF5 path to mesh data
    :param str collocation_option: String specifying "center" to collocate to element center or "ip" for integration points
    :param str element_type: Abaqus element type

    :returns: New ``results`` dictinoary and ``times`` array including reference increment
    '''

    #times = old_times.insert(0, 0.0)
    times = np.insert(old_times, 0, 0.0)
    results = {}
    results[0] = {}

    # Copy old_results to results for all frames after 0
    for key in old_results.keys():
        results[key+1] = old_results[key].copy()

    # Calculate reference coordinates and volume from mesh
    with h5py.File(input_file, 'r') as file:
        node = np.array(file[mesh_path]['node'])
        node_location = np.array(file[mesh_path]['node_location'])
        mesh = np.array(file[mesh_path][f'CPS3_mesh'])
        if collocation_option == 'ip':
            print(f'WARNING! Interpolation for CPS3 has not been implemented!')
        elif collocation_option == 'center':
            COORD = interpolate_to_center(node_location, mesh)
            volume = collect_tri_areas(node_location, mesh)
        else:
            print('Specify valid collocation options')
        results[0]['COORDx'] = COORD[:,0]
        results[0]['COORDy'] = COORD[:,1]
        results[0]['IVOL'] = volume
        
    # Setup initial frame for all other fields
    results[0]['time'] = 0.0
    zero_fields = ['S11', 'S22', 'S33', 'S12', 'U1', 'U2', 'V1', 'V2', 'A1', 'A2']
    for field in zero_fields:
        results[0][field] = np.zeros_like(old_results[0][field])

    return results, times


def new_XDMF_writer(results, output_file, times, ref_density, sets=None):
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
    ## TODO: init_ref is currently implemented because ODBExtract does not play well with multi-step Abaqus jobs and a reference stat is needed. The COORDs are extracted from frame data not in the current reference configuration. Find out how big of a deal this will be --> mainly for the pressurized neper grains simulations
    reference_positions = []
    for i, x in enumerate(results[0]['COORDx']):
        y = results[0]['COORDy'][i]
        reference_positions.append([x,y])
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
            unique_positions.append([x, y])
        unique_positions = np.array(unique_positions)
        print('unique positions', np.shape(unique_positions))
        xdmf.addData(grid, "coord", unique_positions, "Node", dtype='d')

        # get the displacement
        #other_displacements = unique_positions - reference_positions
        unique_displacements = []
        for i, x, in enumerate(results[step_name]['U1']):
            y = results[step_name]['U2'][i]
            unique_displacements.append([x, y])
        unique_displacements = np.array(unique_displacements)
        #print(f'interpolation error = {np.mean(abs(unique_displacements - other_displacements),axis=0)}')
        print('unique displacements', np.shape(unique_displacements))
        xdmf.addData(grid, "disp", unique_displacements, "Node", dtype='d')

        # get the velocity
        unique_velocities = []
        for i, x, in enumerate(results[step_name]['V1']):
            y = results[step_name]['V2'][i]
            unique_velocities.append([x, y])
        unique_velocities = np.array(unique_velocities)
        print(f"shape of unique velocities = {np.shape(unique_velocities)}")
        xdmf.addData(grid, "vel", unique_velocities, "Node", dtype='d')

        # get the acceleration
        unique_accelerations = []
        for i, x, in enumerate(results[step_name]['A1']):
            y = results[step_name]['A2'][i]
            unique_accelerations.append([x, y])
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
            Sxy = results[step_name]['S12'][i]
            data.append([Sxx, Sxy,
                         Sxy, Syy])
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

        # Sets
        if (step_name == step_names[0]) and sets:
            print("Adding element sets to XDMF file!")
            nodeset_mask = np.zeros((ndata, 1), dtype=int)
            for i, set in enumerate(sets.keys()):
                #xdmf.addDomain(grid, set, np.array(sets[set]))
                set_ids = np.array([int(i)-1 for i in sets[set]])
                print(f"shape of set{set} = {np.shape(set_ids.reshape((-1,1)))}")
                # Write the field of points to a Point Array for the Filter to interpret
                xdmf.addData(grid, set, set_ids.reshape((-1,1)), "Node", dtype='i')
                # Write the nodeset for visualization in Paraview
                nodeset_mask[set_ids] = i + 1
                print(f"shape of  nodeset {set} = {np.shape(nodeset_mask)}")
            xdmf.addData(grid, f"grains", nodeset_mask, "Node", dtype='i')

    xdmf.write()
    print("XDMF file written!")

    return 0


def ODBextract_to_XDMF(input_file, output_file, elem_path, node_path, mesh_path, collocation_option, ref_density, element_type,
                       velocities=False, accelerations=False, specific_frames=None, init_ref=None,
                       sets_file=None, num_steps=None):
    '''Convert Abaqus DNS results of 32D bonded grains to XDMF format

    :param str input_file: HDF5 file of Abaqus results extracted using the
        ODB_extract module of WAVES.
    :param str output_file: Name for XDMF file pair output for the Micromorphic
        Filter.
    :param str elem_path: HDF5 path to element data
    :param str node_path: HDF5 path to node data
    :param str mesh_path: HDF5 path to mesh data
    :param str collocation option: String specifying "center" to collocate to element center or "ip" for integration points
    :param str element_type: Abaqus element type
    :param float ref_density: The reference density of the material in g/cm^3
    :param bool velocities: Boolean whether or not to collect DNS velocity data
    :param bool accelerations: Boolean whether or not to collect DNS accelerations data
    :param list specific_frames: An optional list of frame numbers for converting XDMF data
    '''

    # print input arguments and values
    print(f'input_file = {input_file}')
    print(f'output_file = {output_file}')
    print(f'collocation_option = {collocation_option}')
    if specific_frames:
        print(f'specific_frames = {specific_frames}')

    # parse field output of frames from hdf5 file
    results, frames, times = parse_input(input_file, elem_path, node_path, mesh_path, collocation_option, element_type, velocities=velocities, accelerations=accelerations, specific_frames=specific_frames)

    # initialize reference frame if requested
    if init_ref:
        results, times = initialize_reference(results, times, input_file, mesh_path, collocation_option)

    # Open sets dictionary if requested
    sets = None
    if sets_file:
        stream = open(sets_file, 'r')
        sets = yaml.load(stream, Loader=yaml.FullLoader)
        stream.close()

    # output contents to XDMF file pair
    times = [results[f]['time'] for f in results.keys()]
    if num_steps is not None:
        time_0 = times[0]
        time_last = times[-1]
        total_num_times = len(times)
        if total_num_times <= num_steps:
            print(f'Total number of times in simulation is {total_num_times} and less than the\
                    number requested of num_steps={num_steps}. Using all available time steps.')
            times = times
        else:
            print(f'Adjusting times based on num_steps={num_steps}')
            if num_steps > 1:
                times_to_insert = []
                for i in range(1, num_steps):
                    new_time = times[int(np.floor(i*(total_num_times/num_steps)))]
                    times_to_insert.append(float(new_time))
                times = [time_0] + times_to_insert + [time_last]
            else:
                times = [time_0, time_last]
    print(f"times = {times}")
    print(results[0].keys())
    new_XDMF_writer(results, output_file, times, ref_density, sets)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Convert Abaqus DNS results of 2D bonded grains to XDMF format"
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
    parser.add_argument('--element-type', type=str, required=False, default='C3D8',
        help='Abaqus element type')
    parser.add_argument('--init-ref', type=str, required=False, default=None,
        help='A flag (any string) to specify if the reference configuration will be initialized manually')
    parser.add_argument('--sets-file', type=str, required=False, default=None,
        help='A yaml file containing element set information')
    parser.add_argument('--num-steps', type=int, required=False, default=None,
        help='Option to specify how many total timesteps should be written to the XDMF file excluding\
              the reference state. For 1, the final state is used. For num_steps > 1, the final state\
              is written and nearest evenly steps are written.')

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
                element_type=args.element_type,
                init_ref=args.init_ref,
                sets_file=args.sets_file,
                num_steps=args.num_steps,
                ))
