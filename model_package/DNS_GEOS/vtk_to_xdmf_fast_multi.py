#!python
import argparse
import pathlib
import sys
import time

import xml.etree.ElementTree as ET
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy

import file_io.xdmf


def parse_VTP_file(input_file, file_root):
    '''Parse time and subfile information from a VTP file into a dictionary

    :param str input_file: The VTP file to parse
    :param str file_root: The directory within peta_data_copy containing DNS results

    :returns: dictionary containing timestamps and vtm files, list of increment numbers
    '''

    file_dict = {}

    tree=ET.parse(f"{file_root}/{input_file}")
    root=tree.getroot()
    block=root[0]
    incs = [i for i in range(numpy.shape(block)[0])]

    # collect timestep and subfiles from pvd file
    for dataset, inc in zip(block, incs):
        sub_dict = {}
        sub_dict['timestep'] = dataset.attrib['timestep']
        sub_dict['file'] = f"{file_root}/{dataset.attrib['file']}"
        file_dict[inc] = sub_dict

    return file_dict, incs


def multi_blocks_to_array(field, particle_region_block, all_fields_dict, num_partitions=1000, stack='v', id_sort=None):
    '''Extract data for a particular field from a vtk multi block dataset

    :param str field: The name of the field
    :param vtkMultiBlockDataSet particle_region_block: The vtk dataset containing information to extract
    :param dict all_fields_dict: Dictionary for mapping the requested field name to an index in the multi block dataset
    :param int num_partitions: The number of partitions the full data field is split between corresponding to the number of ranks for the GEOS simulations
    :param str stack: Option to stack data vertically or horizontally
    :param array id_sort: Optional ID array for sorting output array by particle ID

    :returns: Array of data
    '''

    data_out = []
    for key in all_fields_dict:
        if field in all_fields_dict[key].keys():
            id = all_fields_dict[key][field]
            for j in range(num_partitions):
                data = vtk_to_numpy(particle_region_block.GetBlock(key).GetBlock(j).GetCellData().GetArray(id))
                if len(data) > 0:
                    data_out.append(data)
    if stack == 'v':
        array_out = numpy.vstack(data_out)
    elif stack == 'h':
        array_out = numpy.hstack(data_out)
    else:
        print('Invalid stacking operation requested!')

    # Sort array based on ids if provided
    if id_sort is not None:
        remap = numpy.argsort(id_sort)
        array_out = array_out[remap]

    return array_out


def create_annulus(coord, rad, x0, y0, annulus_ratio):
    '''Remove points except those within some fraction of the outer radius

    :param array-like coord: Reference coordinates of points
    :param float rad: The radius of the cylinder
    :param float x0: The x-position of the cylinder axis
    :param float y0: The y-position of the cylinder axis
    :param float annulus_ratio: Fraction of the radius of points to keep in the final geometry

    :returns: array of coordinates for points kept, mask of original array to produce final array
    '''

    print('removing points')
    print(f'original number of points = {numpy.shape(coord)[0]}')
    r_min = (1. - annulus_ratio)*rad
    r = numpy.sqrt((coord[:,0] - x0)**2 + (coord[:,1] - y0)**2)
    mask = r >= r_min
    new_coord = coord[mask]
    new_indices = numpy.where(mask)[0]
    print(mask)
    print(f'new number of points = {numpy.shape(new_coord)[0]}')

    return new_coord, new_indices


def collect_and_convert_to_XDMF(vtm_file_dict, increments, output_file,
                                dist_factor, stress_factor, density_factor, annulus_ratio,
                                upscale_damage=None, num_ranks=1000,
                                grain_particle_key=1, binder_particle_key=2):
    '''Write XDMF file of collected GEOS DNS results for Micromorphic Filter

    :param dict results: dictionary of results
    :param dict reference_positions: dictionary of reference particle positions
    :param output_file: Name for XDMF file pair output for the Micromorphic Filter
    :param float dist_factor: Optional argument to scale DNS displacements and coordinates, default=1
    :param float stress_factor: Optional argument to scale DNS stresses, default=1
    :param float density_factor: Optional factor to scale current density (if provided in the DNS results\
                                 to Mg/tonne^3, default=1
    :param float annulus_ratio: Optional fraction of the radius of points to keep in the final geometry
    :param str upscale_damage: Option to specify if damage will be upscaled
    :param int num_ranks: The number of ranks to collect data from
    :param int grain_particle_key: An integer specifying the particle key for grains
    :param int binder_particle_key: An integer specifying the particle key for binder

    :returns: ``{output_file}.xdmf`` and ``{outptu_file}.h5``
    '''

    data_filename=output_file
    xdmf = file_io.xdmf.XDMF(output_filename=data_filename)

    point_name = 'points'
    conn_name = 'connectivity'

    # get steps names
    step_names = [f'timestep_{i}' for i in increments]

    for step_name, i in zip(step_names, increments):
        print(f'step = {step_name}')

        # Grab current time
        t = vtm_file_dict[i]['timestep']

        # vtm file to unpack
        input_file = vtm_file_dict[i]['file']
        reader = vtk.vtkXMLMultiBlockDataReader()
        reader.SetFileName(input_file)
        reader.Update()
        multiblock_data = reader.GetOutput()
        particle_region_block = multiblock_data.GetBlock(0).GetBlock(0).GetBlock(0)

        # get field names and the reference positions
        if i == 0:
            print('collecting reference information')
            # get particle field name and id map
            n_p_blocks = particle_region_block.GetNumberOfBlocks()
            all_fields = {}
            for b in range(n_p_blocks):
                particle_region_block
                n_p_fields = particle_region_block.GetBlock(b).GetBlock(0).GetCellData().GetNumberOfArrays()
                fields = {}
                for j in range(n_p_fields):
                    fields[particle_region_block.GetBlock(b).GetBlock(0).GetCellData().GetArrayName(j)] = j
                print(fields)
                all_fields[b] = fields

            # Parse the relevant particle fields for the required keys
            grain_fields = {grain_particle_key: all_fields[grain_particle_key]}
            binder_fields = {binder_particle_key: all_fields[binder_particle_key]}
            particle_fields = grain_fields | binder_fields
            particle_keys = list(grain_fields.keys()) + list(binder_fields.keys())

            # Ids for sorting
            idox_ids = multi_blocks_to_array(
                'particleID', particle_region_block,
                grain_fields, num_ranks, 'h')
            binder_ids = multi_blocks_to_array(
                'particleID', particle_region_block,
                binder_fields, num_ranks, 'h')

            # ref pos
            idox_reference_positions = dist_factor*multi_blocks_to_array(
                'particleReferencePosition', particle_region_block,
                grain_fields, num_ranks, id_sort=idox_ids)
            binder_reference_positions = dist_factor*multi_blocks_to_array(
                'particleReferencePosition', particle_region_block,
                binder_fields, num_ranks, id_sort=binder_ids)
            reference_positions = numpy.vstack([idox_reference_positions, binder_reference_positions])
            ndata = numpy.shape(reference_positions)[0]

            # process reference positions
            if annulus_ratio is not None:
                xmin, xmax = numpy.min(reference_positions[:,0]), numpy.max(reference_positions[:,0])
                ymin, ymax = numpy.min(reference_positions[:,1]), numpy.max(reference_positions[:,1])
                x0, y0 = 0.5*(xmax - xmin) + xmin, 0.5*(ymax - ymin) + ymin
                rad = numpy.mean([0.5*(xmax - xmin), 0.5*(ymax - ymin)])
                reference_positions, mask = create_annulus(reference_positions, rad, x0, y0, annulus_ratio)
        else:
            # Get IDs for sorting
            idox_ids = multi_blocks_to_array(
                'particleID', particle_region_block,
                grain_fields, num_ranks, 'h')
            binder_ids = multi_blocks_to_array(
                'particleID', particle_region_block,
                binder_fields, num_ranks, 'h')

        ## initialization stuff
        grid = xdmf.addGrid(xdmf.output_timegrid, {})
        xdmf.addTime(grid, t)

        # get the unique displacements
        if i == 0:
            unique_positions = reference_positions
            unique_displacements = numpy.zeros_like(reference_positions)
        else:
            # Get reference positions again since particles may have been reordered
            idox_reference_positions = dist_factor*multi_blocks_to_array(
                'particleReferencePosition', particle_region_block,
                grain_fields, num_ranks, id_sort=idox_ids)
            binder_reference_positions = dist_factor*multi_blocks_to_array(
                'particleReferencePosition', particle_region_block,
                binder_fields, num_ranks, id_sort=binder_ids)
            reference_positions = numpy.vstack([idox_reference_positions, binder_reference_positions])
            # Current positions
            idox_positions = dist_factor*multi_blocks_to_array(
                'particleCenter', particle_region_block,
                grain_fields, num_ranks, id_sort=idox_ids)
            binder_positions = dist_factor*multi_blocks_to_array(
                'particleCenter', particle_region_block,
                binder_fields, num_ranks, id_sort=binder_ids)
            unique_positions = numpy.vstack([idox_positions, binder_positions])
            if annulus_ratio is not None:
                unique_positions = unique_positions[mask]
            unique_displacements = unique_positions - reference_positions
        print(f"shape of unique displacements = {numpy.shape(unique_displacements)}")
        # Now add reference positions
        #xdmf.addPoints(grid, reference_positions, duplicate=point_name)
        xdmf.addPoints(grid, reference_positions)
        xdmf.addConnectivity(grid, "POLYVERTEX", numpy.array([v for v in range(numpy.shape(reference_positions)[0])]).reshape((-1,1)))
        xdmf.addData(grid, "disp", unique_displacements, "Node", dtype='d')

        # get the unique positions
        print(f"shape of unique positions = {numpy.shape(unique_positions)}")
        xdmf.addData(grid, "coord", unique_positions, "Node", dtype='d')

        # get the velocity
        if 'particleVelocity' in particle_keys:
            idox_velocities = dist_factor*multi_blocks_to_array(
                'particleVelocity', particle_region_block,
                grain_fields, num_ranks, id_sort=idox_ids)
            binder_velocities = dist_factor*multi_blocks_to_array(
                'particleVelocity', particle_region_block,
                binder_fields, num_ranks, id_sort=binder_ids)
            unique_velocities = numpy.vstack([idox_velocities, binder_velocities])
            if annulus_ratio is not None:
                unique_velocities = unique_velocities[mask]
        else:
            unique_velocities = numpy.zeros(numpy.shape(unique_positions))
        print(f"shape of unique velocities = {numpy.shape(unique_velocities)}")
        xdmf.addData(grid, "vel", unique_velocities, "Node", dtype='d')

        # get the acceleration, TODO: replace 'particleVelocity' with 'particleAcceleration' Once accelerations have been fixed
        if 'particleAcceleration' in particle_keys:
            idox_accelerations = dist_factor*multi_blocks_to_array(
                'particleVelocity', particle_region_block,
                grain_fields, num_ranks, id_sort=idox_ids)
            binder_accelerations = dist_factor*multi_blocks_to_array(
                'particleVelocity', particle_region_block,
                binder_fields, num_ranks, id_sort=binder_ids)
            unique_accelerations = numpy.vstack([idox_accelerations, binder_accelerations])
            if annulus_ratio is not None:
                unique_accelerations = unique_accelerations[mask]
        else:
            unique_accelerations = numpy.zeros(numpy.shape(unique_positions))
        print(f"shape of unique accelerations = {numpy.shape(unique_accelerations)}")
        xdmf.addData(grid, "acc", unique_accelerations, "Node", dtype='d')

        # get the stresses
        idox_stresses = stress_factor*multi_blocks_to_array(
            'Idox_stress', particle_region_block,
            grain_fields, num_ranks, id_sort=idox_ids)
        binder_stresses = stress_factor*multi_blocks_to_array(
            'EstaneMatrix_stress', particle_region_block,
            binder_fields, num_ranks, id_sort=binder_ids)
        stresses = numpy.vstack([idox_stresses, binder_stresses])
        unique_stresses = numpy.array([stresses[:,0], stresses[:,5], stresses[:,4], #xx, xy, xz
                                    stresses[:,5], stresses[:,1], stresses[:,3], #yx=xy, yy, yz
                                    stresses[:,4], stresses[:,3], stresses[:,2]])#zx=xz, zy=yz, zz         
        unique_stresses = unique_stresses.T
        if annulus_ratio is not None:
            unique_stresses = unique_stresses[mask]
        print(f"shape of stresses = {numpy.shape(unique_stresses)}")
        xdmf.addData(grid, "stress", unique_stresses, "Node", dtype='d')

        #Get the volumes --> update when volumes are added!
        vol_factor = dist_factor*dist_factor*dist_factor
        try:
            idox_volumes = vol_factor*multi_blocks_to_array(
                'particleVolume', particle_region_block,
                grain_fields, num_ranks, id_sort=idox_ids)
            binder_volumes = vol_factor*multi_blocks_to_array(
                'particleVolume', particle_region_block,
                binder_fields, num_ranks, id_sort=binder_ids)
            unique_volumes = numpy.vstack([idox_volumes, binder_volumes])
        except:
            total_volume = 79.44851544702128 #mm^2
            vol = total_volume / ndata
            unique_volumes = vol*numpy.ones(ndata)
        unique_volumes = unique_volumes.reshape((-1,1))
        print(f"total volume = {numpy.sum(unique_volumes)}")
        if annulus_ratio is not None:
            unique_volumes = unique_volumes[mask]
            print(f"reduced volume = {numpy.sum(unique_volumes)}")
        print(f"shape of vol = {numpy.shape(unique_volumes)}")
        xdmf.addData(grid, "volume", unique_volumes, "Node", dtype='d')

        # Get the densities
        idox_densities = density_factor*multi_blocks_to_array(
            'Idox_density', particle_region_block,
            grain_fields, num_ranks, stack='h', id_sort=idox_ids)
        binder_densities = density_factor*multi_blocks_to_array(
            'EstaneMatrix_density', particle_region_block,
            binder_fields, num_ranks, stack='h', id_sort=binder_ids)
        unique_densities = numpy.hstack([idox_densities, binder_densities])
        unique_densities = unique_densities.reshape((-1,1))
        if annulus_ratio is not None:
            unique_densities = unique_densities[mask]
        print(f"shape of density = {numpy.shape(unique_densities)}")
        xdmf.addData(grid, "density", unique_densities, "Node", dtype='d')

        # Get the damages
        if upscale_damage is not None:
            idox_damage = multi_blocks_to_array(
                'Idox_damage', particle_region_block,
                grain_fields, num_ranks, stack='h', id_sort=idox_ids)
            binder_damage = multi_blocks_to_array(
                'EstaneMatrix_damage', particle_region_block,
                binder_fields, num_ranks, stack='h', id_sort=binder_ids)
            unique_damage = numpy.hstack([idox_damage, binder_damage])
            unique_damage = unique_damage.reshape((-1,1))
            print(f"shape of damage = {numpy.shape(unique_damage)}")
            xdmf.addData(grid, "damage", unique_damage, "Node", dtype='d')

    xdmf.write()
    print("XDMF file written!")

    return 0


def convert_VTK_to_XDMF(input_file, file_root, output_file, dist_factor=1, stress_factor=1, density_factor=1,
                        annulus_ratio=None, upscale_damage=None, num_ranks=1000,
                        grain_particle_key=1, binder_particle_key=2):
    '''Driving function to call functions for parsing GEOS VTK results and writing XDMF output

    :param str input_file: The main VTK PVD file containing GEOS DNS results
    :param str file_root: The root directory containing DNS results
    :param str output_file: The output filename for the h5 + XDMF file pair
    :param float dist_factor: Optional argument to scale DNS displacements and coordinates, default=1
    :param float stress_factor: Optional argument to scale DNS stresses, default=1
    :param float density_factor: Optional factor to scale current density (if provided in the DNS results\
                                 to Mg/tonne^3, default=1
    :param float annulus_ratio: Optional fraction of the radius of points to keep in the final geometry
    :param str upscale_damage: Option to specify if damage will be upscaled
    :param int num_ranks: The number of ranks to collect data from
    :param int grain_particle_key: An integer specifying the particle key for grains
    :param int binder_particle_key: An integer specifying the particle key for binder
    '''

    start_time = time.time()

    # parse vtp file
    vtm_file_dict, increments = parse_VTP_file(input_file, file_root)

    # collect and convert to XDMF
    collect_and_convert_to_XDMF(vtm_file_dict, increments, output_file,
                                dist_factor, stress_factor, density_factor,
                                annulus_ratio, upscale_damage, num_ranks,
                                grain_particle_key, binder_particle_key)

    end_time = time.time()
    print(f'executed in {end_time - start_time} seconds')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Convert GEOS DNS results to XDMF format"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-i', '--input-file', type=str, required=True,
        help='Specify the main VTK PVD file containing GEOS DNS results')
    parser.add_argument('--file-root', type=str, required=True,
        help='The root directory containing DNS results')
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='Specify the output filename for the h5 + XDMF file pair')
    parser.add_argument('--dist-factor', type=float, required=False, default=1,
        help='Optional argument to scale DNS displacements and coordinates')
    parser.add_argument('--stress-factor', type=float, required=False, default=1,
        help='Optional argument to scale DNS stresses')
    parser.add_argument('--density-factor', type=float, required=False, default=1,
         help='Optional factor to scale current density (if provided in the DNS results\
               to Mg/tonne^3')
    parser.add_argument('--annulus-ratio', type=float, required=False, default=None,
         help='Optional fraction of the radius of points to keep in the final geometry')
    parser.add_argument('--upscale-damage', type=str, required=False, default=None,
         help='Option to specify if damage will be upscaled')
    parser.add_argument('--num-ranks', type=int, required=False, default=1000,
         help='The number of ranks to collect data from')
    parser.add_argument('--grain-particle-key', type=int, required=False, default=1,
         help='An integer specifying the particle key for grains')
    parser.add_argument('--binder-particle-key', type=int, required=False, default=2,
         help='An integer specifying the particle key for binder')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(convert_VTK_to_XDMF(input_file=args.input_file,
                                 file_root=args.file_root,
                                 output_file=args.output_file,
                                 dist_factor=args.dist_factor,
                                 stress_factor=args.stress_factor,
                                 density_factor=args.density_factor,
                                 annulus_ratio=args.annulus_ratio,
                                 upscale_damage=args.upscale_damage,
                                 num_ranks=args.num_ranks,
                                 grain_particle_key=args.grain_particle_key,
                                 binder_particle_key=args.binder_particle_key,
                                 ))