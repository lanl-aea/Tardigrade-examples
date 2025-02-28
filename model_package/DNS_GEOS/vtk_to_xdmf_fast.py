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


def blocks_to_array(field, main_block, n, fields_dict):
    id = fields_dict[field]
    data = numpy.vstack([vtk_to_numpy(main_block.GetBlock(i).GetCellData().GetArray(id)) for i in range(n)])
    return data


def collect_and_convert_to_XDMF(vtm_file_dict, increments, output_file, dist_factor, stress_factor, density_factor):
    '''Write XDMF file of collected GEOS DNS results for Micromorphic Filter

    :param dict results: dictionary of results
    :param dict reference_positions: dictionary of reference particle positions
    :param output_file: Name for XDMF file pair output for the Micromorphic Filter
    :param float dist_factor: Optional argument to scale DNS displacements and coordinates, default=1
    :param float stress_factor: Optional argument to scale DNS stresses, default=1
    :param float density_factor: Optional factor to scale current density (if provided in the DNS results\
                                 to Mg/tonne^3, default=1

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
        main_block = multiblock_data.GetBlock(0).GetBlock(0).GetBlock(0).GetBlock(0)

        # get field names and the reference positions and
        if i == 0:
            print('collecting reference information')
            # get particle field name and id map
            number_of_blocks = main_block.GetNumberOfBlocks()
            number_of_fields = main_block.GetBlock(0).GetCellData().GetNumberOfArrays()
            fields = {}
            for i in range(number_of_fields):
                fields[main_block.GetBlock(0).GetCellData().GetArrayName(i)] = i
            print(fields)
            # reference positions
            reference_positions = dist_factor*blocks_to_array('particleReferencePosition', main_block, number_of_blocks, fields)
            ndata = reference_positions.shape[0]

        ## initialization stuff
        grid = xdmf.addGrid(xdmf.output_timegrid, {})
        xdmf.addTime(grid, t)
        xdmf.addPoints(grid, reference_positions, duplicate=point_name)
        xdmf.addConnectivity(grid, "POLYVERTEX", numpy.array([v for v in range(ndata)]).reshape((-1,1)), duplicate=conn_name)

        # get the unique displacements
        if i == 0:
            unique_positions = reference_positions
            unique_displacements = numpy.zeros_like(reference_positions)
        else:
            unique_positions = dist_factor*blocks_to_array('particleCenter', main_block, number_of_blocks, fields)
            unique_displacements = unique_positions - reference_positions
        print(f"shape of unique displacements = {numpy.shape(unique_displacements)}")
        xdmf.addData(grid, "disp", unique_displacements, "Node", dtype='d')

        # get the unique positions
        print(f"shape of unique positions = {numpy.shape(unique_positions)}")
        xdmf.addData(grid, "coord", unique_positions, "Node", dtype='d')

        # get the velocity
        if 'particleVelocity' in fields.keys():
            unique_velocities = dist_factor*blocks_to_array('particleVelocity', main_block, number_of_blocks, fields)
        else:
            unique_velocities = numpy.zeros(numpy.shape(unique_positions))
        print(f"shape of unique velocities = {numpy.shape(unique_velocities)}")
        xdmf.addData(grid, "vel", unique_velocities, "Node", dtype='d')

        # get the acceleration
        if 'particleAcceleration' in fields.keys():
            unique_accelerations = dist_factor*blocks_to_array('particleAcceleration', main_block, number_of_blocks, fields)
        else:
            unique_accelerations = numpy.zeros(numpy.shape(unique_positions))
        print(f"shape of unique accelerations = {numpy.shape(unique_accelerations)}")
        xdmf.addData(grid, "acc", unique_accelerations, "Node", dtype='d')

        # get the stresses
        stresses = blocks_to_array('particleStress', main_block, number_of_blocks, fields)
        unique_stresses = stress_factor*numpy.array([stresses[:,0], stresses[:,5], stresses[:,4], #xx, xy, xz
                                                     stresses[:,5], stresses[:,1], stresses[:,3], #yx=xy, yy, yz
                                                     stresses[:,4], stresses[:,3], stresses[:,2]])#zx=xz, zy=yz, zz         
        unique_stresses = unique_stresses.T
        print(f"shape of stresses = {numpy.shape(unique_stresses)}")
        xdmf.addData(grid, "stress", unique_stresses, "Node", dtype='d')

        # Get the volumes
        vol_factor = dist_factor*dist_factor*dist_factor
        unique_volumes = vol_factor*blocks_to_array('particleVolume', main_block, number_of_blocks, fields)
        unique_volumes = unique_volumes.reshape((-1,1))
        print(f"shape of vol = {numpy.shape(unique_volumes)}")
        print(f"total volume = {numpy.sum(unique_volumes)}")
        xdmf.addData(grid, "volume", unique_volumes, "Node", dtype='d')

        # Get the densities
        unique_densities = density_factor*blocks_to_array('particleDensity', main_block, number_of_blocks, fields)
        unique_densities = unique_densities.reshape((-1,1))
        print(f"shape of density = {numpy.shape(unique_densities)}")
        xdmf.addData(grid, "density", unique_densities, "Node", dtype='d')

    xdmf.write()
    print("XDMF file written!")

    return 0


def convert_VTK_to_XDMF(input_file, file_root, output_file, dist_factor=1, stress_factor=1, density_factor=1):
    '''Driving function to call functions for parsing GEOS VTK results and writing XDMF output

    :param str input_file: The main VTK PVD file containing GEOS DNS results
    :param str file_root: The root directory containing DNS results
    :param str output_file: The output filename for the h5 + XDMF file pair
    :param float dist_factor: Optional argument to scale DNS displacements and coordinates, default=1
    :param float stress_factor: Optional argument to scale DNS stresses, default=1
    :param float density_factor: Optional factor to scale current density (if provided in the DNS results\
                                 to Mg/tonne^3, default=1
    '''
    start_time = time.time()

    # parse vtp file
    vtm_file_dict, increments = parse_VTP_file(input_file, file_root)

    # collect and convert to XDMF
    collect_and_convert_to_XDMF(vtm_file_dict, increments, output_file, dist_factor, stress_factor, density_factor)

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
                                 ))