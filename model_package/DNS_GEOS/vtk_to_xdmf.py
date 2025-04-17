#!python
import argparse
import pathlib
import sys
import time

import os
import xml.etree.ElementTree as ET
import numpy
import matplotlib.pyplot

import file_io.xdmf


def get_data(child):
    '''Collect field data from an XML DataArray

    :param node child: The xml child node containing GEOS field data

    :returns: Array containing field values
    '''

    if (child.tag != "DataArray"):
        raise ValueError("This thing needs to be a DataArray!")
    value = None
    # get data type
    if ("Float" in child.attrib['type']):
        dtype = float
    else:
        dtype = int
    if (child.attrib['Name'] == 'ParticleFields/Density'):
        dtype = float
    dtype = float
    value = numpy.hstack([numpy.array(line.split()).astype(dtype) for line in child.text.split("\n")])#.reshape(shape)

    return value


def collect_VTK_output(input_file, file_root):
    '''Parse the GEOS DNS VTK output into a results dictionary

    :param str input_file: The main VTK PVD file containing GEOS DNS results
    :param str file_root: The root directory containing DNS results

    :returns: dictionary of results, array of reference positions
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

    # for each timesteps, find sub file paths for background grid and particle information
    for inc in file_dict.keys():
        sub_dict = {}
        file = file_dict[inc]['file']
        sub_tree = ET.parse(file)
        sub_root = sub_tree.getroot()
        sub_block = sub_root[0] # this is where "backgroundGrid" or "ParticleRegion" live
        for new_block in sub_block:
            if 'particles' in new_block.attrib['name']:
                for dataset in new_block[0][0][0]:
                    rank = dataset.attrib['name']
                    file = dataset.attrib['file']
                    sub_dict[f'particles{rank}'] = file
            file_dict[inc]['results'] = sub_dict

    # Keys to collect data
    good_particle_keys = ['particleDensity', 'particleID', 'particleStress', 'particleVelocity',
                          'particleVolume', 'particleReferencePosition']
    all_keys = ['LESample_density', 'LESample_stress', 'particleSphF', 'particleStress', 'particleBodyForce',
                'particleSPHJacobian', 'particleReferencePorosity', 'particleOverlap', 'particleDensity',
                'particlePlasticStrain', 'particleKineticEnergy', 'particleDamageGradient', 'particleInternalEnergy',
                'particleDamage', 'particleGroup', 'ghostRank', 'particleReferencePosition', 'particlePorosity',
                'particleSurfaceFlag', 'particleSurfaceNormal', 'particleArtificialViscosity', 'particleTemperature',
                'particleRank', 'particleStrengthScale', 'particleWavespeed', 'particleCenter', 'particleVelocity',
                'particleInitialSurfaceNormal', 'particleID', 'particleMass', 'particleInitialMaterialDirection',
                'particleDeformationGradient', 'particleMaterialDirection', 'particleVolume', 'particleHeatCapacity',
                'particleCrystalHealFlag', 'Points', 'connectivity', 'offsets', 'types']

    # parse particle information
    particles_dict = {}
    for inc in file_dict.keys():
        small_dict = {}
        path = 'vtkOutput'
        for results in file_dict[inc]['results']:
            if 'particles' in results:
                # combine filepaths
                sub_tree = ET.parse(f"{file_root}vtkOutput/{file_dict[inc]['results'][results]}")
                sub_root = sub_tree.getroot()
                for domain in sub_root:
                    for collection in domain:
                        for grid in collection:
                            for child in grid:
                                name = child.attrib['Name']
                                #if name in all_keys:
                                data = get_data(child)
                                if name not in small_dict.keys():
                                    small_dict[name] = []
                                small_dict[name] = numpy.hstack([small_dict[name], data])
        particles_dict[inc] = small_dict

    results = {}
    # reorder vaious arrays
    num_particles = numpy.shape(particles_dict[0]['particleDensity'])[0]
    for inc in particles_dict.keys():
        small_dict = {}
        # time
        small_dict['time'] = file_dict[inc]['timestep']
        # volume
        small_dict['volume'] = particles_dict[inc]['particleVolume']
        # density
        small_dict['density'] = particles_dict[inc]['particleDensity']
        # stresses
        stresses = particles_dict[inc]['particleStress'].reshape((num_particles,6))
        small_dict['stresses'] = stresses
        # displacements
        if inc == 0:
            reference_positions = particles_dict[0]['particleReferencePosition'].reshape((num_particles,3))
            small_dict['displacement'] = numpy.zeros_like(reference_positions)
        else:
            new_pos = particles_dict[inc]['particleCenter'].reshape((num_particles,3))
            small_dict['displacement'] = new_pos - reference_positions
        # velocity
        small_dict['velocity'] = particles_dict[inc]['particleVelocity'].reshape((num_particles,3))
        # acceleration

        # store
        results[inc] = small_dict

    return results, reference_positions


def convert_to_XDMF(results, reference_positions, output_file, dist_factor, stress_factor, density_factor):
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

    ndata = reference_positions.shape[0]

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
        xdmf.addConnectivity(grid, "POLYVERTEX", numpy.array([v for v in range(ndata)]).reshape((-1,1)), duplicate=conn_name)

        # get the unique displacements
        unique_displacements = dist_factor*results[step_name]['displacement']
        print(f"shape of unique displacements = {numpy.shape(unique_displacements)}")
        xdmf.addData(grid, "disp", unique_displacements, "Node", dtype='d')

        # get the unique positions
        unique_positions = unique_displacements + reference_positions
        print(f"shape of unique positions = {numpy.shape(unique_positions)}")
        xdmf.addData(grid, "coord", unique_positions, "Node", dtype='d')

        # get the velocity <-- fix for dynamic DNS! And apply conversion!
        unique_velocities = numpy.zeros(numpy.shape(unique_positions))
        print(f"shape of unique velocities = {numpy.shape(unique_velocities)}")
        xdmf.addData(grid, "vel", unique_velocities, "Node", dtype='d')

        # get the acceleration <-- fix for dynamic DNS! And apply conversion!
        unique_accelerations = numpy.zeros(numpy.shape(unique_positions))
        print(f"shape of unique accelerations = {numpy.shape(unique_accelerations)}")
        xdmf.addData(grid, "acc", unique_accelerations, "Node", dtype='d')

        # get the stresses
        unique_stresses = stress_factor*numpy.array([results[step_name]['stresses'][:,0],  #xx
                                                     results[step_name]['stresses'][:,5],  #xy
                                                     results[step_name]['stresses'][:,4],  #xz
                                                     results[step_name]['stresses'][:,5],  #yx=xy
                                                     results[step_name]['stresses'][:,1],  #yy
                                                     results[step_name]['stresses'][:,3],  #yz
                                                     results[step_name]['stresses'][:,4],  #zx=xz
                                                     results[step_name]['stresses'][:,3],  #zy=yz
                                                     results[step_name]['stresses'][:,2]]) #zz
        unique_stresses = unique_stresses.T
        print(f"shape of stresses = {numpy.shape(unique_stresses)}")
        xdmf.addData(grid, "stress", unique_stresses, "Node", dtype='d')

        # Get the volumes
        vol_factor = dist_factor*dist_factor*dist_factor
        unique_volumes = vol_factor*results[step_name]['volume']
        unique_volumes = unique_volumes.reshape((-1,1))
        print(f"shape of vol = {numpy.shape(unique_volumes)}")
        print(f"total volume = {numpy.sum(unique_volumes)}")
        xdmf.addData(grid, "volume", unique_volumes, "Node", dtype='d')

        # Get the densities
        unique_densities = density_factor*results[step_name]['density']
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

    # parse vtk file(s)
    results, reference_positions = collect_VTK_output(input_file, file_root)

    # convert to XDMF
    convert_to_XDMF(results, reference_positions, output_file, dist_factor, stress_factor, density_factor)

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