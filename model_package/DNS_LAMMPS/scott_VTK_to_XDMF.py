# Imports
from pdb import set_trace as dbg
import os.path as osp
import sys
import collections
import re
import os
import xml
import lxml.etree
import xml.etree.ElementTree as ET
import h5py
import inspect
import numpy as np
import argparse
import matplotlib.pyplot as plt
import xdmf_builder_tools as XBT
from itertools import chain

    
def collect_output(path, base_files, timesteps):
    
    # parse contents only for selected timesteps
    results_dict = {}
    frame_number = 0
    
    # map time from 0 to 1 over all timesteps
    num_timesteps = len(timesteps)
    times = np.linspace(0, 1, num_timesteps)
    
    for timestep, time in zip(timesteps, times):
        small_dict = {}
        small_dict['time'] = time
        small_dict['frame'] = frame_number
        files = [f'{path}{base}_{timestep}.vtk' for base in base_files]
        
        for file, base in zip(files, base_files):
            with open(file, 'r') as input:
                all_input=input.readlines()
                data_only=all_input[11:]
            
            # convert to floats
            data = np.array([float(item.split('\n')[0]) for item in data_only])
            
            # stack into small_dict
            small_dict[base] = data
            
        # stack into main results
        results_dict[f'timestep_{timestep}'] = small_dict
        frame_number = frame_number + 1    

    return(results_dict)


def generate_sets(results_dict, timesteps):
    set_dict = {}
        
    # get all
    set_dict['ALL'] = get_sets(results_dict, timesteps)
    
    return(set_dict)


def get_sets(results_dict, timesteps):

    step_name = f'timestep_{timesteps[0]}'
    set_ids = []
    for i, x in enumerate(results_dict[step_name]['xGrid']):
        set_ids.append(i)
        
    return(set_ids)


def convert_to_XDMF(results, output_file, set_dict):

    # create the h5 output
    h5output = h5py.File(output_file + '.h5', 'w')

    # initialize the XML structure
    root, domain, collection = XBT.initialize_xdmf_domain()

    # loop through frames
    frames = [results[key]['frame'] for key in results.keys()]
    frame_names = XBT.sort_string_names(frames)
    print(f"frame names = {frame_names}")
    step_names = list(results.keys())
    print(f"step_names = {step_names}")



    initial_grid = True
    for fn in frame_names:
        grid_data = {'initial_grid':initial_grid}
        initial_grid = False
        
        frame = int(fn)
        step_name = step_names[frame]
        print(f'step = {step_name}')
        
        grid_data.update({'frame':frame})
        grid_data.update({'time':results[step_name]['time']})
        print(f"time = {results[step_name]['time']}")

        # get the unique positions
        unique_positions = []
        factor = 1/1000 # from um to mm
        #if grid_data['initial_grid']:
        for i, x in enumerate(results[step_name]['xGrid']):
            y = results[step_name]['yGrid'][i]
            z = results[step_name]['zGrid'][i]
            unique_positions.append([factor*x,factor*y,factor*z])
        unique_positions = np.array(unique_positions)

        if grid_data['initial_grid']:
            reference_positions = np.copy(unique_positions)

        # get the unique displacements
        unique_displacements = []
        factor = 1/1000 # from um to mm
        for i, x in enumerate(results[step_name]['dxGrid']):
            y = results[step_name]['dxGrid'][i]
            z = results[step_name]['dxGrid'][i]
            unique_displacements.append([factor*x,factor*y,factor*z])
        unique_displacements = np.array(unique_displacements)

        # Add U, V, A, and COORD to grid_data
        print(f"shape of unique positions = {np.shape(unique_positions)}")
        print(f"shape of unique displacements = {np.shape(unique_displacements)}")

        # get the velocity <-- fix for dynamic DNS!
        unique_velocities = np.zeros(np.shape(unique_positions))
        print(f"shape of unique velocities = {np.shape(unique_velocities)}")
        
        # get the acceleration <-- fix for dynamic DNS!
        unique_accelerations = np.zeros(np.shape(unique_positions))
        print(f"shape of unique accelerations = {np.shape(unique_accelerations)}")
        
        pt_ids = np.array([v for v in range(np.shape(unique_positions)[0])])
        grid_data.update({'ids':pt_ids})
        print(f"shape of pt_ids = {np.shape(pt_ids)}")


        grid_data.update({'coordinates':{'values':reference_positions,
                                        'name':'coordinates_0'}})

        #print(grid_data['coordinates']['values'])
        
        grid_data.update({'attributes':{}})
        
        grid_data['attributes'].update({'displacements':{'name':'U',
                                                         'values':unique_displacements,
                                                         'ordering':('x', 'y', 'z')}})

        grid_data['attributes'].update({'velocities':{'name':'V',
                                                      'values':unique_velocities,
                                                      'ordering':('x', 'y', 'z')}})

        grid_data['attributes'].update({'accelerations':{'name':'A',
                                                         'values':unique_accelerations,
                                                         'ordering':('x', 'y', 'z')}})
                                                         
        grid_data['attributes'].update({'COORD':{'name':'COORD',
                                                 'values':unique_positions,
                                                 'ordering':('x', 'y', 'z')}})

        # EVOL
        ## skip for now!

        # IVOL
        factor = 1e-9 # from um^3 to mm^3
        attribute_dict = {'IVOL':{'name':'IVOL'}}
        data = np.array([ivol*factor for ivol in results[step_name]['vol']])
        attribute_dict['IVOL'].update({'values':data})
        attribute_dict['IVOL'].update({'ordering':['']})
        grid_data['attributes'].update(attribute_dict)
        print(f"shape of IVOL = {np.shape(data)}")

        # S
        #grid_data['attributes'].update(attribute_dict)
        attribute_dict = {'S':{'name':'S'}}
        data = []
        for i, s in enumerate(results[step_name]['sig11Grid']):
            factor = 1/1000 # from kPa to MPa
            Sxx = s*factor
            Syy = results[step_name]['sig22Grid'][i]*factor
            Szz = results[step_name]['sig33Grid'][i]*factor
            Syz = results[step_name]['sig23Grid'][i]*factor
            Sxz = results[step_name]['sig13Grid'][i]*factor
            Sxy = results[step_name]['sig12Grid'][i]*factor
            data.append([Sxx, Syy, Szz, Syz, Sxz, Sxy])
        data2 = np.array(data)
        attribute_dict['S'].update({'values':data2})
        attribute_dict['S'].update({'ordering':('xx', 'yy', 'zz', 'yz', 'xz', 'xy')})
        grid_data['attributes'].update(attribute_dict)
        print(f"shape of stresses = {np.shape(data2)}")

        # assign the topology
        grid_data.update({'topology':{'name':'topology_0'}})

        # Detect the sets
        num_pts = len(results[step_names[0]]['xGrid'])
        if frame == 0:
            set_configurations = []
            for set_name in set_dict.keys():
                set_configurations.append({'name':set_name, 'ids':set_dict[set_name]})
            

        grid_data.update({'sets':set_configurations})

        # Density
        factor = 1e-12 # from kg/m^3 to Mg/mm^3
        attribute_dict = {"DENSITY":{'name':'DENSITY'}}
        data = np.array([density*factor for density in results[step_name]['rho']])
        attribute_dict["DENSITY"].update({'values':data})
        attribute_dict["DENSITY"].update({'ordering':['']})
        grid_data['attributes'].update(attribute_dict)
        print(f"shape of density = {np.shape(data)}")
        print(f"min, max, and average densities = {np.min(data)}, {np.max(data)}, {np.average(data)}")

        # build the XML grid
        XBT.construct_grid(collection, grid_data, h5output)

    # Dump the file
    with open(output_file + '.xdmf', 'w') as xml_file:
        xml_file.write(ET.tostring(root, encoding='utf8').decode('utf-8'))
        
    # Pretty-print using lxml
    pretty_xml = lxml.etree.tostring(lxml.etree.parse(output_file + ".xdmf").getroot(),
        pretty_print=True, encoding="unicode")
    with open(output_file + '.xdmf', 'w') as xml_file:
        xml_file.write(pretty_xml)
        
    return 0 


#====================================================================== MAIN ===
    
def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Prep Abaqus DNS for filter"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('--path', type=str, required=True,
        help='Root path to VTK version of LAMMPS results')
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help='Specify the output filename for the h5 + XDMF file pair')
    parser.add_argument('-t', '--timesteps', nargs="+", required=True,
        help='Specify a list of timesteps to extract data from')
    # parser.add_argument('-b', '--block-name', type=str,
    #     help='Specify the part instance "block name" to collect results')
    # parser.add_argument('-n', '--num-steps', type=int,
    #     help='Specify the number of time steps in output file')
        
    return parser


def main(path, output_file, timesteps):

    # convert strings to ints for timesteps
    #timesteps = [int(i) for i in timesteps]
    
    base_files = ['xGrid',
                  'yGrid',
                  'zGrid',
                  'dxGrid',
                  'dyGrid',
                  'dzGrid',
                  'vol',
                  'rho',
                  'sig11Grid',
                  'sig12Grid',
                  'sig13Grid',
                  'sig22Grid',
                  'sig23Grid',
                  'sig33Grid',
                  ]
                  

    # parse output
    results_dict = collect_output(path, base_files, timesteps)
    print(results_dict.keys())
    
    # create set dictionary
    print('Getting sets!')
    set_dict = generate_sets(results_dict, timesteps)
    ## check outputs
    print(set_dict.keys())
    for key in set_dict.keys():
        print(f'set {key} shape = {np.shape(set_dict[key])}')

    # #print(results)
    convert_to_XDMF(results_dict, output_file, set_dict)

    return


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(main(path=args.path,
                  output_file=args.output_file,
                  timesteps=args.timesteps,
                  ))