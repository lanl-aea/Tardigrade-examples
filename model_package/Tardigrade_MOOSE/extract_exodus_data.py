#! /usr/bin/env python

import argparse
import pathlib
import sys

import matplotlib.pyplot
import netCDF4
import numpy
import xarray


def extract_cell_data(exofile, num_times, cell_variable_keys, num_elements=None):
    """Extract cell variables from exodus file

    :params netCDF4-Dataset exofile: Open netCDF4 dataset containing exodus data
    :params int num_times: The number of time steps in the Exodus results file
    :params dict cell_variable_keys: dictionary containing integer keys mapping alphabetically sorted cell data variable names
    :params int num_elements: Optional number of elements to help element extraction for refined meshes

    :returns: dictionary containing element field arrays and integer specifying the number of elements
    """

    equal_block_flag = False
    num_element_blocks = len([item for item in exofile.variables if 'connect' in item])
    if num_elements is None:
        num_elements = num_element_blocks
        equal_block_flag = True
    else:
        multiplier = int(num_elements / num_element_blocks)
    num_variables = len(list(cell_variable_keys.keys()))

    data_dict = {}
    # store each variable as a key in the dictionary
    # for each variable, store the timeseries in the column with a column for each element
    for v in range(1,num_variables+1):
        out_array = numpy.zeros([num_times, num_elements])
        if equal_block_flag == True:
            for e in range(1,num_element_blocks+1):
                string = f'vals_elem_var{v}eb{e}'
                print(string)
                print(cell_variable_keys[v])
                print(exofile.variables[string])
                print(numpy.shape(exofile.variables[string][:].data))
                block_data = exofile.variables[string][:].data
                out_array[:,e-1] = block_data.flatten()
        else:
            for b in range(1,num_element_blocks+1):
                string = f'vals_elem_var{v}eb{b}'
                block_data = exofile.variables[string][:].data
                for j in range(0, multiplier):
                    e = (b - 1)*multiplier + j
                    out_array[:,e] = block_data[:,j].flatten()
        data_dict[cell_variable_keys[v]] = out_array

    return data_dict, num_elements


def extract_node_data(exofile, x, y, z):
    """Extract node variables from exodus file

    :params netCDF4-Dataset exofile: Open netCDF4 dataset containing exodus data
    :params array-like x: The reference x-coordinates for the mesh
    :params array-like y: The reference y-coordinates for the mesh
    :params array-like z: The reference z-coordinates for the mesh

    :returns: dictionary containing nodal field arrays
    """

    variable_keys = {
        1: 'disp_x', 2: 'disp_y', 3: 'disp_z',
        4: 'force_x', 5: 'force_y', 6: 'force_z',
        7: 'phi_xx', 8: 'phi_xy', 9: 'phi_xz',
        10: 'phi_yx', 11: 'phi_yy', 12: 'phi_yz',
        13: 'phi_zx', 14: 'phi_zy', 15: 'phi_zz'}
    num_variables = 15

    data_dict = {}
    # store each variable as a key in the dictionary
    # for each variable, store in the timeseries in a column with a column for each node
    for v in range(1,num_variables+1):
        string = f'vals_nod_var{v}'
        data_dict[variable_keys[v]] = exofile.variables[string][:].data

    return data_dict


def plot_cell_data_over_time(output_file, cell_data, data_key, times):
    """Plot cell data over time for all elements

    :param str output_file: The output plot file name
    :param dict cell_data: A dictionary containing all cell data arrays
    :param str data_key: The dictionary key corresponding to the array to be plotted
    :param array-like times: Array of simulation times

    :returns: Write ``{output_file}`` plot
    """

    data_array = cell_data[data_key]
    num_elements = numpy.shape(data_array)[1]

    matplotlib.pyplot.figure()
    for e in range(num_elements):
        matplotlib.pyplot.plot(times, data_array[:,e])
    matplotlib.pyplot.xlabel('Time (s)')
    matplotlib.pyplot.ylabel(data_key)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_file)

    return 0


def plot_delta_t(times, output_file):
    """Plot timestep history

    :param array-like times: Array of simulation times
    :param str output_file: The output plot file name

    :returns: Write ``{output_file}`` plot
    """

    dts = numpy.diff(times)
    dts = numpy.insert(dts, 0,0)

    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(times, dts)
    matplotlib.pyplot.xlabel('Time (s)')
    matplotlib.pyplot.ylabel('delta t (s)')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(output_file)

    return 0


def get_variable_keys_from_MOOSE_input(input_file):
    """Extract and sort the names of auxilliary variables from a MOOSE input deck corresponding to cell data

    :params str input_file: MOOSE input file used for extracting cell field variable names from AuxVariables
    
    :returns: sorted list of auxiliary variables
    """

    aux_vars = []
    in_aux_block = False

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Detect start and end of [AuxVariables] block
            if line.startswith('[AuxVariables]'):
                in_aux_block = True
                continue
            elif line.startswith('[]') and in_aux_block:
                in_aux_block = False
                continue

            # If inside AuxVariables block, look for variable declarations
            if in_aux_block:
                if line.startswith('[') and line.endswith(']'):
                    name = line.strip('[]').lstrip('./')
                    if name:  # Ignore empty brackets
                        aux_vars.append(name)

    # remove forces since those are nodal quantities
    aux_vars = [var for var in aux_vars if 'force' not in var]
    aux_vars.sort()

    print(f'AuxVariables detected in {input_file}:')
    print(aux_vars)

    return aux_vars


def extract_exodus_data(exodus_file, MOOSE_input_file=None,
                        output_cell_data=None, output_node_data=None,
                        output_plot_base_name=None, output_dt_plot_base_name=None,
                        simulation_type='plastic', num_elements=None):
    """Process results from a MOOSE exodus simulation results file

    :params str exodus_file: The MOOSE exodus simulation results file
    :params str MOOSE_input_file: Optional MOOSE input file used for extracting cell field variable names from AuxVariables
    :params str output_cell_data: Optional output netcdf file containing xarray of collected cell data
    :params str output_node_data: Optional output netcdf file containing xarray of collected node data
    :params str output_plot_base_name: Optional basename for field output plots
    :params str output_dt_plot_base_name: Optional basename for dt history plots
    :params str simulation_type: The simulation type, currently only 'elastic' or 'plastic' is supported

    :returns: Write ``{output_cell_data}`` and ``{output_node_data}``
    """

    exofile = netCDF4.Dataset(exodus_file, mode='r')

    x = exofile.variables['coordx'][:].data
    y = exofile.variables['coordy'][:].data
    z = exofile.variables['coordz'][:].data
    num_nodes = len(x)

    times = exofile.variables['time_whole'][:]
    num_times = len(times.data)

    # Extract cell data variables if MOOSE_input_file provided
    if MOOSE_input_file is not None:
        aux_variables = get_variable_keys_from_MOOSE_input(MOOSE_input_file)
        cell_variable_keys = {}
        for i, var in enumerate(aux_variables):
            cell_variable_keys[i+1] = var
    elif simulation_type == 'plastic':
        cell_variable_keys = {
            1: "macro_gamma", 2: "macro_isv", 3: "micro_gamma",
            4: "micro_gradient_gamma_1",
            5: "micro_gradient_gamma_2",
            6: "micro_gradient_gamma_3",
            7: "micro_gradient_isv_1",
            8: "micro_gradient_isv_2",
            9: "micro_gradient_isv_3",
            10: "micro_isv",
            11: "pk2_11", 12: "pk2_22", 13: "pk2_33",
            14: "sigma_11", 15: "sigma_22", 16: "sigma_33"}
    elif simulation_type == 'elastic':
        cell_variable_keys = {
            1: "pk2_11", 2: "pk2_22", 3: "pk2_33",
            4: "sigma_11", 5: "sigma_22", 6: "sigma_33"}
    else:
        raise NotImplementedError("Combination of MOOSE_input_file=None and simulation_type!= 'elastic'\
                                   or 'plastic' has not bee implemented!")

    # Cell Data
    cell_data, num_elements = extract_cell_data(exofile, num_times, cell_variable_keys, num_elements)
    element_ids = list(range(1,num_elements+1))
    if output_cell_data:
        cell_dataset = xarray.Dataset({key: (('time', 'element'), value) for key, value in cell_data.items()},
                                    coords={'time': times, 'element': element_ids})
        cell_dataset.to_netcdf(output_cell_data)

    # Node Data
    node_data = extract_node_data(exofile, x, y, z)
    node_ids = list(range(1, num_nodes+1))
    if output_node_data:
        node_dataset = xarray.Dataset({key: (('time', 'node'), value) for key, value in node_data.items()},
                                    coords={'time': times, 'node': node_ids})
        node_dataset.to_netcdf(output_node_data)

    # Cell quantity plots
    if output_plot_base_name:
        plot_quantities = ["macro_gamma", "macro_isv", "micro_gamma", "micro_isv",
                           "pk2_33", "sigma_33"]
        for quantity in plot_quantities:
            output_plot_file = f'{output_plot_base_name}_{quantity}.png'
            plot_cell_data_over_time(output_plot_file, cell_data, quantity, times)

    # dt plots
    if output_dt_plot_base_name:
        plot_delta_t(times, f'{output_dt_plot_base_name}_delta_ts.png')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Process results from a MOOSE exodus simulation results file"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)

    parser.add_argument('--exodus-file', type=str, required=True,
        help="The MOOSE exodus simulation results file")
    parser.add_argument('--MOOSE-input-file', type=str, required=True,
        help="Optional MOOSE input file used for extracting cell field variable names from AuxVariables")
    parser.add_argument('--output-cell-data', type=str, required=False, default=None,
        help="Optional output netcdf file containing xarray of collected cell data")
    parser.add_argument('--output-node-data', type=str, required=False, default=None,
        help="Optional output netcdf file containing xarray of collected node data")
    parser.add_argument('--output-plot-base-name', type=str, required=False, default=None,
        help="Optional basename for field output plots")
    parser.add_argument('--output-dt-plot-base-name', type=str, required=False, default=None,
        help="Optional basename for dt history plots")
    parser.add_argument('--simulation-type', type=str, required=False, default='plastic',
        help="The simulation type, currently only 'elastic' or 'plastic' is supported")
    parser.add_argument('--num-elements', type=int, required=False, default=None,
        help="Optional number of elements to help element extraction for refined meshes")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(extract_exodus_data(exodus_file=args.exodus_file,
                                 MOOSE_input_file=args.MOOSE_input_file,
                                 output_cell_data=args.output_cell_data,
                                 output_node_data=args.output_node_data,
                                 output_plot_base_name=args.output_plot_base_name,
                                 output_dt_plot_base_name=args.output_dt_plot_base_name,
                                 simulation_type=args.simulation_type,
                                 num_elements=args.num_elements,
                                 ))
