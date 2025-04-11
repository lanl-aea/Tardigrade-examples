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

    :returns: dictionary containing element field arrays and integer specifying the number of elements
    """

    equal_block_flag = False
    num_element_blocks = len([item for item in exofile.variables if 'connect' in item])
    if num_elements is None:
        num_elements = num_element_blocks
        equal_block_flag = True
    else:
        multiplier = int(num_elements / num_element_blocks)

    data_dict = {}
    # store each variable as a key in the dictionary
    # for each variable, store the timeseries in the column with a column for each element
    for v in cell_variable_keys.keys():
        out_array = numpy.zeros([num_times, num_elements])
        if equal_block_flag == True:
            for e in range(1,num_element_blocks+1):
                string = f'vals_elem_var{v}eb{e}'
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


def extract_node_data(exofile, node_variable_keys, x, y, z):
    """Extract node variables from exodus file

    :params netCDF4-Dataset exofile: Open netCDF4 dataset containing exodus data
    :params array-like x: The reference x-coordinates for the mesh
    :params array-like y: The reference y-coordinates for the mesh
    :params array-like z: The reference z-coordinates for the mesh

    :returns: dictionary containing nodal field arrays
    """

    data_dict = {}

    # store each variable as a key in the dictionary
    # for each variable, store in the timeseries in a column with a column for each node
    for v in node_variable_keys.keys():
        string = f'vals_nod_var{v}'
        data_dict[node_variable_keys[v]] = exofile.variables[string][:].data

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


def decode_chunk(chunk):
    """Convert the strange array of characters to a regular string for variable names from netCDF4

    :params array chunk: An array of encoded bytes representing a variable name

    :returns: properly formatted string
    """
    return ''.join(byte.decode('utf-8') for byte in chunk if byte != b'')


def extract_exodus_data(exodus_file, output_cell_data=None, output_node_data=None,
                        output_plot_base_name=None, output_dt_plot_base_name=None):
    """Process results from a MOOSE exodus simulation results file

    :params str exodus_file: The MOOSE exodus simulation results file
    :params str output_cell_data: Optional output netcdf file containing xarray of collected cell data
    :params str output_node_data: Optional output netcdf file containing xarray of collected node data
    :params str output_plot_base_name: Optional basename for field output plots
    :params str output_dt_plot_base_name: Optional basename for dt history plots

    :returns: Write ``{output_cell_data}`` and ``{output_node_data}``
    """

    exofile = netCDF4.Dataset(exodus_file, mode='r')

    x = exofile.variables['coordx'][:].data
    y = exofile.variables['coordy'][:].data
    z = exofile.variables['coordz'][:].data
    num_nodes = exofile.dimensions['num_nodes'].size
    num_elements = exofile.dimensions['num_elem'].size

    times = exofile.variables['time_whole'][:]
    num_times = len(times.data)

    cell_vars = [decode_chunk(chunk) for chunk in exofile.variables['name_elem_var'][:].data]
    node_vars = [decode_chunk(chunk) for chunk in exofile.variables['name_nod_var'][:].data]

    cell_variable_keys = {i+1: name for i, name in enumerate(cell_vars)}
    node_variable_keys = {i+1: name for i, name in enumerate(node_vars)}

    # Cell Data
    if (output_cell_data is not None) or (output_plot_base_name is not None):
        cell_data, num_elements = extract_cell_data(exofile, num_times, cell_variable_keys, num_elements)
        element_ids = list(range(1,num_elements+1))
        if output_cell_data:
            cell_dataset = xarray.Dataset({key: (('time', 'element'), value) for key, value in cell_data.items()},
                                        coords={'time': times, 'element': element_ids})
            cell_dataset.to_netcdf(output_cell_data)

    # Node Data
    if output_node_data:
        node_data = extract_node_data(exofile, node_variable_keys, x, y, z)
        node_ids = list(range(1, num_nodes+1))
        node_dataset = xarray.Dataset({key: (('time', 'node'), value) for key, value in node_data.items()},
                                    coords={'time': times, 'node': node_ids})
        node_dataset.to_netcdf(output_node_data)

    # Cell quantity plots
    if output_plot_base_name:
        for quantity in cell_vars:
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
    parser.add_argument('--output-cell-data', type=str, required=False, default=None,
        help="Optional output netcdf file containing xarray of collected cell data")
    parser.add_argument('--output-node-data', type=str, required=False, default=None,
        help="Optional output netcdf file containing xarray of collected node data")
    parser.add_argument('--output-plot-base-name', type=str, required=False, default=None,
        help="Optional basename for field output plots")
    parser.add_argument('--output-dt-plot-base-name', type=str, required=False, default=None,
        help="Optional basename for dt history plots")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(extract_exodus_data(exodus_file=args.exodus_file,
                                 output_cell_data=args.output_cell_data,
                                 output_node_data=args.output_node_data,
                                 output_plot_base_name=args.output_plot_base_name,
                                 output_dt_plot_base_name=args.output_dt_plot_base_name,
                                 ))
