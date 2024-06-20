# Imports
from pdb import set_trace as dbg
import os.path as osp
import sys
import collections
import re
import os
import xml
import lxml.etree
import h5py
import inspect
import numpy as np
#import filter_xdmf_writer as fxw
import xml.etree.ElementTree as ET

import argparse
import matplotlib.pyplot as plt

#==================================================== XDMF UTILITY FUNCTIONS ===
def build_xdmf_root():
    """
    Build the root of the XDMF file
    
    """
    root = ET.Element('Xdmf')
    root.set("Version", "3.0")
    root.set("xmlns:xi",r"http://www.w3.org/2001/XInclude")
    return root
    
def add_information(node, **kwargs):
    """
    Add information object

    :param ET.SubElement node: The parent node
    :param dict kwargs: The information to be added in keyword, value form
    """

    information = ET.SubElement(node, "Information")
    for key in kwargs.keys():
        information.set(key, kwargs[key])

    return information
    
def build_domain(root, name="DOMAIN", value="Default domain name"):
    """
    Build a domain in the parent domain

    :param ET.Element root: The root XML object
    :param str name: The keyword name of the domain
    :param str value: Descriptive information for the domain
    """
    domain = ET.SubElement(root, "Domain")
    domain_information = add_information(domain, Name=name, Value=value)

    return domain
    
def build_grid(domain, grid_info, **kwargs):
    """
    Build a grid in the parent domain

    :param ET.SubElement domain: The domain to construct the grid in
    :param dict grid_info: The information for the grid
    :param dict kwargs: Additional information in keword, value form
    """

    grid = ET.SubElement(domain, "Grid")
    for key in grid_info.keys():
        grid.set(key, grid_info[key])
        
    if len(kwargs)>0:
        grid_info = add_information(grid, **kwargs)

    return grid
    
def add_time(grid, value):
    """
    Add a time object to a grid
    
    :param ET.SubElement grid: The XML representation of the grid object
    :param float value: The floating point representation of the time
    """

    time = ET.SubElement(grid, "Time")
    time.set("Value", str(value))

    return time

def add_geometry(grid, geometry_info, **kwargs):
    """
    Add a geometry object

    :param ET.SubElement grid: The XML grid object
    :param dict geometry_info: The geometry values to be added in keyword, value form
    :param dict kwargs: The descriptive information for the geometry to be added in
        keyword, value form
    """

    geometry = ET.SubElement(grid, "Geometry")

    for key in geometry_info:
        geometry.set(key, geometry_info[key])

    if len(kwargs)>0:
        geometry_information = add_information(geometry, **kwargs)

    return geometry

def add_topology(grid, topology_info, **kwargs):
    """
    Add a topology object

    :param ET.SubElement grid: The XML grid object
    :param dict topology_info: The topology values to be added in keyword, value form
    :param dict kwargs: The descriptive information for the topology to be added in
        keyword, value form
    """

    topology = ET.SubElement(grid, "Topology")

    for key in topology_info.keys():
        topology.set(key, topology_info[key])

    if len(kwargs)>0:
        topology_information = add_information(topology, **kwargs)

    return topology

def add_attribute(grid, attribute_info, additional_info=None, h5_info=None, filename=None, dataname=None):
    """
    Add an attribute object to the XDMF file

    :param ET.SubElement grid: The grid to add the attribute to
    :param dict attribute_info: The attribute values to be added in keyword, value form
    :param dict additional_info: Additional information to be added for the attribute
    :param dict h5_info: The data to be added to the h5 file in keyword, value form
    :param str filename: The h5 filename
    :param str dataname: The name of the data in the h5 file
    """

    attribute = ET.SubElement(grid, "Attribute")

    for key in attribute_info:
        attribute.set(key, attribute_info[key])

    if (not (additional_info is None)):
        information = add_information(attribute, **additional_info)

    if (not (filename is None)):
        if (dataname is None):
            raise ValueError("dataname is None")
        if (h5_info is None):
            raise ValueError("h5_info is None")

        data_item = add_h5_data(attribute, filename, dataname, h5_info)
        
    return

def add_h5_data(node, filename, dataname, data_info, h5file=None, data=None, dtype=None, **kwargs):
    """
    Add data to the h5 file

    :param ET.SubElement node: The XML node that is attached to the h5 data
    :param str filename: The h5 filename
    :param dict data_info: The data values to be added in keyword, value form
    :param h5.File h5file: The h5 file
    :param np.array data: The data array
    :param str dtype: The datatype
    :param dict kwargs: Not currently used
    """

    if not ((h5file is None) and (data is None) and (dtype is None)):
        h5file.create_dataset(dataname, data=data, dtype=dtype)

    data_item = ET.SubElement(node, "DataItem")

    for key in data_info.keys():
        data_item.set(key, data_info[key])

    data_item.text = filename + ":" + dataname

    return data_item

def add_set(grid, set_info, h5_info=None, filename=None, dataname=None):
    """
    Add a set object to the XDMF file

    :param ET.SubElement grid: The XML node that is attached to the grid data
    :param dict set_info: The information for the set in keyword, value form
    :param dict h5_info: The information for the h5 file in keyword, value form
    :param str filename: The h5 filename
    :param str dataname: The name of the data in the h5 file
    """

    set_object = ET.SubElement(grid, "Set")

    for key in set_info:
        set_object.set(key, set_info[key])

    if (not (filename is None)):
        if (dataname is None):
            raise ValueError("dataname is None")
        if (h5_info is None):
            raise ValueError("h5_info is None")

        data_item = add_h5_data(set_object, filename, dataname, h5_info)
        
    return
    
def initialize_xdmf_domain():
    """
    Initialize the XDMF domain
    
    """

    # Initialize the XDMF file structure
    root = build_xdmf_root()
    
    # Build the root domain
    domain = build_domain(root)
    
    # Build the grid collection
    collection_options = {"CollectionType":"Temporal",
                          "GridType":"Collection",
                          "Name":"Collection"}
    
    collection_info = {"Name":"Main_Temporal_Collection",
                       "Value":"The main temporal ( or iteration ) collection"}
    
    collection = build_grid(domain, collection_options, **collection_info)

    return root, domain, collection
    
def construct_grid(collection, grid_data, h5file):
    """
    Construct the grid object

    :param ET.SubElement collection: The collection XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param h5py.File h5file: The hdf5 file to write to
    """

    grid = build_grid(collection, {"Name":"Grid"})

    # Set the time
    time_node = add_time(grid, grid_data['time'])

    # Set the geometry
    geometry = construct_geometry(grid, grid_data, h5file)

    # Set the topology
    topology = construct_topology(grid, grid_data, h5file)

    # Set the sets
    sets = [construct_set(grid, grid_data, set_n['name'], set_n['ids'], h5file) for set_n in grid_data['sets']]

    # Set the node ids
    node_attributes = construct_node_ids(grid, grid_data, h5file)

    # Set the element ids
    element_attributes = construct_element_ids(grid, grid_data, h5file)

    # Add the attributes
    #[construct_attribute(grid, grid_data, attribute_name, h5file) for attribute_name in grid_data['attributes'].keys()]
    for attribute_name in grid_data['attributes'].keys():
        print(attribute_name)
        construct_attribute(grid, grid_data, attribute_name, h5file)
        
    return
    
def construct_geometry(grid, grid_data, h5file):
    """
    Construct the geometry for a grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param h5py.File h5file: The hdf5 file to write to
    """

    # Set the coordinates
    coordinates_name = grid_data['coordinates']['name']
    coordinates      = grid_data['coordinates']['values']
    if (grid_data['initial_grid']):

        # Write the data to the h5 file

        h5file.create_dataset(coordinates_name, data=coordinates.flatten(), dtype='d')

    # Create the XML link
    geometry_info = {"Origin":"",
                     "Type":"XYZ"}
    additional_info = {"Name":"Coordinates",
                       "Value":"Coordinates of the nodes in x1, y1, z1, x2, ... format"}
    geometry_data_info = {"DataType":"Float",
                          "Dimensions":f"{coordinates.size}",
                          "Format":"HDF",
                          "Precision":"8"}

    geometry = add_geometry(grid, geometry_info, **additional_info)

    add_h5_data(geometry, h5file.filename, coordinates_name, geometry_data_info)

    return geometry

def construct_topology(grid, grid_data, h5file):
    """
    Construct the topology for a grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param h5py.File h5file: The hdf5 file to write to
    """

    # Set the topology
    topology_name = grid_data['topology']['name']
    if (grid_data['initial_grid']):

        h5file.create_dataset(topology_name, data=grid_data['ids'], dtype='i')

    # Create the XML link
    topology_info = {"TopologyType":"Polyvertex",
                     "NumberOfElements":str(grid_data['ids'].size)}

    topology = add_topology(grid, topology_info)

    topology_data_info = {"DataType":"UInt",
                          "Format":"HDF",
                          "Dimensions":str(grid_data['ids'].size)}

    add_h5_data(topology, h5file.filename, topology_name, topology_data_info)

    return topology

def construct_set(grid, grid_data, name, data, h5file):
    """
    Construct a set in the grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param str name: The set name
    :param np.ndarray data: The data for the set
    :param h5py.File h5file: The hdf5 file to write to
    """

    if (grid_data['initial_grid']):

        h5file.create_dataset(name, data=data, dtype='i')

    set_info = {'Name':name,
                'Type':'Node'}

    data_size = np.shape(data)

    set_data_info = {'DataType':'UInt',
                     'Dimensions':str(data_size),
                     'Format':'HDF',
                     'Precision':'4'}

    return add_set(grid, set_info, h5_info=set_data_info, filename=h5file.filename, dataname=name)

def construct_node_ids(grid, grid_data, h5file):
    """
    Construct the node ids for a grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param h5py.File h5file: The hdf5 file to write to
    """

    # Set the topology
    nodeid_name = 'NODEID'
    if (grid_data['initial_grid']):

        h5file.create_dataset(nodeid_name, data=grid_data['ids'], dtype='i')

    # Create the XML link
    attribute_info = {"Center":"Node",
                      "ElementCell":"",
                      "ElementDegree":"0",
                      "ElementFamily":"",
                      "ItemType":"",
                      "Name":nodeid_name,
                      "Type":"Scalar"}
    additional_info = {"Name":"ID",
                       "Value":"The nodal IDs"}

    node_id_h5info = {"DataType":"UInt",
                      "Dimensions":f"{grid_data['ids'].size}",
                      "Format":"HDF",
                      "Precision":"4"}

    return add_attribute(grid, attribute_info, additional_info=additional_info,
                         h5_info=node_id_h5info, filename=h5file.filename,
                         dataname=nodeid_name)

def construct_element_ids(grid, grid_data, h5file):
    """
    Construct the element ids for a grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param h5py.File h5file: The hdf5 file to write to
    """

    # Set the topology
    elemid_name = 'ELEMID'
    if (grid_data['initial_grid']):

        h5file.create_dataset(elemid_name, data=grid_data['ids'], dtype='i')

    # Create the XML link
    attribute_info = {"Center":"Cell",
                      "ElementCell":"",
                      "ElementDegree":"0",
                      "ElementFamily":"",
                      "ItemType":"",
                      "Name":elemid_name,
                      "Type":"Scalar"}
    additional_info = {"Name":"ID",
                       "Value":"The element IDs"}

    elem_id_h5info = {"DataType":"UInt",
                      "Dimensions":f"{grid_data['ids'].size}",
                      "Format":"HDF",
                      "Precision":"4"}

    return add_attribute(grid, attribute_info, additional_info=additional_info,
                         h5_info=elem_id_h5info, filename=h5file.filename,
                         dataname=elemid_name)
                         

def construct_displacements(grid, grid_data, h5file):
    """
    Construct the displacements for a grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param h5py.File h5file: The hdf5 file to write to
    """

    for i, direction in enumerate(('x', 'y', 'z')):

        name = grid_data['displacements']['name'] + direction  + "_" + grid_data['frame']

        h5file.create_dataset(name, data=grid_data['displacements']['values'][:,i], dtype='d')

        attribute_info = {"Center":"Node",
                          "ElementCell":"",
                          "ElementDegree":"0",
                          "ElementFamily":"",
                          "ItemType":"",
                          "Name":grid_data['displacements']['name'] + direction,
                          "Type":"Scalar"}

        additional_info = {"Name":grid_data['displacements']['name'] + direction,
                           "Value":"The displacement in the {direction} direction"}

        disp_h5info = {"DataType":"Float",
                       "Dimensions":f"{grid_data['displacements']['values'].shape[0]}",
                       "Format":"HDF",
                       "Precision":"8"}
        add_attribute(grid, attribute_info, additional_info=additional_info,
                      h5_info=disp_h5info, filename=h5file.filename,
                      dataname=name)
                      
    return

def construct_velocities(grid, grid_data, h5file):
    """
    Construct the velocities for a grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param h5py.File h5file: The hdf5 file to write to
    """

    for i, direction in enumerate(('x', 'y', 'z')):

        name = grid_data['velocities']['name'] + direction  + "_" + grid_data['frame']

        h5file.create_dataset(name, data=grid_data['velocities']['values'][:,i], dtype='d')

        attribute_info = {"Center":"Node",
                          "ElementCell":"",
                          "ElementDegree":"0",
                          "ElementFamily":"",
                          "ItemType":"",
                          "Name":grid_data['velocities']['name'] + direction,
                          "Type":"Scalar"}

        additional_info = {"Name":grid_data['velocities']['name'] + direction,
                           "Value":"The velocity in the {direction} direction"}

        disp_h5info = {"DataType":"Float",
                       "Dimensions":f"{grid_data['velocities']['values'].shape[0]}",
                       "Format":"HDF",
                       "Precision":"8"}
        add_attribute(grid, attribute_info, additional_info=additional_info,
                      h5_info=disp_h5info, filename=h5file.filename,
                      dataname=name)
    return

def construct_attribute(grid, grid_data, attribute_name, h5file):
    """
    Construct the attribute for the grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param hdf5.File h5file: The hdf5 file to write to
    """

    components = grid_data['attributes'][attribute_name]['ordering']
    values = grid_data['attributes'][attribute_name]['values'].reshape((-1,len(components)))

    for i, component in enumerate(components):

        # fix some of the strings in the "component" variable - Thomas
        #print(f'old component = {component}')
        if "b'" in str(component):
            #print(f'old component: {component}')
            comp0 = str(component)
            comp1 = comp0.replace("'","")
            comp2 = comp1.replace("b","")
            component = str(comp2)
            #print(f'new component: {component}')
        #print(f'new component = {component}')
        #print(np.shape(component))
        
        #name = f"{grid_data['attributes'][attribute_name]['name']}{str(component)}_{grid_data['frame']}"
            
        #print('part 1', grid_data['attributes'][attribute_name]['name'])
        #print('part 2', component + "_" + grid_data['frame'])
        name = grid_data['attributes'][attribute_name]['name'] + str(component) + "_" + str(grid_data['frame'])
        value = values[:,i]

        h5file.create_dataset(name, data=value)

        attribute_info = {"Center":"Node",
                          "ElementCell":"",
                          "ElementDegree":"0",
                          "ElementFamily":"",
                          "ItemType":"",
                          "Name":grid_data['attributes'][attribute_name]['name'] + component,
                          "Type":"Scalar"}

        additional_info = {"Name":grid_data['attributes'][attribute_name]['name'] + component,
                           "Value":"The " + attribute_name + " in the {direction} direction"}

        disp_h5info = {"DataType":"Float",
                       "Dimensions":f"{value.shape[0]}",
                       "Format":"HDF",
                       "Precision":"8"}
        add_attribute(grid, attribute_info, additional_info=additional_info,
                      h5_info=disp_h5info, filename=h5file.filename,
                      dataname=name)
    return

def construct_accelerations(grid, grid_data, h5file):
    """
    Construct the accelerations for a grid
    
    :param ET.SubElement grid: The grid XML node
    :param dict grid_data: The data required for the grid to be built
        in keyword, value pairs
    :param h5py.File h5file: The hdf5 file to write to
    """

    for i, direction in enumerate(('x', 'y', 'z')):

        name = grid_data['accelerations']['name'] + direction  + "_" + grid_data['frame']

        h5file.create_dataset(name, data=grid_data['accelerations']['values'][:,i], dtype='d')

        attribute_info = {"Center":"Node",
                          "ElementCell":"",
                          "ElementDegree":"0",
                          "ElementFamily":"",
                          "ItemType":"",
                          "Name":grid_data['accelerations']['name'] + direction,
                          "Type":"Scalar"}

        additional_info = {"Name":grid_data['accelerations']['name'] + direction,
                           "Value":"The acceleration in the {direction} direction"}

        disp_h5info = {"DataType":"Float",
                       "Dimensions":f"{grid_data['accelerations']['values'].shape[0]}",
                       "Format":"HDF",
                       "Precision":"8"}
        add_attribute(grid, attribute_info, additional_info=additional_info,
                      h5_info=disp_h5info, filename=h5file.filename,
                      dataname=name)
                      
    return
    
def sort_block_names(names):
    block_nums = [n.replace('block','') for n in names]
    return ['block' + n for n in sort_string_names(block_nums)]

def sort_string_names(names):
    return [str(v) for v in np.sort([int(n) for n in names])]