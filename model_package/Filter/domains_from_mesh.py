# read in xdmf mesh, identify nodes and cells

import inspect
import sys
import os
import argparse
import pathlib

import meshio
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import h5py
import shutil
import lxml.etree
import xml.etree.ElementTree as ET

import xdmf_reader_tools as XRT

def load_mesh(filename):
    #filename = 'mesh_36_elements.py'
    mesh=meshio.read(filename)

    nodes = mesh.points
    elements = mesh.cells[0].data

    num_nodes = numpy.shape(nodes)[0]
    num_elems = numpy.shape(elements)[0]

    node_ids = list(range(0, num_nodes))
    elem_ids = list(range(0, num_elems))

    node_dict, elem_dict = {}, {}

    # Make dictionary of nodal coordinates
    for n in node_ids:
        node_dict[n] = nodes[n, :]

    # Make dictionary of element connectivity
    for e in elem_ids:
        elem_dict[e] = elements[e, :]


    # # make sure that we can query element nodes on the fly
    # for e in elem_dict.keys():
        # print(e)
        # for n in elem_dict[e]:
            # print(f'\tnode{n} coords: {node_dict[n]}')
            
    return(node_dict, elem_dict, numpy.array(nodes))


def build_major_planes(macro_node_dict, macro_elem_dict):
    """
    get 6 planes that define the element bounds.
    
    ax + by + cz + d = 0
    
    # list of nodes in each plane CCW
    plane 1 (+y): 0, 4, 5, 1
    plane 2 (-x): 1, 5, 6, 2
    plane 3 (-y): 2, 6, 7, 3
    plane 4 (+y): 3, 7, 4, 0
    plane 5 (+z): 0, 1, 2, 3
    plane 6 (-z): 4, 7, 6, 5
    
    # define planes such that vectors point into the hex
    
    for each plane:
    1. get vector A for first and second node
    2. get vector B for first and third node
    3. cross product to get normal --> (a, b, c)
    4. dot product between normal and first node --> d
    
    
    """
    
    element_ids = list(macro_elem_dict.keys())
    major_plane_dict = {}
    planes = [f"plane_{i}" for i in range(1, 7)]
    plane_tup = [(0, 4, 5),
                 (1, 5, 6),
                 (2, 6, 7),
                 (3, 7, 4),
                 (0, 1, 2),
                 (4, 7, 6)]
    # plane_tup = [(0, 4, 5),
                 # (1, 2, 6),
                 # (2, 3, 7),
                 # (3, 7, 4),
                 # (0, 1, 2),
                 # (4, 5, 6)]
    # plane_tup = [(0, 4, 5, 1),
                 # (1, 5, 6, -1),
                 # (2, 6, 7, -1),
                 # (3, 7, 4, 1),
                 # (0, 1, 2, 1),
                 # (4, 7, 6, -1)]
    
    for eid in element_ids:
        elem_nodes = macro_elem_dict[eid]
        center = numpy.mean(numpy.array([macro_node_dict[n] for n in elem_nodes]), axis=0)
        temp = {}
        for plane, tup in zip(planes, plane_tup):
            # get coors for each point of the plane tuple
            p1 = macro_node_dict[elem_nodes[tup[0]]]
            p2 = macro_node_dict[elem_nodes[tup[1]]]
            p3 = macro_node_dict[elem_nodes[tup[2]]]
            # vector calculations
            vecA = p1 - p2
            vecB = p3 - p2
            n = numpy.cross(vecB, vecA)
            d = -1* numpy.dot(n, p2)
            # get sign
            #s = numpy.sign(numpy.dot(n, center) + d)
            #s = s*tup[3]
            # store
            temp[plane] = [n[0], n[1], n[2], d]
        major_plane_dict[eid] = temp

    return major_plane_dict


def build_minor_planes(macro_node_dict, macro_elem_dict):
    """
    get 6 planes that define the element bounds.
    
    ax + by + cz + d = 0
    
    # list of nodes in each plane CCW
    plane 1 (+y): 0, 4, 5, 1
    plane 2 (-x): 1, 5, 6, 2
    plane 3 (-y): 2, 6, 7, 3
    plane 4 (+y): 3, 7, 4, 0
    plane 5 (+z): 0, 1, 2, 3
    plane 6 (-z): 4, 7, 6, 5
    plane 7 (+y): 9, 11, 15, 13 <-- mid-plane
    plane 8 (+x): 8, 10, 14, 12 <-- mid-plane
    plane 9 (+z): 16, 17, 18, 19 <-- mid-plane
    
    # define planes such that vectors point into the hex
    
    for each plane:
    1. get vector A for first and second node
    2. get vector B for first and third node
    3. cross product to get normal --> (a, b, c)
    4. dot product between normal and first node --> d
    
    
    """
    
    element_ids = list(macro_elem_dict.keys())
    minor_plane_dict = {}
    planes = [f"plane_{i}" for i in ['y', 'x', 'z']]
    plane_tup = [
                 # (0, 4, 5),
                 # (1, 5, 6),
                 # (2, 6, 7),
                 # (3, 7, 4),
                 # (0, 1, 2),
                 # (4, 7, 6),
                 (9, 11, 15),
                 (8, 10, 14),
                 (16, 17, 18)]

    # This time we'll create an element node dictionary that explicitly stores the nodal values
    # Then we'll add in the midpoints 
    # This is because the midpoints don't exist in the node_dict already!
    elem_coords = {}
    for eid in element_ids:
        elem_nodes = macro_elem_dict[eid]
        elem_coords[eid] = [numpy.array(macro_node_dict[n]) for n in elem_nodes]

        # Calculate midpoint of each edge
        c8 = numpy.mean([elem_coords[eid][0], elem_coords[eid][1]], axis=0)
        c9 = numpy.mean([elem_coords[eid][1], elem_coords[eid][2]], axis=0)
        c10 = numpy.mean([elem_coords[eid][2], elem_coords[eid][3]], axis=0)
        c11 = numpy.mean([elem_coords[eid][3], elem_coords[eid][0]], axis=0)
        c12 = numpy.mean([elem_coords[eid][4], elem_coords[eid][5]], axis=0)
        c13 = numpy.mean([elem_coords[eid][5], elem_coords[eid][6]], axis=0)
        c14 = numpy.mean([elem_coords[eid][6], elem_coords[eid][7]], axis=0)
        c15 = numpy.mean([elem_coords[eid][7], elem_coords[eid][0]], axis=0)
        c16 = numpy.mean([elem_coords[eid][0], elem_coords[eid][4]], axis=0)
        c17 = numpy.mean([elem_coords[eid][1], elem_coords[eid][5]], axis=0)
        c18 = numpy.mean([elem_coords[eid][2], elem_coords[eid][6]], axis=0)
        c19 = numpy.mean([elem_coords[eid][3], elem_coords[eid][7]], axis=0)
        midpoints = [c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19]
        #midpoint_ids = list(range(8,20))
        for m in midpoints:
            elem_coords[eid].append(m)

    for eid in element_ids:
        elem_nodes = elem_coords[eid]
        #center = numpy.mean(numpy.array([macro_node_dict[n] for n in elem_nodes]), axis=0)
        temp = {}
        for plane, tup in zip(planes, plane_tup):
            # get coors for each point of the plane tuple
            p1 = elem_nodes[tup[0]]
            p2 = elem_nodes[tup[1]]
            p3 = elem_nodes[tup[2]]
            # vector calculations
            vecA = p1 - p2
            vecB = p3 - p2
            n = numpy.cross(vecB, vecA)
            d = -1* numpy.dot(n, p2)
            # get sign
            #s = numpy.sign(numpy.dot(n, center) + d)
            #s = s*tup[3]
            # store
            temp[plane] = [n[0], n[1], n[2], d]
        minor_plane_dict[eid] = temp

    return minor_plane_dict


def find_major_sets(node_dict, major_plane_dict, option):
    #print(node_dict)
    tol = 1.e-8
    sets = [f'SET_{i}' for i in major_plane_dict.keys()]
    major_set_dict = {}
    for set, i in zip(sets, major_plane_dict.keys()):
        #a, b, c, d = major_plane_dict[i]
        set_ids = []
        for n in node_dict.keys():
            coords = node_dict[n]
            flag = True
            # general algorithm for looping through planes
            if option == 1:
                for plane in major_plane_dict[i].keys():
                    p = major_plane_dict[i][plane]
                    a, b, c, d = p[0], p[1], p[2], p[3]
                    arr = numpy.array([a, b, c])
                    #test = s*(numpy.dot(arr, coords) + d)
                    #test = numpy.sign(numpy.dot(arr, coords) + d)
                    test = numpy.dot(arr, coords) + d
                    if test <= 0.0:
                        flag = False
            # brute force
            elif option == 2:
                # plane 1
                a, b, c, d = major_plane_dict[i]['plane_1']
                arr = numpy.array([a, b, c])
                test = numpy.dot(arr, coords) + d
                if test < 0.0:
                    flag = False
                
                # plane 2
                a, b, c, d = major_plane_dict[i]['plane_2']
                arr = numpy.array([a, b, c])
                test = numpy.dot(arr, coords) + d
                if test <= 0.0:
                    flag = False
                
                # plane 3
                a, b, c, d = major_plane_dict[i]['plane_3']
                arr = numpy.array([a, b, c])
                test = numpy.dot(arr, coords) + d
                if test <= 0.0:
                    flag = False
                
                # plane 4
                a, b, c, d = major_plane_dict[i]['plane_4']
                arr = numpy.array([a, b, c])
                test = numpy.dot(arr, coords) + d
                if test < 0.0:
                    flag = False
                
                # plane 5
                a, b, c, d = major_plane_dict[i]['plane_5']
                arr = numpy.array([a, b, c])
                test = numpy.dot(arr, coords) + d
                if test < 0.0:
                    flag = False
                
                # plane 6
                a, b, c, d = major_plane_dict[i]['plane_6']
                arr = numpy.array([a, b, c])
                test = numpy.dot(arr, coords) + d
                if test <= 0.0:
                    flag = False
                    
            if flag == True:
                #print('Pass!')
                set_ids.append(n)
        major_set_dict[set] = set_ids
            # plane 1
            
    return major_set_dict

def find_minor_sets(node_dict, minor_plane_dict, major_set_dict):

    octants = ['a','b','c','d','e','f','g','h']
    sets = [f'Set_{i}{j}' for i in minor_plane_dict.keys() for j in octants]
    minor_set_dict = {}
    
    # loop over major_set_dict to get relevant nodes
    for major_set_name in major_set_dict.keys():
        major_set = int(major_set_name.split('_')[-1])
        # create empty lists for each octant of macro element
        for oct in octants:
            minor_set_dict[f'Set_{major_set}{oct}'] = []

        # unpack planes
        plane_x = minor_plane_dict[major_set]['plane_x']
        plane_y = minor_plane_dict[major_set]['plane_y']
        plane_z = minor_plane_dict[major_set]['plane_z']
        
        # unpack plane components
        ax, bx, cx, dx = plane_x[0], plane_x[1], plane_x[2], plane_x[3]
        ay, by, cy, dy = plane_y[0], plane_y[1], plane_y[2], plane_y[3]
        az, bz, cz, dz = plane_z[0], plane_z[1], plane_z[2], plane_z[3]
        arrx = numpy.array([ax, bx, cx])
        arry = numpy.array([ay, by, cy])
        arrz = numpy.array([az, bz, cz])
        
        # loop over points in major set
        major_pts = major_set_dict[major_set_name]
        for n in major_pts:
            coords = node_dict[n]
            
            # test which side of x plane
            test_x = numpy.dot(arrx, coords) + dx
            
            # test which side of y plane
            test_y = numpy.dot(arry, coords) + dy
            
            # test which side of z plane
            test_z = numpy.dot(arrz, coords) + dz
            
            # set a - octant 0
            if (test_x <= 0.0) and (test_y <= 0.0) and (test_z <= 0.0):
                minor_set_dict[f'Set_{major_set}a'].append(n)
                
            # set b - octant 1
            if (test_x > 0.0) and (test_y <= 0.0) and (test_z <= 0.0):
                minor_set_dict[f'Set_{major_set}b'].append(n)
            
            # set c - octant 2
            if (test_x > 0.0) and (test_y > 0.0) and (test_z <= 0.0):
                minor_set_dict[f'Set_{major_set}c'].append(n)
            
            # set d - octant 3
            if (test_x <= 0.0) and (test_y > 0.0) and (test_z <= 0.0):
                minor_set_dict[f'Set_{major_set}d'].append(n)
            
            # set e - octant 4
            if (test_x <= 0.0) and (test_y <= 0.0) and (test_z > 0.0):
                minor_set_dict[f'Set_{major_set}e'].append(n)
                
            # set f - octant 5
            if (test_x > 0.0) and (test_y <= 0.0) and (test_z > 0.0):
                minor_set_dict[f'Set_{major_set}f'].append(n)
            
            # set g - octant 6
            if (test_x > 0.0) and (test_y > 0.0) and (test_z > 0.0):
                minor_set_dict[f'Set_{major_set}g'].append(n)
            
            # set h - octant 7
            if (test_x <= 0.0) and (test_y > 0.0) and (test_z > 0.0):
                minor_set_dict[f'Set_{major_set}h'].append(n)
            
    return minor_set_dict


def run_test_major(major_plane_dict, macro_node_dict, macro_nodes, macro_elem_dict, plot=False):

    # find the extent of the nodal coordinates
    xmin = numpy.unique(numpy.min(macro_nodes[:,0]))[0]
    xmax = numpy.unique(numpy.max(macro_nodes[:,0]))[0]
    ymin = numpy.unique(numpy.min(macro_nodes[:,1]))[0]
    ymax = numpy.unique(numpy.max(macro_nodes[:,1]))[0]
    zmin = numpy.unique(numpy.min(macro_nodes[:,2]))[0]
    zmax = numpy.unique(numpy.max(macro_nodes[:,2]))[0]
    
    # create random points within each range
    num = 500
    offset = 0.2
    xs = numpy.random.uniform(low=xmin-offset, high=xmax+offset, size=(num, ))
    ys = numpy.random.uniform(low=ymin-offset, high=ymax+offset, size=(num, ))
    zs = numpy.random.uniform(low=zmin-offset, high=zmax+offset, size=(num, ))
    data = numpy.transpose(numpy.vstack([xs, ys, zs]))
    test_node_dict = {}
    for i, n in enumerate(data):
        test_node_dict[i] = n

    # find the set
    test_set_dict = find_major_sets(test_node_dict, major_plane_dict, 1)
    
    # plot the results for selected elements
    if plot:
        plot_outputs(test_node_dict, macro_node_dict, macro_elem_dict, test_set_dict, [0], 'test_cube_major', bounding=True)
    
    return test_set_dict

def run_test_minor(minor_plane_dict, macro_node_dict, macro_nodes, macro_elem_dict, major_plane_dict, plot=False):

    # find the extent of the nodal coordinates
    xmin = numpy.unique(numpy.min(macro_nodes[:,0]))[0]
    xmax = numpy.unique(numpy.max(macro_nodes[:,0]))[0]
    ymin = numpy.unique(numpy.min(macro_nodes[:,1]))[0]
    ymax = numpy.unique(numpy.max(macro_nodes[:,1]))[0]
    zmin = numpy.unique(numpy.min(macro_nodes[:,2]))[0]
    zmax = numpy.unique(numpy.max(macro_nodes[:,2]))[0]
    
    # create random points within each range
    num = 500
    offset = 0.2
    xs = numpy.random.uniform(low=xmin-offset, high=xmax+offset, size=(num, ))
    ys = numpy.random.uniform(low=ymin-offset, high=ymax+offset, size=(num, ))
    zs = numpy.random.uniform(low=zmin-offset, high=zmax+offset, size=(num, ))
    xs = numpy.linspace(0.0, 1.0, 20)
    ys = numpy.linspace(0.0, 1.0, 20)
    zs = numpy.linspace(0.0, 1.0, 20)
    data = numpy.transpose(numpy.vstack([xs, ys, zs]))
    test_node_dict = {}
    for i, n in enumerate(data):
        test_node_dict[i] = n

    # find the sets
    major_set_dict = find_major_sets(test_node_dict, major_plane_dict, 2)
    test_set_dict = find_minor_sets(test_node_dict, minor_plane_dict, major_set_dict)
    print('minor test set dict')
    print(test_set_dict)
    
    # plot the results for selected elements
    if plot:
        plot_list = [key for key in test_set_dict.keys()]
        plot_outputs(test_node_dict, macro_node_dict, macro_elem_dict, test_set_dict, plot_list, 'test_cube_minor')
    
    return test_set_dict


def plot_outputs(node_dict, macro_node_dict, macro_elem_dict, set_dict, selection, output_name, bounding=False):

    # plot
    ## make the box
    for select in selection:
        print(f'Plotting Set {select}!')
        set_ids = set_dict[select]
        # create bounding box if requested
        if bounding:
            vertices = [macro_node_dict[n] for n in macro_elem_dict[select]]
            #print(f'vertices = {vertices}')
            faces = [[0, 1, 5, 4],
                     [1, 2, 6, 5],
                     [2, 3, 7, 6],
                     [3, 0, 4, 7],
                     [0, 1, 2, 3],
                     [4, 5, 6, 7]]
            poly3d = [[vertices[faces[ix][iy]] for iy in range(len(faces[0]))] for ix in range(len(faces))]
        
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d', proj_type='ortho')
        for i in node_dict.keys():
            n = node_dict[i]
            x, y, z = n[0], n[1], n[2]
            if i in set_ids:
                color = 'blue'
                alpha = 1.0
            else:
                color = 'red'
                alpha = 0.1
            ax.scatter(x, y, z, color=color, alpha=alpha)
        # plot bounding box if requested
        if bounding:
            ax.add_collection3d(Poly3DCollection(poly3d, facecolors='g', alpha=0.1))
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        views = [[0, 0, 'x_plane'], [90, 0, 'y_plane'], [0, 90, 'z_plane'], [45, 45, 'iso']]
        for view in views:
            ax.view_init(elev=view[1], azim=view[0])
            plt.savefig(f'{output_name}_{select}_{view[2]}.png', dpi=300)

    return 0

def create_DNS_node_dict(DNS_file):

    data, geometry, topology = XRT.parse_xdmf_output(DNS_file)
    
    coordsx = data['COORDx'][0]
    coordsy = data['COORDy'][0]
    coordsz = data['COORDz'][0]
    nodeid = data['NODEID'][0]
    
    DNS_node_dict = {}
    
    for n, x, y, z in zip(nodeid, coordsx, coordsy, coordsz):
        DNS_node_dict[n] = numpy.array([x, y, z])

    return DNS_node_dict


def write_out_xdmf(DNS_file, output_xdmf, minor_set_dict, output_hdf5):

    # copy original xdmf file which will include newly detected sets
    shutil.copy(DNS_file, output_xdmf)
    
    # modify the new xdmf file to include sets with reference to new hdf5 file
    tree=ET.parse(output_xdmf)
    root=tree.getroot()
    for grid in root[0][1][1:]:
        for set_name in minor_set_dict.keys():
            set_info = {'Name':set_name, 
                        'Type':'Node'}
            data_size = numpy.shape(minor_set_dict[set_name])
            set_data_info = {'DataType':'UInt',
                             'Dimensions':str(data_size),
                             'Format':'HDF',
                             'Precision':'4'}
            set_object = ET.SubElement(grid, "Set")
            for key in set_info:
                set_object.set(key, set_info[key])
                data_item = ET.SubElement(set_object, "DataItem")
                for info_key in set_data_info.keys():
                    data_item.set(info_key, set_data_info[info_key])
                data_item.text = f'{output_hdf5}:{set_name}'
    
    # Dump the file
    with open(output_xdmf, 'w') as xml_file:
        xml_file.write(ET.tostring(root, encoding='utf8').decode('utf-8'))
        
    # Pretty-print using lxml
    pretty_xml = lxml.etree.tostring(lxml.etree.parse(output_xdmf).getroot(),
        pretty_print=True, encoding="unicode")
    with open(output_xdmf, 'w') as xml_file:
        xml_file.write(pretty_xml)

    return 0


def main(mesh_file, DNS_file=None, output_hdf5=None, output_xdmf=None, test_major=None, test_minor=None, downsample=None, plot_selection=None):

    # get nodes and element connectivity
    macro_node_dict, macro_elem_dict, macro_nodes = load_mesh(mesh_file)

    # create nested dictionary of planes for each element
    print('Finding Plane information')
    major_plane_dict = build_major_planes(macro_node_dict, macro_elem_dict)
    #print(major_plane_dict)

    # test if required
    if test_major:
        print('Creating major test sets')
        test_set_dict = run_test_major(major_plane_dict, macro_node_dict, macro_nodes, macro_elem_dict, plot=True)
    # test if required
    if test_minor:
        # print('Creating major test sets')
        # test_major_set_dict = run_test_major(major_plane_dict, macro_node_dict, macro_nodes, macro_elem_dict, plot=False)
        print('Finding minor plane information!')
        minor_plane_dict = build_minor_planes(macro_node_dict, macro_elem_dict)
        print(minor_plane_dict)
        print('Creating minor test sets')
        test_minor_set_dict = run_test_minor(minor_plane_dict, macro_node_dict, macro_nodes, macro_elem_dict, major_plane_dict, plot=True)

    # read in xdmf file and create DNS_node_dict
    if DNS_file:
        print('Loading DNS information')
        DNS_node_dict = create_DNS_node_dict(DNS_file)

        # find points within each element
        print('Finding major sets!')
        major_set_dict = find_major_sets(DNS_node_dict, major_plane_dict, 2)
        
        # Create nested dictionary of minor planes for each element
        print('Finding minor plane information!')
        minor_plane_dict = build_minor_planes(macro_node_dict, macro_elem_dict)
        
        # Find  octant sets of each element
        print('Finding minor sets!')
        minor_set_dict = find_minor_sets(DNS_node_dict, minor_plane_dict, major_set_dict)
        
        # plot if selection provided:
        if plot_selection:
            if downsample:
                plot_node_dict = {}
                for n in DNS_node_dict.keys():
                    if n % 10 == 0:
                        plot_node_dict[n] = DNS_node_dict[n]
            else:
                plot_node_dict = DNS_node_dict
            #print(major_minor_set_dict)
            selection = [int(i) for i in plot_selection]
            plot_suffix = output_xdmf.split('.')[0]
            plot_outputs(plot_node_dict, macro_node_dict, macro_elem_dict, major_set_dict, selection, f'test_{plot_suffix}', bounding=True)
        
        # Write out hdf5 file
        if output_hdf5:
            print('Writing output hdf5 file!')
            with h5py.File(output_hdf5, 'w') as hdf_out:
                for set_name in minor_set_dict.keys():
                    hdf_out.create_dataset(set_name, data=minor_set_dict[set_name])
            
            # Write out xdmf file
            if output_xdmf:
                print('Writing output xdmf file!')
                write_out_xdmf(DNS_file, output_xdmf, minor_set_dict, output_hdf5)

    return 0


def get_parser():

    basename = pathlib.Path(__file__).name
    cli_description = "Prep Abaqus DNS for filter"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=basename)
    parser.add_argument('-m', '--mesh-file', type=str, required=True,
        help='The name of the input mesh xdmf file')
    parser.add_argument('-d', '--dns-file', type=str, required=False,
        help='The name of the input xdmf file of DNS results')
    parser.add_argument('--output-hdf5', type=str, required=False,
        help='The name of the output hdf5 file containing the set ids of the newly detected sets')
    parser.add_argument('--output-xdmf', type=str, required=False,
        help='The name of the output xdmf file pointing to two hdf5 files, original copied from ``dns-file``')
    parser.add_argument('--test-major', type=str, required=False,
        help='An optional argument for testing the detection of major plane sets')
    parser.add_argument('--test-minor', type=str, required=False,
        help='An optional argument for testing the detection of minor plane sets')
    parser.add_argument('--downsample', type=str, required=False,
        help='An optional argument for downsampling plotting output')
    parser.add_argument('-s', '--plot-selection', nargs="+", required=False,
        help='A sequence of integers indicating which macro elements to plot')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(main(mesh_file=args.mesh_file,
                  DNS_file=args.dns_file,
                  output_hdf5=args.output_hdf5,
                  output_xdmf=args.output_xdmf,
                  test_major=args.test_major,
                  test_minor=args.test_minor,
                  downsample=args.downsample,
                  plot_selection=args.plot_selection,
                  ))