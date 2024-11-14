#!python
import sys
import os
import argparse
import pathlib

import cubit
import numpy
import meshio
import pandas
import scipy


def adjust_centroids(centroids, method='absolute'):

    if method == 'origin':
        # findest closest point to (0,0,0)
        distances = numpy.linalg.norm(centroids, axis=1)
        minx, miny, minz = centroids[numpy.argmin(distances)]
    elif method == 'absolute':
        minx = numpy.min(centroids[:,0])
        miny = numpy.min(centroids[:,1])
        minz = numpy.min(centroids[:,2])

    print(f'(minx, miny, minz) = ({minx}, {miny}, {minz})')

    centroids = pandas.DataFrame(centroids)
    centroids[0] = centroids[0] - minx
    centroids[1] = centroids[1] - miny
    centroids[2] = centroids[2] - minz

    return centroids


def get_element_centroids_from_exouds_mesh(exodus_mesh_map):

    mesh=meshio.read(exodus_mesh_map)
    node_points = mesh.points

    num_blocks = len(mesh.cells)
    connect = []
    for i in range(num_blocks):
        connect.append(mesh.cells[i].data)
    connect = numpy.array(connect)
    tup = numpy.shape(connect)
    connectivity = connect.reshape((tup[0]*tup[1],tup[2]))

    centroids = numpy.array([numpy.mean([node_points[j] for j in connectivity[i]],axis=0)
                             for i in range(0, numpy.shape(connectivity)[0])])

    centroids = adjust_centroids(centroids)

    return centroids


def add_element_blocks_to_mesh(input_mesh, output_mesh, elements, number_existing_blocks=0, exodus_mesh_map=None):
    '''Take an existing exodus mesh, add element blocks for each element, save with new name

    :param str input_mesh: The input exodus mesh file to modify
    :param str output_mesh: The output exodus mesh file with block names defined
    :param int elements: The number of elements in the mesh for which to define a block name

    :returns: Write ``output_mesh``
    '''

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')
    cubit.cmd(f'import mesh geometry "{input_mesh}"')

    # Delete any existing blocks
    if number_existing_blocks == 0:
        cubit.cmd('reset block')
        j = 1
    else:
        cubit.cmd('delete block 1')
        j = number_existing_blocks + 1

    element_list = list(range(0, elements))
    # Match centroids from filter macroscale mesh using KDTree
    if exodus_mesh_map:
        # get macroscale centroids
        macro_centroids = get_element_centroids_from_exouds_mesh(exodus_mesh_map)
        print(f'macro_centroids = \n {macro_centroids}')
        # Collect hex centroids
        hex_centroids = []
        for i in element_list:
            text = cubit.get_center_point('hex', i+1)
            hex_centroids.append(list(text))
        hex_centroids=numpy.array(hex_centroids)
        hex_centroids = adjust_centroids(hex_centroids)
        print(f'hex_centroids = \n {hex_centroids}')
        # Use macro centroids as the KDTree
        centroids_tree = scipy.spatial.KDTree(macro_centroids)
        dist, ids = centroids_tree.query(hex_centroids)
        sorted_ids = numpy.argsort(ids)
        print(f'dist = \n\t {dist}\nids = \n\t {ids}')
        if len(numpy.unique(sorted_ids)) != len(sorted_ids):
            print('Uh oh! There is not a unique map between meshes!')
        # assign blocks to identified elements
        for i, id in zip(element_list, sorted_ids):
            cubit.cmd(f'block {i+j} add hex {id+1}')
            cubit.cmd(f'block {i+j} name "element_{i}"')
    # Assign directly
    else:
        for i in element_list:
            cubit.cmd(f'block {i+j} add hex {i+1}')
            cubit.cmd(f'block {i+j} name "element_{i}"')

    cubit.cmd(f'export mesh "{output_mesh}"  overwrite')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "

    cli_description = "Take an existing exodus mesh, add element blocks for each element, save with new name"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=prog)
    parser.add_argument('--input-mesh', type=str, required=True,
        help="The input exodus mesh file to modify")
    parser.add_argument('--output-mesh', type=str, required=True,
        help="The output exodus mesh file with block names defined")
    parser.add_argument('--elements', type=int, required=True,
        help="The number of elements in the mesh for which to define a block name")
    parser.add_argument('--number-existing-blocks', type=int, required=False, default=0,
        help="The number of existing mesh blocks to keep in final mesh")
    parser.add_argument('--exodus-mesh-map', type=str, required=False, default=None,
        help="An existing macroscale to optionally map centroids for element block numbering")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    sys.exit(add_element_blocks_to_mesh(input_mesh=args.input_mesh,
                                        output_mesh=args.output_mesh,
                                        elements=args.elements,
                                        number_existing_blocks=args.number_existing_blocks,
                                        exodus_mesh_map=args.exodus_mesh_map,
                                        ))
