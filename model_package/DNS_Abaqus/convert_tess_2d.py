#! /usr/bin/env python

import sys
import argparse
import pathlib
import re

sys.path.append('/apps/Cubit-16.12/bin')

import numpy
import cubit
import pandas


def parse_contents(input_file):
    """Parse content of Neper output

    :params str input_file: Input tesselation (.tess) file

    :returns: DataFrame containing rows of vertex coordinates,
        DataFrame containing rows of edge-to-vertex connectivity, and
        dictionary containing ids of faces with lists of face-to-edge connectivity
    """

    vertex_index = []
    edge_index = []
    face_index = []

    with open(input_file, 'r') as f:
        input_contents = numpy.asarray(f.read().splitlines())

    for index, text in enumerate(input_contents):
        if '**vertex' in text:
            vertex_index.append(index)
        if '**edge' in text:
            edge_index.append(index)
        if '**face' in text:
            face_index.append(index)

    faces = {}

    # collect vertices
    vertex_start = vertex_index[0]+2
    num_vertices = int(input_contents[vertex_index[0]+1])
    vertices = pandas.read_csv(input_file, skiprows=vertex_start, nrows=num_vertices,names=['index','x','y','z','blank'],delimiter=r"\s+", index_col=0)

    # collect edges
    edge_start = edge_index[0]+2
    num_edges = int(input_contents[edge_index[0]+1])
    edges = pandas.read_csv(input_file, skiprows=edge_start, nrows=num_edges,names=['index','a','b','blank'],delimiter=r"\s+", index_col=0)

    # collect faces
    face_start = face_index[0]+2
    num_faces = int(input_contents[face_index[0]+1])
    desired_lines = list(input_contents[face_start:face_start+num_faces*4])
    line_to_grab = 1
    for face_id in range(0,num_faces):
        line = desired_lines[face_id*4 + line_to_grab]
        items = re.split(r'\s+', line)
        faces[face_id+1] = [str(abs(int(i))) for i in items[2:]]

    return vertices, edges, faces


def create_geometry(vertices, edges, faces, stl_file):
    """Create geometry from Neper output

    :params DataFrame vertices: Pandas DataFrame containing rows of vertex coordinates
    :params DataFrame edges: Pandas DataFrame containing rows of edge-to-vertex connectivity
    :params dict faces: Dictionary containing ids of faces with lists of face-to-edge connectivity
    :params str stl_file: Optional filename to save STL geometry without extension

    :returns: Write ``{stl_file}.cub``
    """

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')

    # create all the points
    for index, row in vertices.iterrows():
        cubit.cmd(f'create vertex {row["x"]} {row["y"]} {row["z"]}')

    # create all the edges
    for index, row in edges.iterrows():
        cubit.cmd(f'create curve vertex {row["a"]} {row["b"]}')

    # create all the faces
    for face in faces.keys():
        cubit.cmd(f'create surface curve {" ".join(faces[face])}')

    # define 

    cubit.cmd(f'save as "{stl_file}.cub" overwrite')

    #cubit.cmd(f'export stl "{stl_file}.stl" volume {" ".join(str(cell+num_faces) for cell in polyhedra.keys())} fast overwrite')

    return 0


def create_mesh(stl_file, mesh_file, faces, sidesets, seed_size):
    """Create mesh from Neper output

    :params str stl_file: Optional filename to save STL geometry without extension
    :params str mesh_file: Optional filename to create an Abaqus mesh without extension
    :params dict sidesets: Dictionary containing names of surfaces and geometric search strings
    :params float seed_size: The approximate mesh size

    :returns: Write ``{mesh_file}.cub`` and ``{mesh_file}.inp``
    """

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')

    cubit.cmd(f'open "{stl_file}.cub"')
    #cubit.cmd('set developer commands on')
    cubit.cmd('imprint all')
    cubit.cmd('surface all scheme trimesh')
    if seed_size < 1.0:
        cubit.cmd(f'surface all size {seed_size}')
    cubit.cmd('mesh surface all')

    for cell in faces.keys():
        cubit.cmd(f'block {cell} surface {cell}')
    cubit.cmd('block all element type TRI3')

    # new using min-max
    sideset_num = 1
    for sideset in sidesets.keys():
        cubit.cmd(f'nodeset {sideset_num} add curve with {sidesets[sideset]}')
        cubit.cmd(f'nodeset {sideset_num} name "{sideset}"')
        cubit.cmd(f'sideset {sideset_num} add curve with {sidesets[sideset]}')
        cubit.cmd(f'sideset {sideset_num} name "{sideset}"')
        sideset_num = sideset_num + 1

    # Corner
    cubit.cmd(f'nodeset {sideset_num} add vertex with x_max < 0.01 with y_max < 0.01')
    cubit.cmd(f'nodeset {sideset_num} name "corner"')

    # Output
    cubit.cmd(f'save as "{mesh_file}.cub" overwrite')
    blocks = [str(cell) for cell in faces.keys()]
    cubit.cmd(f'export abaqus "{mesh_file}.inp" block {" ".join(blocks)} overwrite partial dimension 2')
    #cubit.cmd(f'export abaqus "{mesh_file}.inp" overwrite partial dimension 2')

    return 0


def convert_tess(input_file, stl_file=None, mesh_file=None, seed_size=1.0):
    """Convert a 2D tesselation file output by Neper to STL and create an Abaqus mesh

    :params str input_file: Input tesselation (.tess) file
    :params str stl_file: Optional filename to save STL geometry without extension
    :params str mesh_file: Optional filename to create an Abaqus mesh without extension
    :params float seed_size: The approximate mesh size
    """

    vertices, edges, faces = parse_contents(input_file)
    print(vertices)
    print(edges)
    print('face')
    print(faces)

    sidesets = {'side1': 'x_min > 0.99',
                'side2': 'y_min > 0.99',
                'side3': 'x_max < 0.01',
                'side4': 'y_max < 0.01',}

    if stl_file:
        create_geometry(vertices, edges, faces, stl_file)
        if mesh_file:
            create_mesh(stl_file, mesh_file, faces, sidesets, seed_size)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Convert a 2D tesselation file output by Neper to STL and create an Abaqus mesh"

    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument("--input-file", type=str, required=True,
        help="Input tesselation (.tess) file")
    parser.add_argument("--stl-file", type=str, required=False, default=None,
        help="Optional filename to save STL of geometry without extension")
    parser.add_argument("--mesh-file", type=str, required=False, default=None,
        help="Optional filename to create an Abaqus mesh without extension")
    parser.add_argument('--seed-size', type=float, required=False, default=1.0,
        help='The approximate mesh size')

    return parser


if __name__ == "__main__":
    parser = get_parser()
    args, unknown = parser.parse_known_args()
    sys.exit(convert_tess(input_file=args.input_file,
                          stl_file=args.stl_file,
                          mesh_file=args.mesh_file,
                          seed_size=args.seed_size,
                          ))
