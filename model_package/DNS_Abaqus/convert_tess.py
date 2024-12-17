#! /usr/bin/env python

import sys
import argparse
import pathlib
import re

import numpy
import cubit
import pandas


def parse_contents(input_file):
    
    lines = []
    vertex_index = []
    edge_index = []
    face_index = []
    bound_index = []
    polyhedron_index = []

    with open(input_file, 'r') as f:
        input_contents = numpy.asarray(f.read().splitlines())

    for index, text in enumerate(input_contents):
        if '**vertex' in text:
            vertex_index.append(index)
        if '**edge' in text:
            edge_index.append(index)
        if '**face' in text:
            face_index.append(index)
        if ' *face' in text:
            bound_index.append(index)
        if '**polyhedron' in text:
            polyhedron_index.append(index)

    faces, bounds, polyhedra = {}, {}, {}

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

    # collect polyhedra
    poly_start = polyhedron_index[0]+2
    num_poly = int(input_contents[polyhedron_index[0]+1])
    for line in input_contents[poly_start:poly_start+num_poly]:
        items = re.split(r'\s+', line)
        polyhedra[int(items[1])] = [str(abs(int(i))) for i in items[3:]]

    return vertices, edges, faces, bounds, polyhedra


def create_geometry(vertices, edges, faces, bounds, polyhedra, stl_file):

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

    # create all the polyhedra
    num_faces = len(faces.keys())
    for cell in polyhedra.keys():
        cubit.cmd(f'create volume surface {" ".join(polyhedra[cell])} heal keep')
        cubit.cmd(f'volume {cell+num_faces} rename "cell {cell}"')

    # define 

    cubit.cmd(f'save as "{stl_file}.cub" overwrite')

    cubit.cmd(f'export stl "{stl_file}.stl" volume {" ".join(str(cell+num_faces) for cell in polyhedra.keys())} fast overwrite')

    return 0


def create_mesh(stl_file, mesh_file, polyhedra, sidesets):

    cubit.init(['cubit', '-noecho', '-nojournal', '-nographics', '-batch'])
    cubit.cmd('new')
    cubit.cmd('reset')

    cubit.cmd(f'import stl "{stl_file}.stl" feature_angle 135.0 merge')
    #cubit.cmd('set developer commands on')
    cubit.cmd('imprint all')
    #cubit.cmd('merge all')
    #cubit.cmd('volume all scheme Polyhedron')
    cubit.cmd('volume all scheme tetmesh')
    cubit.cmd('mesh volume all')
    for cell in polyhedra.keys():
        cubit.cmd(f'block {cell} volume {cell}')

    # new using min-max
    sideset_num = 1
    for sideset in sidesets.keys():
        cubit.cmd(f'nodeset {sideset_num} add surface with {sidesets[sideset]}')
        cubit.cmd(f'nodeset {sideset_num} name "{sideset}"')
        cubit.cmd(f'sideset {sideset_num} add surface with {sidesets[sideset]}')
        cubit.cmd(f'sideset {sideset_num} name "{sideset}"')
        sideset_num = sideset_num + 1

    cubit.cmd(f'save as "{mesh_file}.cub" overwrite')
    blocks = [str(cell) for cell in polyhedra.keys()]
    cubit.cmd(f'export abaqus "{mesh_file}.inp" block {" ".join(blocks)} overwrite partial')

    return 0


def convert_tess(input_file, stl_file=None, mesh_file=None, seed_size=1.0):

    vertices, edges, faces, bounds, polyhedra = parse_contents(input_file)
    print(vertices)
    print(edges)
    print('face')
    print(faces)
    print('bounds')
    print(bounds)
    print('polyhedra')
    print(polyhedra)

    sidesets = {'top': 'z_min > 0.99',
                'bottom': 'z_max < 0.01',
                'side1': 'x_min > 0.99',
                'side2': 'y_min > 0.99',
                'side3': 'x_max < 0.01',
                'side4': 'y_max < 0.01',}

    if stl_file:
        create_geometry(vertices, edges, faces, bounds, polyhedra, stl_file)
        if mesh_file:
            create_mesh(stl_file, mesh_file, polyhedra, sidesets)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)
    prog = f"python {script_name.name} "
    cli_description = "Convert a tesslation file output by Neper to STL and create Abaqus mesh"

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