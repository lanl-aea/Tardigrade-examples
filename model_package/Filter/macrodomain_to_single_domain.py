import inspect
import sys
import os
import argparse

import meshio
import numpy
import subprocess
import shutil


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
            
    return(node_dict, elem_dict, numpy.array(nodes))


def create_new_mesh(node_dict, elem_dict, element, out_name):

    element_nodes = elem_dict[element]
    
    output_nodes = []
    for node in element_nodes:
        output_nodes.append(node_dict[node])
        
    output_nodes = numpy.array(output_nodes)
    print('testing output nodes')
    print(output_nodes)
    
    out_elem = list(range(0,8))
    cells = [('hexahedron', [out_elem])]
    
    out_mesh = meshio.Mesh(output_nodes, cells)
    out_mesh.write(out_name)
    #meshio.write(f'{output}_{element}.xdmf', out_mesh, binary=False)
    
    # convert mesh to ascii using subprocess
    subprocess.run(['ls'])
    subprocess.run([f'meshio ascii {out_name}'], shell=True)
    
    return out_name


def modify_mesh_for_macro(new_mesh_file, new_output_name):

    shutil.copy(new_mesh_file, new_output_name)

    string_1 = '<Xdmf Version="3.0">'\
                '<Domain>'\
                '<Information Name="Domain" Value="Primary data structure from a MOOSE simulation"/>'\
                '<Grid CollectionType="Temporal" GridType="Collection" Name="Collection">'\
                '<Information Name="Collection_0" Value="The collection of temporal grids"/>'\
                '<Grid Name="Grid">'\
                '<Time Value="0">'\
                '<Information Name="Time" Value="This is the current value of the timestep"/>'\
                '</Time>'\
                '<Geometry GeometryType="XYZ">'\
                '<DataItem DataType="Float" Dimensions="8 3" Format="XML" Precision="8">'
    
    string_2 = '</DataItem>'\
               '</Topology>'\
               '<Set Name="macro_domain" Type="Node">'\
               '<DataItem DataType="UInt" Dimensions="8" Format="XML" Precision="4">0 1 2 3 4 5 6 7</DataItem>'\
               '</Set>'\
               '<Attribute Center="Node" ElementCell="" ElementDegree="0" ElementFamily="" ItemType="" Name="NODEID" Type="GlobalID">'\
               '<Information Name="ID" Value="The nodal IDs"/>'\
               '<DataItem DataType="UInt" Dimensions="8" Format="XML" Precision="4">0 1 2 3 4 5 6 7</DataItem>'\
               '</Attribute>'\
               '<Attribute Center="Cell" ElementCell="" ElementDegree="0" ElementFamily="" ItemType="" Name="ELEMID" Type="GlobalId">'\
               '<Information Name="ID" Value="The element IDs"/>'\
               '<DataItem DataType="UInt" Dimensions="1" Format="XML" Precision="4">0</DataItem>'\
               '</Attribute>'\
               '</Grid>'\
               '</Grid>'\
               '</Domain>'\
               '</Xdmf>'\

    with open(new_mesh_file, 'r') as xdmf_file:
        all_input = xdmf_file.readlines()
        #all_input[0].replace('<Grid Name="Grid">', string_1)
        all_input[0] = string_1
        all_input[-1] = string_2

    with open(new_output_name, 'w') as output:
        output.writelines(all_input)
    
    return 0
        

def main(mesh_file, output, element):

    node_dict, elem_dict, nodes = load_mesh(mesh_file)
    
    # write new mesh
    output_name = f'{output}_temp.xdmf'
    base_mesh_file = create_new_mesh(node_dict, elem_dict, element, output_name)
    
    # modify mesh for micromorphic filter
    new_output_name = f'{output}.xdmf'
    modify_mesh_for_macro(base_mesh_file, new_output_name)

    # delete temporary file
    subprocess.run([f'rm {output_name}'], shell=True)
    subprocess.run([f'rm {output_name.split(".")[0]}.h5'], shell=True)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Prep Abaqus DNS for filter"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-m', '--mesh-file', type=str, required=True,
        help='The name of the input mesh xdmf file')
    parser.add_argument('--output', type=str, required=True,
        help='The basename of the output xdmf files')
    parser.add_argument('--element', type=int, required=True,
        help='The element to build a domain from')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(main(mesh_file=args.mesh_file,
                  output=args.output,
                  element=args.element,
                  ))