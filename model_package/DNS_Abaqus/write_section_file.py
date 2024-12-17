import os
import sys
import argparse
import inspect


def write_section_file(output_file, number_grains, material_name):
    '''Write an Abaqus input file for the section definition of grains
    
    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param int number_grains: The number of sections to create corresponding to unqieu grains
    :param str material_name: The name of the material to assign to sections

    :returns: ``output_file``
    '''

    # Write input file
    with open(output_file, 'w') as f:
        for i in range(1, number_grains+1):
            f.write(f'*SOLID SECTION, ELSET=EB{i}, MATERIAL={material_name}\n')

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Write an Abaqus input file for the section definition of grains"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="The name of the Abaqus input file to write")
    parser.add_argument('--number-grains', type=int, required=True,
        help="The number of sections to create corresponding to unqieu grains")
    parser.add_argument('--material-name', type=str, required=True,
        help="The name of the material to assign to sections")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(write_section_file(output_file=args.output_file,
                                number_grains=args.number_grains,
                                material_name=args.material_name,
                                ))