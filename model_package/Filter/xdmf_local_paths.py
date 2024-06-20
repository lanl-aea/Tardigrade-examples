import os
import sys
import argparse
import inspect


def replace_paths(input_file, output_file, oldpath, newpath):
    '''Create a copy of an XDMF file with absolute H5 paths replaced with relative paths

    :param str input_file: The XDMF file output by the Micromorphic Filter with absolute H5 paths
    :param str output_file: The new XDMF file with relative H5 paths
    :param str oldpath: The absolute path to be replaced by ``newpath``
    :param str newpath: The relative path to replace ``oldpath``

    :return: Write ``output_file``
    '''

    # Collect XDMF info and insert line breaks
    data = []
    with open(input_file) as file:
        for line in file:
            data.append(str(line).replace('><','>\n<'))

    # Write out temporary file with line breaks
    temp_file = f"{input_file.split('.')[0]}_temp.xdmf"
    temp = open(temp_file, 'wt')
    for line in data:
        temp.write(line)
    temp.close()

    # Read in new file and replace absolute path with relative path
    out_data = []
    with open(temp_file) as temp:
        for line in temp:
            out_data.append(str(line).replace(oldpath, newpath))

    # Write newly formatted input file
    out = open(output_file, 'wt')
    for line in out_data:
        out.write(line)
    out.close()

    # Delete temporary file
    os.remove(temp_file)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Create a copy of an XDMF file with absolute H5 paths replaced with relative paths"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str, required=True,
                        help="The XDMF file output by the Micromorphic Filter with absolute H5 paths")
    parser.add_argument('-o', '--output-file', type=str, required=True,
                        help="The new XDMF file with relative H5 paths")
    parser.add_argument('--oldpath', type=str, required=True,
                        help="The absolute path to be replaced by ``--newpath``")
    parser.add_argument('--newpath', type=str, required=True,
                        help="The relative path to replace ``--oldpath``")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(replace_paths(input_file=args.input_file,
                           output_file=args.output_file,
                           oldpath=args.oldpath,
                           newpath=args.newpath,
                           ))