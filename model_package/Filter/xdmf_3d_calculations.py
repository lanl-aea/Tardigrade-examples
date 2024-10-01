# Imports
import sys
import os
import inspect
import argparse

import numpy

sys.path.append(r'/projects/tea/PSAAP/tardigrade_filter/src/python')

import file_io.xdmf


def p_q_calculation(xdmf_file, qp_field_name, inc, num_elements):

    p_output = numpy.zeros([1, 1, num_elements, 1])
    q_output = numpy.zeros([1, 1, num_elements, 1])

    data = xdmf_file.getIncrementData(inc, qp_field_name)

    for e in range(num_elements):
        p = (-1/3)*numpy.trace(data[0][0][inc,:])
        dev = numpy.trace(data[0][0][inc,:] + p
        q = numpy.sqrt(numpy.tensordot(dev, dev))
        p_output[0,0,e,:] = p
        q_output[0,0,e,:] = q


    return p_output, q_output

def xdmf_3d_calculations(input_file, num_elements):


    xdmf_file = file_io.xdmf.XDMF(input_file)
    xdmf_file.open()

    num_increments = xdmf_file.getNumIncrements()
    incs = list(range(0, num_increments))
    times = [xdmf_file.getIncrementTime(i)[0] for i in incs]

    # calculate p and q
    for i in incs:
        p_out, q_out = p_q_calculations(xdmf_file, f'cauchy_stress_0', num_elements)


    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Convert Ratel DNS results to XDMF format"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-o', '--input-file', type=str,
        help='Specify the input filename for the h5 + XDMF file pair')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(xdmf_3d_calculations(input_file=args.input_file,
                                 ))