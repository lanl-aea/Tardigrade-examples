# Imports
import sys
import os
import inspect
import argparse

import numpy

sys.path.append(r'/projects/tea/PSAAP/tardigrade_filter/src/python')

import file_io.xdmf


def p_q_calculation(xdmf_file, qp_field_name, inc, num_elements):
    p_output = numpy.zeros([num_elements, 1])
    q_output = numpy.zeros([num_elements, 1])
    data = xdmf_file.getIncrementData(inc, qp_field_name)
    for e in range(num_elements):
        chunk = data[0][0][e,:].reshape(3,3)
        p = (-1/3)*numpy.trace(chunk)
        dev = chunk + p*numpy.eye(3)
        q = numpy.sqrt(numpy.tensordot(dev, dev))
        p_output[e,:] = p
        q_output[e,:] = q
    return p_output, q_output


def average_over_quadrature(xdmf_file_in, qp_field_name, inc, num_elements):

    p_all = numpy.zeros([8, num_elements, 1])
    q_all = numpy.zeros([8, num_elements, 1])
    for qp in range(8):
        ps, qs = p_q_calculation(xdmf_file_in, f'{qp_field_name}_{qp}', inc, num_elements)
        p_all[qp, :, :] = ps
        q_all[qp, :, :] = qs

    return numpy.mean(p_all, axis=0), numpy.mean(q_all, axis=0)


def xdmf_3d_calculations(input_file, output_file, num_elements):


    xdmf_file_in = file_io.xdmf.XDMF(input_file)
    xdmf_file_in.open()

    xdmf_file_out = file_io.xdmf.XDMF(output_filename=output_file)

    num_increments = xdmf_file_in.getNumIncrements()
    incs = list(range(0, num_increments))
    times = [xdmf_file_in.getIncrementTime(i)[0] for i in incs]

    reference_positions = xdmf_file_in.getIncrementReferenceNodePositions(0)[0][0]
    connectivity = xdmf_file_in.getIncrementConnectivity(0)[0][0].reshape((1,-1))

    # calculate p and q
    for i in incs:
        grid = xdmf_file_out.addGrid(xdmf_file_out.output_timegrid, {})
        xdmf_file_out.addTime(grid, times[i])
        xdmf_file_out.addPoints(grid, reference_positions,
                                duplicate='filter_reference_positions')
        xdmf_file_out.addConnectivity(grid, "HEXAHEDRON", connectivity,
                                duplicate='filter_connectivity')

        p_sig, q_sig = average_over_quadrature(
            xdmf_file_in, f'cauchy_stress', i, num_elements)
        xdmf_file_out.addData(grid, f"test_p_sig", p_sig, center='Cell')
        xdmf_file_out.addData(grid, f"test_q_sig", q_sig, center='Cell')
        p_sym, q_sym = average_over_quadrature(
            xdmf_file_in, f'symmetric_micro_stress', i, num_elements)
        xdmf_file_out.addData(grid, f"test_p_sym", p_sym, center='Cell')
        xdmf_file_out.addData(grid, f"test_q_sym", q_sym, center='Cell')

    xdmf_file_out.write()
    print("XDMF file written!")

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Convert Ratel DNS results to XDMF format"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
        help='Specify the input filename for the h5 + XDMF file pair')
    parser.add_argument('-o', '--output-file', type=str,
        help='Specify the output filenmae for the h5 + XDMF file pair')
    parser.add_argument('--num-elements', type=int,
        help='The number of macroscale elements')


    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(xdmf_3d_calculations(input_file=args.input_file,
                                  output_file=args.output_file,
                                  num_elements=args.num_elements,
                                  ))