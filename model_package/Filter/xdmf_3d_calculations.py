# Imports
import sys
import os
import inspect
import argparse

import numpy
import pandas
sys.path.append('/projects/tea/damage_tardigrade_filter/tardigrade_filter/src/python')
import file_io.xdmf


def two_pressure(xdmf_file_in, qp_field_name, inc, num_elements):
    '''Calculate the pressure of a second order tensor

    :params object xdmf_file_in: The XDMF file to calculate quantities from
    :params str qp_field_name: The stress field (defined at the quadrature points) to consider
    :params list inc: The current time increment
    :params int num_elements: The number of macroscale elements
    :params int k: The 3rd index of the third order tensor
    '''

    data = xdmf_file_in.getIncrementData(inc, qp_field_name)
    p_output = numpy.zeros([num_elements, 1])

    for e in range(num_elements):
        chunk = data[0][0][e,:].reshape(3,3)
        p_output[e,:] = (-1/3)*numpy.trace(chunk)

    return p_output


def two_devnorm(xdmf_file_in, qp_field_name, inc, num_elements):
    '''Calculate the norm of the deviatoric part of a second order tensor

    :params object xdmf_file_in: The XDMF file to calculate quantities from
    :params str qp_field_name: The stress field (defined at the quadrature points) to consider
    :params list inc: The current time increment
    :params int num_elements: The number of macroscale elements
    '''

    data = xdmf_file_in.getIncrementData(inc, qp_field_name)
    q_output = numpy.zeros([num_elements, 1])

    for e in range(num_elements):
        chunk = data[0][0][e,:].reshape(3,3)
        p = (-1/3)*numpy.trace(chunk)
        dev = chunk + p*numpy.eye(3)
        q_output[e,:] = numpy.sqrt(numpy.tensordot(dev, dev))

    return q_output


def three_pressure(xdmf_file_in, qp_field_name, inc, num_elements, k):
    '''Calculate the pressure of a third order tensor on 3rd index "k"

    :params object xdmf_file_in: The XDMF file to calculate quantities from
    :params str qp_field_name: The stress field (defined at the quadrature points) to consider
    :params list inc: The current time increment
    :params int num_elements: The number of macroscale elements
    :params int k: The 3rd index of the third order tensor
    '''

    data = xdmf_file_in.getIncrementData(inc, qp_field_name)
    p_output = numpy.zeros([num_elements, 1])

    for e in range(num_elements):
        chunk = data[0][0][e,:].reshape(3,3,3)
        p_output[e,:] = (-1/3)*numpy.trace(chunk[:,:,k])

    return p_output


def three_devnorm(xdmf_file_in, qp_field_name, inc, num_elements, k):
    '''Calculate the norm of the deviatoric part of a third order tensor on 3rd index "k"

    :params object xdmf_file_in: The XDMF file to calculate quantities from
    :params str qp_field_name: The stress field (defined at the quadrature points) to consider
    :params list inc: The current time increment
    :params int num_elements: The number of macroscale elements
    :params int k: The 3rd index of the third order tensor
    '''

    data = xdmf_file_in.getIncrementData(inc, qp_field_name)
    q_output = numpy.zeros([num_elements, 1])

    for e in range(num_elements):
        chunk = data[0][0][e,:].reshape(3,3,3)
        p = (-1/3)*numpy.trace(chunk[:,:,k])
        dev = chunk[:,:,k] + p*numpy.eye(3)
        q_output[e,:] = numpy.sqrt(numpy.tensordot(dev, dev))

    return q_output


def average_over_quadrature(xdmf_file_in, type, qp_field_name, inc, num_elements, k=0):
    '''Call relevant stress calculation function for each quadrature point and return average

    :params object xdmf_file_in: The XDMF file to calculate quantities from
    :params str type: The type of stress quantity to calculate.
                      Specify "p" for pressure of a second order tensor.
                      Specify "q" for norm of the deviatoric part of a second order tensor.
                      Specify "p3" for pressure of a third order tensor on 3rd index "k".
                      Specify "q3" for norm of the deviatoric part of a third order tensor on 3rd index "k".
    :params str qp_field_name: The stress field (defined at the quadrature points) to consider
    :params list inc: The current time increment
    :params int num_elements: The number of macroscale elements
    :params int k: The 3rd index of a third order tensor if type is "p3" or "q3"
    '''

    out_all = numpy.zeros([8, num_elements,1])

    for qp in range(8):
        if type == 'p':
            out_all[qp, :, :] = two_pressure(xdmf_file_in, f'{qp_field_name}_{qp}', inc, num_elements)
        elif type == 'q':
            out_all[qp, :, :] = two_devnorm(xdmf_file_in, f'{qp_field_name}_{qp}', inc, num_elements)
        elif type == 'v':
            out_all[qp, :, :] = numpy.sqrt(3/2)*two_devnorm(xdmf_file_in, f'{qp_field_name}_{qp}', inc, num_elements)
        elif type == 'p3':
            out_all[qp, :, :] = three_pressure(xdmf_file_in, f'{qp_field_name}_{qp}', inc, num_elements, k)
        elif type == 'q3':
            out_all[qp, :, :] = three_devnorm(xdmf_file_in, f'{qp_field_name}_{qp}', inc, num_elements, k)
        else:
            print('Specify a valid field calculation type!')

    return numpy.mean(out_all, axis=0)


def filter_stress_measures(xdmf_file_in, xdmf_file_out, incs, times, num_elements, reference_positions, connectivity):
    '''Calculate a variety of stress norms averaged over the quadrature points of trilinear hexahedral elements and write to an XDMF file

    :params object xdmf_file_in: The XDMF file to calculate quantities from
    :params object xdmf_file_out: The XDMF file to write to
    :params list incs: The list containing time increment indices
    :params list times: The list of unique time stamps
    :params int num_elements: The number of macroscale elements
    :params array-like reference_positions: The reference positions of nodes in an existing XDMF file
    :params array-like connectivity: The element to node connectivity in an existing XDMF file

    :returns: Write data to ``xdmf_file_out``
    '''

    fields = (('p', 'cauchy_stress'),
              ('q', 'cauchy_stress'),
              ('v', 'cauchy_stress'),
              ('v', 'symmetric_micro_stress'),
              ('p', 'symmetric_micro_stress'),
              ('q', 'symmetric_micro_stress'),
              ('p3', 'higher_order_stress'),
              ('q3', 'higher_order_stress'),)

    # calculate p and q
    for i in incs:
        grid = xdmf_file_out.addGrid(xdmf_file_out.output_timegrid, {})
        xdmf_file_out.addTime(grid, times[i])
        xdmf_file_out.addPoints(grid, reference_positions,
                                duplicate='filter_reference_positions')
        xdmf_file_out.addConnectivity(grid, "HEXAHEDRON", connectivity,
                                duplicate='filter_connectivity')
        # Calculate quantities
        for type, field in fields:
            if (type == 'p3') or (type == 'q3'):
                for k in range(0,3):
                    average_over_quadrature(xdmf_file_in, type, field, i, num_elements, k)
                    out_field_name = f'{type}_k={k+1}_{field}'
                    xdmf_file_out.addData(grid, out_field_name, out_field, center='Cell')
            else:
                out_field = average_over_quadrature(xdmf_file_in, type, field, i, num_elements)
                out_field_name = f'{type}_{field}'
                xdmf_file_out.addData(grid, out_field_name, out_field, center='Cell')

    return 0


def assign_calibration_results(xdmf_file_out, calibration_map_file, reference_positions, connectivity):
    '''Write calibrated results onto unique elements of an XDMF file

    :params object xdmf_file_out: The XDMF file to write to
    :params str calibration_map_file: A csv file containing previously calibrated parameters.
    :params array-like reference_positions: The reference positions of nodes in an existing XDMF file
    :params array-like connectivity: The element to node connectivity in an existing XDMF file

    :returns: Write data to ``xdmf_file_out``
    '''

    # setup XDMF output
    grid = xdmf_file_out.addGrid(xdmf_file_out.output_timegrid, {})
    xdmf_file_out.addTime(grid, 0.0)
    xdmf_file_out.addPoints(grid, reference_positions,
                            duplicate='filter_reference_positions')
    xdmf_file_out.addConnectivity(grid, "HEXAHEDRON", connectivity,
                            duplicate='filter_connectivity')

    # output parameters
    df = pandas.read_csv(calibration_map_file)
    df = df.sort_values(by='element')
    parameters = list(df.columns)
    for param in parameters:
        array_out = df[param].values.reshape((1,-1))
        xdmf_file_out.addData(grid, param, array_out, center='Cell')

    return 0


def xdmf_3d_calculations(input_file, output_file, write_type, num_elements=None, calibration_map_file=None):
    '''Create an XDMF file containing a variety of derived quantities

    :params str input_file: Specify the input filename for the h5 + XDMF file pair
    :params str output_file: Specify the output filenmae for the h5 + XDMF file pair
    :params str write_type: The type of quantities to write to XDMF. Choose "filter_stress_measures"
                            to calculate stress invariants directly on filter results.
                            Choose "calibration" results to display calibrations on static mesh.
    :params int num_elements: The number of macroscale elements
    :param str calibration_map_file: A csv file containing previously calibrated parameters.
                                     Required if "--write-type calibration"

    :returns: Call either filter_stress_measures or assign_calibration_results depending on "write_type"
    '''

    xdmf_file_in = file_io.xdmf.XDMF(input_file)
    xdmf_file_in.open()

    xdmf_file_out = file_io.xdmf.XDMF(output_filename=output_file)

    num_increments = xdmf_file_in.getNumIncrements()
    incs = list(range(0, num_increments))


    reference_positions = xdmf_file_in.getIncrementReferenceNodePositions(0)[0][0]
    connectivity = xdmf_file_in.getIncrementConnectivity(0)[0][0].reshape((1,-1))

    if write_type == 'filter_stress_measures':
        assert num_elements != None, "num_elements must be specified!"
        times = [xdmf_file_in.getIncrementTime(i)[0] for i in incs]
        filter_stress_measures(xdmf_file_in, xdmf_file_out, incs, times, num_elements,
                               reference_positions, connectivity)
    elif write_type == 'calibration':
        assert calibration_map_file != None, "input_parameters must be provided!"
        assign_calibration_results(xdmf_file_out, calibration_map_file, reference_positions, connectivity)
    else:
        print('Specify a valid write_type')

    xdmf_file_out.write()
    print("XDMF file written!")

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Create an XDMF file containing a variety of derived quantities"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
        help='Specify the input filename for the h5 + XDMF file pair (no suffix)')
    parser.add_argument('-o', '--output-file', type=str,
        help='Specify the output filenmae for the h5 + XDMF file pair (no suffix)')
    parser.add_argument('--num-elements', type=int, required=False, default=None,
        help='The number of macroscale elements')
    parser.add_argument('--write-type', type=str, required=True,
        help='The type of quantities to write to XDMF. Choose "filter_stress_measures"\
              to calculate stress invariants directly on filter results.\
              Choose "calibration" results to display calibrations on static mesh.')
    parser.add_argument('--calibration-map-file', type=str, required=False, default=None,
        help='A csv file containing previously calibrated parameters.\
              Required if "--write-type calibration"')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    args, unknown = parser.parse_known_args()
    sys.exit(xdmf_3d_calculations(input_file=args.input_file,
                                  output_file=args.output_file,
                                  num_elements=args.num_elements,
                                  write_type=args.write_type,
                                  calibration_map_file=args.calibration_map_file,
                                  ))