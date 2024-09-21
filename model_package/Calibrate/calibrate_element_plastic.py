import sys
import os
import argparse
import inspect
import yaml

import numpy
import matplotlib.pyplot
import scipy
import pandas
from itertools import compress

sys.path.append(r'/projects/tea/tardigrade_plastic/tardigrade_micromorphic_element/src/python')
sys.path.append(r'/projects/tea/tardigrade-examples/model_package')
sys.path.append(f'/projects/tea/tardigrade_plastic/tardigrade_micromorphic_linear_elasticity/src/python')

import micromorphic
import xdmf_reader_tools as XRT
import calibrate_element_temp as calibrate_element


def parse_input_parameters(UI):
    # elastic parameters case
    if len(UI.keys()) <= 5:
        e_params = numpy.hstack([[2], [float(i) for i in UI['line 1'].split(' ')[1:]],
                                 [5], [float(i) for i in UI['line 2'].split(' ')[1:]],
                                 [11],[float(i) for i in UI['line 3'].split(' ')[1:]],
                                 [2], [float(i) for i in UI['line 4'].split(' ')[1:]]])
        return e_params, None
    else:
        e_params = numpy.hstack([[2], [float(i) for i in UI['line 10'].split(' ')[1:]],
                                 [5], [float(i) for i in UI['line 11'].split(' ')[1:]],
                                 [11],[float(i) for i in UI['line 12'].split(' ')[1:]],
                                 [2], [float(i) for i in UI['line 13'].split(' ')[1:]]])
        f_params = numpy.hstack([[float(i) for i in UI['line 01'].split(' ')[1:]],
                                 [float(i) for i in UI['line 02'].split(' ')[1:]],
                                 [float(i) for i in UI['line 03'].split(' ')[1:]]])
        return e_params, f_params


def stack_parameters(params):
    '''Map the elastic and plastic parameters to the fparams vector for use in the Tardigrade-MOOSE micromorphic elasto-plastic model

    :param numpy.ndarray parameters: The parameters vector
        cu0, Hu, cchi0, Hchi, cGchi0, HGchi, all elastic fparams

    :returns: array of fparams
    '''
    offset = 6
    plastic_fparams = numpy.hstack([[2], params[0:2],\
                                    [2], params[2:4],\
                                    [2], params[4:6],\
                                    [2, 0., 0.],\
                                    [2, 0., 0.],\
                                    [2, 0., 0.],\
                                    [2, 0., 0.],\
                                    [2, 0., 0.],\
                                    [2, 0., 0.],\
                                    params[offset:],\
                                    [0.5, 0.5, 0.5, 1e-9, 1e-9],\
                                    #[0.0, 0.0, 0.0, 1e-8, 1e-8],\
                                    ])
    return plastic_fparams

# Objective function evaluation lists
Xstore = []
Lstore = []

def objective(x0, Y, inputs, cal_norm, case, element, increment=None, stresses_to_include=['S','SIGMA','M']):

    model_name=r'LinearElasticityDruckerPragerPlasticity'
    # stack parameters
    XX = x0
    # Evaluate stresses from DNS strain inputs
    if increment:
        max_inc = int(increment[-1])
    else:
        max_inc = len(inputs[5])
    # try to evaluate model
    try:
        PK2_sim, SIGMA_sim, M_sim, SDVS_sim = calibrate_element.evaluate_model(inputs, XX, model_name, stack_parameters, 55, element, max_inc)
    except:
        return numpy.inf

    displacement, grad_u, phi, grad_phi = inputs[1], inputs[2], inputs[3], inputs[4]
    # Parse out stresses from DNS stress data Y
    PK2, SIGMA, M = Y[0], Y[1], Y[2]

    # Number of time steps and elements
    steps = PK2[0].shape[0]
    num_elem = PK2[0].shape[1]

    # Initialize errors and objective
    PK2_error   = []
    SIGMA_error = []
    M_error     = []
    obj = 0

   # define time steps to calibrate against
    if increment and (len(increment) == 1):
        time_steps = [int(increment)]
    elif increment and (len(increment) > 1):
        time_steps = [int(i) for i in increment]
    else:
        time_steps = range(steps)

    # Accumulate errors
    elem = element
    e = 0
    for t in time_steps:
        PK2_error = numpy.hstack([PK2_error, PK2[0][t,0,:,:].flatten() - PK2_sim[0][t,0,:]])
        SIGMA_error = numpy.hstack([SIGMA_error, SIGMA[0][t,0,:,:].flatten() - SIGMA_sim[0][t,0,:]])
        M_error = numpy.hstack([M_error, M[0][t,0,:,:,:].flatten() - M_sim[0][t,0,:]])

    # collect errors
    errors    = {'S':PK2_error, 'SIGMA':SIGMA_error, 'M':M_error}

    sparsity_control = 0.1
    # calculate residual, L1 norm
    for stress in stresses_to_include:
        if cal_norm == 'L1':
            obj += numpy.abs(errors[stress]).sum()
        elif cal_norm == 'L2':
            obj += numpy.dot(errors[stress], errors[stress])
        elif cal_norm == 'L1_sparse':
            obj += numpy.abs(errors[stress]).sum() + sparsity_control*SOMETHING
            print('NOT READY YET!')
        elif cal_norm == 'L2_sparse':
            obj += numpy.dot(errors[stress], errors[stress]) + sparsity_control*SOMETHING
            print('NOT READY YET!')
        elif cal_norm == 'L1-L2':
            obj += 0.5*(numpy.abs(errors[stress]).sum()) + 0.5*(numpy.dot(errors[stress], errors[stress]))
        else:
            print('Specify valid objective!')

    Xstore.append(numpy.copy(XX))
    Lstore.append(obj)

    print(f'obj = {obj}')
    return obj


def opti_options_1(X, Y, inputs, e_params, cal_norm, case, element, calibrate=True, increment=None):
    '''Calibrate macro-plasticity initial cohesion parameter

    '''
    others = [1.e-8,       # macro hardening
              1.e8, 1.e-8, # micro terms
              1.e8, 1.e-8, # micro gradient terms
              ]

    XX = numpy.hstack([X, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_2(X, cohesion, Y, inputs, e_params, cal_norm, case, element, calibrate=True, increment=None):
    '''Calibrate macro-plasticity hardening parameter

    '''
    others = [
              1.e8, 1.e-8, # micro terms
              1.e8, 1.e-8, # micro gradient terms
              ]

    XX = numpy.hstack([cohesion, X, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_3(X, Y, inputs, e_params, cal_norm, case, element, calibrate=True, increment=None):
    '''Calibrate micro-plasticity initial cohesion parameter

    '''
    others = [1.e8, 1.e-8]

    XX = numpy.hstack([others, X, 1.e-8, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_4(X, Y, inputs, e_params, cal_norm, case, element, calibrate=True, increment=None):
    '''Calibrate macro-plasticity initial cohesion and hardening parameters

    '''
    others = [1.e8, 1.e-8]

    XX = numpy.hstack([X, others, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_5(X, Y, inputs, e_params, cal_norm, case, element, calibrate=True, increment=None):
    '''Calibrate micro-plasticity initial cohesion and hardening parameters

    '''
    others = [1.e8, 1.e-8]

    XX = numpy.hstack([others, X, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_6(X, Y, inputs, e_params, cal_norm, case, element, calibrate=True, increment=None):
    '''Calibrate macro-plasticity and micro-plasticity initial cohesion and hardening parameters

    '''
    others = [1.e8, 1.e-8]

    XX = numpy.hstack([X, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def calibrate_plasticity(input_file, output_file, case, input_parameters, element=0, increment=None, plot_file=None, average=True, UQ_file=None):


    PK2_sdstore = []
    SIGMA_sdstore = []
    M_sdstore = []
 
    # Read in the data
    filename = input_file
    data, geometry, topology = XRT.parse_xdmf_output(filename)
    nqp = 8
    ninc = numpy.shape(data['time'])[0]
    nel = numpy.shape(data['cauchy_stress_0_0'])[0]
 
    # Read in the position information
    displacement, gradu, phi, gradphi = XRT.construct_degrees_of_freedom(data, nqp, nel)

    # Read in the stress information
    cauchy, symm, m = XRT.collect_stresses(data, nqp, nel)
    PK2, SIGMA, M = XRT.get_reference_configuration_stresses(data, nqp, nel)

    # Read in the strain information
    E, Ecal, Gamma, F, chi, grad_chi, estrain, h = XRT.compute_deformations(data, nqp, nel)

    # Get times
    times = numpy.unique(data['time'])
    ninc = len(times)

    # always average fields, but only for selected element
    cauchy      = calibrate_element.average_quantities(cauchy, '3x3', element)
    symm        = calibrate_element.average_quantities(symm, '3x3', element)
    PK2         = calibrate_element.average_quantities(PK2, '3x3', element)
    SIGMA       = calibrate_element.average_quantities(SIGMA, '3x3', element)
    E           = calibrate_element.average_quantities(E, '3x3', element)
    displacement    = calibrate_element.average_quantities(displacement, '3', element)
    gradu       = calibrate_element.average_quantities(gradu, '3x3', element)
    phi         = calibrate_element.average_quantities(phi, '3x3', element)
    estrain     = calibrate_element.average_quantities(estrain, '3x3', element)
    h           = calibrate_element.average_quantities(h, '3x3', element)
    gradphi     = calibrate_element.average_quantities(gradphi, '3x3x3', element)

    # store data for calibration
    Y = [PK2, SIGMA, M]
    inputs = [E, displacement, gradu, phi, gradphi, times]

    # Unpack elastic parameters
    stream = open(input_parameters, 'r')
    UI = yaml.load(stream, Loader=yaml.FullLoader)
    stream.close()
    e_params, f_params = parse_input_parameters(UI)

    # Assume we are just calibrating macro plasticity for now!

    # calibrate!
    cal_norm = 'L1'
    maxit = 2000
    # Case 1 calibrated just initial cohesion
    if case == 1:
        parameter_bounds = [[0.0, 10.]]
        param_est = [3.0]
        res = scipy.optimize.differential_evolution(func=opti_options_1,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_1(list(res.x), Y, inputs, e_params, cal_norm, case, element, calibrate=False)
    elif case == 2:
        cohesion = f_params[0]
        parameter_bounds = [[-1000., 1000.]]
        param_est = [-1.0]
        res = scipy.optimize.differential_evolution(func=opti_options_2,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    args=(Y, cohesion, inputs, cal_norm, case, element, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_2(list(res.x), Y, inputs, e_params, cal_norm, case, element, calibrate=False)
    elif case == 3:
        parameter_bounds = [[0.0, 10.]]
        param_est = [3.0]
        res = scipy.optimize.differential_evolution(func=opti_options_3,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_3(list(res.x), Y, inputs, e_params, cal_norm, case, element, calibrate=False)
    elif case == 4:
        parameter_bounds = [[1.0, 100.], [-100., 100]]
        param_est = [3.0, 1.e-8]
        res = scipy.optimize.differential_evolution(func=opti_options_4,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_4(list(res.x), Y, inputs, e_params, cal_norm, case, element, calibrate=False)
    elif case == 5:
        parameter_bounds = [[0.0, 10.], [-100., 10]]
        param_est = [3.0, 1.e-8]
        res = scipy.optimize.differential_evolution(func=opti_options_5,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_5(list(res.x), Y, inputs, e_params, cal_norm, case, element, calibrate=False)
    elif case == 6:
        parameter_bounds = [[0.0, 10.], [-100., 10], [0.0, 10.], [-100., 10]]
        param_est = [3.0, 1.e-8, 3.0, 1.e-8]
        res = scipy.optimize.differential_evolution(func=opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_6(list(res.x), Y, inputs, e_params, cal_norm, case, element, calibrate=False)
    else:
        print('Select valid calibration case!')

    # plot resulting calibration
    if plot_file:
        print('plotting...')
        print(f'parameters = {params}')
        model_name=r'LinearElasticityDruckerPragerPlasticity'
        PK2_sim, SIGMA_sim, M_sim, SDVS_sim = calibrate_element.evaluate_model(inputs, params, model_name, stack_parameters, 55, element)
        #PK2_sim, SIGMA_sim, M_sim, SDVS_sim = calibrate_element.evaluate_model(inputs, XX, model_name, stack_parameters, 55, element, max_inc)
        PK2_sim = XRT.map_sim(PK2_sim, ninc)
        SIGMA_sim = XRT.map_sim(SIGMA_sim, ninc)
        cauchy_sim, symm_sim = XRT.get_current_configuration_stresses(PK2_sim, SIGMA_sim, inputs[2], inputs[3])

        calibrate_element.plot_stresses(estrain, cauchy, cauchy_sim, f'{plot_file}_cauchy_fit_case_{case}.PNG', element)
        calibrate_element.plot_stresses(estrain, symm, symm_sim, f'{plot_file}_symm_fit_case_{case}.PNG', element)

    # output parameters!
    output_filename = output_file
    output_dict = {}
    p = params
    e = e_params
    output_dict['line 01'] = f'2 {p[0]} {p[1]}'
    output_dict['line 02'] = f'2 {p[2]} {p[3]}'
    output_dict['line 03'] = f'2 {p[4]} {p[5]}'
    output_dict['line 04'] = f'2 0. 0.'
    output_dict['line 05'] = f'2 0. 0.'
    output_dict['line 06'] = f'2 0. 0.'
    output_dict['line 07'] = f'2 0. 0.'
    output_dict['line 08'] = f'2 0. 0.'
    output_dict['line 09'] = f'2 0. 0.'
    output_dict['line 10'] = f"2 {e[1]} {e[2]}"
    output_dict['line 11'] = f"5 {e[4]} {e[5]} {e[6]} {e[7]} {e[8]}"
    output_dict['line 12'] = f"11 {e[10]} {e[11]} {e[12]} {e[13]} {e[14]} {e[15]} {e[16]} {e[17]} {e[18]} {e[19]} {e[20]}"
    output_dict['line 13'] = f"2 {e[5]} {e[8]}"
    output_dict['obj'] = f"{res.fun}"
    with open(output_filename, 'w') as f:
        yaml.dump(output_dict, f)
    return
    return


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Calibrate micromorphic linear elasticity for averaged output on a single filter domain (i.e. macroscale element)"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
        help="The homogenized XDMF file output by the Micromorphic Filter")
    parser.add_argument('-o', '--output-file', type=str,
        help="The resulting list of parameters stored in a yaml file")
    parser.add_argument('--element', type=int, default=0,
        help="The macro (filter) element to calibrate")
    parser.add_argument('--increment', nargs="+", required=False, default=None,
        help="An optional list of one or more increments to perform calibration")
    parser.add_argument('--case', type=int, required=True,
        help="Specify the calibration 'case'. 1: two parameter, 2: 7 parameter,\
              3: 7 parameter plus tau7, 4: all 18 parameters")
    parser.add_argument('--input-parameters', type=str, required=True,
        help="A yaml file containing previously calibrated parameters")
    parser.add_argument('--plot-file', type=str, required=False, default=None,
        help="Optional root filename to plot Cauchy and symmetric micro stress\
              comparison between DNS and calibration results")
    parser.add_argument('--average', type=str, required=False, default=True,
        help='Boolean whether or not homogenized DNS results will be averaged')
    parser.add_argument('--UQ-file', type=str, required=False,
        help='Optional csv filename to store function evaluations and parameter sets for UQ')

    return parser


if __name__ == '__main__':
    parser = get_parser()
   
    args, unknown = parser.parse_known_args()
    sys.exit(calibrate_plasticity(
                                  input_file=args.input_file,
                                  output_file=args.output_file,
                                  element=args.element,
                                  increment=args.increment,
                                  case=args.case,
                                  input_parameters=args.input_parameters,
                                  plot_file=args.plot_file,
                                  average=calibrate_element.str2bool(args.average),
                                  UQ_file=args.UQ_file,
                                  ))