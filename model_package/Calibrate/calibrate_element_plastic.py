#!python
import argparse
import pathlib
import sys
import time
import yaml

from itertools import compress
import numpy
import matplotlib.pyplot
import pandas
import scipy

import calibration_tools
import micromorphic
import xdmf_reader_tools as XRT
import calibrate_element


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
                                    #[0.5, 0.5, 0.5, 1e-9, 1e-9],\
                                    [0.0, 0.0, 0.0, 1e-8, 1e-8],\
                                    ])
    return plastic_fparams


# Objective function evaluation lists
Xstore = []
Lstore = []

def objective(x0, Y, inputs, cal_norm, case, element, nqp, increment=None, stresses_to_include=['S','SIGMA','M']):
    '''Primary objective function for calibrating micromorphic elastoplasticity constitutive model against homogenized DNS data

    :param array-like x0: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param int increment: An optional list of one or more increments to perform calibration
    :param list stresses_to_include: Which reference configuration stresses to calculate an error for the objective function, default=['S', 'SIGMA', 'M']

    :returns: the objective function evaluation
    '''

    model_name=r'LinearElasticityDruckerPragerPlasticity'
    # stack parameters
    XX = x0
    # Evaluate stresses from DNS strain inputs
    if increment:
        max_inc = int(increment[-1])
    else:
        max_inc = len(inputs[5])

    # Enforce constraints on plastic parameters for certain cases
    penalty = 0.
    cu0, Hu, cchi0, Hchi = XX[0], XX[1], XX[2], XX[3]
    if case == 10:
        # make sure one of the yield surfaces has softening
        if Hu*Hchi >= 0.:
            penalty += 1.e6 * (abs(Hu) + abs(Hchi) + 1.e-6)
        # make sure the softening surface yields 2nd
        if Hu < 0.:
            # macro-plasticity is softening, make sure cu0 > cchi0
            if cu0 <= cchi0:
                penalty += 1.e6 * (cchi0 - cu0 + 1.e-6)
        if Hchi < 0.:
            # micro-plasticity is softening, make sure cchi0 > cu0
            if cchi0 <= cu0:
                penalty += 1.e6 * (cu0 - cchi0 + 1.e-6)
    # Softening micro-plasticity
    if (case == 8) or (case == 11):
        if cu0 <= cchi0:
            penalty += 1.e6 * (cchi0 - cu0 + 1.e-6)
    # Softening macro-plasticity
    if (case == 9) or (case == 12):
        if cchi0 <= cu0:
            penalty += 1.e6 * (cu0 - cchi0 + 1.e-6)

    # Evaluate stresses from DNS strain inputs
    ## If the solve fails, return objective of infinity
    try:
        PK2_sim, SIGMA_sim, M_sim, SDVS_sim = calibration_tools.evaluate_model(inputs, XX, model_name, stack_parameters, 55, element, nqp)
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
    e = 0
    PK2_dev_error, SIGMA_dev_error, M_dev_error = [], [], []
    ninc = len(inputs[5])
    PK2_sim_mapped = XRT.map_sim(PK2_sim, ninc)
    SIGMA_sim_mapped = XRT.map_sim(SIGMA_sim, ninc)
    M_sim_mapped = XRT.map_sim(M_sim, ninc, third_order=True)
    for t in time_steps:
        error1, error2, error3 = calibration_tools.collect_deviatoric_norm_errors(nqp, t, e,
                                                                            PK2, PK2_sim_mapped,
                                                                            SIGMA, SIGMA_sim_mapped,
                                                                            M, M_sim_mapped)
        PK2_dev_error = numpy.hstack([PK2_dev_error, error1.flatten()])
        SIGMA_dev_error = numpy.hstack([SIGMA_dev_error, error2.flatten()])
        M_dev_error = numpy.hstack([M_dev_error, error3.flatten()])
    PK2_error = numpy.hstack([PK2_error, PK2_dev_error])
    SIGMA_error = numpy.hstack([SIGMA_error, SIGMA_dev_error])
    M_error = numpy.hstack([M_error, M_dev_error])

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

    obj = obj + penalty
    print(f'obj = {obj}')
    return obj


def opti_options_1(X, Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=True, increment=None):
    '''Calibrate macro-plasticity initial cohesion parameter. For case 1.

    :param array-like X: Array of micromorphic plasticity parameters to calibrate
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param array e_params: The elastic fparams
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case".
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    others = [1.e-4,       # macro hardening
              1.e8, 1.e-4, # micro terms
              1.e8, 1.e-4, # micro gradient terms
              ]

    XX = numpy.hstack([X, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, nqp, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_2(X, cohesion, Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=True, increment=None):
    '''Calibrate macro plasticity hardening using an initial estimate/calibration for cohesion. For case 2.

    :param array-like X: Array of micromorphic plasticity parameters to calibrate
    :param float cohesion: The value of the macro-scale initial cohesion parameter
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param array e_params: The elastic fparams
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case".
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    others = [
              1.e8, 1.e-8, # micro terms
              1.e8, 1.e-8, # micro gradient terms
              ]

    XX = numpy.hstack([cohesion, X, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, nqp, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_3(X, Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=True, increment=None):
    '''Calibrate micro-plasticity initial cohesion parameter. For case 3.

    :param array-like X: Array of micromorphic plasticity parameters to calibrate
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param array e_params: The elastic fparams
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case".
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    others = [1.e8, 1.e-8]

    XX = numpy.hstack([others, X, 1.e-8, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, nqp, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_4(X, Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=True, increment=None):
    '''Calibrate macro-plasticity initial cohesion and hardening parameters. For case 4.

    :param array-like X: Array of micromorphic plasticity parameters to calibrate
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param array e_params: The elastic fparams
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case".
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    others = [1.e8, 1.e-8]

    XX = numpy.hstack([X, others, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, nqp, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_5(X, Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=True, increment=None):
    '''Calibrate micro-plasticity initial cohesion and hardening parameters. For case 5.

    :param array-like X: Array of micromorphic plasticity parameters to calibrate
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param array e_params: The elastic fparams
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case".
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    others = [1.e8, 1.e-8]

    XX = numpy.hstack([others, X, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, nqp, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_6(X, Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=True, increment=None):
    '''Calibrate macro-plasticity and micro-plasticity initial cohesion and hardening parameters. For case 6, 8, 9, and 10.

    :param array-like X: Array of micromorphic plasticity parameters to calibrate
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param array e_params: The elastic fparams
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case".
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    others = [1.e8, 1.e-8]

    XX = numpy.hstack([X, others, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, nqp, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return XX


def opti_options_7(X, Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=True, increment=None):
    '''Calibrate macro, micro, and micro gradient plasticity initial cohesion and hardening parameters. For case 7.

    :param array-like X: Array of micromorphic plasticity parameters to calibrate
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param array e_params: The elastic fparams
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case".
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''

    XX = numpy.hstack([X, e_params])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, case, element, nqp, increment=increment, stresses_to_include=['S','SIGMA', 'M']))
    else:
        return XX


def calibrate_plasticity(input_file, output_file, case, input_parameters, element=0, increment=None, plot_file=None, average=False, UQ_file=None, cal_norm='L1'):
    ''' Unpack DNS data and run plasticity calibration routine

    :param str input_file: The homogenized XDMF file output by the Micromorphic Filter
    :param str output_file:  The resulting list of parameters stored in a yaml file
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param str input_parameters: Yaml file containing previously calibrated elastic parameter
    :param int element: The macro (filter) element to calibration, default is zero
    :param int increment: An optional list of one or more increments to perform calibration
    :param str plot_file: Optional root filename to for plotting results
    :param bool average: Boolean whether or not homogenized DNS results will be averaged
    :param str UQ_file: Optional csv filename to store function evaluations and parameter sets for UQ
    :param str cal_norm: The type of norm to use for calibration ("L1", "L2", or "L1-L2")

    :returns: calibrated parameters by minimizing a specified objective function
    '''

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
    start_time = time.time()
    if average == True:
        # average and isolate the quantities needed for calibration
        cauchy = calibration_tools.average_quantities(cauchy, '3x3', element)
        symm = calibration_tools.average_quantities(symm, '3x3', element)
        PK2 = calibration_tools.average_quantities(PK2, '3x3', element)
        SIGMA = calibration_tools.average_quantities(SIGMA, '3x3', element)
        #F = calibration_tools.average_quantities(F, '3x3')
        E = calibration_tools.average_quantities(E, '3x3', element)
        displacement = calibration_tools.average_quantities(displacement, '3', element)
        gradu = calibration_tools.average_quantities(gradu, '3x3', element)
        phi = calibration_tools.average_quantities(phi, '3x3', element)
        estrain = calibration_tools.average_quantities(estrain, '3x3', element)
        h = calibration_tools.average_quantities(h, '3x3', element)
        gradphi = calibration_tools.average_quantities(gradphi, '3x3x3', element)
        nqp = 1
    else:
        # just isolate the quantities needed for calibration
        cauchy = calibration_tools.isolate_element(cauchy, '3x3', element)
        symm = calibration_tools.isolate_element(symm, '3x3', element)
        PK2 = calibration_tools.isolate_element(PK2, '3x3', element)
        SIGMA = calibration_tools.isolate_element(SIGMA, '3x3', element)
        #F = calibration_tools.isolate_element(F, '3x3')
        E = calibration_tools.isolate_element(E, '3x3', element)
        displacement = calibration_tools.isolate_element(displacement, '3', element)
        gradu = calibration_tools.isolate_element(gradu, '3x3', element)
        phi = calibration_tools.isolate_element(phi, '3x3', element)
        estrain = calibration_tools.isolate_element(estrain, '3x3', element)
        h = calibration_tools.isolate_element(h, '3x3', element)
        gradphi = calibration_tools.isolate_element(gradphi, '3x3x3', element)

    # store data for calibration
    Y = [PK2, SIGMA, M]
    inputs = [E, displacement, gradu, phi, gradphi, times]

    # Unpack elastic parameters
    e_params, f_params = calibration_tools.parse_input_parameters(input_parameters)

    # calibrate!
    cal_norm = 'L1'
    maxit = 2000
    num_workers = 8
    if case == 1:
        # Case 1 - Calibrate macro-plasticity initial cohesion parameter
        parameter_bounds = [[1.0, 10.]]
        param_est = [3.0]
        workers = 1
        res = scipy.optimize.differential_evolution(func=opti_options_1,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_1(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 2:
        # Case 2 - Calibrate macro plasticity hardening using an initial estimate/calibration for cohesion
        cohesion = f_params[0]
        parameter_bounds = [[-1000., 1000.]]
        param_est = [-1.0]
        workers = 1
        res = scipy.optimize.differential_evolution(func=opti_options_2,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, cohesion, inputs, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_2(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 3:
        # Case 3 - Calibrate micro-plasticity initial cohesion parameter
        parameter_bounds = [[1.0, 10.]]
        param_est = [3.0]
        workers = 1
        res = scipy.optimize.differential_evolution(func=opti_options_3,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_3(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 4:
        # Case 4 - Calibrate macro-plasticity initial cohesion and hardening parameters
        parameter_bounds = [[1.0, 20.], [1.e-8, 500.]]
        param_est = [3.0, 1.e-8]
        workers = 1
        res = scipy.optimize.differential_evolution(func=opti_options_4,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_4(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 5:
        # Case 5 - Calibrate micro-plasticity initial cohesion and hardening parameters
        parameter_bounds = [[1.0, 20.], [1.e-8, 500.]]
        param_est = [3.0, 1.e-8]
        workers = 1
        res = scipy.optimize.differential_evolution(func=opti_options_5,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_5(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 6:
        # Case 6 - Calibrate macro-plasticity and micro-plasticity initial cohesion and hardening parameters
        parameter_bounds = [[1.0, 20.], [1.e-8, 500.], [1.0, 20.], [1.e-8, 500.]]
        param_est = [3.0, 1.e-4, 3.0, 1.e-4]
        workers = 2
        res = scipy.optimize.differential_evolution(func=opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_6(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 7:
        # Case 7 - Calibrate macro-plasticity, micro-plasticity, and micro-gradient-plasticity initial cohesion and hardening parameters
        parameter_bounds = [[1.0, 20.], [1.e-8, 500.], [1.0, 20.], [1.e-8, 500.], [1.e-3, 20.], [1.e-8, 500.]]
        param_est = [3.0, 1.e-4, 3.0, 1.e-4, 0.1, 1.]
        workers = 4
        res = scipy.optimize.differential_evolution(func=opti_options_7,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_7(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 8:
        # Case 8 - Calibrate macro-plasticity and micro-plasticity initial cohesion and hardening parameters
        ## Ensure macro-plasticity is only hardening and micro-plasticity is only softening
        parameter_bounds = [[1.0, 3.], [1.e-8, 500.], [3.5, 10.], [-500., -1.e-8]]
        param_est = [2.5, 1.e-4, 4.0, -1.e-4]
        workers = 2
        res = scipy.optimize.differential_evolution(func=opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_6(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 9:
        # Case 9 - Calibrate macro-plasticity and micro-plasticity initial cohesion and hardening parameters
        ## Ensure macro-plasticity is only softening and micro-plasticity is only hardening
        #parameter_bounds = [[1.0, 20.], [-500., -1.e-8], [1.0, 20.], [1.e-8, 500.]]
        #param_est = [3.0, -1.e-4, 3.0, 1.e-4]
        parameter_bounds = [[3.5, 20.], [-500., -1.e-8], [1.0, 5.], [1.e-8, 500.]]
        param_est = [4.0, -1.e-4, 2.5, 1.e-4]
        workers = 2
        res = scipy.optimize.differential_evolution(func=opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_6(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 10:
        # Case 10 - Calibrate macro-plasticity and micro-plasticity initial cohesion and hardening parameters
        ## Allow macro-plasticity and micro-plasticity to flip-flop between softening and hardening
        parameter_bounds = [[1.0, 20.], [-500., 500.], [1.0, 20.], [-500., 500.]]
        param_est = [3.0, 1.e-8, 3.1, -1.e-8]
        workers = 2
        res = scipy.optimize.differential_evolution(func=opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_6(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 11:
        # Case 11 - Calibrate macro-plasticity and micro-plasticity initial cohesion and hardening parameters
        ## Ensure macro-plasticity is only hardening and micro-plasticity is only softening
        ## Similar to case 8, used with average=True
        parameter_bounds = [[1.0, 10.], [1.e-4, 500.], [3., 10.], [-500., -1.]]
        param_est = [2.5, 1.e-3, 4.0, -2.]
        workers = 2
        res = scipy.optimize.differential_evolution(func=opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_6(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    elif case == 12:
        # Case 12 - Calibrate macro-plasticity and micro-plasticity initial cohesion and hardening parameters
        ## Ensure macro-plasticity is only softening and micro-plasticity is only hardening
        ## Similar to case 9, used with average=True
        parameter_bounds = [[3.5, 20.], [-500., -1.e-8], [1.0, 10.], [1.e-8, 500.]]
        param_est = [4.0, -1.e-4, 2.5, 1.e-4]
        workers = 2
        res = scipy.optimize.differential_evolution(func=opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=workers,
                                                    args=(Y, inputs, e_params, cal_norm, case, element, nqp, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_6(list(res.x), Y, inputs, e_params, cal_norm, case, element, nqp, calibrate=False)
    else:
        print('Select valid calibration case!')
    end_time = time.time()
    print(f'Finished calibration in {end_time - start_time} seconds!')

    # plot resulting calibration
    if plot_file:
        print('plotting...')
        print(f'parameters = {params}')
        model_name=r'LinearElasticityDruckerPragerPlasticity'
        PK2_sim, SIGMA_sim, M_sim, SDVS_sim = calibration_tools.evaluate_model(inputs, params, model_name, stack_parameters, 55, element, nqp)
        PK2_sim = XRT.map_sim(PK2_sim, ninc)
        SIGMA_sim = XRT.map_sim(SIGMA_sim, ninc)
        M_sim = XRT.map_sim(M_sim, ninc, third_order=True)
        cauchy_sim, symm_sim = XRT.get_current_configuration_stresses(PK2_sim, SIGMA_sim, inputs[2], inputs[3])
        calibration_tools.plot_stresses(E, PK2, PK2_sim, f'{plot_file}_PK2_fit_case_{case}.PNG',
                                        element, nqp, 'E', 'S', increment=increment)
        calibration_tools.plot_stresses(Ecal, SIGMA, SIGMA_sim, f'{plot_file}_SIGMA_fit_case_{case}.PNG',
                                        element, nqp, '\mathcal{E}', '\Sigma', increment=increment)
        calibration_tools.plot_higher_order_stresses(Gamma, M, M_sim, f'{plot_file}_M_fit_case_{case}.PNG',
                                                     element, nqp, increment=increment)
        calibration_tools.plot_stress_norm_calibration_comparison(PK2, PK2_sim, SIGMA, SIGMA_sim, M, M_sim, E, Ecal, Gamma,
                                                                  f'{plot_file}_norms_case_{case}.PNG', nqp, increment=increment)
        # plot the full range of predictions if specific increments were speficied
        if increment:
            calibration_tools.plot_stresses(E, PK2, PK2_sim, f'{plot_file}_PK2_fit_case_{case}_ALL.PNG',
                                            element, nqp, 'E', 'S')
            calibration_tools.plot_stresses(Ecal, SIGMA, SIGMA_sim, f'{plot_file}_SIGMA_fit_case_{case}_ALL.PNG',
                                            element, nqp, '\mathcal{E}', '\Sigma')
            calibration_tools.plot_higher_order_stresses(Gamma, M, M_sim, f'{plot_file}_M_fit_case_{case}_ALL.PNG',
                                                         element, nqp)
            calibration_tools.plot_stress_norm_calibration_comparison(PK2, PK2_sim, SIGMA, SIGMA_sim, M, M_sim, E, Ecal, Gamma,
                                                                      f'{plot_file}_norms_case_{case}_ALL.PNG', nqp)
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
    output_dict['obj_func_value'] = f"{res.fun}"
    with open(output_filename, 'w') as f:
        yaml.dump(output_dict, f)

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Calibrate micromorphic elastoplasticity on a single filter domain (i.e. macroscale element)"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
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
        help="Optional root filename to for plotting results")
    parser.add_argument('--average', type=str, required=False, default=False,
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