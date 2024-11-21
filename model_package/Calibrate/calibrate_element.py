#!python
import argparse
import pathlib
import sys
import yaml

from itertools import compress
import numpy
import matplotlib.pyplot
import pandas
import scipy

import calibration_tools
import micromorphic
import xdmf_reader_tools as XRT
import linear_elastic_parameter_constraint_equations as constraints


elastic_parameter_ordering = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',\
                              'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',\
                              'tau8', 'tau9', 'tau10', 'tau11']

def str2bool(v):
    '''Function for converting string to Boolean. Borrowed from: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse

    :param str/bool v: A string or boolean indicating a True or False value

    :returns: True or False
    '''

    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def parameters_to_fparams(parameters):
    '''Map the elastic parameters to the fparams vector for use in the Tardigrade-MOOSE micromorphic linear elastic material model

    :param numpy.ndarray parameters: The parameters vector
        lambda, mu, eta, tau, kappa, nu, sigma, tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11

    :returns: array of fparams
    '''

    fparams = numpy.hstack([[2], parameters[:2], [5], parameters[2:7], [11], parameters[7:18], 2, [parameters[3], parameters[6]]])

    return fparams


def objective(x0, Y, inputs, cal_norm, nu_targ, case, element, nqp, increment=None, stresses_to_include=['S','SIGMA','M'], dev_norm_errors=False):
    '''Primary objective function for calibrating micromorphic linear elasticity constitutive model against homogenized DNS data

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

    parameter_names = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',
                   'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',
                   'tau8', 'tau9', 'tau10', 'tau11']

    model_name=r'LinearElasticity'
    XX = x0

    consts = calibration_tools.evaluate_constraints(XX)
    cvals = numpy.array([c[0] for c in consts])

    # enforce positivity condition
    if (cvals.min() < 0):
        return numpy.inf
        #return 1.e16

    # enforce Poisson ratio constraint within 2%
    if (case == 1) and (nu_targ > 0.0):
        E = inputs[0]
        poisson = XX[0]/(2.*(XX[0] + XX[1]))
        #poisson = numpy.average([E[0][-1,0,0,0],E[0][-1,0,1,1]])
        if (poisson >= 1.01*nu_targ) or (poisson <= .99*nu_targ):
            return numpy.inf

    # Evaluate stresses from DNS strain inputs
    PK2_sim, SIGMA_sim, M_sim, SDVS_sim = calibration_tools.evaluate_model(inputs, XX, model_name, parameters_to_fparams, 0, element, nqp)
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
    if case == 1:
        for t in time_steps:
            PK2_diff = numpy.array([PK2[q][t,e,2,2].flatten() - PK2_sim[q][t,0,-1] for q in range(nqp)]).flatten()
            SIGMA_diff = numpy.array([SIGMA[q][t,e,2,2].flatten() - SIGMA_sim[q][t,0,-1] for q in range(nqp)]).flatten()
            PK2_error = numpy.hstack([PK2_error, PK2_diff])
            SIGMA_error = numpy.hstack([SIGMA_error, SIGMA_diff])
            M_error = [0.]
    elif case == 4 or case == 5:
        for t in time_steps:
            PK2_diff = numpy.array([PK2[q][t,e,:,:].flatten() - PK2_sim[q][t,0,:] for q in range(nqp)]).flatten()
            SIGMA_diff = numpy.array([SIGMA[q][t,e,:,:].flatten() - SIGMA_sim[q][t,0,:] for q in range(nqp)]).flatten()
            M_diff = numpy.array([M[q][t,e,:,:,:].flatten() - M_sim[q][t,0,:] for q in range(nqp)]).flatten()
            PK2_error = numpy.hstack([PK2_error, PK2_diff])
            SIGMA_error = numpy.hstack([SIGMA_error, SIGMA_diff])
            M_error = numpy.hstack([M_error, M_diff])
    elif case == 6:
        for t in time_steps:
            PK2_error =  [0.]
            SIGMA_error =  [0.]
            M_diff = numpy.array([M[q][t,e,:,:,:].flatten() - M_sim[q][t,0,:] for q in range(nqp)]).flatten()
            M_error = numpy.hstack([M_error, M_diff])
    else:
        for t in time_steps:
            PK2_diff = numpy.array([PK2[q][t,e,:,:].flatten() - PK2_sim[q][t,0,:] for q in range(nqp)]).flatten()
            SIGMA_diff = numpy.array([SIGMA[q][t,e,:,:].flatten() - SIGMA_sim[q][t,0,:] for q in range(nqp)]).flatten()
            PK2_error = numpy.hstack([PK2_error, PK2_diff])
            SIGMA_error = numpy.hstack([SIGMA_error, SIGMA_diff])
            M_error = [0.]

    if dev_norm_errors == True:
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

    print(f'obj = {obj}')
    return obj


def opti_options_1(X, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=True, increment=None):
    '''Objective function number 1 used for calibrating first 2 parameters of micromorphic linear elasticity. For case 1.

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    others = [0., 0., 0., 0., 0.,
              0., 0., 0., 0., 0, 0, 1.e-3, 0., 0, 0., 0.]
    if numpy.isclose(nu_targ, 0.0, atol=1e-4):
        XX = numpy.hstack([0.0, X, others])
    else:
        X1 = X[:2]
        XX = numpy.hstack([X1, others])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, nqp, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return(XX)


def opti_options_2(X, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=True, increment=None, dev_norm_errors=False):
    '''Objective function number 2 used for calibrating 7 parameters of micromorphic linear elasticity. For case 2.

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    X1 = X[:7]
    others = [0., 0., 0., 0., 0, 0, 1.e-3, 0., 0, 0., 0.]
    XX = numpy.hstack([X1, others])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, nqp,
                         increment=increment, stresses_to_include=['S', 'SIGMA'], dev_norm_errors=dev_norm_errors))
    else:
        return(XX)


def opti_options_3(X, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=True, increment=None, dev_norm_errors=False):
    '''Objective function number 3 used for calibrating 8 parameters of micromorphic linear elasticity. For cases 3 and 5.

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''

    X1, X2 = X[:7], X[-1]
    tau1to6 = [0., 0., 0., 0., 0, 0]
    tau8to11 = [0., 0, 0., 0.]
    XX = numpy.hstack([X1, tau1to6, X2, tau8to11])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, nqp,
                         increment=increment, stresses_to_include=['S', 'SIGMA', 'M'], dev_norm_errors=dev_norm_errors))
    else:
        return(XX)


def opti_options_4(X, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=True, increment=None, dev_norm_errors=False):
    '''Objective function number 4 used for calibrating all 18 parameters of micromorphic linear elasticity. For cases 4 and 8.

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''

    XX = X
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, nqp,
                         increment=increment, stresses_to_include=['S', 'SIGMA', 'M'], dev_norm_errors=dev_norm_errors))
    else:
        return(XX)


def opti_options_6(X, Y, inputs, cal_norm, nu_targ, case, element, nqp, second_order_params, calibrate=True, increment=None, dev_norm_errors=False):
    '''Objective function number 6 used for calibrating 11 higher order "tau" parameters of micromorphci linear elasticity. For case 6.

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param array second_order_params: Initial guess values for 7 second order parameters
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''

    XX = numpy.hstack([second_order_params, X])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, nqp,
                         increment=increment, stresses_to_include=['M'], dev_norm_errors=dev_norm_errors))
    else:
        return(XX)


def opti_options_7(X, Y, inputs, cal_norm, nu_targ, case, element, nqp, third_order_parameters, calibrate=True, increment=None, dev_norm_errors=False):
    '''Objective function number 7 used for calibrating 7 parameters of micromorphic linear elasticity with fixed higher order parameters. For case 7.

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param array third_order_parameters: Previously calibrated values for 11 higher order parameters
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional list of one or more increments to perform calibration

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''

    XX = numpy.hstack([X, third_order_parameters])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, nqp,
                         increment=increment, stresses_to_include=['S', 'SIGMA'], dev_norm_errors=dev_norm_errors))
    else:
        return(XX)


def handle_output_for_UQ(Xstore, Lstore, case):
    '''TODO: not currently used

    '''

    UQ_dict = {
        'obj':[],
        'lamb':[], 'mu':[], 'eta':[], 'tau':[], 'kappa':[], 'nu':[], 'sigma':[],
        'tau1':[], 'tau2':[], 'tau3':[], 'tau4':[], 'tau5':[], 'tau6':[], 'tau7':[],
        'tau8':[], 'tau9':[], 'tau10':[], 'tau11':[]}
    # # Store results into dictionary
    Xstore = numpy.array(Xstore).T
    UQ_dict['obj'] = Lstore
    UQ_dict['lamb'] = Xstore[0].flatten()
    UQ_dict['mu'] = Xstore[1].flatten()
    UQ_dict['eta'] = Xstore[2].flatten()
    UQ_dict['tau'] = Xstore[3].flatten()
    UQ_dict['kappa'] = Xstore[4].flatten()
    UQ_dict['nu'] = Xstore[5].flatten()
    UQ_dict['sigma'] = Xstore[6].flatten()
    UQ_dict['tau1'] = Xstore[7].flatten()
    UQ_dict['tau2'] = Xstore[8].flatten()
    UQ_dict['tau3'] = Xstore[9].flatten()
    UQ_dict['tau4'] = Xstore[10].flatten()
    UQ_dict['tau5'] = Xstore[11].flatten()
    UQ_dict['tau6'] = Xstore[12].flatten()
    UQ_dict['tau7'] = Xstore[13].flatten()
    UQ_dict['tau8'] = Xstore[14].flatten()
    UQ_dict['tau9'] = Xstore[15].flatten()
    UQ_dict['tau10'] = Xstore[16].flatten()
    UQ_dict['tau11'] = Xstore[17].flatten()

    # remove zero entries depending on case
    # remove = []
    # if case == 1:
        # remove = ['eta','tau','kappa','nu','sigma','tau1','tau2','tau3',
                  # 'tau4','tau5','tau6','tau7','tau8','tau9','tau10','tau11']
    # elif case == 2:
        # remove = ['tau1','tau2','tau3','tau4','tau5','tau6','tau7','tau8',
                  # 'tau9','tau10','tau11']
    # elif case == 3:
        # remove = ['tau1','tau2','tau3','tau4','tau5','tau6','tau8','tau9',
                  # 'tau10','tau11']
    # for item in remove:
        # UQ_dict.pop(item)

    return(UQ_dict)


def calibrate(input_file, output_file, case, Emod, nu, L, element=0, increment=None, plot_file=None,
              average=False, UQ_file=None, cal_norm='L1', bound_half_width=1.e5, dev_norm_errors=False, input_elastic_parameters=None):
    ''' Unpack DNS data and run calibration routine

    :param str input_file: The homogenized XDMF file output by the Micromorphic Filter
    :param str output_file:  The resulting list of parameters stored in a yaml file
    :param int case: The calibration "case".
                     1: two parameter,
                     2: 7 parameter,
                     3: 7 parameter plus tau7 without error for M,
                     4: all 18 parameters,
                     5: 7 parameter plus tau7 with error for M,
                     6: 11 higher order parameters,
                     7: 7 parameters using fixed higher order parameters determined from case 6,
                     8: 7 parameters using initial guess and tighter bounds for higher order parameters determined from case 6
    :param float Emod: Estimate of a homogenized elastic modulus, used for initial parameter estimation
    :param float nu: Estimate of a homogenized Poisson ratio, used for initial parameter estimation
    :param float L: DNS max dimension (width, height, depth, etc.), used for initial parameter estimation
    :param int element: The macro (filter) element to calibration, default is zero
    :param int increment: An optional list of one or more increments to perform calibration
    :param str plot_file: Optional root filename to for plotting results
    :param bool average: Boolean whether or not homogenized DNS results will be averaged
    :param str UQ_file: Optional csv filename to store function evaluations and parameter sets for UQ
    :param str cal_norm: The type of norm to use for calibration ("L1", "L2", or "L1-L2")
    :param float bound_half_width: The uniform parameter bound "half-width" to apply for all parameters to be calibrated. Bounds for lambda will be [0., bound_half_width]. All other parameter bounds will be [-1*bound_half_width, bound_half_width]
    :param bool dev_norm_errors: Boolean whether to inclue deviatoric stress norms during calibration
    :param str input_elastic_parameters: Yaml file containing previously calibrated elastic parameters

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
    
    if average == True:
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
    
    # get target nu from E
    # define time steps to calibrate against
    if increment and (len(increment) == 1):
        nu_inc = int(increment)
    elif increment and (len(increment) > 1):
        nu_inc = int(increment[-1])
    else:
        nu_inc = -1
    nu_targ = numpy.average([(-1*numpy.average([E[q][nu_inc,0,0,0],
                                                E[q][nu_inc,0,1,1]])/E[q][nu_inc,0,2,2]) for q in range(0,nqp)])

    # Estimate initial parameters
    param_est = calibration_tools.Isbuga_micrormorphic_elasticity_parameters(Emod, nu, L)
 
    # Define the elastic bounds
    upper = bound_half_width
    lower = -1*upper
    parameter_bounds = [[0.0, upper]] + [[lower, upper] for _ in range(6)] + [[lower, upper] for _ in range(11)]
    if len(elastic_parameter_ordering) != len(parameter_bounds):
        raise ValueError(f"The parameter bounds and the parameter names do not have the same length {len(parmaeter_bounds)} vs. {len(elastic_parameter_ordering)}")

    # calibrate!
    maxit = 2000
    num_workers = 8
    if case == 1:
        # Case 1 - calibrate just lambda and mu
        print(f'Target Poisson ratio = {nu_targ}')
        if numpy.isclose(nu_targ, 0.0, atol=1e-4):
            print(f'nu_targ is too close to zero to calibrate lambda')
            nu_targ = 0
            lamb_cal = False
        else:
            lamb_cal = True
        param_mask = [lamb_cal, True, False, False, False, False, False,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        #param_est = [59.25, 70.395]
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_1,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=1,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_1(list(res.x), Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 2:
        # Case 2 - calibrate first 7 parameters
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_2,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_2(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 3:
        # Case 3 - calibrate first 7 parameters and tau 7, do not use error for M in objective function
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      True, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_3,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_3(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 4:
        # Case 4 - calibrate all parameters simultaneously
        param_mask = [True, True, True, True, True, True, True,
                      True, True, True, True, True, True, True, True, True, True, True]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_4,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_4(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 5:
        # Case 5 - Same as case 3, but include M error in objective function
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      True, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_3,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_3(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    elif case == 6:
        # Case 6 - Calibrate only the 11 higher order tau parmaeters
        param_mask = [False, False, False, False, False, False, False,
                      True, True, True, True, True, True, True, True, True, True, True]
        second_order_params = param_est[0:7]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_6,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, second_order_params, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_6(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, second_order_params, calibrate=False)
    elif case == 7:
        # Case 7 - Calibrate 7 second order parmaeters and read in existing higher order parameters which are fixed
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        # unpack previously calibrated parameters
        assert input_elastic_parameters != None, "input_elastic_parameters must be provided!"
        input_parameters, _ = calibration_tools.parse_input_parameters(input_elastic_parameters)
        third_order_params = input_parameters[10:21]
        print(f'input third_order_parameters = {third_order_params}')
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_7,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, third_order_params, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_7(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, third_order_params, calibrate=False)
    elif case == 8:
        # Case 8 - Calibrate all 18 parameters and read in existing higher order parameters which as an initial estimate with adjusted bounds
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        print('initial parameter estimation:')
        print(param_est)
        # unpack previously calibrated parameters
        assert input_elastic_parameters != None, "input_elastic_parameters must be provided!"
        input_parameters, _ = calibration_tools.parse_input_parameters(input_elastic_parameters)
        third_order_params = input_parameters[10:21]
        print(f'input third_order_parameters = {third_order_params}')
        param_est = numpy.hstack([param_est, third_order_params])
        third_order_bounds = [sorted([0.1*p, 10.*p]) for p in third_order_params]
        parameter_bounds = parameter_bounds + third_order_bounds
        res = scipy.optimize.differential_evolution(func=opti_options_4,
                                                    bounds=parameter_bounds,
                                                    maxiter=maxit,
                                                    x0=param_est,
                                                    workers=num_workers,
                                                    args=(Y, inputs, cal_norm, nu_targ, case, element, nqp, True, increment, dev_norm_errors))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_4(res.x, Y, inputs, cal_norm, nu_targ, case, element, nqp, calibrate=False)
    else:
        print('Select valid calibration case!')

    # Make a csv of all function evaluations and parameter sets
    if UQ_file:
        print(f'shape of Xstore = {numpy.shape(Xstore)}')
        print(f'shape of Lstore = {numpy.shape(Lstore)}')
        UQ_dict = handle_output_for_UQ(Xstore, Lstore, case)
        df = pandas.DataFrame(UQ_dict)
        df.to_csv(UQ_file, header=True, sep=',', index=False)

    # look at population energy info
    #population = res.population
    #energies = res.population_energies
    #print(f'size of population = {numpy.shape(population)}')
    #print(f'size of population_energies = {numpy.shape(energies)}')
    #print(f'population = \n {population}\n')
    #print(f'energies = \n {energies}\n')

    # plot resulting calibration
    if plot_file:
        print('plotting...')
        model_name=r'LinearElasticity'
        PK2_sim, SIGMA_sim, M_sim, SDVS_sim = calibration_tools.evaluate_model(inputs, params, model_name, parameters_to_fparams, 0, element, nqp)
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
        calibration_tools.plot_stresses(E, PK2, PK2_sim, f'{plot_file}_PK2_fit_case_{case}_bounded.PNG',
                                        element, nqp, 'E', 'S',increment=increment, find_bounds=True)
        calibration_tools.plot_stresses(Ecal, SIGMA, SIGMA_sim, f'{plot_file}_SIGMA_fit_case_{case}_bounded.PNG',
                                        element, nqp, '\mathcal{E}', '\Sigma', increment=increment, find_bounds=True)
        calibration_tools.plot_higher_order_stresses(Gamma, M, M_sim, f'{plot_file}_M_fit_case_{case}_bounded.PNG',
                                                     element, nqp, increment=increment, find_bounds=True)
        # plot the full range of predictions if specific increments were speficied
        if increment:
            calibration_tools.plot_stresses(E, PK2, PK2_sim, f'{plot_file}_PK2_fit_case_{case}_ALL.PNG',
                                            element, nqp, 'E', 'S')
            calibration_tools.plot_stresses(Ecal, SIGMA, SIGMA_sim, f'{plot_file}_SIGMA_fit_case_{case}_ALL.PNG',
                                            element, nqp, '\mathcal{E}', '\Sigma')
            calibration_tools.plot_higher_order_stresses(Gamma, M, M_sim, f'{plot_file}_M_fit_case_{case}_ALL.PNG',
                                                         element, nqp)

    # output parameters
    output_filename = output_file
    output_dict = {}
    p = params
    output_dict['line 1'] = f"2 {p[0]} {p[1]}"
    output_dict['line 2'] = f"5 {p[2]} {p[3]} {p[4]} {p[5]} {p[6]}"
    output_dict['line 3'] = f"11 {p[7]} {p[8]} {p[9]} {p[10]} {p[11]} {p[12]} {p[13]} {p[14]} {p[15]} {p[16]} {p[17]}"
    output_dict['line 4'] = f"2 {p[3]} {p[6]}"
    output_dict['obj'] = f"{res.fun}"
    with open(output_filename, 'w') as f:
        yaml.dump(output_dict, f)
    return


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Calibrate micromorphic linear elasticity on a single filter domain (i.e. macroscale element)"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-i', '--input-file', type=str,
        help="The homogenized XDMF file output by the Micromorphic Filter")
    parser.add_argument('-o', '--output-file', type=str,
        help="The resulting list of parameters stored in a yaml file")
    parser.add_argument('--Emod', type=float,
        help="DNS elastic modulus, used for initial parameter estimation.")
    parser.add_argument('--nu', type=float,
        help="DNS Poisson's ratio, used for initial parameter estimation.")
    parser.add_argument('--L', type=float,
        help="DNS max dimension (width, height, depth, etc.), used for initial parameter estimation.")
    parser.add_argument('--element', type=int, default=0,
        help="The macro (filter) element to calibrate")
    parser.add_argument('--increment', nargs="+", required=False, default=None,
        help="An optional list of one or more increments to perform calibration")
    parser.add_argument('--case', type=int, required=True,
        help="The calibration 'case'.\
              1: two parameter,\
              2: 7 parameter,\
              3: 7 parameter plus tau7 without error for M,\
              4: all 18 parameters,\
              5: 7 parameter plus tau7 with error for M,\
              6: 11 higher order parameters,\
              7: 7 parameters using fixed higher order parameters determined from case 6,\
              8: 7 parameters using initial guess and tighter bounds for higher order parameters determined from case 6")
    parser.add_argument('--plot-file', type=str, required=False, default=None,
        help="Optional root filename to for plotting results")
    parser.add_argument('--average', type=str, required=False, default=False,
        help='Boolean whether or not homogenized DNS results will be averaged')
    parser.add_argument('--UQ-file', type=str, required=False,
        help='Optional csv filename to store function evaluations and parameter sets for UQ')
    parser.add_argument('--cal-norm', type=str, required=False, default='L1',
        help='The type of norm to use for calibration ("L1", "L2", or "L1-L2")')
    parser.add_argument('--bound-half-width', type=float, required=False, default=1.e5,
        help='The uniform parameter bound "half-width" to apply for all parameters to be calibrated.\
              Bounds for lambda will be [0., bound_half_width].\
              All other parameter bounds will be [-1*bound_half_width, bound_half_width]')
    parser.add_argument('--dev-norm-errors', type=str, required=False, default=False,
        help='Boolean whether to inclue deviatoric stress norms during calibration')
    parser.add_argument('--input-elastic-parameters', type=str, required=False, default=None,
        help='Yaml file containing previously calibrated elastic parameters')

    return parser


if __name__ == '__main__':
    parser = get_parser()

    # Objective function evaluation lists
    Xstore = []
    Lstore = []

    args, unknown = parser.parse_known_args()
    sys.exit(calibrate(input_file=args.input_file,
                       output_file=args.output_file,
                       Emod=args.Emod,
                       nu=args.nu,
                       L=args.L,
                       element=args.element,
                       increment=args.increment,
                       case=args.case,
                       plot_file=args.plot_file,
                       average=str2bool(args.average),
                       UQ_file=args.UQ_file,
                       cal_norm=args.cal_norm,
                       bound_half_width=args.bound_half_width,
                       dev_norm_errors=str2bool(args.dev_norm_errors),
                       input_elastic_parameters=args.input_elastic_parameters,
                       ))