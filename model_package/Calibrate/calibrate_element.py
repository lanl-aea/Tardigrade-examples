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


def average_quantities(quantities, type, elem):
    '''Average tensor quantites over 8 quadrature points

    :param dict quantities: A 2nd or 3rd order tensor dictionary with keys for quadrature points and values storing an array where indices correspond to time, element number, and tensor components
    :param str type: A string specifying the type of tensor to average. Use "3" for a vector. Use "3x3" for a regular second order tensor. Use "9" for a flattened second order tensor. Use "3x3x3" for a third order tensor.
    :param int elem: The macro (filter) element to calibrate

    :returns: ``output`` dict with same indices as ``quantities`` and a single key
    '''

    output = {}
    shapes = numpy.shape(quantities[0])
    
    if type == '9':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2]))
        for k in range(9):
            mean_field = []
            for qp in quantities.keys():
                mean_field.append(quantities[qp][:,elem,k])
            means = numpy.mean(mean_field, axis=0)
            output[0][:,elem,k] = means
    elif type == '3x3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2], shapes[3]))
        for i in range(3):
            for j in range(3):
                mean_field = []
                for qp in quantities.keys():
                    mean_field.append(quantities[qp][:,elem,i,j])
                means = numpy.mean(mean_field, axis=0)
                output[0][:,0,i,j] = means
    elif type == '3x3x3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2], shapes[3], shapes[4]))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    mean_field = []
                    for qp in quantities.keys():
                        mean_field.append(quantities[qp][:,elem,i,j,k])
                    means = numpy.mean(mean_field, axis=0)
                    output[0][:,0,i,j,k] = means
    elif type == '3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2]))
        for i in range(3):
            mean_field = []
            for qp in quantities.keys():
                mean_field.append(quantities[qp][:,elem,i])
            means = numpy.mean(mean_field, axis=0)
            output[0][:,0,i] = means

    return(output)


def plot_stresses(estrain, stress, stress_sim, output_name, element):
    '''Plot comparison of stress vs strain (in the current configuration) between homogenized DNS results against calibrated model predictions

    :param dict estrain: The quantities dict storing Euler-Almansi strain
    :param dict stress: The quantities dict storing either Cauchy or Symmetric micro stress of the homogenized DNS results
    :param dict stress: The quantities dict storing either Cauchy or Symmetric micro stress of the calibrated model predictions
    :param str output_name: The output plot name
    :param int element: The macro (filter) element considered for calibration

    :returns: ``output_name``
    '''

    name = output_name.replace('.PNG','')
    fig1 = matplotlib.pyplot.figure(name, figsize=(11,9))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    ybounds = [-1, 1]
 
    colors = matplotlib.pyplot.rcParams['axes.prop_cycle'].by_key()['color']
    k = 0
    e = 0
    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
            if 'cauchy' in output_name:
                plot_label = r"$\sigma_{" + str(i+1) + str(j+1) + "}$ (MPa)"
            if 'symm' in output_name:
                plot_label = r"$s_{" + str(i+1) + str(j+1) + "}$ (MPa)"

            ax1.plot(estrain[0][:,e,i,j], stress[0][:,e,i,j], '-')
            ax1.plot(estrain[0][:,e,i,j], stress_sim[0][:,e,i,j], 'o')
 
            ax1.set_xlabel(r"$e_{" + str(i+1) + str(j+1) + "}$", fontsize=14)
            ax1.set_ylabel(plot_label, fontsize=14)
            matplotlib.pyplot.ticklabel_format(style='sci', axis='x')
            matplotlib.pyplot.ticklabel_format(style='sci', axis='y')
            if (i == 2) and (j ==2):
                matplotlib.pyplot.xticks(rotation=45)
            k = k + 1
 
    ax1 = axes1[1][2]
    matplotlib.pyplot.figure(name)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(f'{output_name}')
 
    return 0


def initial_estimate(Emod, nu, L):
    '''Calculate initial estimate of 18 parameter micromorphic linear elasticity model parameters using method defined in https://doi.org/10.1016/j.ijengsci.2011.04.006

    :param float Emod: An estimate of homogenized elastic modulus
    :param float nu: An estimate of the homogenized Poisson ratio
    :param float L: An estimate of the length scale parameter

    :returns: array of estimated micromorphic linear elasticity parameters
    '''

    print(f"E, nu, L = {Emod}, {nu}, {L}")
 
    # calculate "classic" lame parameters
    lame_lambda = Emod*nu/((1.+nu)*(1.-2*nu))
    lame_mu     = Emod/(2*(1.+nu)) #shear modulus, K
   
    # estimate characteristic length
    Lc = numpy.sqrt(3*(L**2))
   
    # estimate micromorphic parameters
    lamb = lame_lambda
    mu = lame_mu
    eta = 1.53*lame_lambda
    tau = 0.256*lame_lambda
    kappa = 0.833*lame_mu
    nu_new = 0.667*lame_mu
    sigma = 0.4167*lame_mu
   
    tau_1 = 0.111*(lame_lambda*Lc*Lc)
    tau_2 = 0.185*(lame_lambda*Lc*Lc)
    tau_3 = 0.185*(lame_lambda*Lc*Lc)
    tau_4 = 0.204*(lame_lambda*Lc*Lc)
    tau_5 = 0.1*(lame_lambda*Lc*Lc)
    tau_6 = 0.256*(lame_lambda*Lc*Lc)
    tau_7 = 0.670*(lame_mu*Lc*Lc)
    tau_8 = 0.495*(lame_mu*Lc*Lc)
    tau_9 = 0.495*(lame_mu*Lc*Lc)
    tau_10 = 0.408*(lame_mu*Lc*Lc)
    tau_11 = 0.495*(lame_mu*Lc*Lc)
   
    # collect
    parameters = numpy.array([lamb, mu, eta, tau, kappa, nu_new, sigma,
                           tau_1, tau_2, tau_3, tau_4, tau_5, tau_6,
                           tau_7, tau_8, tau_9, tau_10, tau_11])

    return(parameters)


def evaluate_constraints(parameters, svals=None):
    '''Evaluate Smith conditions by calling tardigrade_micromorphic_linear_elasticity/src/python/linear_elastic_parameter_constraint_equations

    :param array-like parameters: an array of 18 micromorphic linear elasticity parameters
    :param array-like svals: TODO figure out what this is for

    :returns: a dictionary of constants from evaluating the Smith conditions
    '''
 
    elastic_parameter_ordering = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',\
                                  'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',\
                                  'tau8', 'tau9', 'tau10', 'tau11']
 
    parameter_dictionary = dict(zip(elastic_parameter_ordering, parameters[:18]))
 
    consts = [constraints.evaluate_g1, 
              constraints.evaluate_g2, 
              constraints.evaluate_g3,
              constraints.evaluate_g4,
              constraints.evaluate_g5,
              constraints.evaluate_g6, 
              constraints.evaluate_g7, 
              constraints.evaluate_g8,
              constraints.evaluate_g9,
              constraints.evaluate_g10,
              constraints.evaluate_g11,
              constraints.evaluate_g12,
              constraints.evaluate_g13]

    if (svals is None):
        svals = dict([(f's{i+1}', 0) for i in range(len(consts))])
 
    parameter_dictionary.update(svals)
 
    return [const(**parameter_dictionary) for const in consts]


def evaluate_model(inputs, parameters, model_name, parameters_to_fparams, nsdvs, element, maxinc=None, dim=3, maxsubiter=5):
    """Evaluate the model given the parameters. Copied from overlap_coupling/src/python/read_xdmf_output.py.
   
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param numpy.ndarray parameters: The array of parameters
    :param str model_name: The name of the model
    :param func parameters_to_fparams: A function that converts the parameters vector to the fparams vector required
        for the function
    :param int nsdvs: The number of solution dependant state variables
    :param int element: The macro (filter) element to calibration
    :param int maxinc: The maximum increment to evaluate
    :param int dim: The spatial dimension of the problem, default=3
    :param int maxsubiter: The maximum number of sub iterations, default=5

    :returns: evaluated micromorphic simulation quantities for PK2, SIGMA, M, and SDVS
    """

    E, displacement, grad_u, phi, grad_phi, time = inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5]
    ninc = E[0].shape[0]

    nel = 1
    nqp = 1

    if maxinc is None:
        maxinc = ninc-1

    PK2_sim   = dict([(qp,numpy.zeros((maxinc+1,nel,dim * dim))) for qp in range(nqp)])
    SIGMA_sim = dict([(qp,numpy.zeros((maxinc+1,nel,dim * dim))) for qp in range(nqp)])
    M_sim     = dict([(qp,numpy.zeros((maxinc+1,nel,dim * dim * dim))) for qp in range(nqp)])
    SDVS_sim  = dict([(qp,numpy.zeros((maxinc+1,nel,nsdvs))) for qp in range(nqp)])
 
    keys = ['errorCode', 'PK2', 'SIGMA', 'M', 'SDVS',\
            'DPK2Dgrad_u', 'DPK2Dphi', 'DPK2Dgrad_phi',\
            'DSIGMADgrad_u', 'DSIGMADphi', 'DSIGMADgrad_phi',\
            'DMDgrad_u', 'DMDphi', 'DMDgrad_phi',\
            'ADD_TERMS', 'ADD_JACOBIANS', 'output_message']

    tp = 0

    nsubiter = 0

    elem = element
    e = 0
    for qp in range(nqp):
        for i in range(maxinc+1):
            #print("increment: ", i)
            # Map the parameters vector to the function parameters
            fparams = parameters_to_fparams(parameters)

            sp = 0
            ds = 1.

            if (i == 0):
                previous_SDVS_s = numpy.zeros(nsdvs)
            else:
                previous_SDVS_s = numpy.copy(SDVS_sim[qp][i-1,e,:])
       
            while (sp < 1.0):

                s = sp + ds

                time_1     = time[i]
                grad_u_1   = grad_u[qp][i, e, :, :]
                phi_1      = phi[qp][i, e, :, :]
                grad_phi_1 = grad_phi[qp][i, e, :, :, :]

                if (i == 0):
                    time_0     = 0
                    grad_u_0   = numpy.zeros((3,3))
                    phi_0      = numpy.zeros((3,3))
                    grad_phi_0 = numpy.zeros((3,3,3))

                else:
                    time_0     = time[i-1]
                    grad_u_0   = grad_u[qp][i-1, e, :, :]
                    phi_0      = phi[qp][i-1, e, :, :]
                    grad_phi_0 = grad_phi[qp][i-1, e, :, :]

                t                = (time_1 - time_0) * s + time_0
                current_grad_u   = (grad_u_1 - grad_u_0) * s + grad_u_0
                current_phi      = (phi_1 - phi_0) * s + phi_0
                current_grad_phi = (grad_phi_1 - grad_phi_0) * s + grad_phi_0

                tp                = (time_1 - time_0) * sp + time_0
                previous_grad_u   = (grad_u_1 - grad_u_0) * sp + grad_u_0
                previous_phi      = (phi_1 - phi_0) * sp + phi_0
                previous_grad_phi = (grad_phi_1 - grad_phi_0) * sp + grad_phi_0

                current_phi = current_phi.flatten()
                previous_phi = previous_phi.flatten()
               
                current_grad_phi = current_grad_phi.reshape((dim * dim, dim))
                previous_grad_phi = previous_grad_phi.reshape((dim * dim, dim))

                #TODO: add dof and add grad dof not currently used
                current_ADD_DOF = numpy.zeros((1))
                current_ADD_grad_DOF = numpy.zeros((1,3))

                previous_ADD_DOF = numpy.zeros((1))
                previous_ADD_grad_DOF = numpy.zeros((1,3))

                # Evaluate the model
                values = micromorphic.evaluate_model(model_name, numpy.array([t, t - tp]), fparams,
                                                     current_grad_u, current_phi, current_grad_phi,
                                                     previous_grad_u, previous_phi, previous_grad_phi,
                                                     previous_SDVS_s,
                                                     current_ADD_DOF, current_ADD_grad_DOF,
                                                     previous_ADD_DOF, previous_ADD_grad_DOF)

                results = dict(zip(keys, values))

                if (results['errorCode'] == 1):
                    #print("error")
                    ds = 0.5 * ds
                    nsubiter += 1

                    if (nsubiter > maxsubiter):
                        break

                elif (results['errorCode'] == 2):
                    errormessage = f"evaluate_model return error code {results['errorCode']}\n\n"
                    errormessage += results['output_message'].decode("utf-8")
                    raise IOError(errormessage)

                else:
                    sp += ds
                    nsubiter = 0

                    if numpy.isclose(sp, 1):
                        ds = 1
                    else:
                        ds = 1 - sp

                    previous_SDVS_s = numpy.copy(results['SDVS'])

            if (results['errorCode'] != 0):
                errormessage = f"evaluate_model returned error code {results['errorCode']}\n\n"
                errormessage += results['output_message'].decode('utf-8')
                print(parameters, 'fail')

                return numpy.nan

            PK2_sim[qp][i,e,:]   = results['PK2']
            SIGMA_sim[qp][i,e,:] = results['SIGMA']
            M_sim[qp][i,e,:]     = results['M']
            SDVS                 = results['SDVS']
            SDVS_sim[qp][i,e,:]  = results['SDVS']

    return PK2_sim, SIGMA_sim, M_sim, SDVS_sim


def parameters_to_fparams(parameters):
    '''Map the elastic parameters to the fparams vector for use in the Tardigrade-MOOSE micromorphic linear elastic material model

    :param numpy.ndarray parameters: The parameters vector
        lambda, mu, eta, tau, kappa, nu, sigma, tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11

    :returns: array of fparams
    '''

    fparams = numpy.hstack([[2], parameters[:2], [5], parameters[2:7], [11], parameters[7:18], 2, [parameters[3], parameters[6]]])

    return fparams

# Objective function evaluation lists
Xstore = []
Lstore = []


def objective(x0, Y, inputs, cal_norm, nu_targ, case, element, increment=None, stresses_to_include=['S','SIGMA','M']):
    '''Primary objective function for calibrating micromorphic linear elasticity constitutive model against homogenized DNS data

    :param array-like x0: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param int increment: An optional argumet to calibrate only a specific increment, default=None
    :param list stresses_to_include: Which reference configuration stresses to calculate an error for the objective function, default=['S', 'SIGMA', 'M']

    :returns: the objective function evaluation
    '''

    parameter_names = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',
                   'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',
                   'tau8', 'tau9', 'tau10', 'tau11']

    model_name=r'LinearElasticity'
    XX = x0
    # #consts = evaluate_constraints(XX, parameter_names)
    # try:
        # consts = evaluate_constraints(XX)
        # cvals = numpy.array([c[0] for c in consts])
    # except:
        # #return numpy.inf
        # return 1.e16
    consts = evaluate_constraints(XX)
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
    PK2_sim, SIGMA_sim, M_sim, SDVS_sim = evaluate_model(inputs, XX, model_name, parameters_to_fparams, 0, element)
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
    if increment:
        print(f'increment = {increment}')
        time_steps = [increment]
    else:
        time_steps = range(steps)

    # Accumulate errors
    elem = element
    e = 0
    if case == 1:
        for t in time_steps:
            PK2_error = numpy.hstack([PK2_error, PK2[0][t,0,2,2] - PK2_sim[0][t,0,-1]])
            SIGMA_error = numpy.hstack([PK2_error, PK2[0][t,0,2,2] - PK2_sim[0][t,0,-1]])
            M_error = []
    elif case == 4:
        for t in time_steps:
            PK2_error = numpy.hstack([PK2_error, PK2[0][t,0,:,:].flatten() - PK2_sim[0][t,0,:]])
            SIGMA_error = numpy.hstack([SIGMA_error, SIGMA[0][t,0,:,:].flatten() - SIGMA_sim[0][t,0,:]])
            M_error = numpy.hstack([M_error, M[0][t,0,:,:,:].flatten() - M_sim[0][t,0,:]])
    else:
        for t in time_steps:
            PK2_error = numpy.hstack([PK2_error, PK2[0][t,0,:,:].flatten() - PK2_sim[0][t,0,:]])
            SIGMA_error = numpy.hstack([SIGMA_error, SIGMA[0][t,0,:,:].flatten() - SIGMA_sim[0][t,0,:]])
            M_error = []

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


def opti_options_1(X, Y, inputs, cal_norm, nu_targ, case, element, calibrate=True, increment=None):
    '''Objective function number 1 used for calibrating first 2 parameters of micromorphic linear elasticity

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional argumet to calibrate only a specific increment, default=None

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    others = [0e-3, 0e-3, 0e-3, 0e-3, 0e-3,
              0e-3, 0e-3, 0e-3, 0e-3, 0, 0, 1e-3, 0e-3, 0, 0e-3, 0e-3]
    if numpy.isclose(nu_targ, 0.0, atol=1e-4):
        XX = numpy.hstack([0.0, X, others])
    else:
        X1 = X[:2]
        XX = numpy.hstack([X1, others])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, increment=increment, stresses_to_include=['S','SIGMA']))
    else:
        return(XX)


def opti_options_2(X, Y, inputs, cal_norm, nu_targ, case, element, calibrate=True, increment=None):
    '''Objective function number 2 used for calibrating 7 parameters of  micromorphic linear elasticity

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional argumet to calibrate only a specific increment, default=None

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''
    X1 = X[:7]
    others = [0e-3, 0e-3, 0e-3, 0e-3, 0, 0, 1e-3, 0e-3, 0, 0e-3, 0e-3]
    XX = numpy.hstack([X1, others])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, increment=increment, stresses_to_include=['S', 'SIGMA']))
    else:
        return(numpy.hstack([X1, others]))


def opti_options_3(X, Y, inputs, cal_norm, nu_targ, case, element, calibrate=True, increment=None):
    '''Objective function number 3 used for calibrating 8 parameters of  micromorphic linear elasticity

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional argumet to calibrate only a specific increment, default=None

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''

    X1, X2 = X[:7], X[-1]
    tau1to6 = [0e-3, 0e-3, 0e-3, 0e-3, 0, 0]
    tau8to11 = [0e-3, 0, 0e-3, 0e-3]
    XX = numpy.hstack([X1, tau1to6, X2, tau8to11])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, increment=increment, stresses_to_include=['S', 'SIGMA']))
    else:
        return(numpy.hstack([X1, tau1to6, X2, tau8to11]))


def opti_options_4(X, Y, inputs, cal_norm, nu_targ, case, element, calibrate=True, increment=None):
    '''Objective function number 4 used for calibrating all 18 parameters of  micromorphic linear elasticity

    :param array-like X: Array of micromorphic linear elasticity parameters
    :param list Y: List storing dictionaries of DNS quantities for PK2, SIGMA, and M
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param str cal_norm: The form of the norm for the residual, use "L1" or "L2"
    :param float nu_targ: The targeted Poisson ratio if calibrating 2 parameter elasticity
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param int element: The macro (filter) element to calibration
    :param bool calibrate: A flag specifying whether to perform calibration for "True" or to return the stacked list of parameters for "False"
    :param int increment: An optional argumet to calibrate only a specific increment, default=None

    :returns: objective function evaluation by calling primary objective function if calibrate=True or return stacked list of parameters if calibrate=False
    '''

    XX = X
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, case, element, increment=increment, stresses_to_include=['S', 'SIGMA', 'M']))
    else:
        return(XX)

def handle_output_for_UQ(Xstore, Lstore, case):

    UQ_dict = {
        'obj':[],
        'lamb':[], 'mu':[], 'eta':[], 'tau':[], 'kappa':[], 'nu':[], 'sigma':[],
        'tau1':[], 'tau2':[], 'tau3':[], 'tau4':[], 'tau5':[], 'tau6':[], 'tau7':[],
        'tau8':[], 'tau9':[], 'tau10':[], 'tau11':[]}
    # Store results into dictionary
    for X, L in zip(Xstore, Lstore):

        UQ_dict['obj'].append(L)
        UQ_dict['lamb'].append(X[0])
        UQ_dict['mu'].append(X[1])
        UQ_dict['eta'].append(X[2])
        UQ_dict['tau'].append(X[3])
        UQ_dict['kappa'].append(X[4])
        UQ_dict['nu'].append(X[5])
        UQ_dict['sigma'].append(X[6])
        UQ_dict['tau1'].append(X[7])
        UQ_dict['tau2'].append(X[8])
        UQ_dict['tau3'].append(X[9])
        UQ_dict['tau4'].append(X[10])
        UQ_dict['tau5'].append(X[11])
        UQ_dict['tau6'].append(X[12])
        UQ_dict['tau7'].append(X[13])
        UQ_dict['tau8'].append(X[14])
        UQ_dict['tau9'].append(X[15])
        UQ_dict['tau10'].append(X[16])
        UQ_dict['tau11'].append(X[17])

    # remove zero entries depending on case
    remove = []
    if case == 1:
        remove = ['eta','tau','kappa','nu','sigma','tau1','tau2','tau3',
                  'tau4','tau5','tau6','tau7','tau8','tau9','tau10','tau11']
    elif case == 2:
        remove = ['tau1','tau2','tau3','tau4','tau5','tau6','tau7','tau8',
                  'tau9','tau10','tau11']
    elif case == 3:
        remove = ['tau1','tau2','tau3','tau4','tau5','tau6','tau8','tau9',
                  'tau10','tau11']
    for item in remove:
        UQ_dict.pop(item)

    return(UQ_dict)


def calibrate(input_file, output_file, case, Emod, nu, L, element=0, increment=None, plot_file=None, average=True, UQ_file=None):
    ''' Unpack DNS data and run calibration routine

    :param str input_file: The homogenized XDMF file output by the Micromorphic Filter
    :param str output_file:  The resulting list of parameters stored in a yaml file
    :param int case: The calibration "case". 1: two parameter, 2: 7 parameter, 3: 7 parameter plus tau7, 4: all 18 parameters
    :param float Emod: Estimate of a homogenized elastic modulus, used for initial parameter estimation
    :param float nu: Estimate of a homogenized Poisson ratio, used for initial parameter estimation
    :param float L: DNS max dimension (width, height, depth, etc.), used for initial parameter estimation
    :param int element: The macro (filter) element to calibration, default is zero
    :param int increment: An optional argument to calibrate only a specific increment

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
    
    # always average fields, but only for selected element
    cauchy      = average_quantities(cauchy, '3x3', element)
    symm        = average_quantities(symm, '3x3', element)
    PK2         = average_quantities(PK2, '3x3', element)
    SIGMA       = average_quantities(SIGMA, '3x3', element)
    #F           = average_quantities(F, '3x3')
    E           = average_quantities(E, '3x3', element)
    displacement    = average_quantities(displacement, '3', element)
    gradu       = average_quantities(gradu, '3x3', element)
    phi         = average_quantities(phi, '3x3', element)
    estrain     = average_quantities(estrain, '3x3', element)
    h           = average_quantities(h, '3x3', element)
    gradphi     = average_quantities(gradphi, '3x3x3', element)
 
    # store data for calibration
    Y = [PK2, SIGMA, M]
    inputs = [E, displacement, gradu, phi, gradphi, times]
    
    # get target nu from E
    nu_targ = numpy.average(-1*numpy.average([E[0][-1,0,0,0],E[0][-1,0,1,1]])/E[0][-1,0,2,2])
 
    # Estimate initial parameters
    param_est = initial_estimate(Emod, nu, L)
 
    # Define the elastic bounds
    parameter_bounds = [[0.0, 1.e5]] + [[-1.e5, 1.e5] for _ in range(6)] + [[-1.e5, 1.e5] for _ in range(11)]
    if len(elastic_parameter_ordering) != len(parameter_bounds):
        raise ValueError(f"The parameter bounds and the parameter names do not have the same length {len(parmaeter_bounds)} vs. {len(elastic_parameter_ordering)}")

    # calibrate!
    cal_norm = 'L1'
    maxit = 2000
    # TODO: streamline this workflow, very redundant
    # calibrate just lambda and mu
    if case == 1:
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
        res = scipy.optimize.differential_evolution(func=opti_options_1, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, case, element, True, increment))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_1(list(res.x), Y, inputs, cal_norm, nu_targ, case, element, calibrate=False)
    # calibrate first 7 parameters
    elif case == 2:
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      False, False, False, False, False] 
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_2, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, case, element, increment))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_2(res.x, Y, inputs, cal_norm, nu_targ, case, element, calibrate=False)
    # calibrate first 7 parameters and tau 7
    elif case == 3:
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      True, False, False, False, False] 
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_3, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, case, element, increment))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_3(res.x, Y, inputs, cal_norm, nu_targ, case, element, calibrate=False)
    # calibrate all parameters
    elif case == 4:
        param_mask = [True, True, True, True, True, True, True,
                      True, True, True, True, True, True, True, True, True, True, True]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_4, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, case, element, increment))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_4(res.x, Y, inputs, cal_norm, nu_targ, case, element, calibrate=False)
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

    # Manage Objective evaluation for UQ


    # plot resulting calibration
    if plot_file:
        print('plotting...')
        model_name=r'LinearElasticity'
        PK2_sim, SIGMA_sim, M_sim, SDVS_sim = evaluate_model(inputs, params, model_name, parameters_to_fparams, 0, element)
        PK2_sim = XRT.map_sim(PK2_sim, ninc)
        SIGMA_sim = XRT.map_sim(SIGMA_sim, ninc)
        cauchy_sim, symm_sim = XRT.get_current_configuration_stresses(PK2_sim, SIGMA_sim, inputs[2], inputs[3])

        plot_stresses(estrain, cauchy, cauchy_sim, f'{plot_file}_cauchy_fit_case_{case}.PNG', element)
        plot_stresses(estrain, symm, symm_sim, f'{plot_file}_symm_fit_case_{case}.PNG', element)

    # output parameters!
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
    parser.add_argument('--Emod', type=float,
        help="DNS elastic modulus, used for initial parameter estimation.")
    parser.add_argument('--nu', type=float,
        help="DNS Poisson's ratio, used for initial parameter estimation.")
    parser.add_argument('--L', type=float,
        help="DNS max dimension (width, height, depth, etc.), used for initial parameter estimation.")
    parser.add_argument('--element', type=int, default=0,
        help="The macro (filter) element to calibrate")
    parser.add_argument('--increment', type=int, required=False, default=None,
        help="An optional argument to callibrate only for specific increment")
    parser.add_argument('--case', type=int, required=True,
        help="Specify the calibration 'case'. 1: two parameter, 2: 7 parameter,\
              3: 7 parameter plus tau7, 4: all 18 parameters")
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
    sys.exit(calibrate(
                    input_file=args.input_file,
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
                    ))